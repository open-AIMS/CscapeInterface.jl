# =============================================================================
# MCDA_Parallel_Workflow.jl
#
# Multi-loop MCDA-simulation workflow.
#
# SHARED PHASE (main process, once):
#   - Setup R environment
#   - initialise_simulation + run_years!(SHARED_YEARS)
#   - Save output and serialize MainEnvir to disk as baseline
#
# INITIAL MCDA (main process, once):
#   - Rank sites from shared output → initial rankings used by all workers in loop 1
#
# ITERATION PHASE (one R session per worker):
#   - Each worker loads the saved baseline MainEnvir
#   - Repeats N_LOOPS times:
#     1. Select sites from current rankings
#     2. Build deployment schedule + accumulate into per-worker DataFrame
#     3. Apply deployment via modify_simulation_state!
#     4. Run LOOP_DURATION years
#     5. (if not last loop) Re-rank from live R session — no disk I/O
#   - Save single deployment CSV and single simulation RDS at end of all loops
#
# Usage:
#   include("src/MCDA_Parallel_Workflow.jl")
#   results = main()
# =============================================================================

using Distributed
using RCall
using DataFrames
using ADRIAMCDA
using CSV
using Printf
using TimerOutputs


# =============================================================================
# SHARED PHASE — runs once on main process
# =============================================================================

function run_shared_phase(; scenario_id, fpath, fun_path, shared_years, output_dirname)
    println("\n" * "="^60)
    println("SHARED PHASE: setup + run $(first(shared_years)):$(last(shared_years))")
    println("="^60)

    to = TimerOutput()

    @timeit to "setup_r_environment" CscapeInterface.setup_r_environment(fun_path)

    local env
    @timeit to "initialise_simulation" begin
        env = CscapeInterface.initialise_simulation(scenario_id, fpath; fun_path = fun_path)
        CscapeInterface.print_simulation_state(env)
    end

    local elapsed_run
    @timeit to "run_years" elapsed_run = @elapsed CscapeInterface.run_years!(env, shared_years)
    println("Shared years completed in $(round(elapsed_run; digits=1))s")

    out_dir = joinpath(fpath, "adria", output_dirname)
    mkpath(out_dir)
    y0, y1 = first(shared_years), last(shared_years)
    intermittent_file = joinpath(out_dir,
        "Array_scenario_$(scenario_id)_draw_NA_years$(y0)_$(y1).rds")

    @timeit to "finalise_output" CscapeInterface.finalise_simulation(env;
        filename        = intermittent_file,
        export_adria    = false,
        calc_indicators = false,
        keep_open       = true
    )
    println("Saved intermittent output: $intermittent_file")

    baseline_rds = joinpath(out_dir, "MainEnvir_baseline_scenario_$(scenario_id).rds")
    @rput baseline_rds
    @timeit to "save_baseline_rds" R"""
    saveRDS(MainEnvir, file = baseline_rds)
    cat("Baseline MainEnvir saved to:", baseline_rds, "\n")
    """
    println("Saved baseline MainEnvir: $baseline_rds")

    return env, intermittent_file, baseline_rds, to
end

# =============================================================================
# MCDA RANKING
# =============================================================================

# MCDA preferences — shared by main process and workers
const MCDA_PREFS = Dict(
    :names      => ["kappa", "depth", "waves", "relativecover"],
    :weights    => [1.0, 0.5, 0.8, 0.5],
    :directions => [maximum, maximum, minimum, minimum]
)

# Compute rankings from a loaded CscapeOutput (used on main process)
function compute_rankings(output, prefs::Dict)::DataFrame
    n_sites    = length(output.site_ids)
    n_criteria = length(prefs[:names])
    criteria   = Matrix{Union{Missing, Float64}}(missing, n_sites, n_criteria)

    criteria[:, 1] = Float64.(output.kappa)
    criteria[:, 2] = output.spatial[:, :depth_med]
    criteria[:, 3] = output.spatial[:, :ub_med]

    site1_vals          = output.out_array[:, 1, 3, :, :, 1]
    site1_summed        = vec(sum(site1_vals; dims=(2, 3)))
    last_nonzero_global = findlast(x -> x != 0, site1_summed)

    criteria[:, 4] = map(1:n_sites) do site
        if last_nonzero_global !== nothing
            site_vals = output.out_array[:, site, 3, :, :, 1]
            summed    = vec(sum(site_vals; dims=(2, 3)))
            summed[last_nonzero_global] / output.kappa[site]
        else
            0.0
        end
    end

    return DataFrame(
        site_id = output.site_ids,
        rank    = rank_locations(criteria, prefs),
        score   = rank_scores(criteria, prefs)
    )
end

# Rank sites from a pre-computed combined cover array (used for in-loop re-ranking)
function compute_rankings_from_cover(
    cover::Array{Float64,4},   # [years, sites, n_ft, n_enh] — area-weighted combined cover
    kappa::Vector{Float64},
    spatial::DataFrame,
    site_ids::Vector{String},
    prefs::Dict
)::DataFrame
    n_sites    = length(site_ids)
    n_criteria = length(prefs[:names])
    criteria   = Matrix{Union{Missing, Float64}}(missing, n_sites, n_criteria)

    criteria[:, 1] = kappa
    criteria[:, 2] = spatial[:, :depth_med]
    criteria[:, 3] = spatial[:, :ub_med]

    site1_vals          = cover[:, 1, :, :]
    site1_summed        = vec(sum(site1_vals; dims=(2, 3)))
    last_nonzero_global = findlast(x -> x != 0, site1_summed)

    criteria[:, 4] = map(1:n_sites) do site
        if last_nonzero_global !== nothing
            site_vals = cover[:, site, :, :]
            summed    = vec(sum(site_vals; dims=(2, 3)))
            summed[last_nonzero_global] / kappa[site]
        else
            0.0
        end
    end

    return DataFrame(
        site_id = site_ids,
        rank    = rank_locations(criteria, prefs),
        score   = rank_scores(criteria, prefs)
    )
end


# Load shared output, compute initial rankings, and return static spatial data for workers
function run_mcda_phase(fpath, scenario_id, output_dirname, shared_years)
    println("\n" * "="^60)
    println("MCDA PHASE: ranking sites")
    println("="^60)

    to = TimerOutput()

    y0, y1 = first(shared_years), last(shared_years)

    # Use MainEnvir already in R memory if the shared phase ran in this session
    # (saves the readRDS disk read entirely); otherwise read from the intermittent file.
    from_envir = rcopy(R"exists('MainEnvir') && !is.null(MainEnvir[['out_array']])")
    array_file = joinpath(fpath, "adria", output_dirname,
        "Array_scenario_$(scenario_id)_draw_NA_years$(y0)_$(y1).rds")

    if from_envir
        println("  Cover source: MainEnvir in R memory (no disk read)")
    else
        println("  Cover source: $array_file")
    end

    local cover, area_mat, cached_spatial, cached_site_ids,
          cached_years, cached_fts, cached_kappa, cached_meshpoints
    @timeit to "load_ranking_inputs" begin
        cover, area_mat, cached_spatial, cached_site_ids,
        cached_years, cached_fts, cached_kappa, cached_meshpoints =
            CscapeInterface.load_ranking_inputs(fpath, scenario_id;
                array_file = array_file, from_envir = from_envir)
    end

    local rankings
    @timeit to "compute_rankings" rankings = compute_rankings_from_cover(
        cover, cached_kappa, cached_spatial, cached_site_ids, MCDA_PREFS)
    println("Rankings computed for $(length(cached_site_ids)) sites")

    return rankings, cached_spatial, area_mat, cached_meshpoints,
           cached_years, cached_site_ids, cached_fts, cached_kappa, to
end

# =============================================================================
# BUILD ITERATION CONFIGS — packages shared data + variant params for each worker
# =============================================================================

function build_iteration_configs(;
        fpath, scenario_id, output_dirname, baseline_rds, fun_path,
        n_loops, loop_duration, shared_years_end, rankings, variants,
        cached_spatial, cached_area, cached_meshpoints,
        cached_years, cached_site_ids, cached_fts, cached_kappa)

    out_dir = joinpath(fpath, "adria", output_dirname)
    mkpath(out_dir)

    cfgs = Dict{String, Any}[]
    for (i, variant) in enumerate(variants)
        out_file = joinpath(out_dir, "Array_scenario_$(scenario_id)_iter_$(i).rds")
        push!(cfgs, Dict{String, Any}(
            "id"                => i,
            "scenario_id"       => scenario_id,
            "fpath"             => fpath,
            "fun_path"          => fun_path,
            "baseline_rds"      => baseline_rds,
            "n_loops"           => n_loops,
            "loop_duration"     => loop_duration,
            "shared_years_end"  => shared_years_end,
            "output_file"       => out_file,
            "output_dirname"    => output_dirname,
            "rankings"          => rankings,
            "variant"           => variant,
            "cached_spatial"    => cached_spatial,
            "cached_area"       => cached_area,
            "cached_meshpoints" => cached_meshpoints,
            "cached_years"      => cached_years,
            "cached_site_ids"   => cached_site_ids,
            "cached_fts"        => cached_fts,
            "cached_kappa"      => cached_kappa
        ))
    end

    return cfgs
end

# =============================================================================
# WORKER CODE — defined after addprocs so all workers receive it
# Each worker has its own R session. No shared R state.
# =============================================================================

function define_worker_code!()
    expr = quote
        using CscapeInterface
        using RCall
        using DataFrames
        using ADRIAMCDA
        using CSV
        using TimerOutputs

        function worker_run_iteration(iter_cfg::Dict{String, Any})
            id               = iter_cfg["id"]
            scenario_id      = iter_cfg["scenario_id"]
            fpath            = iter_cfg["fpath"]
            fun_path         = iter_cfg["fun_path"]
            baseline_rds     = iter_cfg["baseline_rds"]
            n_loops          = iter_cfg["n_loops"]
            loop_duration    = iter_cfg["loop_duration"]
            shared_years_end = iter_cfg["shared_years_end"]
            output_file      = iter_cfg["output_file"]
            output_dirname   = iter_cfg["output_dirname"]
            rankings         = iter_cfg["rankings"]
            variant          = iter_cfg["variant"]
            cached_spatial    = iter_cfg["cached_spatial"]
            cached_area       = iter_cfg["cached_area"]
            cached_meshpoints = iter_cfg["cached_meshpoints"]
            cached_site_ids   = iter_cfg["cached_site_ids"]
            cached_kappa      = iter_cfg["cached_kappa"]
            n_sites_select   = get(variant, "n_sites", 0)
            selection        = variant["selection"]
            is_counterfactual = selection == "counterfactual"

            out_dir = dirname(output_file)

            to = TimerOutput()

            # Source C-scape R functions
            @timeit to "R_setup" R"""
            FROM_JULIA <- TRUE
            setwd($fun_path)
            source('main.R')
            source('cscape_sim.R')
            source('intervention_setup.R')
            source('ancillary_functions.R')
            source('annual_site_taxa_loops.R')
            source('ipm_pred.R')

            modular_path <- file.path($fun_path, "modules")
            if (dir.exists(modular_path)) {
                for (f in list.files(modular_path, pattern = "\\.R$", full.names = TRUE)) source(f)
            }
            for (f in list.files(file.path($fun_path, "1_toolbox"), pattern = "\\.R$", full.names = TRUE)) source(f)

            ENABLE_PARALLEL <- FALSE
            .fun_path <- $fun_path
            setup_julia_batch(file.path($fun_path, "julia/src"))
            cat("Worker", Sys.getpid(), "functions loaded\n")
            """

            # Load shared baseline state (SHARED_YEARS already run)
            @timeit to "R_load_baseline" R"""
            MainEnvir <- readRDS($baseline_rds)
            cat("Worker", Sys.getpid(), "baseline state loaded\n")

            if (isTRUE(MainEnvir$use_julia_batch)) {
                BATCH_IPM_CACHE   <<- NULL
                BATCH_IPM_N_SITES <<- length(MainEnvir$site_names)
                create_ipm_original <<- create_ipm
                cat("Worker", Sys.getpid(), "batch globals restored: BATCH_IPM_N_SITES =", BATCH_IPM_N_SITES, "\n")
            }
            """

            # Build Julia env handle spanning all loops
            total_start = shared_years_end + 1
            total_end   = shared_years_end + n_loops * loop_duration
            env = Dict{String, Any}(
                "scenario_id"  => scenario_id,
                "fpath"        => fpath,
                "fun_path"     => fun_path,
                "year_start"   => total_start,
                "year_end"     => total_end,
                "current_year" => total_start,
                "initialised"  => true
            )

            prefs = CscapeInterface.MCDA_PREFS

            # Accumulate deployment rows across all loops (written once at end)
            all_deployments = DataFrame[]

            for loop_idx in 1:n_loops
                loop_start = shared_years_end + (loop_idx - 1) * loop_duration + 1
                loop_end   = shared_years_end + loop_idx * loop_duration
                loop_years = loop_start:loop_end

                println("Worker $(myid()): iteration $id loop $loop_idx/$(n_loops) years $(loop_start):$(loop_end)")

                @timeit to "loop_$loop_idx" begin

                    if is_counterfactual
                        if loop_idx == 1
                            @timeit to "clear_deployment" R"""
                            MainEnvir$CoralIntervention <- data.frame(
                                reef_siteid = character(),
                                Year = integer(),
                                ft = integer(),
                                no_int_corals = numeric(),
                                proportion = numeric(),
                                m2 = numeric(),
                                density = numeric(),
                                meshpt_int_corals = character(),
                                Enhancement = integer(),
                                scenario = integer(),
                                stringsAsFactors = FALSE
                            )
                            MainEnvir$Intervened_sites <- character(0)
                            cat("Worker", Sys.getpid(), "counterfactual deployment cleared\n")
                            """
                        end
                    else
                        # Apply selection rule to get sites for this loop
                        n_total = nrow(rankings)
                        deploy_sites = if selection == "top"
                            rankings.site_id[rankings.rank .<= n_sites_select]
                        elseif selection == "bottom"
                            rankings.site_id[rankings.rank .> (n_total - n_sites_select)]
                        else
                            error("Unknown selection rule '$selection': use \"top\", \"bottom\", or \"counterfactual\"")
                        end
                        println("Worker $(myid()): loop $loop_idx — $selection $(length(deploy_sites)) sites")

                        site_year = [(site, yr) for yr in loop_years for site in deploy_sites]
                        n_deploy  = length(site_year)
                        prop_val  = get(variant, "proportion", 1.0)
                        dens_val  = get(variant, "density", variant["no_int_corals"] / variant["m2"])

                        deployment = Dict{String, Any}(
                            "reef_siteid"       => [p[1] for p in site_year],
                            "Year"              => [p[2] for p in site_year],
                            "ft"                => fill(variant["ft"], n_deploy),
                            "no_int_corals"     => fill(variant["no_int_corals"], n_deploy),
                            "proportion"        => fill(prop_val, n_deploy),
                            "m2"                => fill(variant["m2"], n_deploy),
                            "density"           => fill(dens_val, n_deploy),
                            "meshpt_int_corals" => fill(variant["meshpt_int_corals"], n_deploy),
                            "Enhancement"       => fill(variant["Enhancement"], n_deploy)
                        )

                        df = DataFrame(deployment)
                        df[!, :loop] .= loop_idx
                        push!(all_deployments, df)

                        @timeit to "modify_state" CscapeInterface.modify_simulation_state!(env; coral_deployment=deployment)
                    end

                    local elapsed_run
                    @timeit to "run_years" elapsed_run = @elapsed CscapeInterface.run_years!(env, loop_years)
                    println("Worker $(myid()): loop $loop_idx completed in $(round(elapsed_run; digits=1))s")

                    # Re-rank for next loop: cover-only extraction (no full array transfer)
                    if loop_idx < n_loops && !is_counterfactual
                        @timeit to "rerank" begin
                            local cover
                            @timeit to "get_cover" cover =
                                CscapeInterface.get_combined_cover_for_ranking(env, cached_area)
                            @timeit to "compute_rankings" rankings =
                                compute_rankings_from_cover(cover, cached_kappa, cached_spatial, cached_site_ids, prefs)
                            cover = nothing
                        end
                        GC.gc()
                        println("Worker $(myid()): loop $loop_idx re-ranked $(nrow(rankings)) sites")
                    end

                end # @timeit loop_$loop_idx
            end

            # Save deployment DataFrame — single file covering all loops
            if !is_counterfactual && !isempty(all_deployments)
                deployment_file = joinpath(out_dir, "deployment_iter_$(id).csv")
                @timeit to "save_deployment_csv" CSV.write(deployment_file, vcat(all_deployments...))
                println("Worker $(myid()): saved deployment → $deployment_file")
            end

            # Save simulation output — single file at end of all loops
            @timeit to "finalise_output" CscapeInterface.finalise_simulation(env;
                filename        = output_file,
                export_adria    = true,
                calc_indicators = true
            )
            println("Worker $(myid()): saved simulation → $output_file")

            sim_state = CscapeInterface.get_simulation_state(env)

            # Free MainEnvir from R before this worker process may exit.
            # Without this, rmprocs() triggers RCall's atexit cleanup of a large
            # R session, which causes EXCEPTION_ACCESS_VIOLATION on Windows.
            try
                R"""
                if (exists("MainEnvir", envir = .GlobalEnv)) {
                    rm(list = "MainEnvir", envir = .GlobalEnv)
                }
                invisible(gc()); invisible(gc())
                """
            catch
            end

            return Dict{String, Any}(
                "iteration_id" => id,
                "worker_pid"   => myid(),
                "n_loops"      => n_loops,
                "timer_str"    => sprint(show, to),
                "current_year" => get(sim_state, "current_year", missing),
                "kappa_mean"   => get(sim_state, "kappa_mean", missing),
                "plasticity"   => get(sim_state, "plasticity", missing),
                "dhw_enhance"  => get(sim_state, "dhw_enhance", missing),
                "output_file"  => output_file
            )
        end
    end

    eval(expr)
    if !isempty(workers())
        Distributed.remotecall_eval(Main, workers(), :(using CscapeInterface))
        Distributed.remotecall_eval(Main, workers(), :(using RCall))
        Distributed.remotecall_eval(Main, workers(), :(using DataFrames))
        Distributed.remotecall_eval(Main, workers(), :(using ADRIAMCDA))
        Distributed.remotecall_eval(Main, workers(), :(using CSV))
        Distributed.remotecall_eval(Main, workers(), :(using TimerOutputs))
        Distributed.remotecall_eval(Main, workers(), expr)
    end
end

# =============================================================================
# MAIN
# =============================================================================

function run_dynamic_reranking(cfg::Dict{String, Any})
    scenario_id    = cfg["scenario_id"]
    fpath          = cfg["fpath"]
    fun_path       = cfg["fun_path"]
    output_dirname = cfg["output_dirname"]
    n_workers      = cfg["n_workers"]
    shared_years   = cfg["shared_years"]
    n_loops        = cfg["n_loops"]
    loop_duration  = cfg["loop_duration"]
    variants       = cfg["variants"]
    skip_shared    = get(cfg, "skip_shared", false)

    main_to    = TimerOutput()
    main_start = time()

    # 1. Shared phase — skip if baseline files already exist or skip_shared=true
    out_dir  = joinpath(fpath, "adria", output_dirname)
    y0, y1   = first(shared_years), last(shared_years)
    baseline_rds       = joinpath(out_dir, "MainEnvir_baseline_scenario_$(scenario_id).rds")
    intermittent_file  = joinpath(out_dir, "Array_scenario_$(scenario_id)_draw_NA_years$(y0)_$(y1).rds")

    shared_to = nothing
    if skip_shared || (isfile(baseline_rds) && isfile(intermittent_file))
        println("\n" * "="^60)
        println("SHARED PHASE: skipped — using existing baseline")
        println("  $(basename(baseline_rds))")
        println("  $(basename(intermittent_file))")
        println("="^60)
    else
        @timeit main_to "shared_phase" begin
            _, _, baseline_rds, shared_to = run_shared_phase(;
                scenario_id    = scenario_id,
                fpath          = fpath,
                fun_path       = fun_path,
                shared_years   = shared_years,
                output_dirname = output_dirname
            )
        end
    end

    # 2. Initial MCDA ranking from shared output — also cache static fields for workers
    local rankings, cached_spatial, cached_area, cached_meshpoints,
          cached_years, cached_site_ids, cached_fts, cached_kappa, mcda_to
    @timeit main_to "mcda_phase" begin
        rankings, cached_spatial, cached_area, cached_meshpoints,
        cached_years, cached_site_ids, cached_fts, cached_kappa, mcda_to =
            run_mcda_phase(fpath, scenario_id, output_dirname, shared_years)
    end

    # 3. Start workers and push worker code
    @timeit main_to "worker_setup" begin
        current = nprocs() - 1   # nworkers() returns min 1 even with no workers
        if current < n_workers
            addprocs(n_workers - current)
        end
        define_worker_code!()
    end

    # 4. Build per-iteration configs
    iter_cfgs = build_iteration_configs(;
        fpath             = fpath,
        scenario_id       = scenario_id,
        output_dirname    = output_dirname,
        baseline_rds      = baseline_rds,
        fun_path          = fun_path,
        n_loops           = n_loops,
        loop_duration     = loop_duration,
        shared_years_end  = last(shared_years),
        rankings          = rankings,
        variants          = variants,
        cached_spatial    = cached_spatial,
        cached_area       = cached_area,
        cached_meshpoints = cached_meshpoints,
        cached_years      = cached_years,
        cached_site_ids   = cached_site_ids,
        cached_fts        = cached_fts,
        cached_kappa      = cached_kappa
    )

    # 5. Distribute iterations across workers (process 1 participates alongside workers)
    all_procs = WorkerPool(vcat(1, workers()))
    println("\n" * "="^60)
    println("ITERATION PHASE: $(length(iter_cfgs)) variants × $n_loops loops on $(nprocs()) processes")
    println("="^60)
    local results
    @timeit main_to "iterations_pmap" results = pmap(worker_run_iteration, all_procs, iter_cfgs; batch_size=1)

    total_wall = time() - main_start

    println("\n" * "="^60)
    println("RESULTS")
    println("="^60)
    for r in results
        r_display = filter(kv -> kv.first != "timer_str", r)
        println(r_display)
    end

    # ── Main-process phase timer ──────────────────────────────────────────────
    println("\n" * "="^60)
    println("MAIN PROCESS — phase timings")
    println("="^60)
    show(main_to)
    println()

    if shared_to !== nothing
        println("\n  Shared phase detail:")
        show(shared_to)
        println()
    end

    println("\n  MCDA phase detail:")
    show(mcda_to)
    println()

    # ── Per-worker function-level timers ──────────────────────────────────────
    println("\n" * "="^60)
    println("WORKER TIMERS (function breakdown per variant)")
    println("="^60)
    for r in results
        if haskey(r, "timer_str")
            println("\n  — Variant $(r["iteration_id"]) (worker $(r["worker_pid"])) —")
            println(r["timer_str"])
        end
    end
    println("="^60)
    @printf "\nTotal wall time: %.1f s\n" total_wall

    return results
end
