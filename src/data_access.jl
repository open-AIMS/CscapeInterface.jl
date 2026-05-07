#=
Task 5: Julia access to outputs on yearly basis
=#

"""
    CscapeOutput

Container for C-scape simulation output.
"""
struct CscapeOutput
    out_array::Array{Float64,6}   # [years, sites, intervention, ft, enhancement, metrics]
    years::Vector{Int}
    site_ids::Vector{String}
    fts::Vector{String}
    kappa::Vector{Float64}
    area::Matrix{Float64}          # [sites, intervention]
    meshpoints::Array{Float64,2} # [ft, sizes]
    is_juvenile::Matrix{Bool}      # [ft, sizes]
    spatial::DataFrame      # Spatial data for sites
    metadata::Dict{String,Any}
end

Base.size(o::CscapeOutput) = size(o.out_array)
n_years(o::CscapeOutput) = length(o.years)
n_sites(o::CscapeOutput) = length(o.site_ids)
n_fts(o::CscapeOutput) = length(o.fts)


"""
    _build_area_matrix(spatial, coral_df, site_ids) -> Matrix{Float64}

Build the [n_sites × 3] area matrix where columns are non-intervened, intervened, and total area.
"""
function _build_area_matrix(
    spatial::DataFrame,
    coral_df::DataFrame,
    site_ids::Vector{String}
)::Matrix{Float64}

    total_area = Dict{String,Float64}(
        String(row[:reef_siteid]) => Float64(row[:area])
        for row in eachrow(spatial)
    )

    interv_area = Dict{String,Float64}()
    for row in eachrow(coral_df)
        sid = String(row[:reef_siteid])
        if !haskey(interv_area, sid)
            interv_area[sid] = Float64(row[:m2])
        end
    end

    n_sites  = length(site_ids)
    area_mat = zeros(Float64, n_sites, 3)
    for (i, sid) in enumerate(site_ids)
        tot      = get(total_area, sid, 0.0)
        inv      = get(interv_area, sid, 0.0)
        area_mat[i, 1] = tot - inv   # non-intervened
        area_mat[i, 2] = inv         # intervened
        area_mat[i, 3] = tot         # combined
    end
    return area_mat
end


"""
    _compute_combined_dimension!(out_array, area_matrix)

Fill intervention slot 3 of `out_array` in-place:
- Metrics 2:106 (size distributions): additive sum of slots 1 and 2.
- Metric 1 (cover): area-weighted combination, computed per site.
"""
function _compute_combined_dimension!(
    out_array::Array{Float64,6},
    area_matrix::Matrix{Float64}
)::Nothing

    # Size distributions: combined = non-intervened + intervened
    @views out_array[:, :, 3, :, :, 2:106] .=
        out_array[:, :, 1, :, :, 2:106] .+
        out_array[:, :, 2, :, :, 2:106]

    # Cover (metric 1): area-weighted per site
    n_sites = size(out_array, 2)
    for site in 1:n_sites
        tot      = area_matrix[site, 3]
        prop_non = tot > 0.0 ? area_matrix[site, 1] / tot : 1.0
        prop_int = tot > 0.0 ? area_matrix[site, 2] / tot : 0.0
        @views out_array[:, site, 3, :, :, 1] .=
            prop_non .* out_array[:, site, 1, :, :, 1] .+
            prop_int .* out_array[:, site, 2, :, :, 1]
    end
    return nothing
end


"""
    build_cscape_output(raw_out_array, years, site_ids, fts, spatial, coral_df,
                        meshpoints, scenario_id, draw, source_file) -> CscapeOutput

Construct a `CscapeOutput` from already-loaded Julia data.

`raw_out_array` must have 2 intervention dimensions (the raw form returned by `readRDS`).
All array manipulation (combined dimension, area weighting) is done here in Julia.

Use this when you have already loaded spatial/intervention/meshpoint data in memory
and want to avoid re-reading those files from disk.
"""
function build_cscape_output(
    raw_out_array::AbstractArray,
    years::Vector{Int},
    site_ids::Vector{String},
    fts::Vector{String},
    spatial::DataFrame,
    coral_df::DataFrame,
    meshpoints::Matrix{Float64},
    is_juvenile::AbstractMatrix{Bool},
    scenario_id::Int,
    draw::String,
    source_file::String
)::CscapeOutput

    n_years, n_sites, _, n_ft, n_enh, n_metrics = size(raw_out_array)

    # Allocate expanded 3-intervention array; copy raw slots 1 & 2, replacing missing → 0.0
    out_array = Array{Float64,6}(undef, n_years, n_sites, 3, n_ft, n_enh, n_metrics)
    @views out_array[:, :, 1:2, :, :, :] .= Float64.(coalesce.(raw_out_array, 0.0))

    area_mat = _build_area_matrix(spatial, coral_df, site_ids)
    _compute_combined_dimension!(out_array, area_mat)

    kappa    = spatial[!, :k] ./ 100.0
    metadata = Dict{String,Any}(
        "scenario_id" => scenario_id,
        "draw"        => draw,
        "source"      => source_file
    )

    @info "Built CscapeOutput" scenario=scenario_id size=size(out_array) years="$(years[1])-$(years[end])" sites=length(site_ids)

    return CscapeOutput(out_array, years, site_ids, fts, kappa, area_mat, meshpoints, Matrix{Bool}(is_juvenile), spatial, metadata)
end


"""
    build_cscape_output(raw_out_array, years, site_ids, fts, spatial, area_mat,
                        meshpoints, scenario_id, draw, source_file) -> CscapeOutput

Construct a `CscapeOutput` from already-loaded Julia data when the area matrix has
already been built (avoids re-reading `coral_df` from disk on repeated calls).
"""
function build_cscape_output(
    raw_out_array::AbstractArray,
    years::Vector{Int},
    site_ids::Vector{String},
    fts::Vector{String},
    spatial::DataFrame,
    area_mat::Matrix{Float64},
    meshpoints::Matrix{Float64},
    is_juvenile::AbstractMatrix{Bool},
    scenario_id::Int,
    draw::String,
    source_file::String
)::CscapeOutput

    n_years, n_sites, _, n_ft, n_enh, n_metrics = size(raw_out_array)

    out_array = Array{Float64,6}(undef, n_years, n_sites, 3, n_ft, n_enh, n_metrics)
    @views out_array[:, :, 1:2, :, :, :] .= Float64.(coalesce.(raw_out_array, 0.0))
    _compute_combined_dimension!(out_array, area_mat)

    kappa    = spatial[!, :k] ./ 100.0
    metadata = Dict{String,Any}(
        "scenario_id" => scenario_id,
        "draw"        => draw,
        "source"      => source_file
    )

    @info "Built CscapeOutput" scenario=scenario_id size=size(out_array) years="$(years[1])-$(years[end])" sites=length(site_ids)

    return CscapeOutput(out_array, years, site_ids, fts, kappa, area_mat, meshpoints, Matrix{Bool}(is_juvenile), spatial, metadata)
end


"""
    load_output(fpath::String, scenario_id::Int; draw="NA") -> CscapeOutput

Load C-scape output from simulation.

# Example
```julia
output = load_output("/path/to/data", 1)
```
"""
function load_output(fpath::String, scenario_id::Int;
                     draw::String = "NA",
                     output_path::String = "model_outputs",
                     filename::Union{String,Nothing} = nothing,
                     from_envir::Bool = false)

    scenario_file     = joinpath(fpath, "ScenarioID.xlsx")
    intervention_file = joinpath(fpath, "data", "Interventions$(scenario_id).RData")
    isfile(scenario_file) || error("Scenario file not found: $scenario_file")

    # --- extract output array ---
    local raw_out_array, dimnames_r, source_label
    if from_envir
        R"""
        .arr_env      <- MainEnvir$out_array
        .dimnames_env <- dimnames(MainEnvir$out_array)
        """
        raw_out_array = rcopy(R".arr_env")
        R"rm(.arr_env); invisible(gc())"
        dimnames_r  = rcopy(R".dimnames_env")
        R"rm(.dimnames_env); invisible(gc())"
        source_label = "MainEnvir\$out_array"
    else
        array_file = isnothing(filename) ?
            joinpath(fpath, output_path, "Array_scenario_$(scenario_id)_draw_$(draw).rds") :
            filename
        isfile(array_file) || error("Output file not found: $array_file")
        @rput array_file
        R"""
        out_array_R <- readRDS($array_file)
        dimnames_R  <- dimnames(out_array_R)
        """
        raw_out_array = rcopy(R"out_array_R")
        R"rm(out_array_R); invisible(gc())"
        dimnames_r  = rcopy(R"dimnames_R")
        R"rm(dimnames_R); invisible(gc())"
        source_label = array_file
    end

    # spatial / coral / meshpoints — always loaded from disk
    @rput fpath scenario_id scenario_file intervention_file
    R"""
    library(readxl)

    scenario_tbl   <- read_excel($scenario_file)
    scen_row       <- scenario_tbl[scenario_tbl[["ID"]] == $scenario_id, , drop = FALSE]
    spatial_file_R <- scen_row[["Spatial_file"]][1]
    demog_files_R  <- strsplit(scen_row[["Growth_Surv_file"]][1], "/")[[1]]

    spatial_R    <- readRDS(file.path($fpath, "data", spatial_file_R))
    interv_R     <- readRDS($intervention_file)
    coral_df_R   <- interv_R[["Coral"]]

    mesh_R <- array(0, dim = c(length(demog_files_R), 100))
    is_juvenile_R <- array(FALSE, dim = c(length(demog_files_R), 100))
    for (ft in seq_along(demog_files_R)) {
        demog_ft     <- readRDS(file.path($fpath, "data", demog_files_R[ft]))
        mesh_R[ft, ] <- demog_ft[["meshpoints_diam"]]
        juveniles <- as.vector(which(colSums(demog_ft[["fecundity_eggs"]])==0))
        is_juvenile_R[ft, juveniles] <- TRUE
    }
    """

    spatial    = rcopy(R"spatial_R")
    coral_df   = rcopy(R"coral_df_R")
    meshpoints = rcopy(R"mesh_R")
    is_juvenile = rcopy(R"is_juvenile_R")
    R"rm(spatial_R, interv_R, coral_df_R, mesh_R, is_juvenile_R); invisible(gc())"

    years    = parse.(Int, dimnames_r[:year])
    site_ids = String.(dimnames_r[:reef_siteid])
    fts      = String.(dimnames_r[:ft])

    return build_cscape_output(raw_out_array, years, site_ids, fts,
                                spatial, coral_df, meshpoints, is_juvenile,
                                scenario_id, draw, source_label)
end


"""
    load_ranking_inputs(fpath, scenario_id; array_file, from_envir)
        -> (cover, area_mat, spatial, site_ids, years, fts, kappa, meshpoints)

Load only what is needed for MCDA site ranking: area-weighted combined cover
`[years, sites, n_ft, n_enh]` plus spatial metadata.

Avoids loading the full 106-metric array. The cover slice (metric 1) is extracted
inside R before `rcopy`, so only ~1/106 of the data crosses the R-Julia boundary.

# Keywords
- `array_file`: Path to the `.rds` output file (required when `from_envir=false`)
- `from_envir`: If `true`, extract the cover slice directly from `MainEnvir\$out_array`
  already in the R session — no disk read at all
"""
function load_ranking_inputs(
    fpath::String,
    scenario_id::Int;
    array_file::Union{String,Nothing} = nothing,
    from_envir::Bool = false
)
    scenario_file     = joinpath(fpath, "ScenarioID.xlsx")
    intervention_file = joinpath(fpath, "data", "Interventions$(scenario_id).RData")
    isfile(scenario_file) || error("Scenario file not found: $scenario_file")

    @rput fpath scenario_id scenario_file intervention_file
    R"""
    library(readxl)

    scenario_tbl   <- read_excel($scenario_file)
    scen_row       <- scenario_tbl[scenario_tbl[["ID"]] == $scenario_id, , drop = FALSE]
    spatial_file_R <- scen_row[["Spatial_file"]][1]
    demog_files_R  <- strsplit(scen_row[["Growth_Surv_file"]][1], "/")[[1]]

    spatial_R  <- readRDS(file.path($fpath, "data", spatial_file_R))
    interv_R   <- readRDS($intervention_file)
    coral_df_R <- interv_R[["Coral"]]

    mesh_R <- array(0, dim = c(length(demog_files_R), 100))
    for (ft in seq_along(demog_files_R)) {
        demog_ft     <- readRDS(file.path($fpath, "data", demog_files_R[ft]))
        mesh_R[ft, ] <- demog_ft[["meshpoints_diam"]]
    }
    """

    if from_envir
        R"""
        cover_slice_R <- MainEnvir$out_array[, , , , , 1]
        dimnames_R    <- dimnames(MainEnvir$out_array)
        """
    else
        isnothing(array_file) && error("array_file must be provided when from_envir=false")
        isfile(array_file)    || error("Output file not found: $array_file")
        @rput array_file
        R"""
        arr_tmp       <- readRDS($array_file)
        cover_slice_R <- arr_tmp[, , , , , 1]
        dimnames_R    <- dimnames(arr_tmp)
        rm(arr_tmp); invisible(gc())
        """
    end

    cover_raw  = rcopy(R"cover_slice_R")
    R"rm(cover_slice_R); invisible(gc())"

    dimnames_r = rcopy(R"dimnames_R")
    spatial    = rcopy(R"spatial_R")
    coral_df   = rcopy(R"coral_df_R")
    meshpoints = rcopy(R"mesh_R")
    R"rm(dimnames_R, spatial_R, interv_R, coral_df_R, mesh_R); invisible(gc())"

    years    = parse.(Int, dimnames_r[:year])
    site_ids = String.(dimnames_r[:reef_siteid])
    fts      = String.(dimnames_r[:ft])

    area_mat = _build_area_matrix(spatial, coral_df, site_ids)

    # Area-weighted combined cover: [years, sites, n_ft, n_enh]
    n_years, n_sites, _, n_ft, n_enh = size(cover_raw)
    cover = Array{Float64,4}(undef, n_years, n_sites, n_ft, n_enh)
    for site in 1:n_sites
        tot      = area_mat[site, 3]
        prop_non = tot > 0.0 ? area_mat[site, 1] / tot : 1.0
        prop_int = tot > 0.0 ? area_mat[site, 2] / tot : 0.0
        @views cover[:, site, :, :] .=
            prop_non .* Float64.(coalesce.(cover_raw[:, site, 1, :, :], 0.0)) .+
            prop_int .* Float64.(coalesce.(cover_raw[:, site, 2, :, :], 0.0))
    end

    kappa = spatial[!, :k] ./ 100.0

    @info "Loaded ranking inputs" scenario=scenario_id from_envir=from_envir years="$(years[1])-$(years[end])" sites=length(site_ids)

    return cover, area_mat, spatial, site_ids, years, fts, kappa, meshpoints
end


"""
    get_yearly_data(output::CscapeOutput, year::Int) -> Array

Get all data for a specific year.

# Returns
Array [sites, intervention, ft, enhancement, metrics]
"""
function get_yearly_data(output::CscapeOutput, year::Int)
    idx = findfirst(==(year), output.years)
    if isnothing(idx)
        error("Year $year not found. Available: $(output.years[1])-$(output.years[end])")
    end
    return output.out_array[idx, :, :, :, :, :]
end


"""
    get_site_data(output::CscapeOutput, site_id) -> Array

Get time series for a specific site.

# Returns
Array [years, intervention, ft, enhancement, metrics]
"""
function get_site_data(output::CscapeOutput, site_id::String)
    idx = findfirst(==(site_id), output.site_ids)
    if isnothing(idx)
        error("Site '$site_id' not found")
    end
    return output.out_array[:, idx, :, :, :, :]
end

function get_site_data(output::CscapeOutput, site_idx::Int)
    return output.out_array[:, site_idx, :, :, :, :]
end


"""
    yearly_iterator(output::CscapeOutput)

Iterate through years with their data.

# Example
```julia
for (year, data) in yearly_iterator(output)
    mean_cover = mean(data[:, 1, :, :, 1])
    println("Year \$year: mean cover = \$mean_cover")
end
```
"""
function yearly_iterator(output::CscapeOutput)
    return ((year, output.out_array[idx, :, :, :, :, :]) 
            for (idx, year) in enumerate(output.years))
end


"""
    site_iterator(output::CscapeOutput)

Iterate through sites with their data.
"""
function site_iterator(output::CscapeOutput)
    return ((site_id, output.out_array[:, idx, :, :, :, :])
            for (idx, site_id) in enumerate(output.site_ids))
end


"""
    get_cover_timeseries(output::CscapeOutput; intervention_idx=3, ft_idx=nothing) -> Matrix

Get coral cover time series [years × sites].

# Keywords
- `intervention_idx::Int`: 1 = non-intervened, 2 = intervened, 3 = combined
- `ft_idx`: Functional type index (nothing = sum all)
"""
function get_cover_timeseries(output::CscapeOutput; 
                               intervention_idx::Int = 3, 
                               ft_idx::Union{Int,Nothing} = nothing)
    
    if isnothing(ft_idx)
        # Sum across FTs and enhancement
        cover = dropdims(sum(output.out_array[:, :, intervention_idx, :, :, 1], dims=(3,4)), dims=(3,4))
    else
        # Single FT, sum across enhancement
        cover = dropdims(sum(output.out_array[:, :, intervention_idx, ft_idx, :, 1], dims=3), dims=3)
    end
    
    return cover
end


"""
    get_size_distribution(output::CscapeOutput, year::Int, site_id; kwargs...) -> Vector

Get size class distribution (103 classes).
"""
function get_size_distribution(output::CscapeOutput, year::Int, site_id::String;
                                intervention_idx::Int = 3, ft_idx::Int = 1)
    
    year_idx = findfirst(==(year), output.years)
    site_idx = findfirst(==(site_id), output.site_ids)
    
    if isnothing(year_idx) || isnothing(site_idx)
        error("Year or site not found")
    end
    
    # Size classes are metrics 4:106, sum across enhancement
    return vec(sum(output.out_array[year_idx, site_idx, intervention_idx, ft_idx, :, 4:106], dims=1))
end


"""
    get_summary(output::CscapeOutput) -> DataFrame

Get summary statistics.
"""
function get_summary(output::CscapeOutput)
    cover = get_cover_timeseries(output)
    
    return DataFrame(
        year = output.years,
        mean_cover = vec(mean(cover, dims=2)),
        min_cover = vec(minimum(cover, dims=2)),
        max_cover = vec(maximum(cover, dims=2)),
        std_cover = vec(std(cover, dims=2))
    )
end
