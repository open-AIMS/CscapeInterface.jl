using Pkg
Pkg.activate("c:/Users/vhaller/Documents/GitHub/CscapeAccess")
cd("c:/Users/vhaller/Documents/GitHub/CscapeAccess")
include("c:/Users/vhaller/Documents/GitHub/CscapeAccess/src/CscapeJulia.jl")
using .CscapeJulia
using RCall
using DataFrames
using ADRIAMCDA

#Setup the R environment to a separate directory for Cscape
CscapeJulia.setup_r_environment("C:/users/vhaller/Documents/Github/C_scape")

#Setup data location
fpath = "C:/Users/vhaller/OneDrive - Australian Institute of Marine Science/IPMF/InterventionStudy2024_TestingPascal/ChecksTEDRS"
fun_path ="C:/users/vhaller/Documents/Github/C_scape"

# Initialise (loads data, creates arrays, NO years run yet)
env = CscapeJulia.initialise_simulation(19, fpath,fun_path=fun_path)
# Check state before running
CscapeJulia.print_simulation_state(env)

# Run first half - Multiple years
@elapsed CscapeJulia.run_years!(env, 2008:2024)

#Save intermittent results
custom_folder_name = "Test_Workflow"
intermittent_dir = joinpath(fpath, "adria", custom_folder_name)
mkpath(intermittent_dir)

intermittent_file = joinpath(intermittent_dir, "Array_scenario_$(env["scenario_id"])_draw_NA_years2008_2024.rds")
CscapeJulia.finalise_simulation(env, filename=intermittent_file, export_adria=false, calc_indicators=false, keep_open=true)


#Investigate outputs
output = CscapeJulia.load_output(fpath, env["scenario_id"], draw = "NA_years2008_2024",output_path = "adria/$custom_folder_name")


# MCDA site ranking

# This should be adjusted to reflect the system being assessed
prefs = Dict(
    # Name of each criteria
    :names => ["kappa", "depth", "waves", "relativecover"],

    # The desired relative weight to place (0 - 1, higher means more important)
    # Gets normalized so their sum is 0 - 1
    :weights => [1.0, 0.5, 0.8, 0.5],

    # Desired directionality - to prefer lower or higher criteria values
    :directions => [maximum, maximum, minimum, minimum]
)

n_sites = length(output.site_ids)
n_criteria = length(prefs[:names])
criteria = Matrix{Union{Missing, Float64}}(missing, n_sites, n_criteria)
# Fill by numeric column index
criteria[:, 1] = Float64.(output.kappa)
criteria[:, 2] = output.spatial[:, :depth_med]
criteria[:, 3] = output.spatial[:, :ub_med]
# Find the last year with any non-zero values (use site 1 only for speed)
site1_vals = output.out_array[:, 1, 3, :, :, 1]
site1_summed = vec(sum(site1_vals, dims=(2,3)))
last_nonzero_global = findlast(x -> x != 0, site1_summed)

criteria[:, 4] = map(1:n_sites) do site
    # Use the same last non-zero year for all sites
    if last_nonzero_global !== nothing
        site_vals = output.out_array[:, site, 3, :, :, 1]
        summed_by_year = vec(sum(site_vals, dims=(2,3)))
        summed_by_year[last_nonzero_global] / output.kappa[site]
    else
        0.0  # If all zeros, return 0
    end
end

# DEBUG: Check criteria matrix for NaN or missing values
println("Criteria matrix summary:")
for i in 1:n_criteria
    col_data = criteria[:, i]
    n_nan = count(isnan, col_data)
    n_missing = count(ismissing, col_data)
    valid_count = count(x -> !isnan(x) && !ismissing(x), col_data)
    println("Column $i ($(prefs[:names][i])): $valid_count valid, $n_nan NaN, $n_missing missing")
    if n_nan > 0 || n_missing > 0
        println("  Values: $(col_data)")
    end
end

#Calculate rankings
rankings = DataFrame(
    site_id = output.site_ids,
    rank = rank_locations(criteria, prefs),
    score = rank_scores(criteria, prefs)
)


#Setup intervention for the top 5 ranking sites
deploy_sites = rankings.site_id[rankings.rank .<= 5]   
deploy_years = 2025:2029

# Expand to one row per site-year combination
site_year = [(site, year) for year in deploy_years for site in deploy_sites]
n = length(site_year)
deploy_site_ids = [p[1] for p in site_year]
deploy_year_vals = [p[2] for p in site_year]

rcopy(R"MainEnvir$CoralIntervention")
new_deployment = Dict(
    "reef_siteid"       => deploy_site_ids,
    "Year"              => deploy_year_vals,
    "ft"                => fill(1, n),          # functional type 1
    "no_int_corals"     => fill(1000.0, n),     # 1000 corals per site
    "m2"                => fill(500.0, n),       # 500 m² deployment area
    "meshpt_int_corals" => fill("4", n),         # size class
    "Enhancement"       => fill(3, n)            # enhancement class 3
)

CscapeJulia.modify_simulation_state!(env; coral_deployment=new_deployment)

rcopy(R"MainEnvir$CoralIntervention")

# Run intervention
@elapsed CscapeJulia.run_years!(env, 2025:2030)

#Save intermittent results
custom_folder_name = "Test_Workflow"
intermittent_dir = joinpath(fpath, "adria", custom_folder_name)
mkpath(intermittent_dir)

intermittent_file = joinpath(intermittent_dir, "Array_scenario_$(env["scenario_id"])_draw_NA_years2025_2030.rds")
CscapeJulia.finalise_simulation(env, filename=intermittent_file, export_adria=false, calc_indicators=true)