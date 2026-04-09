using Pkg
Pkg.activate("c:/Users/vhaller/Documents/GitHub/CscapeAccess")
cd("c:/Users/vhaller/Documents/GitHub/CscapeAccess")
include("c:/Users/vhaller/Documents/GitHub/CscapeAccess/src/CscapeJulia.jl")
using .CscapeJulia
using RCall
using DataFrames


#Setup the R environment to a separate directory for Cscape
setup_r_environment("C:/users/vhaller/Documents/Github/C_scape")


fpath = "C:/Users/vhaller/OneDrive - Australian Institute of Marine Science/IPMF/InterventionStudy2024_TestingPascal/ChecksTEDRS"
using BenchmarkTools
@elapsed CscapeJulia.run_cscape(19, fpath, fun_path="C:/users/vhaller/Documents/Github/C_scape")



#Investigate outputs
output = CscapeJulia.load_output(fpath, 19)
cover = CscapeJulia.get_cover_timeseries(output)



# MCDA site ranking

using ADRIAMCDA

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
criteria[:, 4] = sum(output.out_array[19, :, 3, :, :, 1], dims=(2,3))[:, 1] ./ output.kappa


#Calculate rankings
rankings = DataFrame(
    site_id = output.site_ids,
    rank = rank_locations(criteria, prefs),
    score = rank_scores(criteria, prefs)
)
