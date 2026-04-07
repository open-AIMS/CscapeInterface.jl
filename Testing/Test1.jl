using Pkg
Pkg.activate("c:/Users/vhaller/Documents/GitHub/CscapeAccess")
cd("c:/Users/vhaller/Documents/GitHub/CscapeAccess")
include("c:/Users/vhaller/Documents/GitHub/CscapeAccess/src/CscapeJulia.jl")
using .CscapeJulia
using RCall


#Setup the R environment to a separate directory for Cscape
setup_r_environment("C:/users/vhaller/Documents/Github/C_scape")


fpath = "C:/Users/vhaller/OneDrive - Australian Institute of Marine Science/IPMF/InterventionStudy2024_TestingPascal/ChecksTEDRS"
using BenchmarkTools
@elapsed CscapeJulia.run_cscape(19, fpath, fun_path="C:/users/vhaller/Documents/Github/C_scape")