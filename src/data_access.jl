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
    metadata::Dict{String,Any}
end

Base.size(o::CscapeOutput) = size(o.out_array)
n_years(o::CscapeOutput) = length(o.years)
n_sites(o::CscapeOutput) = length(o.site_ids)
n_fts(o::CscapeOutput) = length(o.fts)


"""
    load_output(fpath::String, scenario_id::Int; draw="NA") -> CscapeOutput

Load C-scape output from simulation.

# Example
```julia
output = load_output("/path/to/data", 1)
```
"""
function load_output(fpath::String, scenario_id::Int; draw::String = "NA")
    
    array_file = joinpath(fpath, "model_outputs", 
                          "Array_scenario_$(scenario_id)_draw_$(draw).rds")
    scenario_file = joinpath(fpath, "ScenarioID.xlsx")
    intervention_file = joinpath(fpath, "data", "Interventions$(scenario_id).RData")
    
    if !isfile(array_file)
        error("Output file not found: $array_file")
    end

    if !isfile(scenario_file)
        error("Scenario file not found: $scenario_file")
    end
    
    @rput array_file
    @rput intervention_file
    @rput fpath
    @rput scenario_id
    @rput scenario_file
    R"""
    library(readxl)

    # Load outputs and scenario table
    out_array_R <- readRDS($array_file)
    old_out_array <- out_array_R

    scenario_tbl <- read_excel($scenario_file)
    scenario_id_R <- scenario_tbl[scenario_tbl[["ID"]] == $scenario_id, , drop = FALSE]

    spatial_file_R <- scenario_id_R[["Spatial_file"]][1]
    demog_files_R <- strsplit(scenario_id_R[["Growth_Surv_file"]][1], "/")[[1]]

    # Load spatial & intervention data
    spatial_R <- readRDS(file.path($fpath, "data", spatial_file_R))
    intervention <- readRDS($intervention_file)
    coral_df <- intervention[["Coral"]]

    # Load demographic data for meshpoints
    meshpoints <- array(0, dim = c(length(demog_files_R), 100))
    for (ft in seq_along(demog_files_R)) {
        demog_R <- readRDS(file.path($fpath, "data", demog_files_R[ft]))
        meshpoints[ft, ] <- demog_R[["meshpoints_diam"]] # diameter in cm
    }

    dims_R <- dim(out_array_R)
    dimnames_R <- dimnames(out_array_R)

    # Calculate intervened=3 the combined dimension
    dims_R[3] <- 3
    dimnames_R[["intervened"]][3] <- "combined"
    out_array_R <- array(NA, dim = dims_R)
    out_array_R[, , -3, , , ] <- old_out_array

    idx_sizes <- seq(2, 106)
    out_array_R[, , 3, , , idx_sizes] <- out_array_R[, , 1, , , idx_sizes] + out_array_R[, , 2, , , idx_sizes]

    site_names_R <- unique(spatial_R[["reef_siteid"]])
    intervened_sites_R <- unique(coral_df[["reef_siteid"]])
    for (site in seq_along(site_names_R)) {
        if (site_names_R[site] %in% intervened_sites_R) {
            area <- spatial_R[["area"]][spatial_R[["reef_siteid"]] == site_names_R[site]]
            intervened_area <- coral_df[["m2"]][coral_df[["reef_siteid"]] == site_names_R[site]][1]
            prop_area <- c((area - intervened_area), intervened_area) / area
        } else {
            prop_area <- c(1, 0)
        }
        out_array_R[, , 3, , , 1] <- prop_area[1] * out_array_R[, , 1, , , 1] + prop_area[2] * out_array_R[, , 2, , , 1]
    }

    # Setup area
    area_R <- array(0, dim = c(length(site_names_R), 3))
    area_R[, 1] <- spatial_R[["area"]]
    for (site in intervened_sites_R) {
        area_R[site_names_R == site, 2] <- coral_df[["m2"]][coral_df[["reef_siteid"]] == site][1]
        area_R[site_names_R == site, 1] <- area_R[site_names_R == site, 1] - coral_df[["m2"]][coral_df[["reef_siteid"]] == site][1]
    }
    area_R[, 3] <- area_R[, 1] + area_R[, 2]

    """
    
    out_array = rcopy(R"out_array_R")
    dims = rcopy(R"dims_R")
    dimnames_r = rcopy(R"dimnames_R")
    spatial = rcopy(R"spatial_R")
    meshpoints = rcopy(R"meshpoints")       
    
    # RCall returns Symbol keys
    years = parse.(Int, dimnames_r[:year]) 
    site_ids = String.(dimnames_r[:reef_siteid])
    fts = String.(dimnames_r[:ft])
    
    # Load kappa and area from spatial file
    #kappa = ones(length(site_ids))  # Default, can be updated
    kappa = spatial[!, :k]/100
    area = rcopy(R"area_R")


    metadata = Dict{String,Any}(
        "scenario_id" => scenario_id,
        "draw" => draw,
        "source" => array_file
    )
    
    @info "Loaded output" scenario=scenario_id size=size(out_array) years="$(years[1])-$(years[end])" sites=length(site_ids)
    
    # Replace missing values with 0.0
    out_array_clean = replace(out_array, missing => 0.0)
    
    return CscapeOutput(Float64.(out_array_clean), years, site_ids, fts, kappa, area, meshpoints, metadata)
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
