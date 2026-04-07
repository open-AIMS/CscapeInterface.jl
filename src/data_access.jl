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
    
    if !isfile(array_file)
        error("Output file not found: $array_file")
    end
    
    @rput array_file
    R"""
    out_array_R <- readRDS($array_file)
    dims_R <- dim(out_array_R)
    dimnames_R <- dimnames(out_array_R)
    """
    
    out_array = rcopy(R"out_array_R")
    dims = rcopy(R"dims_R")
    dimnames_r = rcopy(R"dimnames_R")
    
    # RCall returns Symbol keys
    years = parse.(Int, dimnames_r[:year])
    site_ids = String.(dimnames_r[:reef_siteid])
    fts = String.(dimnames_r[:ft])
    
    # Load kappa from spatial file
    kappa = ones(length(site_ids))  # Default, can be updated
    
    metadata = Dict{String,Any}(
        "scenario_id" => scenario_id,
        "draw" => draw,
        "source" => array_file
    )
    
    @info "Loaded output" scenario=scenario_id size=size(out_array) years="$(years[1])-$(years[end])" sites=length(site_ids)
    
    # Replace missing values with 0.0
    out_array_clean = replace(out_array, missing => 0.0)
    
    return CscapeOutput(Float64.(out_array_clean), years, site_ids, fts, kappa, metadata)
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
    get_cover_timeseries(output::CscapeOutput; intervention_idx=1, ft_idx=nothing) -> Matrix

Get coral cover time series [years × sites].

# Keywords
- `intervention_idx::Int`: 1 = non-intervened, 2 = intervened
- `ft_idx`: Functional type index (nothing = sum all)
"""
function get_cover_timeseries(output::CscapeOutput; 
                               intervention_idx::Int = 1, 
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
                                intervention_idx::Int = 1, ft_idx::Int = 1)
    
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
