#=
Task 6: Connect with ADRIAIndicators
=#

using ADRIAIndicators
using JLD2

"""
    to_adria_format(output::CscapeOutput; intervention_idx=1) -> NamedTuple

Convert CscapeOutput to ADRIAIndicators format.

# Returns
NamedTuple with:
- `cover`: Array [timesteps, groups, sizes, locations]
- `habitable_area`: Vector [locations]
- `reef_area`: Vector [locations]
- `area`: Area data
- `meshpoints`: Mesh points data
"""

function to_adria_format(output::CscapeOutput; intervention_idx::Int = 3)
    
    n_years = length(output.years)
    n_sites =length(output.site_ids)
    n_fts = length(output.fts)
    n_sizes = 100  # Dimension level 7:106, size classes 6:105
    
    # Extract absolute coral cover in m2: [years, sites, int, ft, enh, metrics] -> [timesteps, groups, sizes, locations]
    absolutecover = zeros(Float64, n_years, n_fts, n_sizes, n_sites)
    
    for t in 1:n_years, g in 1:n_fts, s in 1:n_sizes, l in 1:n_sites
        # Convert diameter (cm) to area (m²): π * (diameter/200)²
        diameter_cm = output.meshpoints[g, s]
        area_m2 = π * (diameter_cm / 200.0)^2
        absolutecover[t, g, s, l] = area_m2 * sum(output.out_array[t, l, intervention_idx, g, :, 6+s])
    end
    
    # Adjust cover to be relative to habitable area (kappa)
    relativecover = zeros(Float64, n_years, n_fts, n_sizes, n_sites)
    
    for l in 1:n_sites
        if output.area[l, intervention_idx] > 0
            relativecover[:, :, :, l] = absolutecover[:, :, :, l] / (output.area[l, intervention_idx] * output.kappa[l])
        end
    end
    
    clamp!(relativecover, 0.0, 1.0)
    
    return (
        cover = relativecover,
        habitable_area = Float64.(output.kappa*output.area[:, intervention_idx]),
        reef_area = Float64.(output.area[:, intervention_idx]),
        area = output.area[:, intervention_idx],
        meshpoints = output.meshpoints
    )
end


"""
    calculate_indicators(output::CscapeOutput; juvenile_threshold=5) -> Dict

Calculate reef indicators using ADRIAIndicators.

# Returns
Dict with indicator arrays.
"""
function calculate_indicators(output::CscapeOutput; )
    
    adria = to_adria_format(output)
    cover = adria.cover
    habitable_area = adria.habitable_area
    reef_area = adria.reef_area

    colony_mean_diam_cm = mean(output.meshpoints)  # Approximate mean diameter for shelter volume calculation
    #planar_area_params
    #reference
    #is_juvenile

    n_sizes = size(cover, 3)
    
    results = Dict{String, Any}()
    
    # Calculate cover metrics
    results["relative_cover"] = ADRIAIndicators.relative_cover(cover)
    results["relative_loc_taxa_cover"] = ADRIAIndicators.relative_loc_taxa_cover(cover)
    results["relative_taxa_cover"] = ADRIAIndicators.relative_taxa_cover(cover, habitable_area)
    results["ltmp_cover"] = ADRIAIndicators.relative_cover(cover, habitable_area, reef_area)

    
    # Calculate metrics
    results["relative_shelter_volumne"] = ADRIAIndicators.relative_shelter_volume(cover,colony_mean_diam_cm,planar_area_params, habitable_area, reference)
    results["coral_diversity"] = ADRIAIndicators.coral_diversity(results["relative_loc_taxa_cover"])
    results["coral_evenness"] = ADRIAIndicators.coral_evenness(results["relative_loc_taxa_cover"])   
    results["relative_juveniles"] = ADRIAIndicators.relative_juveniles(cover, is_juvenile)
    results["juvenile_indicator"] = ADRIAIndicators.juvenile_indicator(cover, is_juvenile, habitable_area, colony_mean_diam_cm,15)

    #Calculate reef indices
    results["reef_biodiversity_condition_index"] = ADRIAIndicators.reef_biodiversity_condition_index(results["relative_loc_taxa_cover"],results["coral_diversity"], results["relative_shelter_volumne"])
    results["reef_condition_index"] = ADRIAIndicators.reef_condition_index(results["ltmp_cover"], results["relative_shelter_volumne"], results["juvenile_indicator"])
    results["reef_fish_index"] = ADRIAIndicators.reef_fish_index(results["relative_cover"])
    results["reef_tourism_index"] = ADRIAIndicators.reef_tourism_index_no_rubble(results["ltmp_cover"], results["coral_evenness"], results["relative_shelter_volumne"], results["relative_juveniles"])
    
    

    # Metadata
    results["years"] = output.years
    results["site_ids"] = output.site_ids
    results["fts"] = output.fts
    
    @info "Calculated indicators" n_years=length(output.years) n_sites=length(output.site_ids)
    
    return results
end


"""
    indicator_summary(results::Dict; year=nothing)

Print summary of indicators.
"""
function indicator_summary(results::Dict; year::Union{Int,Nothing} = nothing)
    years = results["years"]
    year_idx = isnothing(year) ? length(years) : findfirst(==(year), years)
    
    println("=" ^ 50)
    println("INDICATOR SUMMARY - Year $(years[year_idx])")
    println("=" ^ 50)
    
    for name in ["relative_cover", "relative_juveniles"]
        if haskey(results, name)
            data = results[name]
            year_data = data[year_idx, :]
            println("\n$name:")
            println("  Mean: $(round(mean(year_data), digits=4))")
            println("  Min:  $(round(minimum(year_data), digits=4))")
            println("  Max:  $(round(maximum(year_data), digits=4))")
        end
    end
end


"""
    export_for_adria(output::CscapeOutput, filepath::String)

Save a `CscapeOutput` to a JLD2 file. Load downstream with:

```julia
using JLD2
output = load(filepath, "cscape_output")
adria  = to_adria_format(output)
```
"""
function export_for_adria(output::CscapeOutput, filepath::String)
    mkpath(dirname(filepath))
    jldsave(filepath; cscape_output = output)
    @info "Exported CscapeOutput to $filepath"
end
