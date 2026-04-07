#=
Task 6: Connect with ADRIAIndicators
=#

using ADRIAIndicators

"""
    to_adria_format(output::CscapeOutput; intervention_idx=1) -> NamedTuple

Convert CscapeOutput to ADRIAIndicators format.

# Returns
NamedTuple with:
- `cover`: Array [timesteps, groups, sizes, locations]
- `habitable_area`: Vector [locations]
- `reef_area`: Vector [locations]
"""
function to_adria_format(output::CscapeOutput; intervention_idx::Int = 1)
    
    n_years = length(output.years)
    n_sites = length(output.site_ids)
    n_fts = length(output.fts)
    n_sizes = 103  # Size classes 4:106
    
    # Extract: [years, sites, int, ft, enh, metrics] -> [timesteps, groups, sizes, locations]
    cover = zeros(Float64, n_years, n_fts, n_sizes, n_sites)
    
    for t in 1:n_years, g in 1:n_fts, s in 1:n_sizes, l in 1:n_sites
        cover[t, g, s, l] = sum(output.out_array[t, l, intervention_idx, g, :, 3+s])
    end
    
    # Normalize by kappa
    for l in 1:n_sites
        if output.kappa[l] > 0
            cover[:, :, :, l] ./= output.kappa[l]
        end
    end
    
    clamp!(cover, 0.0, 1.0)
    
    return (
        cover = cover,
        habitable_area = Float64.(output.kappa),
        reef_area = Float64.(output.kappa)
    )
end


"""
    calculate_indicators(output::CscapeOutput; juvenile_threshold=5) -> Dict

Calculate reef indicators using ADRIAIndicators.

# Returns
Dict with indicator arrays.
"""
function calculate_indicators(output::CscapeOutput; juvenile_threshold::Int = 5)
    
    adria = to_adria_format(output)
    cover = adria.cover
    habitable_area = adria.habitable_area
    reef_area = adria.reef_area
    
    n_sizes = size(cover, 3)
    is_juvenile = falses(n_sizes)
    is_juvenile[1:juvenile_threshold] .= true
    
    results = Dict{String, Any}()
    
    # Use ADRIAIndicators functions
    results["relative_cover"] = ADRIAIndicators.relative_cover(cover)
    results["relative_loc_taxa_cover"] = ADRIAIndicators.relative_loc_taxa_cover(cover)
    results["relative_juveniles"] = ADRIAIndicators.relative_juveniles(cover, is_juvenile)
    results["relative_taxa_cover"] = ADRIAIndicators.relative_taxa_cover(cover, habitable_area)
    
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

Export output to file for external use.
"""
function export_for_adria(output::CscapeOutput, filepath::String)
    
    adria = to_adria_format(output)
    
    data = Dict(
        "cover" => adria.cover,
        "habitable_area" => adria.habitable_area,
        "reef_area" => adria.reef_area,
        "years" => output.years,
        "site_names" => output.site_ids,
        "functional_types" => output.fts,
        "n_timesteps" => size(adria.cover, 1),
        "n_groups" => size(adria.cover, 2),
        "n_sizes" => size(adria.cover, 3),
        "n_locations" => size(adria.cover, 4)
    )
    
    @rput data filepath
    R"""
    dir.create(dirname($filepath), showWarnings = FALSE, recursive = TRUE)
    saveRDS(data, $filepath)
    """
    
    @info "Exported to $filepath"
end
