#=
Task 6: Connect with ADRIAIndicators
=#

using ADRIAIndicators
using JLD2

"""
    CscapeIndicators

Container for all reef indicators calculated from a CscapeOutput.
"""
struct CscapeIndicators
    # Cover metrics
    relative_cover::Array{Float64,2}               # [timesteps × locations]
    relative_loc_taxa_cover::Array{Float64,3}      # [timesteps × groups × locations]
    relative_taxa_cover::Array{Float64,2}          # [timesteps × groups]
    ltmp_cover::Array{Float64,2}                   # [timesteps × locations]
    # Metrics
    relative_shelter_volume::Array{Float64,4}      # [timesteps × groups × sizes × locations]
    coral_diversity::Array{Float64,2}              # [timesteps × locations]
    coral_evenness::Array{Float64,2}               # [timesteps × locations]
    relative_juveniles::Array{Float64,2}           # [timesteps × locations]
    juvenile_indicator::Array{Float64,2}           # [timesteps × locations]
    relative_loc_taxa_juveniles::Array{Float64,3}  # [timesteps × groups × locations]
    relative_taxa_juveniles::Array{Float64,2}      # [timesteps × groups]
    # Reef indices
    reef_biodiversity_condition_index::Array{Float64,2}  # [timesteps × locations]
    reef_condition_index::Array{Float64,2}               # [timesteps × locations]
    reef_fish_index::Array{Float64,2}                    # [timesteps × locations]
    reef_tourism_index::Array{Float64,2}                 # [timesteps × locations]
    # Metadata
    years::Vector{Int}
    site_ids::Vector{String}
    fts::Vector{String}
end

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
        habitable_area = Float64.(output.kappa .* output.area[:, intervention_idx]),
        reef_area = Float64.(output.area[:, intervention_idx]),
        area = output.area[:, intervention_idx],
        meshpoints = output.meshpoints
    )
end


"""
    calculate_indicators(output::CscapeOutput) -> CscapeIndicators

Calculate reef indicators using ADRIAIndicators.
"""
function calculate_indicators(output::CscapeOutput)
    adria = to_adria_format(output)
    cover = adria.cover
    habitable_area = adria.habitable_area
    reef_area = adria.reef_area

    colony_mean_diam_cm = output.meshpoints
    n_sizes = size(cover, 3)
    pap_2d = planar_area_params()  # [groups, 2]
    pap = repeat(reshape(pap_2d, size(pap_2d, 1), 1, size(pap_2d, 2)), 1, n_sizes, 1)  # [groups, sizes, 2]
    reference = (mean(colony_mean_diam_cm), -8.79, 3.14)  # Tuple{T, T, T}: (mean diameter, intercept, coefficient)
    is_juvenile = output.is_juvenile

    # Cover metrics
    rel_cover          = ADRIAIndicators.relative_cover(cover)
    rel_loc_taxa_cover = ADRIAIndicators.relative_loc_taxa_cover(cover)
    rel_taxa_cover     = ADRIAIndicators.relative_taxa_cover(cover, habitable_area)
    ltmp_cov           = ADRIAIndicators.ltmp_cover(cover, habitable_area, reef_area)

    # Metrics
    rsv                  = ADRIAIndicators.relative_shelter_volume(cover, colony_mean_diam_cm, pap, habitable_area, reference)
    rsv_2d               = dropdims(sum(rsv, dims=(2, 3)), dims=(2, 3))  # [timesteps × locations] for reef indices
    coral_div            = ADRIAIndicators.coral_diversity(rel_loc_taxa_cover)
    coral_even           = ADRIAIndicators.coral_evenness(rel_loc_taxa_cover)
    rel_juv              = ADRIAIndicators.relative_juveniles(cover, is_juvenile)
    juv_ind              = ADRIAIndicators.juvenile_indicator(cover, is_juvenile, habitable_area, colony_mean_diam_cm, 15.0)
    rel_loc_taxa_juv     = ADRIAIndicators.relative_loc_taxa_juveniles(cover, is_juvenile)
    rel_taxa_juv         = ADRIAIndicators.relative_taxa_juveniles(cover, is_juvenile, habitable_area)

    # Reef indices
    rbci  = ADRIAIndicators.reef_biodiversity_condition_index(rel_cover, coral_div, rsv_2d)
    rci   = ADRIAIndicators.reef_condition_index(ltmp_cov, rsv_2d, juv_ind)
    rfi   = ADRIAIndicators.reef_fish_index(rel_cover)
    rti   = ADRIAIndicators.reef_tourism_index_no_rubble(ltmp_cov, coral_even, rsv_2d, rel_juv)

    @info "Calculated indicators" n_years=length(output.years) n_sites=length(output.site_ids)

    return CscapeIndicators(
        rel_cover, rel_loc_taxa_cover, rel_taxa_cover, ltmp_cov,
        rsv, coral_div, coral_even, rel_juv, juv_ind, rel_loc_taxa_juv, rel_taxa_juv,
        rbci, rci, rfi, rti,
        output.years, output.site_ids, output.fts
    )
end


"""
    indicator_summary(ind::CscapeIndicators; year=nothing)

Print summary of indicators.
"""
function indicator_summary(ind::CscapeIndicators; year::Union{Int,Nothing} = nothing)
    year_idx = isnothing(year) ? length(ind.years) : findfirst(==(year), ind.years)

    println("=" ^ 50)
    println("INDICATOR SUMMARY - Year $(ind.years[year_idx])")
    println("=" ^ 50)

    for (name, data) in [("relative_cover", ind.relative_cover), ("relative_juveniles", ind.relative_juveniles)]
        year_data = data[year_idx, :]
        println("\n$name:")
        println("  Mean: $(round(mean(year_data), digits=4))")
        println("  Min:  $(round(minimum(year_data), digits=4))")
        println("  Max:  $(round(maximum(year_data), digits=4))")
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



"""
    planar_area_params()

Colony planar area parameters (see Fig 2B in [1])
First column is `b`, second column is `a`
log(S) = b + a * log(x)

# References
1. Aston Eoghan A., Duce Stephanie, Hoey Andrew S., Ferrari Renata (2022).
    A Protocol for Extracting Structural Metrics From 3D Reconstructions of Corals.
    Frontiers in Marine Science, 9.
    https://doi.org/10.3389/fmars.2022.854395
"""
function planar_area_params()
    return Array{Float64,2}([
        -8.95 2.80      # Tabular Acropora
        -9.13 2.94      # Corymbose Acropora
        -8.90 2.94      # Corymbose non-Acropora (using branching pocillopora values from fig2B)
        -8.87 2.30      # Small massives
        -8.87 2.30      # Large massives
        -8.90 2.94      # Corymbose non-Acropora Brooders (using branching pocillopora values from fig2B)
    ])
end