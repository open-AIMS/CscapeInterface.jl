# =============================================================================
# Julia Batch Processing for C-scape
# Processes all sites × FTs in one Julia call to avoid R↔Julia overhead
# =============================================================================

module BatchProcess

using Statistics
using LinearAlgebra
import SpecialFunctions: erf

export process_year_batch, create_ipm_jl, ipm_pred_jl, batch_create_ipms_full

# -----------------------------------------------------------------------------
# Helper: Convert diameter to area (cm²)
# -----------------------------------------------------------------------------
function to_area(diameter)
    return π .* (diameter ./ 2).^2
end

# -----------------------------------------------------------------------------
# Helper: Get coral cover from population vector
# -----------------------------------------------------------------------------
function get_cover(nt_pred::AbstractVector, meshpoints::AbstractVector, polygon_area::Float64)
    # nt_pred is 103 elements: [3 early stages, 100 size classes]
    # meshpoints is 100 elements
    areas = to_area(meshpoints)
    coral_area = sum(areas .* nt_pred[4:103])
    cover_pct = (coral_area / polygon_area) * 100
    return cover_pct
end

function get_cover_summed(nt_t::AbstractMatrix, meshpoints::AbstractVector, polygon_area::Float64)
    # nt_t is 103 × n_enhancement, sum across enhancement classes
    areas = to_area(meshpoints)
    # Sum population across enhancement classes, then compute cover
    nt_sum = vec(sum(nt_t, dims=2))  # 103 elements
    coral_area = sum(areas .* nt_sum[4:103])
    cover_pct = (coral_area / polygon_area) * 100
    return cover_pct
end

# -----------------------------------------------------------------------------
# Create IPM Matrix (Julia version)
# -----------------------------------------------------------------------------
function create_ipm_jl(
    pred_growth::Vector{Float64},      # 100 predicted growth values
    pred_surv::Vector{Float64},        # 100 survival probabilities
    sigma_growth::Vector{Float64},     # 100 growth SDs
    log_area_L::Float64,               # Lower bound
    log_area_U::Float64,               # Upper bound
    log_area_mshpts::Vector{Float64},  # 100 log(area) meshpoints
    eggs::Vector{Float64},             # 100 fecundity values
    larv_prob::Float64,                # Larvae probability
    settle_prob::Float64,              # Settlement probability
    sett_to_cont::Vector{Float64},     # Settler to continuous probs
    d_i::Float64,                      # Density dependence
    temp_growth::Float64,              # Temperature growth modifier
    initial_handle::Bool               # Include early life transitions
)
    n_sizes = 100
    
    # Calculate bin edges
    log_area_edges = range(log_area_L, log_area_U, length=n_sizes + 1)
    
    # Clamp bounds
    min_log = minimum(log_area_mshpts)
    max_log = maximum(log_area_mshpts)
    
    # Pre-allocate probability arrays
    prob_array_g = zeros(n_sizes, n_sizes)
    prob_array_g_s = zeros(n_sizes, n_sizes)
    
    # Main loop over size classes
    for n in 1:n_sizes
        mu_size_t1 = pred_growth[n]
        mu_surv = pred_surv[n]
        sigma = sigma_growth[n]
        
        # Clamp to valid range
        mu_size_t1 = clamp(mu_size_t1, min_log, max_log)
        
        # CDF method for growth probabilities (vectorized)
        for m in 1:n_sizes
            cdf_upper = 0.5 * (1 + erf((log_area_edges[m+1] - mu_size_t1) / (sigma * sqrt(2))))
            cdf_lower = 0.5 * (1 + erf((log_area_edges[m] - mu_size_t1) / (sigma * sqrt(2))))
            prob_array_g[m, n] = cdf_upper - cdf_lower
        end
        
        # Normalize column
        col_sum = sum(prob_array_g[:, n])
        if col_sum > 0
            prob_array_g[:, n] ./= col_sum
        end
        
        # Density dependence on lower triangle
        growth_sum = sum(prob_array_g[n:n_sizes, n])
        
        # Reduce below-diagonal by d_i * temp_growth
        for m in (n+1):n_sizes
            prob_array_g[m, n] *= d_i * temp_growth
        end
        
        # Add difference to diagonal
        growth_sum_post = sum(prob_array_g[n:n_sizes, n])
        prob_array_g[n, n] += (growth_sum - growth_sum_post)
        
        # Apply survival (cap at 0.975)
        surv_capped = min(mu_surv, 0.975)
        prob_array_g_s[:, n] = prob_array_g[:, n] .* surv_capped
    end
    
    # Build full IPM matrix (103 × 103)
    Pmat = zeros(n_sizes + 3, n_sizes + 3)
    
    # Fill top rows (fecundity)
    Pmat[1, 4:103] = eggs
    Pmat[2, 4:103] = eggs .* larv_prob
    if initial_handle
        Pmat[3, 4:103] = eggs .* larv_prob .* settle_prob
    end
    
    # Fill side columns (early life transitions)
    Pmat[2, 1] = larv_prob
    if initial_handle
        Pmat[3, 2] = settle_prob
    end
    
    # Settler to continuous transitions
    n_sett = length(sett_to_cont)
    Pmat[4:(3 + n_sett), 3] = sett_to_cont
    
    # Fill main growth/survival matrix
    Pmat[4:103, 4:103] = prob_array_g_s
    
    return Pmat
end

# -----------------------------------------------------------------------------
# DHW Mortality (Julia version)
# -----------------------------------------------------------------------------
function create_dhw_mortality(
    dhw::Float64,
    dhw_adjustment::Vector{Float64},  # Per tolerance class
    depth::Float64,
    fogging::Float64,
    bleaching_suscept::Float64
)
    n_classes = length(dhw_adjustment)
    p_mort = zeros(n_classes)
    
    for i in 1:n_classes
        if dhw <= 3
            p_mort[i] = 0.0
        else
            depth_coef = (0.420 + 0.272 * abs(depth))^(-1)
            m_init = min(
                bleaching_suscept * depth_coef * fogging *
                (exp(0.168 + 0.347 * (dhw - dhw_adjustment[i])) - 1),
                100.0
            ) / 100.0
            p_mort[i] = 1 - (1 - m_init)^6
        end
    end
    
    # Clamp to [0, 1]
    p_mort = clamp.(p_mort, 0.0, 1.0)
    
    return p_mort
end

# -----------------------------------------------------------------------------
# Cyclone Mortality (Julia version)
# Using pre-extracted coefficients instead of R model objects
# -----------------------------------------------------------------------------
function cyclone_mortality_jl(
    cyclone_cat::Int,
    depth::Float64,
    is_sensitive::Bool,
    # Model coefficients: [intercept, slope] for each model
    branching_shallow_coef::Vector{Float64},  # [intercept, slope]
    branching_deep_coef::Vector{Float64},
    massive_coef::Vector{Float64}
)
    if cyclone_cat <= 0
        return 0.0
    end
    
    # Convert category to windspeed
    windspeed = if cyclone_cat == 1
        24.5
    elseif cyclone_cat == 2
        32.5
    elseif cyclone_cat == 3
        44.2
    elseif cyclone_cat == 4
        55.3
    else  # 5+
        65.0
    end
    
    # Select model based on sensitivity and depth
    if is_sensitive
        if abs(depth) < 20
            # Branching shallow model (logistic)
            logit = branching_shallow_coef[1] + branching_shallow_coef[2] * windspeed
            p_mort = 1 / (1 + exp(-logit))
        else
            # Branching deep model (logistic)
            logit = branching_deep_coef[1] + branching_deep_coef[2] * windspeed
            p_mort = 1 / (1 + exp(-logit))
        end
    else
        # Massive model (linear)
        p_mort = massive_coef[1] + massive_coef[2] * windspeed
    end
    
    return max(p_mort, 0.0)
end

# -----------------------------------------------------------------------------
# COTS Mortality (Julia version)
# -----------------------------------------------------------------------------
function cots_mortality_jl(
    nt_pred::Matrix{Float64},         # 103 × n_enhancement
    ft::Int,
    meshpoints::Vector{Float64},
    cots_mort_yearly::Vector{Float64},  # 4 values
    odds::Vector{Float64},            # Per FT
    cots_ti::Vector{Float64},         # 4 COTS size classes
    cover_ft::Vector{Float64}         # Cover per FT
)
    n_enhance = size(nt_pred, 2)
    
    # Calculate total coral area for this FT
    areas = to_area(meshpoints)
    total_coral_sizeclass = areas .* nt_pred[4:103, :]  # 100 × n_enhance
    total_coral = sum(total_coral_sizeclass)
    
    if total_coral <= 0
        return nt_pred[4:103, :]  # Return unchanged
    end
    
    # Calculate COTS consumption
    all_cots_consumption = sum(cots_ti .* cots_mort_yearly)
    
    # Susceptibility ratio
    sr_even = odds ./ sum(odds)
    pr = sr_even ./ (1 / length(odds))
    
    # Odds given cover
    odds_given_cover = cover_ft .* pr
    
    if sum(odds_given_cover) > 0
        sr = odds_given_cover[ft] / sum(odds_given_cover)
    else
        sr = 0.0
    end
    
    species_consumption = sr * all_cots_consumption
    
    if isnan(species_consumption)
        species_consumption = 0.0
    end
    
    # Proportion removed
    p_removed = species_consumption / total_coral
    p_removed = min(p_removed, 0.95)  # Cap at 95%
    
    # Apply mortality
    cm2_pred = (1 - p_removed) .* total_coral_sizeclass
    nt_new = cm2_pred ./ areas
    
    return nt_new
end

# -----------------------------------------------------------------------------
# Apply cover-to-number differential (for DHW mortality)
# -----------------------------------------------------------------------------
function cov2num_diff(
    nt_pred::Matrix{Float64},     # 103 × n_enhancement
    meshpoints::Vector{Float64},
    polygon_area::Float64,
    p_mort::Vector{Float64}       # Per enhancement class
)
    n_enhance = length(p_mort)
    areas = to_area(meshpoints)
    
    result = similar(nt_pred[4:103, :])
    
    for e in 1:n_enhance
        # Apply mortality: reduce by (1 - p_mort)
        result[:, e] = nt_pred[4:103, e] .* (1 - p_mort[e])
    end
    
    return result
end

# -----------------------------------------------------------------------------
# Apply cover-to-number (for cyclone mortality)
# -----------------------------------------------------------------------------
function cov2num(
    nt_pred::Matrix{Float64},
    meshpoints::Vector{Float64},
    polygon_area::Float64,
    p_mort::Float64
)
    # Apply uniform mortality across all enhancement classes
    result = nt_pred[4:103, :] .* (1 - p_mort)
    return result
end

# -----------------------------------------------------------------------------
# Single site × FT processing (equivalent to taxa_loop)
# -----------------------------------------------------------------------------
function process_site_ft(
    # Population state
    nt::Matrix{Float64},              # 105 × n_enhancement (current population)
    
    # IPM parameters
    pred_growth::Vector{Float64},
    pred_surv::Vector{Float64},
    sigma_growth::Vector{Float64},
    log_area_L::Float64,
    log_area_U::Float64,
    log_area_mshpts::Vector{Float64},
    eggs::Vector{Float64},
    larv_prob::Float64,
    settle_prob_ipm::Float64,
    sett_to_cont::Vector{Float64},
    meshpoints::Vector{Float64},
    
    # Density and growth
    d_i::Float64,
    temp_growth::Float64,
    
    # Disturbances
    dhw_ti::Float64,
    cyclone_ti::Int,
    cots_ti::Vector{Float64},
    
    # Site properties
    depth_ti::Float64,
    kappa_i::Float64,
    area_i::Float64,
    polygon_area::Float64,
    
    # Larvae/settlers
    larvae_received::Vector{Float64},   # Per enhancement class
    external_eggs::Vector{Float64},     # Per enhancement class
    settle_prob::Float64,
    ext_larv_surv_sett::Float64,
    rubble::Float64,
    max_juvis::Float64,
    
    # Mortality parameters
    dhw_enhance::Vector{Float64},       # DHW adjustment per class
    fogging::Float64,
    bleaching_suscept::Float64,
    is_cyclone_sensitive::Bool,
    cyclone_coefs::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}},
    cots_mort_yearly::Vector{Float64},
    odds::Vector{Float64},
    cover_ft_all::Vector{Float64},      # Cover for all FTs
    ft::Int,
    
    # Cover calculation
    cover_tot::Float64
)
    n_enhance = size(nt, 2)
    
    # Create IPM matrix
    ipm_mat = create_ipm_jl(
        pred_growth, pred_surv, sigma_growth,
        log_area_L, log_area_U, log_area_mshpts,
        eggs, larv_prob, settle_prob_ipm, sett_to_cont,
        d_i, temp_growth, false  # InitialHandle = false
    )
    
    # Extract nt for continuous stages (columns 3:105 → indices 1:103 in 0-indexed)
    # In R: nt[,3:105], in Julia nt[:, 3:105] but nt is already 105 columns
    # Actually nt comes as 105 × n_enhance, we need columns 3:105
    nt_pred = copy(nt[3:105, :])  # 103 × n_enhance
    
    # Zero early stages (they were handled in previous timestep)
    nt_pred[1:3, :] .= 0.0  # eggs, larvae, settlers
    
    # IPM matrix without diagonal
    ipm_pred_cont = ipm_mat - Diagonal(diag(ipm_mat))
    
    # --- SETTLERS ---
    # Calculate raw settlers
    settlers_i = (larvae_received .* settle_prob .+ external_eggs .* ext_larv_surv_sett) .* (1 - rubble / 100)
    settlers_i_tot = sum(settlers_i)
    
    # Juvenile cap
    cont_prob = sum(ipm_pred_cont[4:103, 3])
    juveniles_resulting = settlers_i_tot * cont_prob
    
    # Calculate available space
    available_space = area_i * ((kappa_i - cover_tot) / 100)
    max_juvis_cc = max_juvis * available_space
    
    # Existing juveniles (simplified)
    juvi_threshold = 5.0  # cm
    juvi_indices = findall(meshpoints .< juvi_threshold)
    
    # Pre-multiply to get future juveniles
    pop_vector = vec(sum(nt_pred, dims=2))  # Sum across enhancement
    added_temp = ipm_pred_cont * pop_vector
    nt_temp = pop_vector .+ added_temp
    
    juvi_count = sum(nt_temp[juvi_indices .+ 3])  # +3 for offset
    if (kappa_i - cover_tot) > 0
        juvi_count = juvi_count / ((kappa_i - cover_tot) / 100)
    end
    
    new_juvis_allowed = max(max_juvis_cc - juvi_count, 0.0)
    
    # Apply juvenile cap
    if juveniles_resulting > new_juvis_allowed && settlers_i_tot > 0
        settlers_allowed = new_juvis_allowed / cont_prob
        d_sett = settlers_allowed / settlers_i_tot
        settlers_i = settlers_i .* d_sett
    end
    
    settlers_i = round.(settlers_i)
    
    # Add settlers to population
    nt_pred[3, :] .+= settlers_i
    
    # --- CYCLONE MORTALITY ---
    if cyclone_ti > 0
        p_mort_cyclone = cyclone_mortality_jl(
            cyclone_ti, depth_ti, is_cyclone_sensitive,
            cyclone_coefs[1], cyclone_coefs[2], cyclone_coefs[3]
        )
        nt_pred[4:103, :] = cov2num(nt_pred, meshpoints, polygon_area, p_mort_cyclone)
    end
    
    # --- DHW MORTALITY ---
    p_mort_dhw = create_dhw_mortality(dhw_ti, dhw_enhance, depth_ti, fogging, bleaching_suscept)
    
    # Apply only to non-juveniles (size >= 4.37 cm)
    juvenile_threshold = 4.37
    juvi_max_idx = findlast(meshpoints .< juvenile_threshold)
    if juvi_max_idx === nothing
        juvi_max_idx = 0
    end
    
    # Apply DHW mortality to adult sizes only
    dhw_result = cov2num_diff(nt_pred, meshpoints, polygon_area, p_mort_dhw)
    nt_pred[(juvi_max_idx + 4):103, :] = dhw_result[(juvi_max_idx + 1):100, :]
    
    # --- COTS MORTALITY ---
    nt_pred[4:103, :] = cots_mortality_jl(
        nt_pred, ft, meshpoints, cots_mort_yearly, odds, cots_ti, cover_ft_all
    )
    
    # --- GROWTH AND SURVIVAL ---
    added = ipm_pred_cont * nt_pred
    
    # Round final population
    nt1_pred = round.(nt_pred .+ added)
    
    # Store settlers in output
    nt1_pred[3, :] = settlers_i
    
    return nt1_pred
end

# -----------------------------------------------------------------------------
# Main batch processor: process entire year
# Called once from R with all data
# -----------------------------------------------------------------------------
function process_year_batch(
    # Dimensions
    n_sites::Int,
    n_fts::Int,
    n_enhance::Int,
    
    # Population state: 4D array [site, ft, enhance, 105]
    nt_all::Array{Float64, 4},
    
    # Per-site data
    dhw::Vector{Float64},             # [n_sites]
    cyclone::Vector{Int},             # [n_sites]
    cots::Matrix{Float64},            # [n_sites, 4]
    depth::Vector{Float64},           # [n_sites]
    kappa::Vector{Float64},           # [n_sites]
    area::Vector{Float64},            # [n_sites] in m²
    polygon_area::Vector{Float64},    # [n_sites] in cm²
    rubble::Vector{Float64},          # [n_sites]
    fogging::Vector{Float64},         # [n_sites]
    temp_growth::Matrix{Float64},     # [n_sites, n_fts]
    
    # Pre-computed larvae dispersal: [n_sites, n_enhance, n_fts]
    larvae_received::Array{Float64, 3},
    external_eggs::Array{Float64, 3},
    
    # Per-FT IPM parameters
    pred_growth_all::Matrix{Float64},     # [100, n_fts] - will be site-specific
    pred_surv_all::Matrix{Float64},       # [100, n_fts]
    sigma_growth_all::Matrix{Float64},    # [100, n_fts]
    log_area_L::Vector{Float64},          # [n_fts]
    log_area_U::Vector{Float64},          # [n_fts]
    log_area_mshpts::Matrix{Float64},     # [100, n_fts]
    eggs_all::Matrix{Float64},            # [100, n_fts]
    larv_prob::Vector{Float64},           # [n_fts]
    settle_prob_ipm::Vector{Float64},     # [n_fts]
    sett_to_cont::Matrix{Float64},        # [max_len, n_fts]
    meshpoints::Matrix{Float64},          # [100, n_fts]
    
    # Mortality parameters
    dhw_enhance::Vector{Float64},         # Enhancement classes
    bleaching_suscept::Vector{Float64},   # [n_fts]
    is_sensitive::Vector{Bool},           # [n_fts]
    cyclone_shallow_coef::Vector{Float64},
    cyclone_deep_coef::Vector{Float64},
    cyclone_massive_coef::Vector{Float64},
    cots_mort_yearly::Vector{Float64},    # [4]
    odds::Vector{Float64},                # [n_fts]
    
    # Settlement parameters
    settle_prob::Vector{Float64},         # [n_fts]
    ext_larv_surv_sett::Vector{Float64},  # [n_fts]
    max_juvis::Vector{Float64}            # [n_fts]
)
    # Output array: [site, ft, enhance, 107]
    # [cover, ext_eggs, larvae_received, 103 pop values, lambda]
    values_out = zeros(n_sites, n_fts, n_enhance, 107)
    
    cyclone_coefs = (cyclone_shallow_coef, cyclone_deep_coef, cyclone_massive_coef)
    
    # Process each site
    Threads.@threads for site in 1:n_sites
        # Calculate total cover across all FTs for this site
        cover_ft = zeros(n_fts)
        for ft in 1:n_fts
            nt_site_ft = nt_all[site, ft, :, 3:105]  # [n_enhance, 103]
            cover_ft[ft] = get_cover_summed(
                permutedims(nt_site_ft, (2, 1)),  # 103 × n_enhance
                meshpoints[:, ft],
                polygon_area[site]
            )
        end
        cover_tot = sum(cover_ft)
        
        # Calculate density dependence d_i (same for all FTs at site)
        d_i = if kappa[site] <= 0
            0.0
        elseif kappa[site] > cover_tot
            (kappa[site] - cover_tot) / (kappa[site] - (kappa[site] * 0.1 / 90))
        else
            0.0
        end
        d_i = max(d_i, 0.0)
        
        # Process each functional type
        for ft in 1:n_fts
            # Extract population for this site/ft
            nt = permutedims(nt_all[site, ft, :, :], (2, 1))  # 105 × n_enhance
            
            # Get site-specific IPM predictions
            # NOTE: In full version, these would be looked up by WaveExposure/DepthCat
            pred_growth = pred_growth_all[:, ft]
            pred_surv = pred_surv_all[:, ft]
            sigma_growth = sigma_growth_all[:, ft]
            
            # Process site × FT
            nt1_pred = process_site_ft(
                nt,
                pred_growth, pred_surv, sigma_growth,
                log_area_L[ft], log_area_U[ft], log_area_mshpts[:, ft],
                eggs_all[:, ft], larv_prob[ft], settle_prob_ipm[ft],
                sett_to_cont[:, ft], meshpoints[:, ft],
                d_i, temp_growth[site, ft],
                dhw[site], cyclone[site], cots[site, :],
                depth[site], kappa[site], area[site], polygon_area[site],
                larvae_received[site, :, ft], external_eggs[site, :, ft],
                settle_prob[ft], ext_larv_surv_sett[ft],
                rubble[site], max_juvis[ft],
                dhw_enhance, fogging[site], bleaching_suscept[ft],
                is_sensitive[ft], cyclone_coefs,
                cots_mort_yearly, odds, cover_ft, ft,
                cover_tot
            )
            
            # Calculate cover from result
            for e in 1:n_enhance
                cover = get_cover(nt1_pred[:, e], meshpoints[:, ft], polygon_area[site])
                
                # Store results
                values_out[site, ft, e, 1] = cover
                values_out[site, ft, e, 2] = external_eggs[site, e, ft]
                values_out[site, ft, e, 3] = larvae_received[site, e, ft]
                values_out[site, ft, e, 4:106] = nt1_pred[:, e]
                values_out[site, ft, e, 107] = NaN  # lambda not computed
            end
        end
    end
    
    return values_out
end

# -----------------------------------------------------------------------------
# Batch create all IPMs for a year
# Called from R via JuliaCall to pre-compute all IPMs
# -----------------------------------------------------------------------------
function batch_create_ipms_full(
    n_sites::Int,
    n_fts::Int,
    pred_growth_all::Array{Float64, 3},   # [100, n_sites, n_fts]
    pred_surv_all::Array{Float64, 3},     # [100, n_sites, n_fts]
    sigma_growth_all::Matrix{Float64},    # [100, n_fts]
    log_area_L::Vector{Float64},          # [n_fts]
    log_area_U::Vector{Float64},          # [n_fts]
    log_area_mshpts::Matrix{Float64},     # [100, n_fts]
    eggs_all::Matrix{Float64},            # [100, n_fts]
    larv_prob::Vector{Float64},           # [n_fts]
    settle_prob_ipm::Vector{Float64},     # [n_fts]
    sett_to_cont_all::Matrix{Float64},    # [max_len, n_fts]
    all_d_i::Vector{Float64},             # [n_sites]
    all_temp_growth::Matrix{Float64}      # [n_sites, n_fts]
)
    results = Array{Matrix{Float64}}(undef, n_sites, n_fts)
    
    for site in 1:n_sites
        d_i = all_d_i[site]
        for ft in 1:n_fts
            temp_growth = all_temp_growth[site, ft]
            results[site, ft] = create_ipm_jl(
                pred_growth_all[:, site, ft],
                pred_surv_all[:, site, ft],
                sigma_growth_all[:, ft],
                log_area_L[ft],
                log_area_U[ft],
                log_area_mshpts[:, ft],
                eggs_all[:, ft],
                larv_prob[ft],
                settle_prob_ipm[ft],
                sett_to_cont_all[:, ft],
                d_i,
                temp_growth,
                false
            )
        end
    end
    return results
end

end  # module BatchProcess