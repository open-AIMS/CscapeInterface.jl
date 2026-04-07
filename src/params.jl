#=
Model parameters using ModelParameters.jl
For ADRIA interoperability

# ============================================================
# USAGE EXAMPLES
# ============================================================

# Load model
model = load_params(1, fpath)

# View parameters
model           # Shows table
model[:val]     # Get values as tuple
model[:bounds]  # Get bounds
collect(model)  # Get values as vector (for Optim.jl)

# ============================================================
# WITH HELPER FUNCTIONS (recommended)
# ============================================================

set_rcp!(model, 2)              # Change RCP scenario
set_year_end!(model, 2050)      # Change end year
set_year_start!(model, 2024)    # Change start year
set_plasticity!(model, 0.5)     # Change plasticity
set_dhw_enhance!(model, 3.0)    # Change DHW enhancement
set_draw!(model, 5)             # Change posterior draw
set_cyclone_rep!(model, 2)      # Change cyclone replicate
set_scenario_id!(model, 1)      # Change scenario ID

run_cscape(model)

# ============================================================
# WITHOUT HELPER FUNCTIONS
# ============================================================

# Parameter order: (scenario_id, Cyclone_rep, rcp, year_start, year_end, Plasticity, DHW_enhance, draw)
# Index:          (     1      ,      2     ,  3 ,     4     ,    5    ,     6     ,      7     ,   8 )

# Method 1: Set ALL values at once
model[:val] = (1, 1, 2, 2024, 2050, 0.5, 3.0, 0)

# Method 2: Modify single value using collect/Tuple
v = collect(model[:val])
v[3] = 2      # rcp at index 3
v[5] = 2050   # year_end at index 5
model[:val] = Tuple(v)

# Method 3: Using splatting (change rcp at index 3)
v = model[:val]
model[:val] = (v[1:2]..., 2, v[4:8]...)

run_cscape(model)

=#

using ModelParameters: Model, Param, stripparams, parent, update!

"""
    CscapeParams

C-scape simulation parameters.
Only numeric optimization parameters use Param.
String/config fields are plain values.
"""
Base.@kwdef struct CscapeParams{A,B,C,D,E,F,G,H}
    # === Numeric parameters (wrapped in Param for optimization) ===
    scenario_id::A = Param(1; val=1, bounds=(1, 1000), description="Scenario ID")
    Cyclone_rep::B = Param(1; val=1, bounds=(1, 100), description="Cyclone replicate")
    rcp::C = Param(1; val=1, bounds=(1, 4), description="RCP scenario")
    year_start::D = Param(2024; val=2024, bounds=(2000, 2100), description="Start year")
    year_end::E = Param(2050; val=2050, bounds=(2000, 2100), description="End year")
    Plasticity::F = Param(0; val=0, bounds=(0, 1), description="Plasticity")
    DHW_enhance::G = Param(0.0; val=0.0, bounds=(0.0, 10.0), description="DHW enhancement")
    draw::H = Param(0; val=0, bounds=(0, 1000), description="Posterior draw (0=mean)")
    
    # === Config fields (not optimized) ===
    rootdir_data::String = ""
    region::String = ""
    simulation_name::String = ""
    Intervention::String = "No"
    sites::Bool = true
    init_cover::String = ""
    fts::Vector{String} = String[]
    HeatTolerance::String = ""
    HeatInit::String = ""
    Heritability::String = ""
    output::Bool = true
    TradeOff::String = ""
    Spatial_file::String = ""
    Connectivity_file::String = ""
    Disturbance_file::String = ""
    Growth_Surv_file::Vector{String} = String[]
    intervened_sites::String = ""
    log_file::String = ""
    use_cached_IPM::Bool = true
    temp_growth_switch::Bool = true
end


"""
    to_dict(model) -> Dict{String, Any}

Convert CscapeParams Model to Dict for R interop.
"""
function to_dict(model::Model)
    params = parent(model)
    
    # Get draw value, convert 0 to nothing for R's NA
    draw_val = stripparams(params.draw)
    draw_out = draw_val == 0 ? nothing : draw_val
    
    Dict{String, Any}(
        "rootdir_data" => params.rootdir_data,
        "scenario_id" => stripparams(params.scenario_id),
        "region" => params.region,
        "simulation_name" => params.simulation_name,
        "draw" => draw_out,
        "Cyclone_rep" => stripparams(params.Cyclone_rep),
        "rcp" => stripparams(params.rcp),
        "year_start" => stripparams(params.year_start),
        "year_end" => stripparams(params.year_end),
        "Intervention" => params.Intervention,
        "sites" => params.sites,
        "init_cover" => params.init_cover,
        "fts" => params.fts,
        "HeatTolerance" => params.HeatTolerance,
        "HeatInit" => params.HeatInit,
        "Heritability" => params.Heritability,
        "Plasticity" => stripparams(params.Plasticity),
        "DHW_enhance" => stripparams(params.DHW_enhance),
        "output" => params.output,
        "TradeOff" => params.TradeOff,
        "Spatial_file" => params.Spatial_file,
        "Connectivity_file" => params.Connectivity_file,
        "Disturbance_file" => params.Disturbance_file,
        "Growth_Surv_file" => params.Growth_Surv_file,
        "intervened_sites" => params.intervened_sites,
        "log_file" => params.log_file,
        "use_cached_IPM" => params.use_cached_IPM,
        "temp_growth_switch" => params.temp_growth_switch
    )
end


"""
    from_dict(d::Dict) -> Model

Convert Dict to CscapeParams Model.
"""
function from_dict(d::Dict)
    # Helper to safely get values
    function get_int(k, default)
        v = get(d, k, default)
        if isnothing(v) || ismissing(v)
            return default
        elseif v isa Number && isnan(v)
            return default
        else
            return Int(v)
        end
    end
    
    function get_float(k, default)
        v = get(d, k, default)
        if isnothing(v) || ismissing(v)
            return default
        elseif v isa Number && isnan(v)
            return default
        elseif v isa String
            return default
        else
            return Float64(v)
        end
    end
    
    function get_str(k, default)
        v = get(d, k, default)
        if isnothing(v) || ismissing(v)
            return default
        elseif v isa Number && isnan(v)
            return default
        else
            return string(v)
        end
    end
    
    function get_bool(k, default)
        v = get(d, k, default)
        if isnothing(v) || ismissing(v)
            return default
        else
            return Bool(v)
        end
    end
    
    # Handle fts - can be String or Vector
    fts_raw = get(d, "fts", String[])
    fts = if fts_raw isa Vector
        String[string(x) for x in fts_raw]
    elseif isnothing(fts_raw) || ismissing(fts_raw)
        String[]
    else
        String[string(fts_raw)]
    end
    
    # Handle Growth_Surv_file
    gsf_raw = get(d, "Growth_Surv_file", String[])
    gsf = if gsf_raw isa Vector
        String[string(x) for x in gsf_raw]
    elseif isnothing(gsf_raw) || ismissing(gsf_raw)
        String[]
    else
        String[string(gsf_raw)]
    end
    
    # Handle draw - convert nothing/NA to 0
    draw_raw = get(d, "draw", nothing)
    draw_val = if isnothing(draw_raw) || ismissing(draw_raw)
        0
    elseif draw_raw isa Number && isnan(draw_raw)
        0
    else
        Int(draw_raw)
    end
    
    # Handle Plasticity - can be Int or Float
    plasticity_raw = get(d, "Plasticity", 0)
    plasticity_val = if isnothing(plasticity_raw) || ismissing(plasticity_raw)
        0
    elseif plasticity_raw isa Number && isnan(plasticity_raw)
        0
    elseif plasticity_raw isa String
        0
    else
        Int(plasticity_raw)
    end
    
    params = CscapeParams(
        # Param-wrapped numeric fields
        scenario_id = Param(get_int("scenario_id", 1); val=get_int("scenario_id", 1), bounds=(1, 1000)),
        Cyclone_rep = Param(get_int("Cyclone_rep", 1); val=get_int("Cyclone_rep", 1), bounds=(1, 100)),
        rcp = Param(get_int("rcp", 1); val=get_int("rcp", 1), bounds=(1, 4)),
        year_start = Param(get_int("year_start", 2024); val=get_int("year_start", 2024), bounds=(2000, 2100)),
        year_end = Param(get_int("year_end", 2050); val=get_int("year_end", 2050), bounds=(2000, 2100)),
        Plasticity = Param(plasticity_val; val=plasticity_val, bounds=(0, 1)),
        DHW_enhance = Param(get_float("DHW_enhance", 0.0); val=get_float("DHW_enhance", 0.0), bounds=(0.0, 10.0)),
        draw = Param(draw_val; val=draw_val, bounds=(0, 1000)),
        
        # Plain config fields
        rootdir_data = get_str("rootdir_data", ""),
        region = get_str("region", ""),
        simulation_name = get_str("simulation_name", ""),
        Intervention = get_str("Intervention", "No"),
        sites = get_bool("sites", true),
        init_cover = get_str("init_cover", ""),
        fts = fts,
        HeatTolerance = get_str("HeatTolerance", ""),
        HeatInit = get_str("HeatInit", ""),
        Heritability = get_str("Heritability", ""),
        output = get_bool("output", true),
        TradeOff = get_str("TradeOff", ""),
        Spatial_file = get_str("Spatial_file", ""),
        Connectivity_file = get_str("Connectivity_file", ""),
        Disturbance_file = get_str("Disturbance_file", ""),
        Growth_Surv_file = gsf,
        intervened_sites = get_str("intervened_sites", ""),
        log_file = get_str("log_file", ""),
        use_cached_IPM = get_bool("use_cached_IPM", true),
        temp_growth_switch = get_bool("temp_growth_switch", true)
    )
    return Model(params)
end


# ============================================================
# Helper functions to modify individual parameters
# ============================================================

"""
    set_param!(model, index, value)

Set parameter at given index. Use `model` to see parameter order.
"""
function set_param!(model::Model, index::Int, value)
    v = collect(model[:val])
    v[index] = value
    model[:val] = Tuple(v)
    return model
end

"""Set scenario_id (index 1)"""
set_scenario_id!(model::Model, val::Int) = set_param!(model, 1, val)

"""Set Cyclone_rep (index 2)"""
set_cyclone_rep!(model::Model, val::Int) = set_param!(model, 2, val)

"""Set rcp (index 3) - values 1-4"""
set_rcp!(model::Model, val::Int) = set_param!(model, 3, val)

"""Set year_start (index 4)"""
set_year_start!(model::Model, val::Int) = set_param!(model, 4, val)

"""Set year_end (index 5)"""
set_year_end!(model::Model, val::Int) = set_param!(model, 5, val)

"""Set Plasticity (index 6)"""
set_plasticity!(model::Model, val) = set_param!(model, 6, val)

"""Set DHW_enhance (index 7)"""
set_dhw_enhance!(model::Model, val::Float64) = set_param!(model, 7, val)

"""Set draw (index 8) - 0 means use mean"""
set_draw!(model::Model, val::Int) = set_param!(model, 8, val)


# =============================================================================
# HELPER FUNCTIONS FOR STRING/BOOL/VECTOR PARAMETERS
# =============================================================================
# 
# These functions modify parameters that are NOT wrapped in Param()
# Since ModelParameters structs are immutable, each helper:
#   1. Converts model to Dict
#   2. Changes the value
#   3. Converts back to new Model
#   4. Returns new Model
#
# Usage: model = set_region!(model, "Moore")  # Must reassign!
# =============================================================================

# -----------------------------------------------------------------------------
# Scenario ID Group
# -----------------------------------------------------------------------------

"""
    set_region!(model, val) -> Model

Set region name. Returns new model (must reassign).

# Example
```julia
model = set_region!(model, "Moore")
```
"""
function set_region!(model::Model, val::String)
    d = to_dict(model)
    d["region"] = val
    return from_dict(d)
end

"""
    set_simulation_name!(model, val) -> Model

Set simulation folder name. Returns new model (must reassign).
"""
function set_simulation_name!(model::Model, val::String)
    d = to_dict(model)
    d["simulation_name"] = val
    return from_dict(d)
end

# -----------------------------------------------------------------------------
# Paths Group
# -----------------------------------------------------------------------------

"""
    set_rootdir_data!(model, val) -> Model

Set data directory path. Returns new model (must reassign).
"""
function set_rootdir_data!(model::Model, val::String)
    d = to_dict(model)
    d["rootdir_data"] = val
    return from_dict(d)
end

# -----------------------------------------------------------------------------
# Data Files Group
# -----------------------------------------------------------------------------

"""
    set_spatial_file!(model, val) -> Model

Set spatial data file name. Returns new model (must reassign).

# Example
```julia
model = set_spatial_file!(model, "New_k_MooreReefCluster")
```
"""
function set_spatial_file!(model::Model, val::String)
    d = to_dict(model)
    d["Spatial_file"] = val
    return from_dict(d)
end

"""
    set_connectivity_file!(model, val) -> Model

Set connectivity file name. Returns new model (must reassign).
"""
function set_connectivity_file!(model::Model, val::String)
    d = to_dict(model)
    d["Connectivity_file"] = val
    return from_dict(d)
end

"""
    set_disturbance_file!(model, val) -> Model

Set disturbance file name. Returns new model (must reassign).
"""
function set_disturbance_file!(model::Model, val::String)
    d = to_dict(model)
    d["Disturbance_file"] = val
    return from_dict(d)
end

"""
    set_growth_surv_file!(model, val) -> Model

Set growth/survival file names. Returns new model (must reassign).

# Example
```julia
model = set_growth_surv_file!(model, ["file1.RData", "file2.RData"])
```
"""
function set_growth_surv_file!(model::Model, val::Vector{String})
    d = to_dict(model)
    d["Growth_Surv_file"] = val
    return from_dict(d)
end

# -----------------------------------------------------------------------------
# Cover/Population Group
# -----------------------------------------------------------------------------

"""
    set_init_cover!(model, val) -> Model

Set initial cover file or values. Returns new model (must reassign).
"""
function set_init_cover!(model::Model, val::String)
    d = to_dict(model)
    d["init_cover"] = val
    return from_dict(d)
end

# -----------------------------------------------------------------------------
# Coral Traits Group
# -----------------------------------------------------------------------------

"""
    set_heat_tolerance!(model, val) -> Model

Set heat tolerance groups. Returns new model (must reassign).

# Example
```julia
model = set_heat_tolerance!(model, "3groups")
```
"""
function set_heat_tolerance!(model::Model, val::String)
    d = to_dict(model)
    d["HeatTolerance"] = val
    return from_dict(d)
end

"""
    set_heat_init!(model, val) -> Model

Set initial heat tolerance distribution. Returns new model (must reassign).
"""
function set_heat_init!(model::Model, val::String)
    d = to_dict(model)
    d["HeatInit"] = val
    return from_dict(d)
end

"""
    set_fts!(model, val) -> Model

Set functional types. Returns new model (must reassign).

# Example
```julia
model = set_fts!(model, ["Acropora", "Pocillopora", "Massive"])
```
"""
function set_fts!(model::Model, val::Vector{String})
    d = to_dict(model)
    d["fts"] = val
    return from_dict(d)
end

# -----------------------------------------------------------------------------
# Model Settings Group
# -----------------------------------------------------------------------------

"""
    set_tradeoff!(model, val) -> Model

Set trade-off setting. Returns new model (must reassign).
"""
function set_tradeoff!(model::Model, val::String)
    d = to_dict(model)
    d["TradeOff"] = val
    return from_dict(d)
end

"""
    set_output!(model, val) -> Model

Set whether to save output files. Returns new model (must reassign).

# Example
```julia
model = set_output!(model, true)
```
"""
function set_output!(model::Model, val::Bool)
    d = to_dict(model)
    d["output"] = val
    return from_dict(d)
end

"""
    set_use_cached_ipm!(model, val) -> Model

Set whether to use cached IPM. Returns new model (must reassign).

# Example
```julia
model = set_use_cached_ipm!(model, false)  # Recompute IPM
```
"""
function set_use_cached_ipm!(model::Model, val::Bool)
    d = to_dict(model)
    d["use_cached_IPM"] = val
    return from_dict(d)
end

"""
    set_temp_growth_switch!(model, val) -> Model

Set temperature-dependent growth switch. Returns new model (must reassign).

# Example
```julia
model = set_temp_growth_switch!(model, true)
```
"""
function set_temp_growth_switch!(model::Model, val::Bool)
    d = to_dict(model)
    d["temp_growth_switch"] = val
    return from_dict(d)
end

# -----------------------------------------------------------------------------
# Interventions Group
# -----------------------------------------------------------------------------

"""
    set_intervention!(model, val) -> Model

Set intervention enabled ("Yes" or "No"). Returns new model (must reassign).

# Example
```julia
model = set_intervention!(model, "Yes")
```
"""
function set_intervention!(model::Model, val::String)
    d = to_dict(model)
    d["Intervention"] = val
    return from_dict(d)
end

"""
    set_intervened_sites!(model, val) -> Model

Set intervened site IDs. Returns new model (must reassign).

# Example
```julia
model = set_intervened_sites!(model, "1,2,3,4,5")
```
"""
function set_intervened_sites!(model::Model, val::String)
    d = to_dict(model)
    d["intervened_sites"] = val
    return from_dict(d)
end

# -----------------------------------------------------------------------------
# Spatial Group
# -----------------------------------------------------------------------------

"""
    set_sites!(model, val) -> Model

Set sites to use (true for all, or vector of site names). Returns new model (must reassign).

# Example
```julia
model = set_sites!(model, true)  # Use all sites
model = set_sites!(model, ["Site1", "Site2", "Site3"])  # Specific sites
```
"""
function set_sites!(model::Model, val::Union{Bool, Vector{String}})
    d = to_dict(model)
    d["sites"] = val
    return from_dict(d)
end


# =============================================================================
# CONVENIENCE FUNCTION: Modify multiple parameters at once
# =============================================================================

"""
    set_params!(model; kwargs...) -> Model

Modify multiple parameters at once. Returns new model (must reassign).

# Example
```julia
model = set_params!(model;
    region = "Moore",
    rcp = 2,
    year_end = 2050,
    plasticity = 0.5,
    use_cached_IPM = false
)
```
"""
function set_params!(model::Model; kwargs...)
    d = to_dict(model)
    
    # Map keyword arguments to dict keys
    key_map = Dict(
        :region => "region",
        :simulation_name => "simulation_name",
        :rootdir_data => "rootdir_data",
        :spatial_file => "Spatial_file",
        :connectivity_file => "Connectivity_file",
        :disturbance_file => "Disturbance_file",
        :growth_surv_file => "Growth_Surv_file",
        :init_cover => "init_cover",
        :heat_tolerance => "HeatTolerance",
        :heat_init => "HeatInit",
        :fts => "fts",
        :tradeoff => "TradeOff",
        :output => "output",
        :use_cached_ipm => "use_cached_IPM",
        :temp_growth_switch => "temp_growth_switch",
        :intervention => "Intervention",
        :intervened_sites => "intervened_sites",
        :sites => "sites",
        # Numeric params (also supported)
        :scenario_id => "scenario_id",
        :cyclone_rep => "Cyclone_rep",
        :rcp => "rcp",
        :year_start => "year_start",
        :year_end => "year_end",
        :plasticity => "Plasticity",
        :heritability => "Heritability",
        :dhw_enhance => "DHW_enhance",
        :draw => "draw"
    )
    
    for (key, val) in kwargs
        dict_key = get(key_map, key, nothing)
        if isnothing(dict_key)
            @warn "Unknown parameter: $key (skipping)"
        else
            d[dict_key] = val
        end
    end
    
    return from_dict(d)
end