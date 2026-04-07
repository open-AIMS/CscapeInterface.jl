# CscapeJulia

Julia wrapper for C-scape with ModelParameters.jl for ADRIA interoperability.

## Setup

```julia
using Pkg
Pkg.activate("C:/path/C_scape/julia")
Pkg.add("ModelParameters")
Pkg.develop(path="C:/path/ADRIAIndicators.jl")
Pkg.instantiate()
```

## Quick Start

```julia
using Pkg
Pkg.activate("C:/path/C_scape/julia")
cd("C:/path/C_scape/julia")
include("src/CscapeJulia.jl")
using .CscapeJulia

setup_r_environment(n_cores=4) #specify num of cores
setup_r_environment(enable_parallel=false) #set parallel to false
setup_r_environment()  #dafult



fpath = "C:/path/TestingInterventions"
using BenchmarkTools
@elapsed run_cscape(1, fpath)

# Run ALL scenarios in ScenarioID.xlsx
run_cscape(fpath)

# Run specific scenarios
run_cscape([1, 2, 3, 4, 5, 56], fpath)
```

---

## Running Multiple Scenarios

```julia
# Run ALL scenarios from ScenarioID.xlsx
run_cscape(fpath)

# Run specific scenario IDs
run_cscape([1, 2, 3, 4, 5, 56], fpath)

# Customise settings and apply to multiple scenarios
input_data = load_input_data(1, fpath)
input_data["scenario_id"] = [1, 2, 3, 4, 5, 56]   # vector of IDs
input_data["year_end"] = 2050
input_data["Plasticity"] = 0.6
run_cscape(input_data)
# Each scenario loads its own data from ScenarioID.xlsx,
# then year_end and Plasticity overrides are applied to all

# Skip indicators for speed during batch runs
failed = run_cscape([1, 2, 3, 4, 5], fpath; calc_indicators=false)
# Returns list of any scenario IDs that failed
```

### run_cscape Summary

| Call | What it does |
|------|-------------|
| `run_cscape(fpath)` | All scenarios in ScenarioID.xlsx |
| `run_cscape(1, fpath)` | Single scenario |
| `run_cscape([1,2,3], fpath)` | Multiple scenarios by ID |
| `run_cscape(input_data)` | Single or multiple (if `scenario_id` is a vector), with custom overrides |

---

## Using ModelParameters

### Load and View Model

```julia
model = CscapeJulia.load_params(1, fpath)

# View parameters table (shows only Param-wrapped numeric fields)
model

# Get all numeric values as tuple
model[:val]

# Get all bounds
model[:bounds]

# Get as vector (for Optim.jl)
collect(model)

# View ALL fields including strings/bools
parent(model)

# Convert to Dict to see everything
to_dict(model)
```

---

## Helper Functions

### Important: Reassignment Rules

| Helper Type | Example | Reassign needed? |
|-------------|---------|------------------|
| **Numeric** | `set_rcp!(model, 2)` | ❌ No |
| **String/Bool/Vector** | `model = set_region!(model, "Moore")` | ✅ Yes |
| **set_params!** | `model = set_params!(model; ...)` | ✅ Yes |

---

### Numeric Helpers (No reassignment needed)

```julia
model = CscapeJulia.load_params(1, fpath)

# These modify in place - no reassignment needed
set_rcp!(model, 2)              # Change RCP scenario
set_year_end!(model, 2050)      # Change end year
set_year_start!(model, 2024)    # Change start year
set_plasticity!(model, 0.5)     # Change plasticity
set_dhw_enhance!(model, 3.0)    # Change DHW enhancement
set_draw!(model, 5)             # Change posterior draw
set_cyclone_rep!(model, 2)      # Change cyclone replicate
set_scenario_id!(model, 1)      # Change scenario ID
set_heritability!(model, 0.3)   # Change heritability

CscapeJulia.run_cscape(model)
```

| Function | Parameter | Bounds |
|----------|-----------|--------|
| `set_scenario_id!(model, val)` | scenario_id | (1, 1000) |
| `set_cyclone_rep!(model, val)` | Cyclone_rep | (1, 100) |
| `set_rcp!(model, val)` | rcp | (1, 4) |
| `set_year_start!(model, val)` | year_start | (2000, 2100) |
| `set_year_end!(model, val)` | year_end | (2000, 2100) |
| `set_plasticity!(model, val)` | Plasticity | (0, 1) |
| `set_dhw_enhance!(model, val)` | DHW_enhance | (0.0, 10.0) |
| `set_draw!(model, val)` | draw | (0, 1000) |
| `set_heritability!(model, val)` | Heritability | (0, 1) |

---

### String/Bool/Vector Helpers (Reassignment required)

```julia
model = CscapeJulia.load_params(1, fpath)

# These return a NEW model - MUST reassign!
model = set_region!(model, "Moore")
model = set_simulation_name!(model, "my_simulation")
model = set_rootdir_data!(model, "/path/to/data")

# Data files
model = set_spatial_file!(model, "New_k_MooreReefCluster")
model = set_connectivity_file!(model, "connectivity_matrix")
model = set_disturbance_file!(model, "disturbance_data")
model = set_growth_surv_file!(model, ["file1.RData", "file2.RData"])

# Coral traits
model = set_heat_tolerance!(model, "3groups")
model = set_heat_init!(model, "equal")
model = set_fts!(model, ["Acropora", "Pocillopora", "Massive"])
model = set_init_cover!(model, "default_cover")

# Model settings
model = set_tradeoff!(model, "growth")
model = set_output!(model, true)
model = set_use_cached_ipm!(model, false)
model = set_temp_growth_switch!(model, true)

# Interventions
model = set_intervention!(model, "Yes")
model = set_intervened_sites!(model, "1,2,3,4,5")

# Spatial
model = set_sites!(model, true)  # Use all sites
model = set_sites!(model, ["Site1", "Site2"])  # Specific sites

CscapeJulia.run_cscape(model)
```

| Function | Parameter | Type |
|----------|-----------|------|
| `set_region!(model, val)` | region | String |
| `set_simulation_name!(model, val)` | simulation_name | String |
| `set_rootdir_data!(model, val)` | rootdir_data | String |
| `set_spatial_file!(model, val)` | Spatial_file | String |
| `set_connectivity_file!(model, val)` | Connectivity_file | String |
| `set_disturbance_file!(model, val)` | Disturbance_file | String |
| `set_growth_surv_file!(model, val)` | Growth_Surv_file | Vector{String} |
| `set_init_cover!(model, val)` | init_cover | String |
| `set_heat_tolerance!(model, val)` | HeatTolerance | String |
| `set_heat_init!(model, val)` | HeatInit | String |
| `set_fts!(model, val)` | fts | Vector{String} |
| `set_tradeoff!(model, val)` | TradeOff | String |
| `set_output!(model, val)` | output | Bool |
| `set_use_cached_ipm!(model, val)` | use_cached_IPM | Bool |
| `set_temp_growth_switch!(model, val)` | temp_growth_switch | Bool |
| `set_intervention!(model, val)` | Intervention | String |
| `set_intervened_sites!(model, val)` | intervened_sites | String |
| `set_sites!(model, val)` | sites | Bool/Vector |

---

### set_params! - Modify Multiple Parameters At Once

Use `set_params!` to modify multiple parameters (numeric, string, bool) in one call. **Always requires reassignment.**

```julia
model = CscapeJulia.load_params(1, fpath)

# Modify multiple parameters at once
model = set_params!(model;
    rcp = 2,
    year_end = 2050,
    plasticity = 0.5,
    region = "Moore",
    spatial_file = "New_k_MooreReefCluster",
    temp_growth_switch = true,
    use_cached_ipm = false
)

CscapeJulia.run_cscape(model)
```

**Available keyword arguments for set_params!:**

| Keyword | Maps to |
|---------|---------|
| `scenario_id` | scenario_id |
| `cyclone_rep` | Cyclone_rep |
| `rcp` | rcp |
| `year_start` | year_start |
| `year_end` | year_end |
| `plasticity` | Plasticity |
| `heritability` | Heritability |
| `dhw_enhance` | DHW_enhance |
| `draw` | draw |
| `region` | region |
| `simulation_name` | simulation_name |
| `rootdir_data` | rootdir_data |
| `spatial_file` | Spatial_file |
| `connectivity_file` | Connectivity_file |
| `disturbance_file` | Disturbance_file |
| `growth_surv_file` | Growth_Surv_file |
| `init_cover` | init_cover |
| `heat_tolerance` | HeatTolerance |
| `heat_init` | HeatInit |
| `fts` | fts |
| `tradeoff` | TradeOff |
| `output` | output |
| `use_cached_ipm` | use_cached_IPM |
| `temp_growth_switch` | temp_growth_switch |
| `intervention` | Intervention |
| `intervened_sites` | intervened_sites |
| `sites` | sites |

---

### WITHOUT Helper Functions

```julia
model = CscapeJulia.load_params(1, fpath)

# Parameter order and index:
# Index:     1            2           3      4           5         6          7            8
# Param: scenario_id, Cyclone_rep, rcp, year_start, year_end, Plasticity, DHW_enhance, draw

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

CscapeJulia.run_cscape(model)
```

### Parameter Index Reference

| Index | Parameter | Bounds | Description |
|-------|-----------|--------|-------------|
| 1 | scenario_id | (1, 1000) | Scenario ID |
| 2 | Cyclone_rep | (1, 100) | Cyclone replicate |
| 3 | rcp | (1, 4) | RCP scenario |
| 4 | year_start | (2000, 2100) | Start year |
| 5 | year_end | (2000, 2100) | End year |
| 6 | Plasticity | (0, 1) | Plasticity |
| 7 | DHW_enhance | (0.0, 10.0) | DHW enhancement |
| 8 | draw | (0, 1000) | Posterior draw (0=mean) |

---

## Example: Sensitivity Analysis

```julia
# Test different RCP scenarios (loop with modifications)
for rcp_val in 1:4
    model = CscapeJulia.load_params(1, fpath)
    set_rcp!(model, rcp_val)
    set_year_end!(model, 2050)
    CscapeJulia.run_cscape(model)
end
```

---

## Using Dict (Simpler Alternative)

```julia
input_data = CscapeJulia.load_input_data(1, fpath)
input_data["year_end"] = 2050
input_data["rcp"] = 2
input_data["region"] = "Moore"
CscapeJulia.run_cscape(input_data)
```

---

## Access Outputs

```julia
output = CscapeJulia.load_output(fpath, 1)
cover = CscapeJulia.get_cover_timeseries(output)

for (year, data) in CscapeJulia.yearly_iterator(output)
    println("Year $year: $(mean(data[:, 1, :, :, 1]))")
end
```

---

## ADRIA Indicators (Automatic)

ADRIA indicators are calculated and saved automatically when a simulation completes, both for full runs (`run_cscape`) and year-by-year runs (`finalise_simulation`). Indicators are saved as an RDS file alongside the simulation output.

**Output file:** `model_outputs/Indicators_scenario_{id}_draw_{draw}.rds`

The saved RDS contains a list with:
- `relative_cover` — Relative coral cover per site over time
- `relative_juveniles` — Relative juvenile abundance per site over time
- `relative_taxa_cover` — Relative cover by functional type
- `relative_loc_taxa_cover` — Relative cover by location and functional type
- `years`, `site_ids`, `fts` — Dimension labels

### Automatic (default)

```julia
# Indicators calculated automatically after simulation
run_cscape(1, fpath)

# Also automatic after year-by-year runs
env = initialise_simulation(1, fpath)
run_years!(env, 2008:2018)
finalise_simulation(env)   # indicators saved here
```

### Skip Indicator Calculation

```julia
run_cscape(1, fpath; calc_indicators=false)
finalise_simulation(env, calc_indicators=false)
```

### Manual Calculation

```julia
output = CscapeJulia.load_output(fpath, 1)
indicators = CscapeJulia.calculate_indicators(output)
CscapeJulia.indicator_summary(indicators)
```

### Load Saved Indicators in R

```r
indicators <- readRDS("model_outputs/Indicators_scenario_1_draw_NA.rds")
names(indicators)
# [1] "relative_cover" "relative_juveniles" "relative_taxa_cover"
# [4] "relative_loc_taxa_cover" "years" "site_ids" "fts"
```

---

## Per-Timestep Control (Year by Year)

Run simulation year by year with ability to inspect/modify state between years.

### Basic Year-by-Year Run

```julia
using .CscapeJulia

setup_r_environment()
fpath = "C:/Users/Jojo/Downloads/TestingInterventions"

# Initialise (loads data, creates arrays, NO years run yet)
env = initialise_simulation(1, fpath)

# Run year by year
for year in 2008:2018
    output = run_single_year!(env, year)
    println("Year $year complete")
end

# Save outputs (indicators calculated automatically)
finalise_simulation(env)
```

### With State Modification Mid-Simulation

```julia
env = initialise_simulation(1, fpath)

# Check state before running
print_simulation_state(env)

# Single year
@elapsed run_single_year!(env, 2008)

# Run first half - Multiple years
@elapsed run_years!(env, 2008:2012)

# Modify state mid-simulation
modify_simulation_state!(env; 
    kappa_scale = 0.9,          # Habitat degradation
    plasticity = 0.6,           # Higher plasticity
    heritability = "0.2_0.01",  # Higher heritability
    connectivity_scale = 0.8    # Reduced connectivity
)

# Run second half
@elapsed run_years!(env, 2013:2018)

# Save with different options
finalise_simulation(env)                                    # Default naming
finalise_simulation(env, filename="my_run")                 # Custom filename
finalise_simulation(env, filename="C:/results/my_run.rds")  # Full path
finalise_simulation(env, export_adria=false)                # Skip ADRIA export
finalise_simulation(env, calc_indicators=false)             # Skip indicators
```

### Mid-Simulation Modifiable Parameters

| Parameter | Julia Argument | Type | What It Does |
|-----------|----------------|------|--------------|
| Kappa scale | `kappa_scale` | Float | Multiply all kappa by factor |
| Kappa values | `kappa_values` | Vector{Float64} | Set kappa per site directly |
| DHW threshold | `dhw_threshold` | Vector{Float64} | Set threshold per site directly |
| DHW enhance | `dhw_enhance` | Float | Set DHW enhancement factor |
| Plasticity | `plasticity` | Float or String | Set plasticity value |
| Heritability | `heritability` | String | Set heritability value |
| Connectivity scale | `connectivity_scale` | Float | Multiply connectivity by factor |
| Fogging reduction | `fogging_reduction` | Float | Set fogging heat reduction |
| Rubble handle | `rubble_handle` | Bool | Toggle rubble dynamics |
| Coral deployment | `coral_deployment` | Dict | Replace coral deployment schedule |
| Fogging schedule | `fogging_schedule` | Dict | Replace fogging schedule |
| COTS mortality | `cots_mortality` | Dict | Modify COTS pressure in temporal data |

### Adaptive Coral Deployment (Julia-only)

Julia can inspect site cover mid-simulation and adaptively redirect coral deployment 
to sites with available capacity.

```julia
env = initialise_simulation(1, fpath)

# Run first half - Multiple years
@elapsed run_years!(env, 2008:2012);

cover_state = get_site_cover(env);
available = findall(i -> cover_state[:proportion_full][i] < 0.8, 
                    1:length(cover_state[:site_names]))
deploy_sites = cover_state[:site_names][available[1:min(5, length(available))]]

n = length(deploy_sites)
new_deployment = Dict(
    "reef_siteid"       => deploy_sites,
    "Year"              => fill(2013, n),
    "ft"                => fill(1, n),          # functional type 1
    "no_int_corals"     => fill(1000.0, n),     # 1000 corals per site
    "m2"                => fill(500.0, n),       # 500 m² deployment area
    "meshpt_int_corals" => fill("4", n),         # size class
    "Enhancement"       => fill(3, n)            # enhancement class 3
)

modify_simulation_state!(env; coral_deployment=new_deployment)


run_years!(env, 2013:2018)
finalise_simulation(env)
```

### Modifying Fogging Schedule Mid-Simulation

```julia
# Expand fogging to new sites based on bleaching risk
modify_simulation_state!(env;
    fogging_schedule = Dict(
        "reef_siteid" => ["Moore_001", "Moore_002", "Moore_003"],
        "year"        => [2015, 2015, 2015],
        "reduction"   => [0.3, 0.3, 0.25]
    )
)
```

### Modifying COTS  Mid-Simulation

```julia
# Scale COTS for a specific year (e.g. simulate COTS control program)
modify_simulation_state!(env;
    cots_mortality = Dict("year" => 2015, "scale" => 0.5)
)

# Or set COTS directly at specific sites
modify_simulation_state!(env;
    cots_mortality = Dict(
        "year"         => 2015,
        "site_indices" => [1, 2, 3],
        "values"       => [0.0 0.0 0.0;    # bin1
                           0.0 0.0 0.0;    # bin2
                           0.0 0.0 0.0;    # bin3
                           0.0 0.0 0.0]    # bin4
    )
)
```

### Inspect State During Simulation

```julia
env = initialise_simulation(1, fpath)

# Check simulation state (rcopy returns Symbol keys)
state = get_simulation_state(env)
println("Sites: $(state[:n_sites])")
println("FTs: $(state[:fts])")

# Check site cover
cover = get_site_cover(env)
println("First site: $(cover[:site_names][1]), full: $(round(cover[:proportion_full][1]*100, digits=1))%")

# Run and get output
run_single_year!(env, 2008)
output = get_simulation_output(env)
```

### Per-Timestep Functions Reference

| Function | Description |
|----------|-------------|
| `initialise_simulation(id, fpath)` | Setup without running any years |
| `run_single_year!(env, year)` | Run single year, returns year output |
| `run_years!(env, 2008:2012)` | Run multiple years |
| `get_simulation_output(env)` | Get full output array |
| `get_simulation_state(env)` | Get current state info |
| `get_site_cover(env)` | Get cover, kappa, remaining capacity per site |
| `modify_simulation_state!(env; ...)` | Modify state between years |
| `finalise_simulation(env)` | Save outputs, calculate indicators, cleanup |

### Compare: Full Run vs Year-by-Year

**Full run (simple):**
```julia
run_cscape(1, fpath)
```

**Year-by-year (flexible):**
```julia
env = initialise_simulation(1, fpath)
for year in 2008:2018
    run_single_year!(env, year)
end
finalise_simulation(env)
```

---

## R Equivalent (Year by Year)

```r
source("modules/run_year.R")
source("modules/save_outputs.R")

MainEnvir <- initialise_simulation(InputData, fun_path)

for (year in 2008:2018) {
  run_single_year(year, MainEnvir)
}

save_outputs(MainEnvir$out_array, MainEnvir$InputData, NA)
```