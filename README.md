# CscapeInterface.jl

Julia interface to the C-scape individual-based coral reef simulation model. Provides R integration via RCall, parallel batch scenario execution, MCDA-guided adaptive management workflows, and interoperability with [ADRIA.jl](https://github.com/open-AIMS/ADRIA.jl) for sensitivity analysis and visualisation.

---

## Installation

```julia
using Pkg
Pkg.activate("path/to/CscapeInterface.jl")
Pkg.instantiate()
```

**Requirements:**
- Julia ≥ 1.9
- The C-scape R package installed locally. Its root path (the folder containing `run_cscape.R`, `ipm_pred.R`, etc.) is passed to `setup_r_environment` at the start of each session.

---

## Required Input Files

All input files live in `fpath/data/`. The `ScenarioID.xlsx` file in `fpath/` records which file each scenario uses.

| File | ScenarioID column | Description |
|------|-------------------|----|
| `<name>.RData` | `Spatial_file` | Reef site geometries, kappa, area, depths |
| `<name>.RData` | `Connectivity_file` | Larval connectivity matrices (one per FT) |
| `<name>.RData` | `Disturbance_file` | Temporal DHW and disturbance time series |
| `demog_inputs_<ft>*.rds` | `Growth_Surv_file` | **Demographic IPM kernel, one file per functional type** |

### Demographic files

Demographic files are prepared independently of the C-scape simulation package using a dedicated IPM fitting workflow — they are not generated automatically. One RDS is needed per functional type; `Growth_Surv_file` stores all filenames for a scenario separated by `/`:

```
demog_inputs_acro_tableALLDAT.rds/demog_inputs_acro_corymALLDAT.rds/demog_inputs_corym_non_acroALLDAT.rds/...
```

Fields read by Julia (from `src/data_access.jl`):

| Field | Type | Description |
|-------|------|-------------|
| `meshpoints_diam` | numeric vector, length 100 | Size-class midpoints in diameter (cm) |
| `fecundity_eggs` | matrix | Fecundity-by-size; columns with `colSums == 0` are flagged as juvenile classes |

The C-scape R engine (`ipm_pred.R`) reads additional fields from the same file (growth kernels, survival vectors, recruitment distributions). See the C-scape R package documentation for the complete IPM input specification.

> **Note:** filenames in `Growth_Surv_file` are case-sensitive — they must match the on-disk filenames exactly.

---

## Output Files

Each completed scenario writes up to three files, all inside `fpath`:

| | File pattern | Folder | Format | Size |
|---|---|---|---|---|
| 1 | `Array_scenario_{id}_draw_{draw}.rds` | `model_outputs/` | R RDS | Large |
| 2 | `adria_scenario_{id}_draw_{draw}.jld2` | `adria_exports/` | Julia JLD2 | Very large |
| 3 | `Indicators_scenario_{id}_draw_{draw}.jld2` | `adria_exports/` | Julia JLD2 | Small |

### Option 1 — C-scape raw output (`model_outputs/`)

The original simulation array written by the C-scape R package. Contains the full 6D array `[years × sites × intervention × functional_types × enhancement × 106 metrics]`. This is the source of truth for all downstream products.

```julia
output = load_output(fpath, 1)             # load scenario 1, draw = NA
output = load_output(fpath, 1; draw="5")   # specific posterior draw
```

### Option 2 — ADRIA full export (`adria_exports/adria_*.jld2`)

The complete `CscapeOutput` Julia struct serialised to JLD2. Preserves every field (out_array, kappa, area, meshpoints, spatial metadata) so Julia-side analysis can be re-run without touching the R output again.

```julia
using JLD2
output = load(joinpath(fpath, "adria_exports", "adria_scenario_1_draw_NA.jld2"), "cscape_output")
```

> **Recommendation for large batch runs:** these files are very large. For hundreds or thousands of scenarios, disable the full export and keep only the indicator files (Option 3) — `CScapeResultSet` analysis only needs those:
> ```julia
> run_cscape(sid, fpath; export_adria=false, calc_indicators=true)
> ```

### Option 3 — ADRIA indicators (`adria_exports/Indicators_*.jld2`)

Pre-computed summary metrics produced by `ADRIAIndicators.jl`. Small and self-contained. This is the file scanned by `load_results(CScapeResultSet, fpath)` to build a `CScapeResultSet` for sensitivity analysis and visualisation.

Each file stores:

| Indicator | Dimensions |
|-----------|------------|
| `relative_cover` | timesteps × locations |
| `relative_taxa_cover` | timesteps × functional groups |
| `relative_juveniles` | timesteps × locations |
| `coral_diversity` | timesteps × locations |
| `coral_evenness` | timesteps × locations |
| `relative_shelter_volume` | timesteps × groups × sizes × locations |
| `reef_biodiversity_condition_index` | timesteps × locations |
| `reef_condition_index` | timesteps × locations |
| spatial metadata | kappa, area, meshpoints, site IDs |

```julia
# load_results scans adria_exports/ for all Indicators_*.jld2 files
rs = load_results(CScapeResultSet, fpath)
```

---

## Workflow 1 — Build Input Data and Run Scenarios

> Source: `sandbox/CscapeInterface_MCDA_PAWN/Setup/DataWorkflow.jl`

### Step 1 — Create C-scape input data (R via Julia)

```julia
using RCall, CscapeInterface

fun_path = "C:/path/to/C_scape"          # root of the C-scape R package
fpath    = "C:/path/to/model_runs/MyRun"  # output directory

# Setup R and source the build helpers
setup_r_environment(fun_path)

@rput fun_path
R"""
setwd(fun_path)
source("simulated_reef/build_spatial_file.R")
source("simulated_reef/build_connectivity.R")
source("simulated_reef/build_temporal.R")
source("simulated_reef/build_scenario_table.R")
library(sf)

output_dir  <- file.path("C:/path/to/model_runs/MyRun", "data")
start_year  <- 2025
end_year    <- 2050

# 1. Create spatial data
spatial_data <- create_spatial_data(...)

# 2. Create connectivity matrix
connectivity_matrix <- create_connectivity_matrix(spatial_file_path = ..., output_path = ...)

# 3. Create temporal/disturbance files
temporal_results <- create_temporal_file(spatial_file_path = ..., ...)

# 4. Build ScenarioID.xlsx
scenario_data <- rbind(scenario_data_int, scenario_data_counter)
scenario_data$ID <- 1:nrow(scenario_data)
create_scenario_table(data = scenario_data,
                      output_path = file.path(output_dir, "ScenarioID.xlsx"),
                      sheet_name  = "ScenarioID")
"""
@info "Input data created"
```

See `DataWorkflow.jl` for a complete example that builds spatial, connectivity, and temporal files for a simulated 50-site reef with 4 disturbance levels and intervention/counterfactual scenarios.

### Step 2 — Run scenarios in parallel

```julia
using Distributed, ProgressMeter, CscapeInterface, RCall

const fpath    = "C:/path/to/model_runs/MyRun"
const fun_path = "C:/path/to/C_scape"

# Read all scenario IDs
@rput fpath
R"library(readxl); .scn <- read_excel(paste0($fpath, '/data/ScenarioID.xlsx'), sheet='ScenarioID')"
all_ids = Int.(rcopy(R"as.integer(.scn$ID)"))
@info "Running $(length(all_ids)) scenarios"

# Spawn workers
n_workers = max(1, Sys.CPU_THREADS - 1)
addprocs(n_workers)

@everywhere begin
    using CscapeInterface
    _fpath    = "C:/path/to/model_runs/MyRun"
    _fun_path = "C:/path/to/C_scape"
end

# Progress bar
p           = Progress(length(all_ids); desc="Scenarios: ", showspeed=true)
progress_ch = RemoteChannel(() -> Channel{Bool}(length(all_ids)))
@async for _ in all_ids; take!(progress_ch); next!(p); end

failed = pmap(all_ids) do sid
    result = try
        setup_r_environment(_fun_path; enable_parallel=false)
        run_cscape(sid, _fpath; export_adria=true, calc_indicators=true)
        nothing
    catch e
        @error "Scenario $sid failed" exception=e
        sid
    end
    put!(progress_ch, true)
    result
end |> x -> filter(!isnothing, x)

finish!(p)
isempty(failed) ? @info("All scenarios complete") : @warn("Failed: $failed")
rmprocs(workers())
```

**Memory per worker:** each worker runs an independent R session with its own copy of all C-scape state (IPM kernels, temporal arrays, output arrays). Observed usage is **~12 GB RAM per worker** for typical runs (50 sites, 6 FTs, 25 years). As a starting point use:

```julia
n_workers = max(1, floor(Int, Sys.total_memory() / (12 * 2^30)) - 1)
```

Using `enable_parallel=false` inside each worker prevents R from spawning additional sub-threads, keeping memory usage predictable.

### Step 3 — Load and save results

```julia
using CscapeInterface

# Scan the results directory and build a CScapeResultSet
rs = load_results(CScapeResultSet, fpath)

# Persist to a JLD2 file for fast future reloads
save_results(rs, joinpath(fpath, "all_results.jld2"))

# Fast reload — skips all indicator file scanning
rs = load_results(CScapeResultSet, joinpath(fpath, "all_results.jld2"))
```

| Function | Description |
|----------|-------------|
| `load_results(CScapeResultSet, fpath)` | Scan directory; load all indicator RDS files |
| `load_results(CScapeResultSet, "file.jld2")` | Fast reload from a saved JLD2 |
| `load_grouped_results(fpath)` | Load separately by spatial/connectivity group |
| `load_grouped_results(fpath; group=1)` | Load a specific group |
| `save_results(rs, "file.jld2")` | Save for fast future reload |
| `list_groups(fpath)` | Preview available groups without loading data |

---

## Workflow 2 — MCDA Dynamic Reranking

> Source: `sandbox/CscapeInterface_MCDA_PAWN/Setup/Workflow_MCDA_Evaluation.jl`

This workflow runs repeated MCDA → simulation cycles. At the end of each cycle, the sites selected for coral deployment are re-evaluated based on updated cover, so the intervention adapts over time.

### Configuration

```julia
using CscapeInterface

cfg = Dict{String, Any}(
    "scenario_id"    => 20,
    "fpath"          => "C:/path/to/model_runs/MyRun",
    "fun_path"       => "C:/path/to/C_scape",
    "output_dirname" => "MCDA_results",
    "n_workers"      => 2,
    "shared_years"   => 2008:2024,   # burn-in years run once for all variants
    "n_loops"        => 3,           # MCDA→simulation cycles
    "loop_duration"  => 5,           # years per cycle  →  3 × 5 = 15 intervention years
    "skip_shared"    => false,       # set true to reuse a previously saved baseline
    "variants"       => [
        # Intervention variant: deploy to top-5 ranked sites
        Dict{String,Any}(
            "n_sites"            => 5,
            "selection"          => "top",
            "ft"                 => 1,
            "no_int_corals"      => 5000.0,
            "proportion"         => 1.0,
            "m2"                 => 3000.0,
            "density"            => 5000/3000,
            "meshpt_int_corals"  => "1.21",
            "Enhancement"        => 3
        ),
        # Counterfactual: no intervention
        Dict{String,Any}(
            "n_sites"   => 0,
            "selection" => "counterfactual",
            "ft"        => 1, "no_int_corals" => 0, "proportion" => 0,
            "m2" => 0, "density" => 0, "meshpt_int_corals" => "1.21", "Enhancement" => 0
        ),
    ]
)
```

### Running

```julia
@elapsed CscapeInterface.run_dynamic_reranking(cfg)
```

### Analysing outputs

```julia
using CscapeInterface, Statistics

_fpath  = cfg["fpath"]
_outdir = cfg["output_dirname"]
_sid    = cfg["scenario_id"]

outdir = joinpath(_fpath, "adria", _outdir)

# Load intervention and counterfactual outputs
iter1_out = load_output(_fpath, _sid; filename = joinpath(outdir, "Array_scenario_$(_sid)_iter_1.rds"))
iter2_out = load_output(_fpath, _sid; filename = joinpath(outdir, "Array_scenario_$(_sid)_iter_2.rds"))

# Mean relative coral cover (cover / kappa) across sites over time
function mean_rel_cover(out::CscapeOutput)
    cover = get_cover_timeseries(out; intervention_idx=3)   # [years × sites]
    return vec(mean(cover ./ out.kappa', dims=2))
end

rel_iter1 = mean_rel_cover(iter1_out)   # intervention
rel_iter2 = mean_rel_cover(iter2_out)   # counterfactual
years     = iter1_out.years
```

### Visualising

```julia
using Plots

# Time-series comparison
p1 = plot(years, rel_iter1; label="Top 5 (intervened)", lw=2, color=:steelblue)
plot!(p1, years, rel_iter2; label="Counterfactual",      lw=2, color=:coral)
xlabel!(p1, "Year"); ylabel!(p1, "Mean cover (%)"); title!(p1, "Coral cover / habitable area")
savefig(p1, joinpath(outdir, "cover_timeseries.png"))

# Difference
p2 = plot(years, rel_iter1 .- rel_iter2;
    label="Top 5 vs counterfactual", lw=2, color=:steelblue)
xlabel!(p2, "Year"); ylabel!(p2, "Δ Mean cover (%)")
savefig(p2, joinpath(outdir, "cover_difference.png"))
```

Per-loop intervention site maps are generated by reading the `deployment_iter_1.csv` output and the reef polygon geometry from `ScenarioID.xlsx`. See `Workflow_MCDA_Evaluation.jl` for the full map-plotting code.

---

## Workflow 3 — Visualisation and Sensitivity Analysis (ADRIA)

> Source: `sandbox/CscapeInterface_MCDA_PAWN/Analysis/Visualisation.jl`

> **Important:** WGLMakie (or GLMakie) must be loaded *before* ADRIA. This triggers ADRIA's `AvizExt` extension which provides all `ADRIA.viz.*` functions.

### Load results

```julia
using WGLMakie, GeoMakie, GraphMakie   # must come before ADRIA
WGLMakie.activate!()
using ADRIA, CscapeInterface, Statistics

fpath = "C:/path/to/model_runs/MyRun"
rs    = load_results(CScapeResultSet, joinpath(fpath, "all_results.jld2"))

# Inspect available content
println("Locations : ", n_locations(rs))
println("Scenarios : ", n_scenarios(rs))
println("Timesteps : ", collect(timesteps(rs)))
println("Outcomes  : ", sort(collect(keys(rs.outcomes))))
println("Inputs    : ", names(rs.inputs))
```

### Scenario grouping

Group scenarios by any column in `rs.inputs` — here by intervention status:

```julia
let interv = uppercase.(strip.(string.(rs.inputs[!, :Intervention])))
    global scen_groups = filter!(
        kv -> any(kv.second),
        Dict{Symbol,BitVector}(
            :counterfactual => interv .== "NO",
            :guided         => interv .== "YES",
        )
    )
end
```

### Scenario time-series plots

```julia
rel_cover = scenario_outcome(rs, :relative_cover)   # area-weighted mean [T × S]

fig = Figure(size=(900, 400))
g   = fig[1, 1] = GridLayout()
ax  = Axis(g[1, 1]; title="Mean relative coral cover", xlabel="Year", ylabel="Relative cover")
ADRIA.viz.scenarios!(g, ax, rel_cover, scen_groups;
    opts=Dict{Symbol,Any}(:by_RCP => true, :histogram => false))
```

### Taxonomy (functional groups)

```julia
taxa_cover = rs.outcomes[:relative_taxa_cover]
fig = Figure(size=(1200, 600))
ADRIA.viz.taxonomy!(fig[1,1] = GridLayout(), taxa_cover, scen_groups)
```

### PAWN global sensitivity analysis

```julia
ts_all   = collect(timesteps(rs))
ts_range = max(1, length(ts_all)-4):length(ts_all)

# Scalar outcome: mean over last 5 timesteps per scenario
y_cover = vec(mean(parent(rel_cover)[ts_range, :]; dims=1))

pawn_si  = ADRIA.sensitivity.pawn(rs, y_cover)
pawn_fig = ADRIA.viz.pawn(pawn_si)
display(pawn_fig)
```

### Temporal Sensitivity Analysis (TSA)

```julia
tsa_si  = ADRIA.sensitivity.tsa(rs, parent(rel_cover))
tsa_fig = ADRIA.viz.tsa(rs, tsa_si)
display(tsa_fig)
```

---

## Working with CScapeResultSet

| Accessor | Returns |
|----------|---------|
| `n_locations(rs)` | Number of reef sites |
| `n_scenarios(rs)` | Number of scenarios |
| `timesteps(rs)` | Year vector |
| `loc_k(rs)` | Carrying capacity per site |
| `rs.outcomes[:relative_cover]` | YAXArray `(timesteps × locations × scenarios)` |
| `scenario_outcome(rs, :relative_cover)` | Area-weighted mean over locations `(timesteps × scenarios)` |
| `rs.inputs` | DataFrame of all input parameters |

---

## Advanced: Per-Timestep Control

Run a simulation year by year, inspecting or modifying state between years.

```julia
using CscapeInterface

setup_r_environment(fun_path)
fpath = "C:/path/to/model_runs/MyRun"

env = initialise_simulation(1, fpath)    # load data, allocate arrays — no years run yet

run_years!(env, 2008:2012)               # run first half

# Inspect current state
print_simulation_state(env)
cover = get_site_cover(env)              # returns site names, kappa, proportion_full

# Modify mid-simulation
modify_simulation_state!(env;
    kappa_scale        = 0.9,            # habitat degradation
    plasticity         = 0.6,
    connectivity_scale = 0.8
)

run_years!(env, 2013:2018)               # run second half

finalise_simulation(env)                 # save output + calculate indicators
```

### Mid-simulation modifiable parameters

| Argument | Type | What it does |
|----------|------|--------------|
| `kappa_scale` | Float | Multiply all kappa by a factor |
| `kappa_values` | Vector{Float64} | Set kappa per site directly |
| `dhw_threshold` | Vector{Float64} | Set bleaching threshold per site |
| `dhw_enhance` | Float | DHW enhancement factor |
| `plasticity` | Float or String | Plasticity value |
| `heritability` | String | Heritability value |
| `connectivity_scale` | Float | Scale connectivity matrix |
| `fogging_reduction` | Float | Fogging heat reduction |
| `coral_deployment` | Dict | Replace deployment schedule |
| `fogging_schedule` | Dict | Replace fogging schedule |
| `cots_mortality` | Dict | Modify COTS pressure in temporal data |

### Per-timestep function reference

| Function | Description |
|----------|-------------|
| `initialise_simulation(id, fpath)` | Set up without running any years |
| `run_single_year!(env, year)` | Run a single year, return year output |
| `run_years!(env, 2008:2018)` | Run a range of years |
| `get_simulation_state(env)` | Current state info (sites, FTs, params) |
| `get_site_cover(env)` | Cover, kappa, remaining capacity per site |
| `modify_simulation_state!(env; ...)` | Modify state between years |
| `finalise_simulation(env)` | Save output, calculate indicators, cleanup |

---

## Parameter Reference

### Numeric helpers — modify in place, no reassignment needed

```julia
model = load_params(1, fpath)
set_rcp!(model, 2)
set_year_end!(model, 2050)
set_plasticity!(model, 0.5)
run_cscape(model)
```

| Function | Parameter | Bounds |
|----------|-----------|--------|
| `set_scenario_id!(model, val)` | scenario_id | (1, 1000) |
| `set_cyclone_rep!(model, val)` | Cyclone_rep | (1, 100) |
| `set_rcp!(model, val)` | rcp | (1, 4) |
| `set_year_start!(model, val)` | year_start | (2000, 2100) |
| `set_year_end!(model, val)` | year_end | (2000, 2100) |
| `set_plasticity!(model, val)` | Plasticity | (0, 1) |
| `set_dhw_enhance!(model, val)` | DHW_enhance | (0, 10) |
| `set_draw!(model, val)` | draw | (0, 1000) |
| `set_heritability!(model, val)` | Heritability | (0, 1) |

### String/Bool/Vector helpers — must reassign

```julia
model = set_region!(model, "Moore")
model = set_intervention!(model, "Yes")
model = set_intervened_sites!(model, "Site1/Site2/Site3")
model = set_growth_surv_file!(model, ["demog_acro.rds", "demog_massive.rds"])
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
| `set_sites!(model, val)` | sites | Bool / Vector |

### `set_params!` — modify multiple parameters at once

Always requires reassignment:

```julia
model = set_params!(model;
    rcp               = 2,
    year_end          = 2050,
    plasticity        = 0.5,
    region            = "Moore",
    spatial_file      = "reef_sites_1_nogeo.RData",
    temp_growth_switch = true
)
```

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
