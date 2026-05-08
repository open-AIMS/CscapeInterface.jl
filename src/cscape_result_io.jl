using ADRIA: ResultSet, EnvLayer, SimConstants
import ADRIA.sensitivity: pawn
using JLD2
using DataFrames
using YAXArrays
using DimensionalData: Dim
using RCall
using XLSX

# ─────────────────────────────────────────────────────────────────────────────
# Struct
# ─────────────────────────────────────────────────────────────────────────────

"""
    CScapeResultSet

Multi-scenario result container for C~scape simulations. Extends ADRIA's `ResultSet`
interface so results can be used with ADRIA analysis tooling.

Loaded from self-contained indicator JLD2 files produced by `_calculate_and_save_indicators`.
Each JLD2 stores the `CscapeIndicators` plus all spatial/domain metadata required to
build this struct — the large output JLD2 is never opened during loading.
"""
struct CScapeResultSet <: ResultSet
    name::String
    RCP::String

    loc_ids                          # Vector{String}
    loc_area::Vector{Float64}
    loc_max_coral_cover::Vector{Float64}
    loc_centroids                    # Matrix{Float64} [n_sites × 2] (lon, lat), or nothing
    env_layer_md::EnvLayer
    connectivity_data                # DataFrame
    loc_data                         # DataFrame: spatial table (from output.spatial)

    inputs::DataFrame                # one row per scenario (scenario_id, draw, ...)
    sim_constants                    # SimConstants()
    model_spec::DataFrame

    outcomes                         # Dict{Symbol, YAXArray}
    coral_size_diameter              # Matrix{Float64} = meshpoints [ft × sizes]

    raw_output_paths::Vector{String} # paths to output JLD2 files (on-demand raw access)
end

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    _parse_indicator_filename(fpath) -> (scenario_id::Int, draw::String)

Extract scenario ID and draw tag from an indicator JLD2 filename.
Standard pattern: `Indicators_scenario_{id}_draw_{draw}.jld2`
MCDA pattern:     `Indicators_Array_scenario_{id}_iter_{n}.jld2`

Falls back to loading the JLD2 metadata when neither pattern matches.
"""
function _parse_indicator_filename(fpath::String)
    fname = basename(fpath)
    # Standard pattern
    m = match(r"[Ii]ndicators_scenario_(\d+)_draw_(.+)\.jld2", fname)
    !isnothing(m) && return parse(Int, m[1]), String(m[2])
    # MCDA output-collocated pattern (Indicators_Array_scenario_{id}_iter_{n}.jld2)
    m2 = match(r"scenario_(\d+)", fname)
    !isnothing(m2) && return parse(Int, m2[1]), "NA"
    # Last resort: open the file and read metadata
    try
        meta = JLD2.load(fpath, "metadata")
        return Int(meta["scenario_id"]), string(get(meta, "draw", "NA"))
    catch
    end
    error("Cannot parse scenario ID from indicator filename: $fname")
end

"""
    _get_result_paths(result_dir) -> Vector{String}

Return sorted paths to all indicator JLD2 files in `result_dir`.
Matches both the standard name (`Indicators_scenario_*.jld2`) and the
output-collocated name produced by the MCDA workflow
(`Indicators_Array_scenario_*.jld2`).
"""
function _get_result_paths(result_dir::String)::Vector{String}
    !isdir(result_dir) && throw(ArgumentError("Directory not found: $result_dir"))
    files = filter(isfile, readdir(result_dir; join=true))
    files = filter(x -> occursin(r"[Ii]ndicators_.*\.jld2", x), files)
    isempty(files) && throw(ArgumentError("No indicator JLD2 files found in: $result_dir"))
    return sort(files)
end

"""
    _infer_output_path(indicator_path) -> String

Derive the corresponding output JLD2 path from an indicator JLD2 path.
"""
function _infer_output_path(indicator_path::String)::String
    dir = dirname(indicator_path)
    fname = basename(indicator_path)
    output_fname = replace(fname, r"[Ii]ndicators_" => "adria_"; count=1)
    return joinpath(dir, output_fname)
end

"""
    _load_scenario_spec(data_dir) -> DataFrame

Read `ScenarioID.xlsx` from `data_dir` into a DataFrame.
"""
function _load_scenario_spec(data_dir::String)::DataFrame
    xlsx_path = joinpath(data_dir, "ScenarioID.xlsx")
    isfile(xlsx_path) || error("ScenarioID.xlsx not found in: $data_dir")
    xf = XLSX.readxlsx(xlsx_path)
    sheet = xf[XLSX.sheetnames(xf)[1]]
    data = XLSX.getdata(sheet)
    headers = Symbol.(string.(data[1, :]))
    rows = [data[i, :] for i in 2:size(data, 1)]
    df = DataFrame(
        [headers[j] => [rows[i][j] for i in eachindex(rows)] for j in eachindex(headers)]...
    )
    return df
end

"""
    _get_connectivity_col(scenario_spec) -> Symbol

Find the connectivity file column in the scenario spec (handles case variation).
"""
function _get_connectivity_col(scenario_spec::DataFrame)::Symbol
    cols = names(scenario_spec)
    idx = findfirst(c -> lowercase(c) == "connectivity_file", lowercase.(cols))
    isnothing(idx) && error("Cannot find Connectivity_file column in ScenarioID.xlsx. Columns: $(cols)")
    return Symbol(cols[idx])
end

"""
    _get_spatial_col(scenario_spec) -> Symbol

Find the spatial file column in the scenario spec.
"""
function _get_spatial_col(scenario_spec::DataFrame)::Symbol
    cols = names(scenario_spec)
    idx = findfirst(c -> lowercase(c) == "spatial_file", lowercase.(cols))
    isnothing(idx) && error("Cannot find Spatial_file column in ScenarioID.xlsx. Columns: $(cols)")
    return Symbol(cols[idx])
end

"""
    _validate_domain_uniformity(jld2_data, scenario_spec)

Error if scenarios span different site lists or connectivity files — each `CScapeResultSet`
represents one domain. Use `load_grouped_results` to handle mixed directories.
"""
function _validate_domain_uniformity(
    jld2_data::Vector{<:Dict}, scenario_spec::DataFrame
)::Nothing
    conn_col = _get_connectivity_col(scenario_spec)
    id_col = first(filter(c -> lowercase(c) == "id", names(scenario_spec)))

    ref_sites = jld2_data[1]["cscape_indicators"].site_ids
    ref_scen_id = jld2_data[1]["metadata"]["scenario_id"]
    ref_row = scenario_spec[scenario_spec[!, Symbol(id_col)] .== ref_scen_id, :]
    ref_conn = isempty(ref_row) ? "" : string(ref_row[1, conn_col])

    for d in jld2_data[2:end]
        sites = d["cscape_indicators"].site_ids
        scen_id = d["metadata"]["scenario_id"]
        row = scenario_spec[scenario_spec[!, Symbol(id_col)] .== scen_id, :]
        conn = isempty(row) ? "" : string(row[1, conn_col])

        if sites != ref_sites
            error(
                "Scenarios have different locations (site count: $(length(ref_sites)) vs " *
                "$(length(sites))). Use `load_grouped_results` to auto-group by domain."
            )
        end
        if conn != ref_conn
            error(
                "Scenarios use different connectivity files ($ref_conn vs $conn). " *
                "Use `load_grouped_results` to auto-group by domain."
            )
        end
    end
    return nothing
end

"""
    _load_connectivity(data_dir, scenario_spec, scenario_id) -> DataFrame

Load the connectivity matrix from the RData/RDS file specified in the scenario spec.
"""
function _load_connectivity(
    data_dir::String, scenario_spec::DataFrame, scenario_id::Int
)
    conn_col = _get_connectivity_col(scenario_spec)
    id_col = first(filter(c -> lowercase(c) == "id", names(scenario_spec)))
    row = scenario_spec[scenario_spec[!, Symbol(id_col)] .== scenario_id, :]
    isempty(row) && error("Scenario $scenario_id not found in ScenarioID.xlsx")

    conn_rel = string(row[1, conn_col])
    isempty(conn_rel) && return DataFrame()

    conn_file = joinpath(data_dir, "data", conn_rel)
    isfile(conn_file) || (conn_file = joinpath(data_dir, conn_rel))
    isfile(conn_file) || (@warn "Connectivity file not found: $conn_file"; return DataFrame())

    @rput conn_file
    R"""
    conn_R <- tryCatch(
        readRDS(conn_file),
        error = function(e) {
            env_tmp <- new.env()
            load(conn_file, envir = env_tmp)
            as.data.frame(env_tmp[[ls(env_tmp)[1]]])
        }
    )
    """
    return rcopy(R"conn_R")
end

"""
    _build_env_layer(data_dir, conn_path, years) -> EnvLayer
"""
function _build_env_layer(
    data_dir::String, conn_rel::String, years::Vector{Int}
)::EnvLayer
    conn_path = isempty(conn_rel) ? "" : joinpath(data_dir, "data", conn_rel)
    return EnvLayer(
        data_dir, "spatial", "reef_siteid", "", "", conn_path, "", "",
        years[1]:years[end]
    )
end

"""
    _build_inputs_df(jld2_data, scenario_spec) -> DataFrame

Build an inputs DataFrame (one row per scenario) from JLD2 metadata + scenario spec.
"""
function _build_inputs_df(
    jld2_data::Vector{<:Dict}, scenario_spec::DataFrame
)::DataFrame
    id_col = first(filter(c -> lowercase(c) == "id", names(scenario_spec)))
    rows = DataFrame[]
    for d in jld2_data
        meta = d["metadata"]
        scenario_id = meta["scenario_id"]
        draw = meta["draw"]

        row_df = DataFrame(scenario_id=[scenario_id], draw=[draw])

        spec_rows = scenario_spec[scenario_spec[!, Symbol(id_col)] .== scenario_id, :]
        if !isempty(spec_rows)
            for col in names(spec_rows)
                lowercase(col) == "id" && continue
                row_df[!, col] = [spec_rows[1, col]]
            end
        end
        push!(rows, row_df)
    end
    return vcat(rows...; cols=:union)
end

"""
    _build_model_spec(inputs) -> DataFrame

Build a minimal model spec DataFrame from the inputs DataFrame columns.
"""
function _build_model_spec(inputs::DataFrame)::DataFrame
    cols = names(inputs)
    return DataFrame(
        component  = fill("CScape", length(cols)),
        fieldname  = Symbol.(cols),
        name       = cols,
        ptype      = fill("continuous", length(cols)),
        is_constant = fill(false, length(cols))
    )
end

"""
    _build_outcomes(indicators_list, years, site_ids, fts) -> Dict{Symbol, YAXArray}

Stack per-scenario `CscapeIndicators` along a new scenarios axis and return a
`Dict{Symbol, YAXArray}` compatible with ADRIA's outcomes interface.
"""
function _build_outcomes(
    indicators_list::Vector{CscapeIndicators},
    years::Vector{Int},
    site_ids::Vector{String},
    fts::Vector{String}
)::Dict{Symbol, YAXArray}
    n_s = length(indicators_list)
    outcomes = Dict{Symbol, YAXArray}()

    ax_t  = Dim{:timesteps}(years)
    ax_l  = Dim{:locations}(site_ids)
    ax_sc = Dim{:scenarios}(1:n_s)
    ax_sp = Dim{:species}(1:length(fts))  # integer labels — ADRIA.viz.taxonomy! requires Int species axis

    # 2D fields [timesteps × locations] → stack → [timesteps, locations, scenarios]
    for field in (
        :relative_cover, :ltmp_cover, :coral_diversity, :coral_evenness,
        :relative_juveniles, :juvenile_indicator,
        :reef_biodiversity_condition_index, :reef_condition_index,
        :reef_fish_index, :reef_tourism_index
    )
        data = cat([getfield(ind, field) for ind in indicators_list]...; dims=3)
        outcomes[field] = YAXArray((ax_t, ax_l, ax_sc), data)
    end

    # 3D fields [timesteps × groups × locations] → [timesteps, species, locations, scenarios]
    for field in (:relative_loc_taxa_cover, :relative_loc_taxa_juveniles)
        data = cat([getfield(ind, field) for ind in indicators_list]...; dims=4)
        outcomes[field] = YAXArray((ax_t, ax_sp, ax_l, ax_sc), data)
    end

    # 2D fields [timesteps × groups] → [timesteps, species, scenarios]
    for field in (:relative_taxa_cover, :relative_taxa_juveniles)
        data = cat([getfield(ind, field) for ind in indicators_list]...; dims=3)
        outcomes[field] = YAXArray((ax_t, ax_sp, ax_sc), data)
    end

    # 4D [timesteps × groups × sizes × locations] → sum groups+sizes → [timesteps, locations, scenarios]
    rsv_data = cat(
        [dropdims(sum(ind.relative_shelter_volume, dims=(2, 3)), dims=(2, 3))
         for ind in indicators_list]...;
        dims=3
    )
    outcomes[:relative_shelter_volume] = YAXArray((ax_t, ax_l, ax_sc), rsv_data)

    return outcomes
end

# ─────────────────────────────────────────────────────────────────────────────
# load_results
# ─────────────────────────────────────────────────────────────────────────────

"""
    load_results(CScapeResultSet, path; scenario_ids) -> CScapeResultSet
    load_results(CScapeResultSet, data_dir, result_dir; scenario_ids) -> CScapeResultSet
    load_results(CScapeResultSet, data_dir, indicator_files) -> CScapeResultSet

Load a multi-scenario `CScapeResultSet` from indicator JLD2 files, or from a previously
saved ResultSet JLD2 (see `save_results`).

When `path` is a `.jld2` file written by `save_results`, it is loaded directly.
When `path` is a directory, indicator files are discovered from `<path>/adria_exports/`.

Only indicator files that **exist on disk** are loaded — the ScenarioID table is used
solely to enrich metadata, not to decide which scenarios to load. To restrict which
scenarios are loaded, pass `scenario_ids` (a vector or range of integer IDs).

# Examples
```julia
rs = load_results(CScapeResultSet, "/path/to/data")                      # all found files
rs = load_results(CScapeResultSet, "/path/to/data"; scenario_ids=1:50)   # scenarios 1–50 only
rs = load_results(CScapeResultSet, "/path/to/data", "/path/to/exports"; scenario_ids=[1,5,20])
rs = load_results(CScapeResultSet, "/path/to/data", ["file1.jld2", "file2.jld2"])
rs = load_results(CScapeResultSet, "/path/to/saved_rs.jld2")             # reload saved ResultSet
```
"""
function load_results(
    ::Type{CScapeResultSet}, path::String;
    scenario_ids=nothing, show_progress::Bool=true
)::CScapeResultSet
    # If path points to a saved ResultSet JLD2, load it directly
    if isfile(path) && endswith(lowercase(path), ".jld2")
        @info "Loading saved CScapeResultSet from $path"
        return JLD2.load(path, "result_set")
    end
    return load_results(
        CScapeResultSet, path,
        joinpath(path, "adria_exports");
        scenario_ids=scenario_ids, show_progress=show_progress
    )
end

function load_results(
    ::Type{CScapeResultSet}, data_dir::String, result_dir::String;
    scenario_ids=nothing, show_progress::Bool=true
)::CScapeResultSet
    files = _get_result_paths(result_dir)
    if !isnothing(scenario_ids)
        id_set = Set{Int}(scenario_ids)
        files = filter(f -> _parse_indicator_filename(f)[1] in id_set, files)
        isempty(files) && error(
            "No indicator files found for scenario_ids=$scenario_ids in $result_dir"
        )
    end
    return load_results(CScapeResultSet, data_dir, files; show_progress=show_progress)
end

function load_results(
    ::Type{CScapeResultSet}, data_dir::String, indicator_files::Vector{String};
    show_progress::Bool=true
)::CScapeResultSet
    isempty(indicator_files) && error("No indicator files provided.")
    !isdir(data_dir) && error("data_dir not found: $data_dir")

    scenario_spec = _load_scenario_spec(data_dir)
    conn_col = _get_connectivity_col(scenario_spec)
    id_col = first(filter(c -> lowercase(c) == "id", names(scenario_spec)))

    # Load all JLD2 files
    @info "Loading $(length(indicator_files)) indicator file(s)"
    jld2_data = [JLD2.load(f) for f in indicator_files]

    _validate_domain_uniformity(jld2_data, scenario_spec)

    # Extract shared domain metadata from first file
    first_d      = jld2_data[1]
    first_ind    = first_d["cscape_indicators"]
    loc_ids      = first_ind.site_ids
    loc_area     = first_d["loc_area"]
    kappa        = first_d["kappa"]
    meshpoints   = first_d["meshpoints"]
    spatial      = first_d["spatial"]
    years        = first_ind.years
    fts          = first_ind.fts

    # Centroids (from long/lat columns if available)
    loc_centroids = if hasproperty(spatial, :long) && hasproperty(spatial, :lat)
        Matrix{Float64}(hcat(spatial[!, :long], spatial[!, :lat]))
    else
        nothing
    end

    # Connectivity
    first_scenario_id = first_d["metadata"]["scenario_id"]
    spec_row = scenario_spec[scenario_spec[!, Symbol(id_col)] .== first_scenario_id, :]
    conn_rel = isempty(spec_row) ? "" : string(spec_row[1, conn_col])
    connectivity = _load_connectivity(data_dir, scenario_spec, first_scenario_id)

    env_layer_md = _build_env_layer(data_dir, conn_rel, years)

    # Inputs and model spec
    inputs     = _build_inputs_df(jld2_data, scenario_spec)
    model_spec = _build_model_spec(inputs)

    # RCP string (unique RCP values, comma-separated)
    rcp_str = if :RCP in propertynames(inputs)
        join(unique(string.(inputs.RCP)), ", ")
    elseif :rcp in propertynames(inputs)
        join(unique(string.(inputs.rcp)), ", ")
    else
        ""
    end

    # Outcomes
    indicators_list = CscapeIndicators[d["cscape_indicators"] for d in jld2_data]
    outcomes = _build_outcomes(indicators_list, years, loc_ids, fts)

    # Raw output paths (optional; may not exist)
    raw_output_paths = _infer_output_path.(indicator_files)

    name = basename(rstrip(data_dir, ['/', '\\']))

    @info "Loaded CScapeResultSet" scenarios=length(indicator_files) locations=length(loc_ids) timesteps="$(years[1]):$(years[end])"

    return CScapeResultSet(
        name,
        rcp_str,
        loc_ids,
        loc_area,
        kappa,
        loc_centroids,
        env_layer_md,
        connectivity,
        spatial,
        inputs,
        SimConstants(),
        model_spec,
        outcomes,
        meshpoints,
        raw_output_paths
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# load_grouped_results
# ─────────────────────────────────────────────────────────────────────────────

"""
    list_groups(data_dir; result_dir, scenario_ids) -> DataFrame

List the domain groups present in `result_dir` without loading any data.
Shows group index, spatial file, connectivity file, and scenario count for each group.
Use the index with `load_grouped_results(...; group=i)` to load a specific group.

# Example
```julia
list_groups("/path/to/data")
# ┌────────┬─────────────────┬───────────────────┬──────────┐
# │ group  │ spatial_file    │ connectivity_file  │ n_scens  │
# ├────────┼─────────────────┼───────────────────┼──────────┤
# │ 1      │ reef_sites_1... │ reef_sites_1co...  │ 16       │
# │ 2      │ reef_sites_2... │ reef_sites_2co...  │ 16       │
# └────────┴─────────────────┴───────────────────┴──────────┘
```
"""
function list_groups(
    data_dir::String;
    result_dir::String = joinpath(data_dir, "adria_exports"),
    scenario_ids=nothing
)::DataFrame
    indicator_files = _get_result_paths(result_dir)
    if !isnothing(scenario_ids)
        id_set = Set{Int}(scenario_ids)
        indicator_files = filter(f -> _parse_indicator_filename(f)[1] in id_set, indicator_files)
    end
    scenario_spec = _load_scenario_spec(data_dir)
    conn_col      = _get_connectivity_col(scenario_spec)
    spat_col      = _get_spatial_col(scenario_spec)
    id_col        = first(filter(c -> lowercase(c) == "id", names(scenario_spec)))

    ordered_keys = Tuple{String,String}[]
    counts       = Dict{Tuple{String,String}, Int}()
    for fpath in indicator_files
        scenario_id, _ = _parse_indicator_filename(fpath)
        row = scenario_spec[scenario_spec[!, Symbol(id_col)] .== scenario_id, :]
        spat = isempty(row) ? "(unknown)" : string(row[1, spat_col])
        conn = isempty(row) ? "(unknown)" : string(row[1, conn_col])
        key = (spat, conn)
        if !haskey(counts, key)
            push!(ordered_keys, key)
            counts[key] = 0
        end
        counts[key] += 1
    end

    return DataFrame(
        group             = 1:length(ordered_keys),
        spatial_file      = [k[1] for k in ordered_keys],
        connectivity_file = [k[2] for k in ordered_keys],
        n_scenarios       = [counts[k] for k in ordered_keys]
    )
end

"""
    load_grouped_results(data_dir; result_dir, scenario_ids, group, show_progress)

Group all indicator JLD2 files in `result_dir` by `(Spatial_file, Connectivity_file)` and
return one `CScapeResultSet` per unique domain combination.

Use `list_groups(data_dir)` first to see which groups are present and their indices.
Pass `group=i` to load only that group (returns a single `CScapeResultSet`).
Omit `group` to load all groups (returns `Vector{CScapeResultSet}`).

# Examples
```julia
# Inspect groups first
list_groups("/path/to/data")

# Load all groups
result_sets = load_grouped_results("/path/to/data")

# Load only group 2
rs = load_grouped_results("/path/to/data"; group=2)

# Load a subset of scenarios from group 1
rs = load_grouped_results("/path/to/data"; scenario_ids=1:16, group=1)
```
"""
function load_grouped_results(
    data_dir::String;
    result_dir::String = joinpath(data_dir, "adria_exports"),
    scenario_ids=nothing,
    group::Union{Int,Nothing} = nothing,
    show_progress::Bool = true
)
    indicator_files = _get_result_paths(result_dir)
    if !isnothing(scenario_ids)
        id_set = Set{Int}(scenario_ids)
        indicator_files = filter(f -> _parse_indicator_filename(f)[1] in id_set, indicator_files)
        isempty(indicator_files) && error(
            "No indicator files found for scenario_ids=$scenario_ids in $result_dir"
        )
    end
    scenario_spec = _load_scenario_spec(data_dir)
    conn_col      = _get_connectivity_col(scenario_spec)
    spat_col      = _get_spatial_col(scenario_spec)
    id_col        = first(filter(c -> lowercase(c) == "id", names(scenario_spec)))

    # Map each file to its domain key preserving insertion order
    ordered_keys = Tuple{String,String}[]
    groups_map   = Dict{Tuple{String,String}, Vector{String}}()
    for fpath in indicator_files
        scenario_id, _ = _parse_indicator_filename(fpath)
        row = scenario_spec[scenario_spec[!, Symbol(id_col)] .== scenario_id, :]
        spat = isempty(row) ? "" : string(row[1, spat_col])
        conn = isempty(row) ? "" : string(row[1, conn_col])
        key = (spat, conn)
        if !haskey(groups_map, key)
            push!(ordered_keys, key)
            groups_map[key] = String[]
        end
        push!(groups_map[key], fpath)
    end

    n_groups = length(ordered_keys)
    if n_groups > 1
        @info "Found $n_groups domain group(s). Use list_groups() to inspect. Pass group=i to load a specific one."
    end

    if !isnothing(group)
        1 <= group <= n_groups || error("group=$group is out of range 1:$n_groups — call list_groups() to see available groups")
        key = ordered_keys[group]
        @info "Loading group $group (spatial=$(key[1]), connectivity=$(key[2]))"
        return load_results(
            CScapeResultSet, data_dir, sort(groups_map[key]); show_progress=show_progress
        )::CScapeResultSet
    end

    return CScapeResultSet[
        load_results(CScapeResultSet, data_dir, sort(groups_map[k]); show_progress=show_progress)
        for k in ordered_keys
    ]
end

# ─────────────────────────────────────────────────────────────────────────────
# Display
# ─────────────────────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", rs::CScapeResultSet)
    tf = rs.env_layer_md.timeframe
    n_scens = size(rs.inputs, 1)
    n_locs  = length(rs.loc_ids)
    n_raw   = count(isfile, rs.raw_output_paths)
    return println(io, """
        CScapeResultSet
          Name:          $(rs.name)
          RCP(s):        $(rs.RCP)
          Scenarios:     $(n_scens)
          Locations:     $(n_locs)
          Timesteps:     $(tf)
          Outcomes:      $(join(string.(keys(rs.outcomes)), ", "))
          Raw JLD2s available: $(n_raw)/$(length(rs.raw_output_paths))
          Results path:  $(rs.env_layer_md.dpkg_path)
    """)
end

# ─────────────────────────────────────────────────────────────────────────────
# save_results / reload
# ─────────────────────────────────────────────────────────────────────────────

"""
    save_results(rs::CScapeResultSet, filepath::String) -> String

Save a `CScapeResultSet` to a JLD2 file so it can be reloaded instantly without
re-scanning indicator files or re-reading the ScenarioID table.

Returns the absolute path written. Adds a `.jld2` extension if not already present.

# Example
```julia
save_results(rs, "/path/to/my_results.jld2")
rs2 = load_results(CScapeResultSet, "/path/to/my_results.jld2")
```
"""
function save_results(rs::CScapeResultSet, filepath::String)::String
    endswith(lowercase(filepath), ".jld2") || (filepath = filepath * ".jld2")
    mkpath(dirname(abspath(filepath)))
    jldsave(filepath; result_set=rs)
    @info "Saved CScapeResultSet to $filepath" scenarios=size(rs.inputs, 1) locations=length(rs.loc_ids)
    return abspath(filepath)
end

"""
    save_results(result_sets::Vector{CScapeResultSet}, filepath::String) -> Vector{String}

Save each `CScapeResultSet` in the vector to a separate JLD2 file.
The group index is inserted before the `.jld2` extension, e.g.
`"my_results.jld2"` → `"my_results_group1.jld2"`, `"my_results_group2.jld2"`, …

Returns the list of absolute paths written.

# Example
```julia
result_sets = load_grouped_results(fpath)
save_results(result_sets, joinpath(fpath, "results.jld2"))
# → ["…/results_group1.jld2", "…/results_group2.jld2"]
```
"""
function save_results(result_sets::Vector{CScapeResultSet}, filepath::String)::Vector{String}
    endswith(lowercase(filepath), ".jld2") || (filepath = filepath * ".jld2")
    base = filepath[1:end-5]  # strip .jld2
    return [save_results(rs, "$(base)_group$(i).jld2") for (i, rs) in enumerate(result_sets)]
end

# ─────────────────────────────────────────────────────────────────────────────
# reformat_cube (kept for downstream compatibility)
# ─────────────────────────────────────────────────────────────────────────────

"""
    reformat_cube(cscape_cube) -> YAXArray

Rename and reorder C~scape dimension names to ADRIA's expected convention:
`:year` → `:timesteps`, `:ft` → `:species`, `:reef_sites` → `:locations`,
`:draws` → `:scenarios`.
"""
function reformat_cube(cscape_cube::YAXArray)::YAXArray
    name_map = Dict(:year => :timesteps, :ft => :species,
                    :reef_sites => :locations, :draws => :scenarios)
    dim_names = YAXArrays.name.(cscape_cube.axes)
    for (old, new) in name_map
        if old in dim_names
            cscape_cube = YAXArrays.renameaxis!(cscape_cube, old => new)
        end
    end
    if haskey(cscape_cube.properties, "units") &&
            cscape_cube.properties["units"] == "percent"
        cscape_cube = cscape_cube ./ 100
    end
    return cscape_cube
end

# ─────────────────────────────────────────────────────────────────────────────
# ADRIA compatibility helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    pawn(rs::CScapeResultSet, y; S) -> YAXArray

Override ADRIA's generic `ResultSet` dispatch so that only numeric columns from
`rs.inputs` are passed to the matrix-based PAWN implementation. The generic path
calls `Matrix(rs.inputs)` which yields `Matrix{Any}` whenever `inputs` contains
non-Real columns (e.g. `draw::String`, connectivity/spatial path strings).
"""
function pawn(
    rs::CScapeResultSet,
    y::Union{AbstractVector{<:Real}, YAXArrays.Cubes.YAXArray};
    S::Int64=10
)::YAXArray
    numeric_cols = filter(c -> eltype(rs.inputs[!, c]) <: Real, names(rs.inputs))
    return pawn(select(rs.inputs, numeric_cols), y; S=S)
end

"""
    scenario_outcome(rs::CScapeResultSet, key::Symbol) -> YAXArray [T × S]

Area-weighted mean of any pre-computed `[T × L × S]` outcome across locations,
returning a `[timesteps × scenarios]` YAXArray suitable for `ADRIA.viz.scenarios!`,
`ADRIA.sensitivity.pawn`, and `ADRIA.sensitivity.tsa`.

All 2D-per-location outcomes in `rs.outcomes` are supported:
    :relative_cover, :relative_juveniles, :reef_fish_index,
    :reef_condition_index, :reef_tourism_index,
    :reef_biodiversity_condition_index, :relative_shelter_volume, etc.

# Examples
    scenario_outcome(rs, :relative_cover)
    scenario_outcome(rs, :reef_fish_index)
"""
function scenario_outcome(rs::CScapeResultSet, key::Symbol)::YAXArray
    outcome = rs.outcomes[key]                              # YAXArray [T × L × S]
    data_3d = Float64.(collect(outcome))                    # Array{Float64,3} [T × L × S]
    data_2d = ADRIAIndicators.scenario_metric(
        data_3d, rs.loc_area, 2                            # dim 2 = locations
    )                                                       # Array{Float64,2} [T × S]
    ax_t = Dim{:timesteps}(collect(rs.env_layer_md.timeframe))
    ax_s = Dim{:scenarios}(1:size(data_2d, 2))
    return YAXArray((ax_t, ax_s), data_2d)
end

"""
    scenario_groups(rs::CScapeResultSet; by=nothing) -> Dict{Symbol,BitVector}

Group scenarios into named sets for use with ADRIA visualisation and sensitivity
analysis. Pass the result to `ADRIA.viz.scenarios!(g, ax, outcomes, scen_groups)`.

`by` may be any column name (as a Symbol) in `rs.inputs`. If `nothing`, all
scenarios are returned as a single `:counterfactual` group (the key ADRIA uses
for baseline/unmanaged scenarios in its standard color scheme).

When using custom `by` keys (e.g. `:disturbance_intensity`), pass
`opts=Dict(:by_RCP => true)` to `ADRIA.viz.scenarios!` so that the non-standard
key names are accepted:
    ADRIA.viz.scenarios!(g, ax, outcome, scen_groups; opts=Dict(:by_RCP => true))

# Examples
    scenario_groups(rs)                            # single :counterfactual group
    scenario_groups(rs; by=:disturbance_intensity) # one group per unique value
"""
function scenario_groups(rs::CScapeResultSet; by::Union{Symbol,Nothing}=nothing)::Dict{Symbol,BitVector}
    isnothing(by) && return Dict(:counterfactual => trues(nrow(rs.inputs)))
    col = rs.inputs[!, by]
    return Dict(Symbol(string(by, "_", v)) => col .== v for v in unique(col))
end
