import FLOWMath

function _aerodyn_uncommented_line(line)
    return strip(first(split(line, '!'; limit = 2)))
end

function _aerodyn_tokens(line)
    clean = _aerodyn_uncommented_line(line)
    return isempty(clean) ? String[] : split(clean)
end

function _aerodyn_keyword_value(line, keyword)
    tokens = _aerodyn_tokens(line)
    length(tokens) >= 2 || return nothing
    return tokens[2] == keyword ? tokens[1] : nothing
end

function _aerodyn_parse_numeric_prefix(line, count)
    tokens = _aerodyn_tokens(line)
    length(tokens) >= count || return nothing
    values = Float64[]
    for token in tokens[1:count]
        value = tryparse(Float64, token)
        value === nothing && return nothing
        push!(values, value)
    end
    return values
end

function _aerodyn_find_keyword_value(lines, keyword)
    for line in lines
        value = _aerodyn_keyword_value(line, keyword)
        value === nothing || return value
    end
    throw(ArgumentError("AeroDyn keyword $(repr(keyword)) was not found"))
end

function _aerodyn_find_optional_keyword_value(lines, keyword)
    for line in lines
        value = _aerodyn_keyword_value(line, keyword)
        value === nothing || return value
    end
    return nothing
end

function _aerodyn_parse_bool(token, keyword)
    token_lower = lowercase(strip(token))
    token_lower == "true" && return true
    token_lower == "false" && return false
    throw(ArgumentError("AeroDyn keyword $(repr(keyword)) must be True or False"))
end

function _aerodyn_exact_integer(value, name)
    rounded = round(Int, value)
    isapprox(value, rounded; atol = 0.0, rtol = 0.0) ||
        throw(ArgumentError("$name must be an integer-valued field"))
    return rounded
end

"""
    readAeroDynBladeFile(filename)

Read an AeroDyn v15 blade-definition file and return a named tuple with the
geometry columns needed by the CCBlade HAWT verification path. Twist and
curve-angle values are returned both in the file's degree units and in radians.

The reader consumes exactly `NumBlNds` blade rows and ignores later notes or
extra rows. It validates strictly increasing span, positive chord, finite
numeric geometry, and integer-valued airfoil IDs before returning data.
"""
function readAeroDynBladeFile(filename)
    lines = readlines(filename)
    num_nodes = parse(Int, _aerodyn_find_keyword_value(lines, "NumBlNds"))

    rows = Vector{Vector{Float64}}()
    num_nodes_line =
        findfirst(line -> _aerodyn_keyword_value(line, "NumBlNds") !== nothing, lines)
    num_nodes_line === nothing &&
        throw(ArgumentError("AeroDyn blade file does not define NumBlNds"))

    for line in lines[(num_nodes_line+1):end]
        values = _aerodyn_parse_numeric_prefix(line, 7)
        values === nothing && continue
        push!(rows, values)
        length(rows) == num_nodes && break
    end

    length(rows) == num_nodes || throw(
        ArgumentError(
            "AeroDyn blade file expected $num_nodes blade rows but found $(length(rows))",
        ),
    )

    span = [row[1] for row in rows]
    curve_ac = [row[2] for row in rows]
    sweep_ac = [row[3] for row in rows]
    curve_angle_deg = [row[4] for row in rows]
    twist_deg = [row[5] for row in rows]
    chord = [row[6] for row in rows]
    airfoil_ids = [_aerodyn_exact_integer(row[7], "BlAFID") for row in rows]

    all(isfinite, span) &&
    all(isfinite, curve_ac) &&
    all(isfinite, sweep_ac) &&
    all(isfinite, curve_angle_deg) &&
    all(isfinite, twist_deg) &&
    all(isfinite, chord) ||
        throw(ArgumentError("AeroDyn blade geometry must contain only finite values"))
    all(diff(span) .> zero(first(span))) ||
        throw(ArgumentError("AeroDyn blade span stations must be strictly increasing"))
    all(>(0.0), chord) ||
        throw(ArgumentError("AeroDyn blade chord values must be positive"))
    all(>(0), airfoil_ids) ||
        throw(ArgumentError("AeroDyn blade airfoil IDs must be positive"))

    return (
        source = String(filename),
        num_nodes = num_nodes,
        span = span,
        curve_ac = curve_ac,
        sweep_ac = sweep_ac,
        curve_angle_deg = curve_angle_deg,
        curve_angle_rad = deg2rad.(curve_angle_deg),
        twist_deg = twist_deg,
        twist_rad = deg2rad.(twist_deg),
        chord = chord,
        airfoil_ids = airfoil_ids,
    )
end

function _aerodyn_airfoil_metadata(lines, start_index)
    metadata = Dict{String,String}()
    line_index = start_index
    while line_index <= length(lines)
        tokens = _aerodyn_tokens(lines[line_index])
        if length(tokens) >= 2
            metadata[tokens[2]] = tokens[1]
            tokens[2] == "NumAlf" && return metadata, parse(Int, tokens[1]), line_index
        end
        line_index += 1
    end
    throw(ArgumentError("AeroDyn airfoil table is missing NumAlf"))
end

function _aerodyn_required_metadata_float(metadata, keyword)
    haskey(metadata, keyword) ||
        throw(ArgumentError("AeroDyn airfoil table is missing $(repr(keyword))"))
    return parse(Float64, metadata[keyword])
end

function _aerodyn_optional_metadata_float(metadata, keyword)
    return haskey(metadata, keyword) ? parse(Float64, metadata[keyword]) : nothing
end

function _aerodyn_required_metadata_bool(metadata, keyword)
    haskey(metadata, keyword) ||
        throw(ArgumentError("AeroDyn airfoil table is missing $(repr(keyword))"))
    return _aerodyn_parse_bool(metadata[keyword], keyword)
end

function _aerodyn_unquote(token)
    stripped = strip(String(token))
    if length(stripped) >= 2 && first(stripped) == '"' && last(stripped) == '"'
        return stripped[2:(end-1)]
    end
    return stripped
end

function _aerodyn_required_int(lines, keyword)
    value = _aerodyn_find_keyword_value(lines, keyword)
    parsed = tryparse(Int, value)
    parsed !== nothing && return parsed
    numeric = tryparse(Float64, value)
    numeric !== nothing && return _aerodyn_exact_integer(numeric, keyword)
    throw(ArgumentError("AeroDyn keyword $(repr(keyword)) must be an integer"))
end

function _aerodyn_required_bool(lines, keyword)
    return _aerodyn_parse_bool(_aerodyn_find_keyword_value(lines, keyword), keyword)
end

function _aerodyn_optional_bool(lines, keyword)
    raw_value = _aerodyn_find_optional_keyword_value(lines, keyword)
    return raw_value === nothing ? nothing : _aerodyn_parse_bool(raw_value, keyword)
end

function _aerodyn_required_float(lines, keyword)
    value = tryparse(Float64, _aerodyn_find_keyword_value(lines, keyword))
    value === nothing &&
        throw(ArgumentError("AeroDyn keyword $(repr(keyword)) must be numeric"))
    return value
end

function _aerodyn_optional_float(lines, keyword)
    raw_value = _aerodyn_find_optional_keyword_value(lines, keyword)
    raw_value === nothing && return nothing
    value = tryparse(Float64, raw_value)
    value === nothing &&
        throw(ArgumentError("AeroDyn keyword $(repr(keyword)) must be numeric"))
    return value
end

function _aerodyn_collect_airfoil_files(lines, num_airfoil_files)
    af_line = findfirst(line -> _aerodyn_keyword_value(line, "AFNames") !== nothing, lines)
    af_line === nothing && throw(ArgumentError("AeroDyn primary file is missing AFNames"))

    files = String[]
    for line in lines[af_line:end]
        tokens = _aerodyn_tokens(line)
        isempty(tokens) && continue
        startswith(strip(tokens[1]), "\"") || break
        push!(files, _aerodyn_unquote(tokens[1]))
        length(files) == num_airfoil_files && break
    end

    length(files) == num_airfoil_files || throw(
        ArgumentError(
            "AeroDyn primary file expected $num_airfoil_files airfoil files but found $(length(files))",
        ),
    )
    return files
end

function _aerodyn_collect_blade_files(lines)
    files = String[]
    for line in lines
        tokens = _aerodyn_tokens(line)
        length(tokens) >= 2 || continue
        startswith(tokens[2], "ADBlFile") || continue
        push!(files, _aerodyn_unquote(tokens[1]))
    end
    isempty(files) &&
        throw(ArgumentError("AeroDyn primary file is missing ADBlFile entries"))
    return files
end

function _aerodyn_required_input_column(lines, keyword)
    column = _aerodyn_required_int(lines, keyword)
    column >= 0 ||
        throw(ArgumentError("AeroDyn keyword $(repr(keyword)) must be nonnegative"))
    return column
end

"""
    readAeroDynPrimaryFile(filename)

Read the AeroDyn primary input fields needed to reproduce a steady BEM HAWT
case with OWENSAero's CCBlade adapter. The returned named tuple includes wake
and airfoil model switches, BEM options, input-column mapping, airfoil-file
paths, blade-file paths, and the `UseBlCm` moment flag.

The reader intentionally does not execute AeroDyn or resolve secondary-file
paths. It preserves file paths exactly as written so callers can choose whether
to resolve them relative to the AeroDyn primary file, the driver file, or the
current working directory.
"""
function readAeroDynPrimaryFile(filename)
    lines = readlines(filename)
    num_airfoil_files = _aerodyn_required_int(lines, "NumAFfiles")
    num_airfoil_files > 0 ||
        throw(ArgumentError("AeroDyn primary file must list at least one airfoil file"))

    return (
        source = String(filename),
        wake_model = _aerodyn_required_int(lines, "WakeMod"),
        airfoil_aero_model = _aerodyn_required_int(lines, "AFAeroMod"),
        skew_model = _aerodyn_required_int(lines, "SkewMod"),
        tip_loss = _aerodyn_required_bool(lines, "TipLoss"),
        hub_loss = _aerodyn_required_bool(lines, "HubLoss"),
        tangential_induction = _aerodyn_required_bool(lines, "TanInd"),
        axial_induction_drag = _aerodyn_required_bool(lines, "AIDrag"),
        tangential_induction_drag = _aerodyn_required_bool(lines, "TIDrag"),
        induction_tolerance = _aerodyn_unquote(
            _aerodyn_find_keyword_value(lines, "IndToler"),
        ),
        max_iterations = _aerodyn_required_int(lines, "MaxIter"),
        dbemt_model = _aerodyn_required_int(lines, "DBEMT_Mod"),
        tau1_const = _aerodyn_required_float(lines, "tau1_const"),
        unsteady_aero_model = _aerodyn_required_int(lines, "UAMod"),
        f_lookup = _aerodyn_required_bool(lines, "FLookup"),
        airfoil_table_model = _aerodyn_required_int(lines, "AFTabMod"),
        input_columns = (
            alpha = _aerodyn_required_input_column(lines, "InCol_Alfa"),
            cl = _aerodyn_required_input_column(lines, "InCol_Cl"),
            cd = _aerodyn_required_input_column(lines, "InCol_Cd"),
            cm = _aerodyn_required_input_column(lines, "InCol_Cm"),
            cpmin = _aerodyn_required_input_column(lines, "InCol_Cpmin"),
        ),
        num_airfoil_files = num_airfoil_files,
        airfoil_files = _aerodyn_collect_airfoil_files(lines, num_airfoil_files),
        use_blade_cm = _aerodyn_required_bool(lines, "UseBlCm"),
        blade_files = _aerodyn_collect_blade_files(lines),
    )
end

"""
    readAeroDynAirfoilInfo(filename; table_index=1)

Read one table from an AeroDyn AirfoilInfo v1 file. The returned named tuple
contains table metadata, unsteady-aero coefficients commonly used for
validation, and static aerodynamic coefficient arrays. Angles are returned in
both degrees and radians; `Cm` is included and defaults to zeros only when the
file table omits a moment column. Unsteady-aero metadata fields are returned as
`nothing` when `InclUAdata = False` and the file omits those fields.
"""
function readAeroDynAirfoilInfo(filename; table_index = 1)
    table_index isa Integer && table_index > 0 ||
        throw(ArgumentError("table_index must be a positive integer"))

    lines = readlines(filename)
    num_tabs = parse(Int, _aerodyn_find_keyword_value(lines, "NumTabs"))
    table_index <= num_tabs || throw(
        ArgumentError(
            "AeroDyn airfoil file has $num_tabs table(s), cannot read table $table_index",
        ),
    )

    re_lines = findall(line -> _aerodyn_keyword_value(line, "Re") !== nothing, lines)
    length(re_lines) >= num_tabs || throw(
        ArgumentError(
            "AeroDyn airfoil file declares $num_tabs table(s) but only $(length(re_lines)) Re line(s) were found",
        ),
    )

    table_start = re_lines[table_index]
    metadata, num_alpha, num_alpha_line = _aerodyn_airfoil_metadata(lines, table_start)

    alpha_deg = Float64[]
    cl = Float64[]
    cd = Float64[]
    cm = Float64[]

    for line in lines[(num_alpha_line+1):end]
        values = _aerodyn_parse_numeric_prefix(line, 3)
        values === nothing && continue
        all(isfinite, values) ||
            throw(ArgumentError("AeroDyn airfoil coefficient rows must be finite"))
        push!(alpha_deg, values[1])
        push!(cl, values[2])
        push!(cd, values[3])

        four_values = _aerodyn_parse_numeric_prefix(line, 4)
        push!(cm, four_values === nothing ? 0.0 : four_values[4])
        length(alpha_deg) == num_alpha && break
    end

    length(alpha_deg) == num_alpha || throw(
        ArgumentError(
            "AeroDyn airfoil table expected $num_alpha rows but found $(length(alpha_deg))",
        ),
    )
    all(isfinite, alpha_deg) &&
    all(isfinite, cl) &&
    all(isfinite, cd) &&
    all(isfinite, cm) ||
        throw(ArgumentError("AeroDyn airfoil data must contain only finite values"))
    all(diff(alpha_deg) .> zero(first(alpha_deg))) ||
        throw(ArgumentError("AeroDyn airfoil alpha values must be strictly increasing"))

    return (
        source = String(filename),
        table_index = Int(table_index),
        num_tables = num_tabs,
        reynolds_million = _aerodyn_required_metadata_float(metadata, "Re"),
        user_property = _aerodyn_required_metadata_float(metadata, "UserProp"),
        includes_unsteady_aero = _aerodyn_required_metadata_bool(metadata, "InclUAdata"),
        alpha0_deg = _aerodyn_optional_metadata_float(metadata, "alpha0"),
        alpha1_deg = _aerodyn_optional_metadata_float(metadata, "alpha1"),
        alpha2_deg = _aerodyn_optional_metadata_float(metadata, "alpha2"),
        cn_alpha = _aerodyn_optional_metadata_float(metadata, "C_nalpha"),
        cd0 = _aerodyn_optional_metadata_float(metadata, "Cd0"),
        cm0 = _aerodyn_optional_metadata_float(metadata, "Cm0"),
        num_alpha = num_alpha,
        alpha_deg = alpha_deg,
        alpha_rad = deg2rad.(alpha_deg),
        cl = cl,
        cd = cd,
        cm = cm,
    )
end

"""
    aeroDynAirfoilFunction(table)

Create an OWENSAero/CCBlade-compatible airfoil lookup from a table returned by
`readAeroDynAirfoilInfo`. The returned function accepts angle of attack in
radians and returns `(cl, cd)` by default or `(cl, cd, cm)` when called with
`return_cm=true`.

Queries outside the tabulated alpha range throw `ArgumentError`; this keeps the
validation path from silently extrapolating beyond the AeroDyn polar.
"""
function aeroDynAirfoilFunction(table)
    alpha = table.alpha_rad
    cl = table.cl
    cd = table.cd
    cm = table.cm
    length(alpha) == length(cl) == length(cd) == length(cm) ||
        throw(ArgumentError("AeroDyn airfoil table arrays must have matching lengths"))
    all(diff(alpha) .> zero(first(alpha))) ||
        throw(ArgumentError("AeroDyn airfoil alpha values must be strictly increasing"))

    cl_spline = FLOWMath.Akima(alpha, cl)
    cd_spline = FLOWMath.Akima(alpha, cd)
    cm_spline = FLOWMath.Akima(alpha, cm)
    alpha_min = first(alpha)
    alpha_max = last(alpha)

    function af(alpha_query, reynolds_number, mach; return_cm = false)
        alpha_query isa Real && isfinite(alpha_query) ||
            throw(ArgumentError("alpha_query must be a finite real value in radians"))
        alpha_min <= alpha_query <= alpha_max || throw(
            ArgumentError(
                "alpha_query=$(alpha_query) is outside the AeroDyn table range [$alpha_min, $alpha_max]",
            ),
        )
        cl_value = cl_spline(alpha_query)
        cd_value = cd_spline(alpha_query)
        cm_value = cm_spline(alpha_query)
        return return_cm ? (cl_value, cd_value, cm_value) : (cl_value, cd_value)
    end

    return af
end

function _aerodyn_uniform_driver_pitch(pitch_deg)
    isempty(pitch_deg) && throw(ArgumentError("AeroDyn driver must define blade pitch"))
    reference = first(pitch_deg)
    all(pitch -> isapprox(pitch, reference; atol = 1e-12, rtol = 0.0), pitch_deg) || throw(
        ArgumentError(
            "AeroDyn driver blade pitches must be identical for a rigid CCBlade solve",
        ),
    )
    return reference
end

function _read_aerodyn_steady_inflow(inflow_file)
    lines = readlines(inflow_file)
    wind_type = _aerodyn_required_int(lines, "WindType")
    wind_type == 1 ||
        throw(ArgumentError("Only InflowWind WindType=1 steady inflow is supported"))
    return (
        source = inflow_file,
        wind_type = wind_type,
        horizontal_wind_speed = _aerodyn_required_float(lines, "HWindSpeed"),
        reference_height = _aerodyn_required_float(lines, "RefHt"),
        power_law_exponent = _aerodyn_required_float(lines, "PLExp"),
    )
end

"""
    readAeroDynDriverFile(filename; base_directory=dirname(filename),
                          turbine_index=1)

Read the operating-point fields needed to replay a steady HAWT AeroDyn-driver
case through OWENSAero's CCBlade adapter. The reader supports the common
`AnalysisType=1` constant-rotor-speed driver layout with either direct steady
wind (`CompInflow=0`) or an InflowWind `WindType=1` steady-wind file
(`CompInflow=1`).

The returned named tuple includes resolved AeroDyn/InflowWind filenames,
environment properties, rotor speed in rpm and rad/s, blade pitch in degrees
and radians, Basic HAWT geometry when present, and driver hub-radius metadata.
For Basic HAWT drivers, `HubRad` is returned as the blade-root offset used by
CCBlade; for advanced drivers, the per-blade `BldHubRad_bl` values are retained.
The driver `VTKHubRad` is retained for traceability but is not treated as a
physical BEM hub radius.
"""
function readAeroDynDriverFile(
    filename;
    base_directory = dirname(filename),
    turbine_index = 1,
)
    turbine_index isa Integer && turbine_index >= 1 ||
        throw(ArgumentError("turbine_index must be a positive integer"))
    lines = readlines(filename)

    analysis_type = _aerodyn_required_int(lines, "AnalysisType")
    analysis_type == 1 ||
        throw(ArgumentError("Only AeroDyn driver AnalysisType=1 is supported"))
    num_turbines = _aerodyn_required_int(lines, "NumTurbines")
    turbine_index <= num_turbines ||
        throw(ArgumentError("turbine_index must select one of the driver turbines"))

    comp_inflow = _aerodyn_required_int(lines, "CompInflow")
    aero_file = _aerodyn_resolve_path(
        base_directory,
        _aerodyn_unquote(_aerodyn_find_keyword_value(lines, "AeroFile")),
    )
    inflow_file =
        comp_inflow == 1 ?
        _aerodyn_resolve_path(
            base_directory,
            _aerodyn_unquote(_aerodyn_find_keyword_value(lines, "InflowFile")),
        ) : nothing
    inflow = if comp_inflow == 0
        (
            source = nothing,
            wind_type = 0,
            horizontal_wind_speed = _aerodyn_required_float(lines, "HWindSpeed"),
            reference_height = _aerodyn_required_float(lines, "RefHt"),
            power_law_exponent = _aerodyn_required_float(lines, "PLExp"),
        )
    elseif comp_inflow == 1
        _read_aerodyn_steady_inflow(inflow_file)
    else
        throw(ArgumentError("Only AeroDyn driver CompInflow=0 or 1 is supported"))
    end

    num_blades = _aerodyn_required_int(lines, "NumBlades($turbine_index)")
    num_blades > 0 || throw(ArgumentError("AeroDyn driver NumBlades must be positive"))

    basic_hawt_format =
        something(_aerodyn_optional_bool(lines, "BasicHAWTFormat($turbine_index)"), false)
    hub_radius = nothing
    hub_height = nothing
    overhang = nothing
    shaft_tilt_deg = nothing
    precone_deg = 0.0
    tower_to_shaft = nothing
    if basic_hawt_format
        hub_radius = _aerodyn_required_float(lines, "HubRad($turbine_index)")
        hub_radius > 0.0 || throw(ArgumentError("AeroDyn driver HubRad must be positive"))
        hub_height = _aerodyn_required_float(lines, "HubHt($turbine_index)")
        overhang = _aerodyn_required_float(lines, "Overhang($turbine_index)")
        shaft_tilt_deg = _aerodyn_required_float(lines, "ShftTilt($turbine_index)")
        precone_deg = _aerodyn_required_float(lines, "Precone($turbine_index)")
        tower_to_shaft = _aerodyn_required_float(lines, "Twr2Shft($turbine_index)")
    end

    pitch_deg = if basic_hawt_format
        fill(_aerodyn_required_float(lines, "BldPitch($turbine_index)"), num_blades)
    else
        [
            _aerodyn_required_float(lines, "BldPitch($(turbine_index)_$blade_index)")
            for blade_index = 1:num_blades
        ]
    end
    uniform_pitch_deg = _aerodyn_uniform_driver_pitch(pitch_deg)
    blade_hub_radius = if basic_hawt_format
        fill(hub_radius, num_blades)
    else
        [
            _aerodyn_required_float(lines, "BldHubRad_bl($(turbine_index)_$blade_index)") for blade_index = 1:num_blades
        ]
    end
    all(>=(0.0), blade_hub_radius) ||
        throw(ArgumentError("AeroDyn driver blade hub radii must be nonnegative"))

    rot_speed_rpm = _aerodyn_required_float(lines, "RotSpeed($turbine_index)")
    fluid_density = _aerodyn_required_float(lines, "FldDens")
    kinematic_viscosity = _aerodyn_required_float(lines, "KinVisc")
    speed_of_sound = _aerodyn_required_float(lines, "SpdSound")
    nac_yaw_deg = _aerodyn_optional_float(lines, "NacYaw($turbine_index)")

    return (
        source = filename,
        base_directory = base_directory,
        aero_file = aero_file,
        inflow_file = inflow_file,
        analysis_type = analysis_type,
        comp_inflow = comp_inflow,
        num_turbines = num_turbines,
        turbine_index = turbine_index,
        basic_hawt_format = basic_hawt_format,
        num_blades = num_blades,
        fluid_density = fluid_density,
        kinematic_viscosity = kinematic_viscosity,
        dynamic_viscosity = fluid_density * kinematic_viscosity,
        speed_of_sound = speed_of_sound,
        inflow = inflow,
        inflow_speed = inflow.horizontal_wind_speed,
        rot_speed_rpm = rot_speed_rpm,
        rotor_speed_rad_per_s = rot_speed_rpm * 2pi / 60,
        blade_pitch_deg = pitch_deg,
        pitch_deg = uniform_pitch_deg,
        pitch_rad = deg2rad(uniform_pitch_deg),
        blade_hub_radius = blade_hub_radius,
        hub_radius = basic_hawt_format ? hub_radius : maximum(blade_hub_radius),
        hub_height = hub_height,
        overhang = overhang,
        shaft_tilt_deg = shaft_tilt_deg,
        shaft_tilt_rad = shaft_tilt_deg === nothing ? nothing : deg2rad(shaft_tilt_deg),
        precone_deg = precone_deg,
        precone_rad = deg2rad(precone_deg),
        tower_to_shaft = tower_to_shaft,
        nac_yaw_deg = nac_yaw_deg,
        nac_yaw_rad = nac_yaw_deg === nothing ? nothing : deg2rad(nac_yaw_deg),
        vtk_hub_radius = _aerodyn_optional_float(lines, "VTKHubRad"),
    )
end

function _aerodyn_resolve_path(base_directory, path)
    isabspath(path) && return normpath(path)
    return normpath(joinpath(base_directory, path))
end

function _aerodyn_hawt_station_indices(span, root_station_policy)
    if root_station_policy == :drop_zero_span
        indices = findall(x -> x > 0.0, span)
        isempty(indices) && throw(
            ArgumentError("root_station_policy=:drop_zero_span removed every station"),
        )
        return indices
    elseif root_station_policy == :strict_positive
        any(x -> x <= 0.0, span) && throw(
            ArgumentError(
                "AeroDyn blade contains nonpositive span stations; use root_station_policy=:drop_zero_span or provide a positive-span blade file",
            ),
        )
        return collect(eachindex(span))
    end
    throw(ArgumentError("root_station_policy must be :drop_zero_span or :strict_positive"))
end

function _ccblade_tip_correction_from_aerodyn(primary)
    if primary.tip_loss && primary.hub_loss
        return CCBlade.PrandtlTipHub()
    elseif primary.tip_loss && !primary.hub_loss
        return CCBlade.PrandtlTip()
    elseif !primary.tip_loss && !primary.hub_loss
        return nothing
    end
    throw(
        ArgumentError(
            "AeroDyn hub-loss without tip-loss is not representable by the current CCBlade adapter; pass an explicit tip_correction override",
        ),
    )
end

function _aerodyn_hawt_comparison_notes(primary, hub_radius)
    notes = String[]
    if primary.axial_induction_drag || primary.tangential_induction_drag
        push!(
            notes,
            "AeroDyn drag-in-induction flags are enabled; confirm the selected CCBlade correction settings before comparing coefficients.",
        )
    else
        push!(
            notes,
            "AeroDyn drag-in-induction flags are disabled, while CCBlade's BEM formulation includes drag terms in the induction solve.",
        )
    end
    if primary.hub_loss && hub_radius == 0.0
        push!(
            notes,
            "AeroDyn HubLoss is enabled but hub_radius=0.0 was supplied to CCBlade, so the root-loss effect is not equivalent to a finite hub radius.",
        )
    end
    return notes
end

"""
    aeroDynHAWTCCBladeInputs(primary_file; blade_index=1,
                             base_directory=dirname(primary_file),
                             root_station_policy=:drop_zero_span,
                             tip_correction=:from_aerodyn,
                             hub_radius=0.0)

Build the geometry, airfoil callables, and metadata needed to run OWENSAero's
CCBlade HAWT adapter from AeroDyn primary, blade, and AirfoilInfo files.
`root_station_policy=:drop_zero_span` removes the common AeroDyn root station at
zero span because CCBlade sections require positive radial positions.
`hub_radius` is added to AeroDyn `BlSpn` values because AeroDyn blade files store
span from the blade root, while CCBlade sections use rotor radii from the
center of rotation.

The returned named tuple does not run the BEM solve. It keeps the parsed
AeroDyn primary options, blade data, polar tables, resolved files, station
indices, selected CCBlade tip correction, and comparison notes visible so
validation scripts can document any convention mismatch.
"""
function aeroDynHAWTCCBladeInputs(
    primary_file;
    blade_index = 1,
    base_directory = dirname(primary_file),
    root_station_policy = :drop_zero_span,
    tip_correction = :from_aerodyn,
    hub_radius = 0.0,
)
    _validate_nonnegative_real_input(hub_radius, "hub_radius")
    primary = readAeroDynPrimaryFile(primary_file)
    primary.wake_model == 1 ||
        throw(ArgumentError("Only AeroDyn WakeMod=1 BEM inputs are supported"))
    primary.airfoil_aero_model == 1 ||
        throw(ArgumentError("Only AeroDyn AFAeroMod=1 steady airfoil inputs are supported"))
    primary.airfoil_table_model == 1 ||
        throw(ArgumentError("Only AeroDyn AFTabMod=1 first-table lookup is supported"))
    primary.tangential_induction || throw(
        ArgumentError("AeroDyn TanInd=false cannot be represented by this CCBlade path"),
    )
    primary.input_columns == (alpha = 1, cl = 2, cd = 3, cm = 4, cpmin = 0) || throw(
        ArgumentError(
            "Only the standard AeroDyn Alpha/Cl/Cd/Cm input-column layout is supported",
        ),
    )
    blade_index isa Integer && 1 <= blade_index <= length(primary.blade_files) ||
        throw(ArgumentError("blade_index must select one of primary.blade_files"))

    airfoil_files =
        [_aerodyn_resolve_path(base_directory, file) for file in primary.airfoil_files]
    airfoil_tables = readAeroDynAirfoilInfo.(airfoil_files)
    airfoils = aeroDynAirfoilFunction.(airfoil_tables)

    blade_file = _aerodyn_resolve_path(base_directory, primary.blade_files[blade_index])
    blade = readAeroDynBladeFile(blade_file)
    maximum(blade.airfoil_ids) <= length(airfoils) || throw(
        ArgumentError(
            "AeroDyn blade references airfoil ID $(maximum(blade.airfoil_ids)) but only $(length(airfoils)) airfoil files were parsed",
        ),
    )

    station_indices = _aerodyn_hawt_station_indices(blade.span, root_station_policy)
    station_airfoils = [airfoils[blade.airfoil_ids[i]] for i in station_indices]
    selected_tip_correction =
        tip_correction == :from_aerodyn ? _ccblade_tip_correction_from_aerodyn(primary) :
        tip_correction
    _validate_hawt_tip_correction(selected_tip_correction)

    return (
        primary = primary,
        blade = blade,
        airfoil_tables = airfoil_tables,
        airfoil_files = airfoil_files,
        blade_file = blade_file,
        station_indices = station_indices,
        radial_positions = blade.span[station_indices] .+ hub_radius,
        blade_span = blade.span[station_indices],
        chord = blade.chord[station_indices],
        twist = blade.twist_rad[station_indices],
        airfoils = station_airfoils,
        num_blades = length(primary.blade_files),
        hub_radius = hub_radius,
        tip_radius = maximum(blade.span) + hub_radius,
        tip_correction = selected_tip_correction,
        root_station_policy = root_station_policy,
        comparison_notes = _aerodyn_hawt_comparison_notes(primary, hub_radius),
    )
end

"""
    ccbladeHAWTSolveFromAeroDyn(primary_file, rotor_speed, inflow_speed, rho; kwargs...)

Read an AeroDyn steady-BEM HAWT input set and run `ccbladeHAWTSolve` with the
parsed blade geometry and polar files. Keyword arguments accepted by
`aeroDynHAWTCCBladeInputs` control file resolution and root-station handling;
remaining operating keywords are forwarded to `ccbladeHAWTSolve`.

This helper is an input-normalization bridge for validation. It does not make
CCBlade numerically identical to AeroDyn; use the returned `comparison_notes`
when deciding which OpenFAST channels can be compared directly.
"""
function ccbladeHAWTSolveFromAeroDyn(
    primary_file,
    rotor_speed,
    inflow_speed,
    rho;
    blade_index = 1,
    base_directory = dirname(primary_file),
    root_station_policy = :drop_zero_span,
    hub_radius = 0.0,
    tip_radius = nothing,
    pitch = 0.0,
    precone = 0.0,
    mu = 1.7894e-5,
    asound = 340.0,
    npts = 10,
    tip_correction = :from_aerodyn,
)
    inputs = aeroDynHAWTCCBladeInputs(
        primary_file;
        blade_index,
        base_directory,
        root_station_policy,
        tip_correction,
        hub_radius,
    )
    solve_tip_radius = tip_radius === nothing ? inputs.tip_radius : tip_radius
    result = ccbladeHAWTSolve(
        inputs.radial_positions,
        inputs.chord,
        inputs.twist,
        inputs.airfoils,
        rotor_speed,
        inflow_speed,
        rho;
        num_blades = inputs.num_blades,
        hub_radius,
        tip_radius = solve_tip_radius,
        pitch,
        precone,
        mu,
        asound,
        npts,
        tip_correction = inputs.tip_correction,
    )

    return (
        result...,
        aerodyn_primary = inputs.primary,
        aerodyn_blade = inputs.blade,
        aerodyn_airfoil_tables = inputs.airfoil_tables,
        aerodyn_airfoil_files = inputs.airfoil_files,
        aerodyn_blade_file = inputs.blade_file,
        aerodyn_station_indices = inputs.station_indices,
        root_station_policy = inputs.root_station_policy,
        comparison_notes = inputs.comparison_notes,
    )
end

"""
    ccbladeHAWTSolveFromAeroDynDriver(driver_file; kwargs...)

Read a steady AeroDyn-driver operating point with `readAeroDynDriverFile`, then
run `ccbladeHAWTSolveFromAeroDyn` using the driver's AeroDyn primary file,
steady inflow, density, viscosity, sound speed, rotor speed, and uniform blade
pitch. Geometry-resolution keywords are forwarded to
`ccbladeHAWTSolveFromAeroDyn`; by default `hub_radius` is the largest
`BldHubRad_bl` value in the driver file.

The returned named tuple includes the CCBlade solve outputs plus
`aerodyn_driver` metadata. This helper is meant to make validation scripts
explicit about operating conditions before comparing against OpenFAST channels.
"""
function ccbladeHAWTSolveFromAeroDynDriver(
    driver_file;
    base_directory = dirname(driver_file),
    turbine_index = 1,
    blade_index = 1,
    root_station_policy = :drop_zero_span,
    hub_radius = nothing,
    tip_radius = nothing,
    precone = nothing,
    npts = 10,
    tip_correction = :from_aerodyn,
)
    driver = readAeroDynDriverFile(driver_file; base_directory, turbine_index)
    solve_hub_radius =
        hub_radius === nothing ? maximum(driver.blade_hub_radius) : hub_radius
    solve_precone = precone === nothing ? driver.precone_rad : precone
    result = ccbladeHAWTSolveFromAeroDyn(
        driver.aero_file,
        driver.rotor_speed_rad_per_s,
        driver.inflow_speed,
        driver.fluid_density;
        blade_index,
        base_directory = dirname(driver.aero_file),
        root_station_policy,
        hub_radius = solve_hub_radius,
        tip_radius,
        pitch = driver.pitch_rad,
        precone = solve_precone,
        mu = driver.dynamic_viscosity,
        asound = driver.speed_of_sound,
        npts,
        tip_correction,
    )

    return (
        result...,
        aerodyn_driver = driver,
        aerodyn_driver_file = driver.source,
        driver_hub_radius_used = solve_hub_radius,
    )
end
