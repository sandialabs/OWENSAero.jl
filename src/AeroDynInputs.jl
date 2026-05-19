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

function _aerodyn_required_float(lines, keyword)
    value = tryparse(Float64, _aerodyn_find_keyword_value(lines, keyword))
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
