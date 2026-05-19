function _validated_hawt_station_vector(values, name; positive = false)
    values isa AbstractVector || throw(ArgumentError("$name must be a vector"))
    isempty(values) && throw(ArgumentError("$name must not be empty"))
    all(x -> x isa Real && isfinite(x), values) ||
        throw(ArgumentError("$name must contain only finite real values"))
    positive &&
        any(x -> x <= zero(x), values) &&
        throw(ArgumentError("$name must contain only finite positive real values"))
    return collect(values)
end

function _validate_hawt_station_inputs(radial_positions, chord, twist)
    r = _validated_hawt_station_vector(
        radial_positions,
        "radial_positions";
        positive = true,
    )
    c = _validated_hawt_station_vector(chord, "chord"; positive = true)
    theta = _validated_hawt_station_vector(twist, "twist")
    length(r) == length(c) == length(theta) ||
        throw(ArgumentError("radial_positions, chord, and twist must have the same length"))
    all(diff(r) .> zero(first(r))) ||
        throw(ArgumentError("radial_positions must be strictly increasing"))
    return r, c, theta
end

function _hawt_airfoil_vector(airfoils, nstations)
    if airfoils isa AbstractVector
        isempty(airfoils) && throw(ArgumentError("airfoils must not be empty"))
        length(airfoils) == nstations || throw(
            ArgumentError("airfoils must be a single callable or one callable per station"),
        )
        return collect(airfoils)
    end
    return fill(airfoils, nstations)
end

function _ccblade_airfoil_adapter(af)
    return function (alpha, reynolds_number, mach)
        cl, cd, _ = _airfoil_coefficients(af, alpha, reynolds_number, mach)
        return cl, cd
    end
end

function _validate_hawt_operating_values(
    inflow_speed,
    rotor_speed,
    rho,
    pitch,
    mu,
    asound,
    precone,
)
    _validate_positive_real_value(inflow_speed, "inflow_speed")
    _validate_finite_real_value(rotor_speed, "rotor_speed")
    rotor_speed == zero(rotor_speed) && throw(ArgumentError("rotor_speed must be nonzero"))
    _validate_positive_real_value(rho, "rho")
    _validate_finite_real_value(pitch, "pitch")
    _validate_positive_real_value(mu, "mu")
    _validate_positive_real_value(asound, "asound")
    _validate_finite_real_value(precone, "precone")
    return nothing
end

function _validate_hawt_rotor_geometry(radial_positions, hub_radius, tip_radius)
    _validate_nonnegative_real_input(hub_radius, "hub_radius")
    _validate_positive_real_value(tip_radius, "tip_radius")
    hub_radius < tip_radius ||
        throw(ArgumentError("hub_radius must be smaller than tip_radius"))
    hub_radius <= first(radial_positions) || throw(
        ArgumentError("hub_radius must be less than or equal to the first radial station"),
    )
    tip_radius >= last(radial_positions) || throw(
        ArgumentError(
            "tip_radius must be greater than or equal to the last radial station",
        ),
    )
    return nothing
end

function _validate_hawt_tip_correction(tip_correction)
    tip_correction === nothing && return nothing
    tip_correction isa CCBlade.TipCorrection ||
        throw(ArgumentError("tip_correction must be nothing or a CCBlade.TipCorrection"))
    return nothing
end

"""
    ccbladeHAWTSections(radial_positions, chord, twist, airfoils)

Construct `CCBlade.Section` objects for a horizontal-axis rotor. `radial_positions`
are dimensional blade-station radii, `chord` values are dimensional chord
lengths, and `twist` uses CCBlade's section-angle convention in radians.
`airfoils` may be a single OWENSAero/CCBlade-style callable or one callable per
station. If a callable returns `(cl, cd, cm)`, only `cl` and `cd` are passed to
CCBlade because CCBlade's BEM solve does not consume moment coefficient.
"""
function ccbladeHAWTSections(radial_positions, chord, twist, airfoils)
    r, c, theta = _validate_hawt_station_inputs(radial_positions, chord, twist)
    station_airfoils = _hawt_airfoil_vector(airfoils, length(r))
    ccblade_airfoils = _ccblade_airfoil_adapter.(station_airfoils)
    return CCBlade.Section.(r, c, theta, ccblade_airfoils)
end

"""
    ccbladeHAWTOperatingPoints(radial_positions, inflow_speed, rotor_speed, rho;
                               pitch=0.0, mu=1.7894e-5, asound=340.0,
                               precone=0.0)

Construct CCBlade operating points for a rigid HAWT rotor at the supplied
radial stations. `inflow_speed` is the axial freestream speed, `rotor_speed` is
the rotor angular speed in rad/s, and `pitch`/`precone` are radians. Negative
`rotor_speed` is accepted for caller-side sign-convention studies, but
zero-speed operation is rejected because the BEM solve and nondimensional
power metrics are singular there.
"""
function ccbladeHAWTOperatingPoints(
    radial_positions,
    inflow_speed,
    rotor_speed,
    rho;
    pitch = 0.0,
    mu = 1.7894e-5,
    asound = 340.0,
    precone = 0.0,
)
    r = _validated_hawt_station_vector(
        radial_positions,
        "radial_positions";
        positive = true,
    )
    all(diff(r) .> zero(first(r))) ||
        throw(ArgumentError("radial_positions must be strictly increasing"))
    _validate_hawt_operating_values(
        inflow_speed,
        rotor_speed,
        rho,
        pitch,
        mu,
        asound,
        precone,
    )
    return CCBlade.simple_op.(inflow_speed, rotor_speed, r, rho; pitch, mu, asound, precone)
end

function _collect_ccblade_outputs(outputs, field)
    hasproperty(outputs, field) && return getproperty(outputs, field)
    return [getfield(output, field) for output in outputs]
end

"""
    ccbladeHAWTSolve(radial_positions, chord, twist, airfoils, rotor_speed,
                     inflow_speed, rho; num_blades=3, hub_radius=0.0,
                     tip_radius=maximum(radial_positions), pitch=0.0,
                     precone=0.0, mu=1.7894e-5, asound=340.0, npts=10,
                     tip_correction=CCBlade.PrandtlTipHub())

Run a rigid-rotor CCBlade HAWT BEM solve and return a named tuple containing
the CCBlade rotor, sections, operating points, per-station outputs, integrated
thrust/torque/power, and nondimensional `CP`, `CT`, and `CQ`.
`tip_correction` is passed directly to `CCBlade.Rotor`; by default this uses
CCBlade's Prandtl tip-plus-hub correction, while `tip_correction = nothing`
disables the correction for comparison studies.

This is an adapter for HAWT aero verification and examples. It does not yet
map HAWT blade loads into OWENS structural dynamics and it does not apply Oye
dynamic inflow by itself; callers can use `oyeDynamicInflowStep` with the
returned quasi-steady induction arrays when building unsteady HAWT studies.
"""
function ccbladeHAWTSolve(
    radial_positions,
    chord,
    twist,
    airfoils,
    rotor_speed,
    inflow_speed,
    rho;
    num_blades = 3,
    hub_radius = zero(first(radial_positions)),
    tip_radius = maximum(radial_positions),
    pitch = 0.0,
    precone = 0.0,
    mu = 1.7894e-5,
    asound = 340.0,
    npts = 10,
    tip_correction = CCBlade.PrandtlTipHub(),
)
    num_blades isa Integer && num_blades > 0 ||
        throw(ArgumentError("num_blades must be a positive integer"))
    npts isa Integer && npts > 0 || throw(ArgumentError("npts must be a positive integer"))

    r, c, theta = _validate_hawt_station_inputs(radial_positions, chord, twist)
    _validate_hawt_rotor_geometry(r, hub_radius, tip_radius)
    _validate_hawt_operating_values(
        inflow_speed,
        rotor_speed,
        rho,
        pitch,
        mu,
        asound,
        precone,
    )
    _validate_hawt_tip_correction(tip_correction)

    rotor = CCBlade.Rotor(
        hub_radius,
        tip_radius,
        num_blades;
        precone,
        turbine = true,
        tip = tip_correction,
    )
    sections = ccbladeHAWTSections(r, c, theta, airfoils)
    operating_points = ccbladeHAWTOperatingPoints(
        r,
        inflow_speed,
        rotor_speed,
        rho;
        pitch,
        mu,
        asound,
        precone,
    )
    outputs = CCBlade.solve.(Ref(rotor), sections, operating_points; npts)
    thrust, torque = CCBlade.thrusttorque(rotor, sections, outputs)
    power = torque * rotor_speed
    swept_area = pi * tip_radius^2
    dynamic_pressure_area = 0.5 * rho * swept_area * inflow_speed^2

    return (
        rotor = rotor,
        sections = sections,
        operating_points = operating_points,
        outputs = outputs,
        thrust = thrust,
        torque = torque,
        power = power,
        CP = power / (dynamic_pressure_area * inflow_speed),
        CT = thrust / dynamic_pressure_area,
        CQ = torque / (dynamic_pressure_area * tip_radius),
        Np = _collect_ccblade_outputs(outputs, :Np),
        Tp = _collect_ccblade_outputs(outputs, :Tp),
        a = _collect_ccblade_outputs(outputs, :a),
        ap = _collect_ccblade_outputs(outputs, :ap),
        alpha = _collect_ccblade_outputs(outputs, :alpha),
        cl = _collect_ccblade_outputs(outputs, :cl),
        cd = _collect_ccblade_outputs(outputs, :cd),
        cn = _collect_ccblade_outputs(outputs, :cn),
        ct = _collect_ccblade_outputs(outputs, :ct),
        phi = _collect_ccblade_outputs(outputs, :phi),
        W = _collect_ccblade_outputs(outputs, :W),
        F = _collect_ccblade_outputs(outputs, :F),
        G = _collect_ccblade_outputs(outputs, :G),
    )
end
