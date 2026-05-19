function _validate_finite_real_value(value, name)
    value isa Real && isfinite(value) ||
        throw(ArgumentError("$name must be a finite real value"))
    return value
end

function _validate_positive_real_value(value, name)
    _validate_finite_real_value(value, name)
    value > zero(value) || throw(ArgumentError("$name must be a finite positive real value"))
    return value
end

function _validate_finite_real_input(value, name)
    if value isa AbstractArray
        isempty(value) && throw(ArgumentError("$name must not be empty"))
        all(x -> x isa Real && isfinite(x), value) ||
            throw(ArgumentError("$name must contain only finite real values"))
    else
        _validate_finite_real_value(value, name)
    end
    return value
end

function _validate_positive_real_input(value, name)
    _validate_finite_real_input(value, name)
    all(x -> x > zero(x), value isa AbstractArray ? value : (value,)) ||
        throw(ArgumentError("$name must contain only finite positive real values"))
    return value
end

function _validate_nonnegative_real_input(value, name)
    _validate_finite_real_input(value, name)
    all(x -> x >= zero(x), value isa AbstractArray ? value : (value,)) ||
        throw(ArgumentError("$name must contain only finite nonnegative real values"))
    return value
end

function _validate_oye_k(k)
    _validate_finite_real_value(k, "k")
    zero(k) <= k <= one(k) || throw(ArgumentError("k must be between 0 and 1"))
    return k
end

"""
    oyeDynamicInflowTimeConstants(rotor_radius, wind_speed, axial_induction, radial_position)

Return Oye dynamic-inflow time constants `(tau1, tau2)` for a HAWT radial
station or vector of radial stations. `rotor_radius`, `wind_speed`, and
`radial_position` use consistent dimensional units; `axial_induction` is the
rotor-average axial induction used in the AeroDyn DBEMT `tau1` expression.
"""
function oyeDynamicInflowTimeConstants(
    rotor_radius,
    wind_speed,
    axial_induction,
    radial_position,
)
    _validate_positive_real_value(rotor_radius, "rotor_radius")
    _validate_positive_real_value(wind_speed, "wind_speed")
    _validate_finite_real_value(axial_induction, "axial_induction")
    _validate_nonnegative_real_input(radial_position, "radial_position")
    all(r -> r <= rotor_radius, radial_position isa AbstractArray ? radial_position : (radial_position,)) ||
        throw(ArgumentError("radial_position must be less than or equal to rotor_radius"))

    tau1 =
        1.1 / (1 - 1.3 * min(axial_induction, one(axial_induction) / 2)) *
        rotor_radius / wind_speed
    tau2 = @. (0.39 - 0.26 * (radial_position / rotor_radius)^2) * tau1
    return tau1, tau2
end

"""
    oyeDynamicInflowDerivative(reduced_induction, dynamic_induction,
                               quasi_steady_induction, tau1, tau2; k=0.6)

Return the continuous-time Oye dynamic-inflow state derivatives
`(d_reduced_induction, d_dynamic_induction)`. The inputs may be scalars or
arrays that broadcast together. The state follows AeroDyn's continuous DBEMT
form where `reduced_induction + k * quasi_steady_induction` is the intermediate
inflow state.
"""
function oyeDynamicInflowDerivative(
    reduced_induction,
    dynamic_induction,
    quasi_steady_induction,
    tau1,
    tau2;
    k = 0.6,
)
    _validate_finite_real_input(reduced_induction, "reduced_induction")
    _validate_finite_real_input(dynamic_induction, "dynamic_induction")
    _validate_finite_real_input(quasi_steady_induction, "quasi_steady_induction")
    _validate_positive_real_input(tau1, "tau1")
    _validate_positive_real_input(tau2, "tau2")
    _validate_oye_k(k)

    reduced_rate = @. ((1 - k) * quasi_steady_induction - reduced_induction) / tau1
    dynamic_rate =
        @. (reduced_induction + k * quasi_steady_induction - dynamic_induction) / tau2
    return reduced_rate, dynamic_rate
end

function _oye_dynamic_inflow_coupling_gain(dt, tau1, tau2)
    if tau1 == tau2
        return dt / tau2 * exp(-dt / tau2)
    end
    return (exp(-dt / tau1) - exp(-dt / tau2)) / (tau2 * (1 / tau2 - 1 / tau1))
end

"""
    oyeDynamicInflowStep(reduced_induction, dynamic_induction,
                         quasi_steady_induction, dt, tau1, tau2; k=0.6)

Advance the continuous-time Oye dynamic-inflow state exactly for a constant
quasi-steady induction over time step `dt`. Returns
`(next_reduced_induction, next_dynamic_induction)`.

This is intended as the unsteady inflow primitive for future CCBlade-based HAWT
work. It is deliberately separate from the DMS wake-speed filter because HAWT
dynamic inflow should track axial and tangential induction states at radial
stations.
"""
function oyeDynamicInflowStep(
    reduced_induction,
    dynamic_induction,
    quasi_steady_induction,
    dt,
    tau1,
    tau2;
    k = 0.6,
)
    _validate_finite_real_input(reduced_induction, "reduced_induction")
    _validate_finite_real_input(dynamic_induction, "dynamic_induction")
    _validate_finite_real_input(quasi_steady_induction, "quasi_steady_induction")
    _validate_nonnegative_real_input(dt, "dt")
    _validate_positive_real_input(tau1, "tau1")
    _validate_positive_real_input(tau2, "tau2")
    _validate_oye_k(k)

    reduced_steady = @. (1 - k) * quasi_steady_induction
    reduced_decay = @. exp(-dt / tau1)
    dynamic_decay = @. exp(-dt / tau2)
    coupling_gain = _oye_dynamic_inflow_coupling_gain.(dt, tau1, tau2)

    next_reduced =
        @. reduced_steady + (reduced_induction - reduced_steady) * reduced_decay
    next_dynamic =
        @. quasi_steady_induction +
           (dynamic_induction - quasi_steady_induction) * dynamic_decay +
           (reduced_induction - reduced_steady) * coupling_gain
    return next_reduced, next_dynamic
end
