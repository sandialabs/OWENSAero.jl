module OWENSAero
import Statistics:mean
import Interpolations
import OWENSOpenFASTWrappers
# Common
export Unsteady_Step
export Turbine, Environment, UnsteadyParams
export cpValidationMetrics
export wholeRevolutionIndexRange, wholeRevolutionMean
export jointDragForce

# Actuator Cylinder
export AC, radialforce, pInt

# DMS
export DMS, streamtube, readaerodyn, readaerodyn_BV

# Dynamic Stall
export Boeing_Vertol
# export Leishman_Beddoes

# Unsteady Method
export Unsteady_Step

# Module Path
const path,_ = splitdir(@__FILE__)

# Common Structs

function _canonical_dynamic_stall_model(model)
    token = lowercase(strip(string(model)))
    if token in ("bv", "boeing-vertol", "boeing_vertol")
        return "BV"
    elseif token in ("none", "no", "off", "nods", "no_ds", "no-ds")
        return "none"
    elseif token == "analytic"
        return "analytic"
    elseif token in ("lb", "leishman-beddoes", "leishman_beddoes")
        throw(ArgumentError("DynamicStallModel = \"LB\" is not implemented; use \"BV\" or \"none\"."))
    end
    throw(ArgumentError("DynamicStallModel must be \"BV\", \"none\", or the internal \"analytic\" test mode; got $(repr(model))."))
end

function _canonical_aero_model(model)
    token = uppercase(strip(string(model)))
    if token == "DMS"
        return "DMS"
    elseif token == "AC"
        return "AC"
    end
    throw(ArgumentError("AeroModel must be \"DMS\" or \"AC\"; got $(repr(model))."))
end

"""
    Turbine(R::TF,r::TAF,z::TF,chord::TAF3,twist::TAF5,delta::TAF,omega::TAF4,B::TI,af::TFN,ntheta::TI,r_delta_influence::TB,centerX::TAF2,centerY::TAF2)
    Turbine(R,r,z,chord,twist,delta,omega,B,af,ntheta,r_delta_infl) = Turbine(R,r,z,chord,twist,delta,omega,B,af,ntheta,r_delta_infl,zeros(Real,size(R)),zeros(Real,size(R)))
    Turbine(R,r,chord,twist,delta,omega,B,af,ntheta,r_delta_infl) = Turbine(R,r,1.0,chord,twist,delta,omega,B,af,ntheta,r_delta_infl,zeros(Real,size(R)),zeros(Real,size(R)))

Contains specications for turbine slice (geometry, location, airfoil)

# Inputs
* `R::TF`: Nominal turbine radius (m)
* `r::TAF`: Array of local radaii corresponding to each azimuthal position for the slice, allows for active blade deformation (m)
* `z::TF`: Vertical location of slice (only used when calling inflow-wind turbulent input)(m)
* `chord::TAF3`: Array of chord corresponding to each azimuthal position, allows for active blade deformation (m)
* `twist::TAF5`: Array of blade twist corresponding to each azimuthal position, allows for active blade deformation (rad)
* `delta::TAF`: Array of blade slope corresponding to each azimuthal position, allows for active blade deformation (rad)
* `omega::TAF4`: Array of rotational rate corresponding to each azimuthal position, allows for active blade deformation (rad/s)
* `B::TI`: Number of blades
* `af::TFN`: Airfoil function - see tests for example of how to create
* `ntheta::TI`: Number of azimuthal discretizations
* `r_delta_influence::TB`: Specification of whether local radius and blade slope are used in the influence coefficients for the actuator cylinder method
* `centerX::TAF2`: Turbine center x location (only used if multiple turbines are modeled)
* `centerY::TAF2`: Turbine center y location (only used if multiple turbines are modeled)

# Outputs:
* `none`:

"""
struct Turbine{TF1,TF2,TI1,TI2,TAF0,TAF1,TAF2,TAF3,TAF4,TAF5,TAF6,TAF7,TFN,TB,TAI}
    R::TF1
    r::TAF1
    z::TF2
    chord::TAF6
    thick::TAF3
    twist::TAF5
    delta::TAF0
    omega::TAF4
    B::TI1
    af::TFN
    ntheta::TI2
    r_delta_influence::TB
    centerX::TAF2
    centerY::TAF2
    helical_offset::TAI
    rhoA::TAF7
end

Turbine(R,r,z,chord,twist,delta,omega,B,af,ntheta,r_delta_infl) = Turbine(R,r,z,chord,0.18,twist,delta,omega,B,af,ntheta,r_delta_infl,zeros(Real,size(R)),zeros(Real,size(R)),zeros(1),zeros(Real,size(R)))
Turbine(R,r,chord,twist,delta,omega,B,af,ntheta,r_delta_infl) = Turbine(R,r,1.0,chord,0.18,twist,delta,omega,B,af,ntheta,r_delta_infl,zeros(Real,size(R)),zeros(Real,size(R)),zeros(1),zeros(Real,size(R)))

"""
Environment(rho::TF,mu::TF,V_x::TAF #Vinf is Vx,V_y::TAF,V_z::TAF,V_twist::TAF,windangle::TF #radians,DynamicStallModel::TS,AeroModel::TS,aw_warm::TVF,steplast::TAI,idx_RPI::TAI,V_wake_old::TVF2,BV_DynamicFlagL::TAI,BV_DynamicFlagD::TAI,alpha_last::TAF2,suction::TB)
Environment(rho,mu,V_x,V_y,V_z,V_twist,windangle,DynamicStallModel,AeroModel,aw_warm) = Environment(rho,mu,V_x,V_y,V_z,V_twist,windangle,DynamicStallModel,AeroModel,aw_warm,zeros(Int,1),zeros(Int,length(V_x)),deepcopy(V_x),zeros(Int,1),zeros(Int,1),zeros(Real,1),false)
Environment(rho,mu,V_x,DynamicStallModel,AeroModel,aw_warm) = Environment(rho,mu,V_x,zeros(Real,size(V_x)),zeros(Real,size(V_x)),zeros(Real,size(V_x)),0.0,DynamicStallModel,AeroModel,aw_warm,zeros(Int,1),zeros(Int,length(V_x)),deepcopy(V_x),zeros(Int,1),zeros(Int,1),zeros(Real,1),false)

Contains specications for turbine slice environment/operating conditions as well as some backend memory for dynamic stall and unsteady calculations

# Inputs
* `rho::TF`: Working fluid density (kg/m^3)
* `mu::TF`: Working fluid viscosity (standard SI units)
* `V_x::TAF` Vinf is Vx for simple simulations (m/s), array corresponding to each azimuthal position
* `V_y::TAF`: y input velocity (m/s), array corresponding to each azimuthal position
* `V_z::TAF`: z input velocity (m/s), array corresponding to each azimuthal position
* `V_twist::TAF`: rotational velocity from active twist (rad/s), array corresponding to each azimuthal position
* `windangle::TF`: angle of mean oncoming wind (rad)
* `DynamicStallModel::TS`: dynamic stall model ("BV" or "none" or "LB" - once it is finished)
* `AeroModel::TS`: aero model used ("DMS" or "AC")
* `aw_warm::TVF`: warm start induction factor array, first half corresponding to u, second half to v
* `steplast::TAI`: prior simulation step index, used for unsteady wake propogation
* `idx_RPI::TAI`: used to specify the azimuthal indices needed for a partial solve (i.e. not every azimuthal index), such as is used in the RPI method
* `V_wake_old::TVF2`: Prior step's mean wake velocity (m/s)
* `BV_DynamicFlagL::TAI`: Boeing-vertol dynamic stall lift flag
* `BV_DynamicFlagD::TAI`: Boeing-vertol dynamic stall drag flag
* `alpha_last::TAF2`: Boeing-vertol dynamic stall prior step's angle of attack
* `suction::TB`: DMS flag for alternate induction model

# Outputs:
* `none`:

"""
struct Environment{TF,TB,TAFx,TAFy,TAF2,TS1,TS2,TVF,TVF2,TAI,TAF3,TAF4,TAF5}
    rho::TF
    mu::TF
    V_x::TAFx #Vinf is Vx
    V_y::TAFy
    V_z::TAF3
    V_twist::TAF3
    windangle::TF #radians
    DynamicStallModel::TS1
    AeroModel::TS2
    Aero_AddedMass_Active::TB
    Aero_Buoyancy_Active::TB
    Aero_RotAccel_Active::TB
    AddedMass_Coeff_Ca::TF
    centrifugal_force_flag::TB
    aw_warm::TVF
    steplast::TAI
    idx_RPI::TAI
    V_wake_old::TVF2
    BV_DynamicFlagL::TAI
    BV_DynamicFlagD::TAI
    alpha_last::TAF2
    suction::TB
    accel_flap::TAF4
    accel_edge::TAF4
    gravity::TAF5
end
Environment(rho,mu,V_x,V_y,V_z,V_twist,windangle,DynamicStallModel,AeroModel,aw_warm) =
    Environment(rho,mu,V_x,V_y,V_z,V_twist,windangle,_canonical_dynamic_stall_model(DynamicStallModel),_canonical_aero_model(AeroModel),false,false,false,1.0,false,aw_warm,zeros(Int,1),zeros(Int,length(V_x)),deepcopy(V_x),zeros(Int,length(V_x)),zeros(Int,length(V_x)),zeros(Real,length(V_x)),false,zeros(Real,length(V_x)),zeros(Real,length(V_x)),[0.0,0.0,-9.81])
Environment(rho,mu,V_x,V_y,V_z,V_twist,windangle,DynamicStallModel,AeroModel,Aero_AddedMass_Active,Aero_Buoyancy_Active,Aero_RotAccel_Active,AddedMass_Coeff_Ca,centrifugal_force_flag,aw_warm) =
    Environment(rho,mu,V_x,V_y,V_z,V_twist,windangle,_canonical_dynamic_stall_model(DynamicStallModel),_canonical_aero_model(AeroModel),Aero_AddedMass_Active,Aero_Buoyancy_Active,Aero_RotAccel_Active,AddedMass_Coeff_Ca,centrifugal_force_flag,aw_warm,zeros(Int,1),zeros(Int,length(V_x)),deepcopy(V_x),zeros(Int,length(V_x)),zeros(Int,length(V_x)),zeros(Real,length(V_x)),false,zeros(Real,length(V_x)),zeros(Real,length(V_x)),[0.0,0.0,-9.81])
Environment(rho,mu,V_x,DynamicStallModel,AeroModel,aw_warm) =
    Environment(rho,mu,V_x,zeros(Real,size(V_x)),zeros(Real,size(V_x)),zeros(Real,size(V_x)),0.0,_canonical_dynamic_stall_model(DynamicStallModel),_canonical_aero_model(AeroModel),false,false,false,1.0,false,aw_warm,zeros(Int,1),zeros(Int,length(V_x)),deepcopy(V_x),zeros(Int,length(V_x)),zeros(Int,length(V_x)),zeros(Real,length(V_x)),false,zeros(Real,length(V_x)),zeros(Real,length(V_x)),[0.0,0.0,-9.81])

"""
UnsteadyParams(RPI::TB,tau::TAF,ifw::TB,IECgust::TB,nominalVinf::TF,G_amp::TF,gustX0::TF,gustT::TF)
UnsteadyParams(RPI,tau,ifw) = UnsteadyParams(RPI,tau,ifw,false,1.0,0.0,1.0,1.0)

Contains specications for turbine slice unsteady inputs

# Inputs
* `RPI::TB`: Flag to specify if RPI is being used
* `tau::TAF`: Unsteady method wake propogation weighting [3.0,0.3]
* `ifw::TB`: Flag to specify if inflow-wind is being used
* `IECgust::TB`: Flag to specify if the simple sin-cos gust profile in the x-direction will be used
* `nominalVinf::TF`: Nominal velocity used to calculate the IEC gust size (m/s)
* `G_amp::TF`: IEC gust amplitude (m/s)
* `gustX0::TF`: IEC gust normalized starting point (x-location divided by reference radius)
* `gustT::TF`: IEC gust duration (s)

# Outputs:
* `none`:

"""
struct UnsteadyParams{TB,TF,TAF}
    RPI::TB
    tau::TAF
    ifw::TB
    IECgust::TB
    nominalVinf::TF
    G_amp::TF
    gustX0::TF
    gustT::TF
end

UnsteadyParams(RPI,tau,ifw) = UnsteadyParams(RPI,tau,ifw,false,1.0,0.0,1.0,1.0)

"""
    added_mass_flap_volume_per_unit_span(chord)

Return the current OWENSAero flap-direction added-mass reference volume per
unit span. The model uses a circular projected section with diameter equal to
the local chord.
"""
added_mass_flap_volume_per_unit_span(chord) = pi * (chord / 2)^2

"""
    added_mass_edge_volume_per_unit_span(thickness)

Return the current OWENSAero edge-direction added-mass reference volume per
unit span. This preserves the legacy `thickness / 10` projected-diameter
convention so the formula is explicit and test-pinned.
"""
added_mass_edge_volume_per_unit_span(thickness) = pi * ((thickness / 10) / 2)^2

"""
    buoyancy_section_area_per_unit_span(chord, thickness)

Return the current triangular-section buoyancy area per unit span used by the
DMS and actuator-cylinder paths.
"""
buoyancy_section_area_per_unit_span(chord, thickness) = chord * thickness / 2

"""
    jointDragForce(rho, velocity, CdA)

Return the lumped bluff-body drag force vector for a joint or other ancillary
body in the same frame as `velocity`. `CdA` is the drag coefficient times
reference area in square meters, and the returned force opposes the supplied
relative velocity. This helper is not coupled into DMS or AC induction.
"""
function jointDragForce(rho, velocity, CdA)
    rho isa Real && isfinite(rho) && rho >= 0 ||
        throw(ArgumentError("rho must be a finite nonnegative real value"))
    CdA isa Real && isfinite(CdA) && CdA >= 0 ||
        throw(ArgumentError("CdA must be a finite nonnegative real value"))
    velocity isa AbstractVector && length(velocity) in (2, 3) ||
        throw(ArgumentError("velocity must be a two- or three-component vector"))
    all(x -> x isa Real && isfinite(x), velocity) ||
        throw(ArgumentError("velocity must contain only finite real values"))

    speed = sqrt(sum(abs2, velocity))
    return @. -0.5 * rho * CdA * speed * velocity
end

function _finite_real_vector(values, name)
    values isa AbstractVector || throw(ArgumentError("$(name) must be a vector"))
    isempty(values) && throw(ArgumentError("$(name) must not be empty"))
    all(x -> x isa Real, values) ||
        throw(ArgumentError("$(name) must contain only real values"))
    vector = Float64.(collect(values))
    all(isfinite, vector) || throw(ArgumentError("$(name) must contain only finite values"))
    return vector
end

function _sorted_curve_points(x, y, label; minimum_points = 1)
    x_vector = _finite_real_vector(x, "$(label)_x")
    y_vector = _finite_real_vector(y, "$(label)_y")
    length(x_vector) == length(y_vector) ||
        throw(ArgumentError("$(label) x and y vectors must have the same length"))
    length(x_vector) >= minimum_points ||
        throw(ArgumentError("$(label) must contain at least $(minimum_points) points"))

    order = sortperm(x_vector)
    x_sorted = x_vector[order]
    y_sorted = y_vector[order]
    all(diff(x_sorted) .> 0.0) ||
        throw(ArgumentError("$(label) x values must be unique"))
    return x_sorted, y_sorted
end

function _linear_interpolate_sorted(x, y, x_query)
    if x_query < first(x) || x_query > last(x)
        throw(ArgumentError("query point $(x_query) is outside the model curve range"))
    end
    idx = searchsortedlast(x, x_query)
    idx == length(x) && return y[end]
    idx == 0 && return y[1]
    x0 = x[idx]
    x1 = x[idx+1]
    y0 = y[idx]
    y1 = y[idx+1]
    fraction = (x_query - x0) / (x1 - x0)
    return y0 + fraction * (y1 - y0)
end

"""
    cpValidationMetrics(model_tsr, model_cp, reference_tsr, reference_cp)

Compare a modeled CP curve against reference CP data on the overlapping
reference TSR points. The model curve is linearly interpolated to each
reference point, and the returned named tuple includes RMSE, bias, mean
absolute error, maximum absolute error, and peak-CP diagnostics.
"""
function cpValidationMetrics(model_tsr, model_cp, reference_tsr, reference_cp)
    model_x, model_y =
        _sorted_curve_points(model_tsr, model_cp, "model"; minimum_points = 2)
    reference_x, reference_y =
        _sorted_curve_points(reference_tsr, reference_cp, "reference")

    overlap = findall(x -> first(model_x) <= x <= last(model_x), reference_x)
    isempty(overlap) && throw(ArgumentError("reference curve does not overlap model curve"))

    reference_x_overlap = reference_x[overlap]
    reference_y_overlap = reference_y[overlap]
    model_on_reference =
        [_linear_interpolate_sorted(model_x, model_y, x) for x in reference_x_overlap]
    error = model_on_reference .- reference_y_overlap
    abs_error = abs.(error)

    model_peak_idx = argmax(model_y)
    reference_peak_idx = argmax(reference_y_overlap)
    return (
        n = length(reference_x_overlap),
        rmse = sqrt(mean(error .^ 2)),
        mean_bias = mean(error),
        mean_abs_error = mean(abs_error),
        max_abs_error = maximum(abs_error),
        model_peak_tsr = model_x[model_peak_idx],
        model_peak_cp = model_y[model_peak_idx],
        reference_peak_tsr = reference_x_overlap[reference_peak_idx],
        reference_peak_cp = reference_y_overlap[reference_peak_idx],
        peak_cp_error = model_y[model_peak_idx] - reference_y_overlap[reference_peak_idx],
        model_cp_on_reference = model_on_reference,
        reference_tsr = reference_x_overlap,
        reference_cp = reference_y_overlap,
    )
end

"""
    wholeRevolutionIndexRange(azimuth; revolutions=nothing, period=2*pi, atol=1e-10, allow_partial=false)

Return a suffix index range that covers complete revolutions in a monotonically
increasing azimuth history. The terminal repeated phase is excluded so arithmetic
means do not double-count the revolution boundary.

When `revolutions` is omitted, all complete revolutions available at the end of
the signal are used. If `revolutions` exceeds the number available, the request
is clamped to the available count.
"""
function wholeRevolutionIndexRange(
    azimuth;
    revolutions = nothing,
    period = 2*pi,
    atol = 1e-10,
    allow_partial = false,
)
    isempty(azimuth) && throw(ArgumentError("azimuth history must not be empty"))
    period <= 0 && throw(ArgumentError("period must be positive"))

    az = collect(float.(azimuth))
    all(isfinite, az) ||
        throw(ArgumentError("azimuth history must contain only finite values"))
    if any(diff(az) .< -atol)
        throw(ArgumentError("azimuth history must be monotonically increasing"))
    end

    available_revolutions = floor(Int, max((az[end] - az[1]) / period, 0.0) + atol)
    if available_revolutions < 1
        allow_partial || throw(
            ArgumentError("azimuth history must span at least one complete revolution"),
        )
        return firstindex(azimuth):lastindex(azimuth)
    end

    requested_revolutions = if isnothing(revolutions)
        available_revolutions
    else
        requested = try
            Int(revolutions)
        catch
            throw(ArgumentError("revolutions must be an integer count"))
        end
        requested == revolutions ||
            throw(ArgumentError("revolutions must be an integer count"))
        requested
    end
    n_revolutions = min(requested_revolutions, available_revolutions)
    n_revolutions < 1 && throw(ArgumentError("revolutions must be at least 1"))

    stop_azimuth = az[end]
    start_azimuth = stop_azimuth - n_revolutions * period
    local_start = findfirst(x -> x >= start_azimuth - atol, az)
    local_stop = findlast(x -> x < stop_azimuth - atol, az)
    isnothing(local_start) &&
        throw(ArgumentError("whole-revolution window has no start sample"))
    isnothing(local_stop) &&
        throw(ArgumentError("whole-revolution window has no stop sample"))

    idx_start = firstindex(azimuth) + local_start - 1
    idx_stop = firstindex(azimuth) + local_stop - 1
    idx_start <= idx_stop || throw(
        ArgumentError("whole-revolution window has no samples before the terminal phase"),
    )
    return idx_start:idx_stop
end

"""
    wholeRevolutionMean(values, azimuth; kwargs...)

Compute the arithmetic mean of `values` over the window returned by
[`wholeRevolutionIndexRange`](@ref). This is intended for unsteady CP, torque,
and load histories where partial revolutions should not bias validation metrics.
"""
function wholeRevolutionMean(values, azimuth; kwargs...)
    length(values) == length(azimuth) ||
        throw(ArgumentError("values and azimuth must have the same length"))
    window = wholeRevolutionIndexRange(azimuth; kwargs...)
    return mean(values[window])
end

"""
    steady(turbine::Turbine, env::Env; w=zeros(Real,2*turbine.ntheta), idx_RPI=1:2*turbine.ntheta,solve=true,ifw=false)

Calculates steady state aerodynamics for a single VAWT slice

# Inputs
* `turbine::Turbine`: Turbine struct, see ?Turbine for details
* `env::Env`: Env struct, see ?Env for details
* `w::Array(<:Real)`: Optional, used if solve=false, induction factor array, first half corresponding to u, second half to v
* `idx_RPI::Array(<:Int)`: Optional, used to specify the azimuthal indices needed for a partial solve (i.e. not every azimuthal index), such as is used in the RPI method
* `solve::Bool`: Optional, False is used when you want the model outputs for a given set of induction factors without resolving them.
* `ifw::Bool`: Optional, used to tell the Vinf lookup to attempt to use the dynamic inflow wind library, requires preprocessing as is shown in the test cases.


# Outputs:
* `CP`: This slice's coefficient of performance
* `Th`: This slice's thrust coefficient
* `Q`: Torque (N0m)
* `Rp`: Radial force per height (N)
* `Tp`: Tangential force per height (N)
* `Zp`: Vertical force per height (N)
* `Vloc`: Local velocity array for each azimuthal position (includes induction) (m/s)
* `CD`: This slice's drag coefficient
* `CT`: This slice's thrust coefficient (should equal drag, but may no depending on usage or solver status)
* `amean`: Mean turbine induction in the streamwise direction
* `astar`: Solved induction factors for each azimuthal location. First half are streamwise (u), second are cross-steam (v)
* `alpha`: Local angle of attack array for each azimuthal position (includes induction) (rad)
* `cl`: Local lift coefficient used for each azimuthal position
* `cd_af`: Local drag coefficient used for each azimuthal position
* `thetavec`: Azimuthal location of each discretization (rad)
* `Re`: Reynolds number for each azimuthal position
"""
function steady(turbine, env; w=zeros(Real,2*turbine.ntheta), idx_RPI=1:2*turbine.ntheta,solve=true,ifw=false)
    if env.AeroModel=="DMS"
        return DMS(turbine, env; w, idx_RPI, solve)
    elseif env.AeroModel=="AC"
        turbines = Array{OWENSAero.Turbine}(undef,1)
        turbines[1] = turbine
        return AC(turbines, env; w, idx_RPI, solve, ifw)
        # return AC_steady(turbines, env)
    else
        error("AeroModel not recognized, choose DMS or AC")
    end
end

@inline function safeakima(x,y,xpt)
    if minimum(xpt)<(minimum(x)-abs(minimum(x))*0.1) || maximum(xpt)>(maximum(x)+abs(maximum(x))*0.1)
        msg="Extrapolating on akima spline results in undefined solutions minimum(xpt)<minimum(x) $(minimum(xpt))<$(minimum(x)) or maximum(xpt)<maximum(x) $(maximum(xpt))>$(maximum(x))"
        throw(OverflowError(msg))
    end
    return FLOWMath.akima(x,y,xpt)
end

_zero_like(x) = zero(x)

function _split_airfoil_coefficients(coefficients)
    cl = coefficients[1]
    cd = coefficients[2]
    cm = length(coefficients) >= 3 ? coefficients[3] : _zero_like(cl)
    return cl, cd, cm
end

function _airfoil_coefficients(af, alpha, Re, mach)
    coefficients = try
        af(alpha, Re, mach; return_cm = true)
    catch err
        if err isa MethodError
            af(alpha, Re, mach)
        else
            rethrow()
        end
    end
    return _split_airfoil_coefficients(coefficients)
end

function _airfoil_coefficients(
    af,
    alpha,
    Re,
    mach,
    env,
    V_twist,
    chord,
    dt,
    U;
    solvestep = false,
    idx = 1,
)
    coefficients = try
        af(
            alpha,
            Re,
            mach,
            env,
            V_twist,
            chord,
            dt,
            U;
            solvestep,
            idx,
            return_cm = true,
        )
    catch err
        if err isa MethodError
            af(alpha, Re, mach, env, V_twist, chord, dt, U; solvestep, idx)
        else
            rethrow()
        end
    end
    return _split_airfoil_coefficients(coefficients)
end

include("DMS.jl")
include("./vawt-ac/src/airfoilread.jl") #TODO: switch for the CCBlade airfoil reading library
include("./vawt-ac/src/acmultiple.jl")
include("Unsteady_Step.jl")
include("Boeing_Vertol.jl")
include("advanceTurbine.jl")

end #module
