module OWENSAero
import Statistics:mean
import Dierckx
import OWENSOpenFASTWrappers
# Common
export Unsteady_Step
export Turbine, Environment, UnsteadyParams

# Actuator Cylinder
export AC_steady, radialforce, pInt

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
struct Turbine{TF1,TF2,TI1,TI2,TAF0,TAF1,TAF2,TAF3,TAF4,TAF5,TFN,TB,TAI}
    R::TF1
    r::TAF1
    z::TF2
    chord::TAF3
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
end

Turbine(R,r,z,chord,twist,delta,omega,B,af,ntheta,r_delta_infl) = Turbine(R,r,z,chord,0.18,twist,delta,omega,B,af,ntheta,r_delta_infl,zeros(Real,size(R)),zeros(Real,size(R)),zeros(1))
Turbine(R,r,chord,twist,delta,omega,B,af,ntheta,r_delta_infl) = Turbine(R,r,1.0,chord,0.18,twist,delta,omega,B,af,ntheta,r_delta_infl,zeros(Real,size(R)),zeros(Real,size(R)),zeros(1))

"""
Environment(rho::TF,mu::TF,V_x::TAF #Vinf is Vx,V_y::TAF,V_z::TAF,V_twist::TAF,windangle::TF #radians,DSModel::TS,AModel::TS,aw_warm::TVF,steplast::TAI,idx_RPI::TAI,V_wake_old::TVF2,BV_DynamicFlagL::TAI,BV_DynamicFlagD::TAI,alpha_last::TAF2,suction::TB)
Environment(rho,mu,V_x,V_y,V_z,V_twist,windangle,DSModel,AModel,aw_warm) = Environment(rho,mu,V_x,V_y,V_z,V_twist,windangle,DSModel,AModel,aw_warm,zeros(Int,1),zeros(Int,length(V_x)),deepcopy(V_x),zeros(Int,1),zeros(Int,1),zeros(Real,1),false)
Environment(rho,mu,V_x,DSModel,AModel,aw_warm) = Environment(rho,mu,V_x,zeros(Real,size(V_x)),zeros(Real,size(V_x)),zeros(Real,size(V_x)),0.0,DSModel,AModel,aw_warm,zeros(Int,1),zeros(Int,length(V_x)),deepcopy(V_x),zeros(Int,1),zeros(Int,1),zeros(Real,1),false)

Contains specications for turbine slice environment/operating conditions as well as some backend memory for dynamic stall and unsteady calculations

# Inputs
* `rho::TF`: Working fluid density (kg/m^3)
* `mu::TF`: Working fluid viscosity (standard SI units)
* `V_x::TAF` Vinf is Vx for simple simulations (m/s), array corresponding to each azimuthal position
* `V_y::TAF`: y input velocity (m/s), array corresponding to each azimuthal position
* `V_z::TAF`: z input velocity (m/s), array corresponding to each azimuthal position
* `V_twist::TAF`: rotational velocity from active twist (rad/s), array corresponding to each azimuthal position
* `windangle::TF`: angle of mean oncoming wind (rad)
* `DSModel::TS`: dynamic stall model ("BV" or "none" or "LB" - once it is finished)
* `AModel::TS`: aero model used ("DMS" or "AC")
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
struct Environment{TF,TB,TAFx,TAFy,TAF2,TS,TVF,TVF2,TAI,TAF3}
    rho::TF
    mu::TF
    V_x::TAFx #Vinf is Vx
    V_y::TAFy
    V_z::TAF3
    V_twist::TAF3
    windangle::TF #radians
    DSModel::TS
    AModel::TS
    AM_flag::TB
    buoy_flag::TB
    rotAccel_flag::TB
    AM_Coeff_Ca::TF
    aw_warm::TVF
    steplast::TAI
    idx_RPI::TAI
    V_wake_old::TVF2
    BV_DynamicFlagL::TAI
    BV_DynamicFlagD::TAI
    alpha_last::TAF2
    suction::TB
end
Environment(rho,mu,V_x,V_y,V_z,V_twist,windangle,DSModel,AModel,aw_warm) = Environment(rho,mu,V_x,V_y,V_z,V_twist,windangle,DSModel,AModel,false,false,false,1.0,aw_warm,zeros(Int,1),zeros(Int,length(V_x)),deepcopy(V_x),zeros(Int,length(V_x)),zeros(Int,length(V_x)),zeros(Real,length(V_x)),false)
Environment(rho,mu,V_x,V_y,V_z,V_twist,windangle,DSModel,AModel,AM_flag,buoy_flag,rotAccel_flag,AM_Coeff_Ca,aw_warm) = Environment(rho,mu,V_x,V_y,V_z,V_twist,windangle,DSModel,AModel,AM_flag,buoy_flag,rotAccel_flag,AM_Coeff_Ca,aw_warm,zeros(Int,1),zeros(Int,length(V_x)),deepcopy(V_x),zeros(Int,length(V_x)),zeros(Int,length(V_x)),zeros(Real,length(V_x)),false)
Environment(rho,mu,V_x,DSModel,AModel,aw_warm) = Environment(rho,mu,V_x,zeros(Real,size(V_x)),zeros(Real,size(V_x)),zeros(Real,size(V_x)),0.0,DSModel,AModel,false,false,false,1.0,aw_warm,zeros(Int,1),zeros(Int,length(V_x)),deepcopy(V_x),zeros(Int,length(V_x)),zeros(Int,length(V_x)),zeros(Real,length(V_x)),false)

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
    if env.AModel=="DMS"
        return DMS(turbine, env; w, idx_RPI, solve)
    elseif env.AModel=="AC"
        turbines = Array{OWENSAero.Turbine}(undef,1)
        turbines[1] = turbine
        return AC(turbines, env; w, idx_RPI, solve, ifw)
        # return AC_steady(turbines, env)
    else
        error("AModel not recognized, choose DMS or AC")
    end
end

@inline function safeakima(x,y,xpt)
    if minimum(xpt)<minimum(x) || maximum(xpt)>maximum(x)
        msg="Extrapolating on akima spline results in undefined solutions minimum(xpt)<minimum(x) $(minimum(xpt))<$(minimum(x)) or maximum(xpt)<maximum(x) $(maximum(xpt))>$(maximum(x))"
        throw(OverflowError(msg))
    end
    return FLOWMath.akima(x,y,xpt)
end

include("DMS.jl")
include("./vawt-ac/src/airfoilread.jl") #TODO: switch for the CCBlade airfoil reading library
include("./vawt-ac/src/acmultiple.jl")
include("Unsteady_Step.jl")
include("Boeing_Vertol.jl")
include("advanceTurbine.jl")

end #module
