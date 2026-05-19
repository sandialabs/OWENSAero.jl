# This module is based on https://doi.org/10.5194/wes-2019-44 A Double Multiple
# Streamtube model for Vertical Axis Wind Turbines of arbitrary rotor loading
# The equation references align to this report
# Some equations also taken from (are explicitely referenced) doi:10.5194/wes-1-327-2016 Actuator cylinder theory for multiple vertical axis wind turbines

import FLOWMath
import Statistics
import Statistics: mean

"""
INTERNAL streamtube(a,theta,turbine,env;output_all=false,Vxwake=nothing,solvestep=false)

Double multiple streamtube individual streamtube calculation

Output:

if output_all
    return Th, Q, Rp, Tp, Zp, Vloc, CD, CT, alpha, cl, cd_af, Re,
           M_addedmass_Np, M_addedmass_Tp, F_addedmass_Np, F_addedmass_Tp,
           F_buoy, cm_af, M25
else
    return CD-CT # Residual, section 2.4
end

"""
function streamtube(
    a,
    theta,
    turbine,
    env;
    output_all = false,
    Vxwake = nothing,
    solvestep = false,
    finite_span_factor = 1.0,
)

    # Unpack Vars
    B = turbine.B
    k = 1.0
    af = turbine.af
    chord = turbine.chord[1]
    thickness = turbine.thick[1]
    ntheta = turbine.ntheta
    suction = false
    rho = env.rho
    mu = env.mu
    rhoA = turbine.rhoA
    gravity = env.gravity
    Aero_AddedMass_Active = env.Aero_AddedMass_Active
    centrifugal_force_flag = env.centrifugal_force_flag
    Aero_Buoyancy_Active = env.Aero_Buoyancy_Active
    Aero_RotAccel_Active = env.Aero_RotAccel_Active
    # winddir = env.winddir

    dtheta = 2 * pi / (ntheta) #Assuming discretization is fixed equidistant (but omega can change between each point)

    if length(turbine.r) > 1
        idx = round(Int, (theta + dtheta / 2) / dtheta)
        r = turbine.r[idx]
        twist = turbine.twist[idx]
        delta = turbine.delta[idx]
        omega = turbine.omega[idx]
        rotation = sign(turbine.omega[idx])
        V_xt = env.V_x[idx] # Vinf is V_xt, t is for turbine direction f.o.r.
        V_yt = env.V_y[idx]
        V_zt = env.V_z[idx]
        V_twist = env.V_twist[idx]
        accel_flap = env.accel_flap[idx]
        accel_edge = env.accel_edge[idx]
        # V_delta = env.V_delta[idx] # Does not apply since the model calculation is centered around the point of rotation
        # V_sweep = env.V_sweep[idx] # Does not apply since the model calculation is centered around the point of rotation

    else
        idx = 1
        r = turbine.r
        twist = turbine.twist
        delta = turbine.delta
        omega = turbine.omega
        rotation = sign(mean(turbine.omega))
        V_xt = env.V_x[1]
        V_yt = env.V_y[1]
        V_zt = env.V_z[1]
        V_twist = env.V_twist[1]
        accel_flap = env.accel_flap[1]
        accel_edge = env.accel_edge[1]
        # V_delta = env.V_delta[1] # Does not apply since the model calculation is centered around the point of rotation
        # V_sweep = env.V_sweep[1] # Does not apply since the model calculation is centered around the point of rotation
    end
    finite_span = _finite_span_factor_at(finite_span_factor, idx)

    if Vxwake != nothing
        V_xt = Vxwake
    end

    if suction
        # Drag Coefficient: Eq. 5
        if a <= 0.7
            CD = 4 / 3 * a * (3 - a) / (1 + a)
        elseif a > 0.7 && a <= 1
            CD = 0.889 - (0.0203 - (a - 0.143)^2) / 0.6427
        else
            error("a > 1 is undefined")
        end
    else
        # Drag Coefficient: Eq. 2
        if a <= 0.4
            CD = 4 * a * (1 - a)
        else
            CD = 0.889 - (0.0203 - (a - 0.143)^2) / 0.6427
        end
    end

    # Cycle-average Thrust Coefficient
    # TSR = abs(omega)*r/sqrt(V_x^2+V_y^2) # Tip Speed Ratio, section 2.3, however since this is being used for relative local velocities, we use the local radius r, and not the nominal radius R
    # Vt = V_x*(1-a)*cos(theta)+V_y*sin(theta)+V_y*(a)*cos(theta)+TSR

    Vn =
        (V_xt * (1 - a) * sin(theta) - V_yt * cos(theta) + V_yt * (a) * sin(theta)) *
        cos.(delta) + V_zt * sin(delta)
    Vt =
        V_xt * ((1 - a) * cos(theta)) +
        V_yt * sin(theta) +
        V_yt * (1 - a) * cos(theta) +
        abs(omega) * r
    Vloc = sqrt(Vn^2 + Vt^2)
    phi = atan(Vn, Vt)
    alpha = phi - twist

    # Original Equations: I broke them out above to be more readable and to know where to put the cos(delta)
    # phi = atan(((1-a)*sin(theta)).*cos.(delta) , ((1-a)*cos(theta)+TSR)) # Eq. 7 but correct for slope from the 3D to 2D frame of reverence, also must use atan2 to determine which quadrant and get the signs correct
    # Vloc = Vinf .* ((1-a)^2 + 2*(1-a)*TSR*cos(theta).*cos.(delta) + TSR^2).^0.5 # Eq. 8

    Re = rho * Vloc * chord / mu
    dt = dtheta / abs.(omega)
    v_sound = 343.0 #m/s #TODO: calculate this using Atmosphere.jl
    mach = Vloc / v_sound
    if env.DynamicStallModel == "BV"
        cl, cd_af, cm_af = OWENSAero._airfoil_coefficients(
            af,
            alpha,
            Re,
            mach,
            env,
            V_twist,
            chord,
            dt,
            Vloc;
            solvestep,
            idx,
        )
    elseif env.DynamicStallModel == "LB"
        error("LB Dynamic Stall Model Not Implemented Yet")
    else
        cl, cd_af, cm_af = OWENSAero._airfoil_coefficients(af, alpha, Re, mach)
    end
    ct = cd_af * cos(phi) - cl * sin(phi) # Eq. 9
    cn = cd_af * sin(phi) + cl * cos(phi) # Eq. 10

    if Aero_AddedMass_Active
        Vol_flap = added_mass_flap_volume_per_unit_span(chord)
        Vol_edge = added_mass_edge_volume_per_unit_span(thickness)

        if Aero_RotAccel_Active
            accel_rot = omega^2 * r
        else
            accel_rot = 0.0
        end

        M_addedmass_flap = rho * env.AddedMass_Coeff_Ca * Vol_flap
        M_addedmass_edge = rho * env.AddedMass_Coeff_Ca * Vol_edge

        F_addedmass_flap = M_addedmass_flap * (accel_flap + accel_rot)
        F_addedmass_edge = M_addedmass_edge * (accel_edge + accel_rot)

        M_addedmass_Np = M_addedmass_flap * cos(twist) + M_addedmass_edge * sin(twist) # Go from the beam frame of reference to the normal and tangential direction #TODO: verify the directions
        M_addedmass_Tp = M_addedmass_edge * cos(twist) - M_addedmass_flap * sin(twist)

        F_addedmass_Np = F_addedmass_flap * cos(twist) + F_addedmass_edge * sin(twist) # Go from the beam frame of reference to the normal and tangential direction #TODO: verify the directions
        F_addedmass_Tp = F_addedmass_edge * cos(twist) - F_addedmass_flap * sin(twist)

    else
        M_addedmass_Np = 0.0
        M_addedmass_Tp = 0.0
        F_addedmass_Np = 0.0
        F_addedmass_Tp = 0.0
    end

    if Aero_Buoyancy_Active
        section_area = buoyancy_section_area_per_unit_span(chord, thickness)
        dcm = [
            cos(theta) -sin(theta) 0
            sin(theta) cos(theta) 0
            0 0 1
        ]
        F_buoy = dcm * -gravity .* (rho * section_area - rhoA) #buoyancy mass minus structural mass since added mass requires moving the gravity here
    else
        F_buoy = [0.0, 0.0, 0.0]
    end

    if centrifugal_force_flag
        f_centrifugal = rhoA * omega^2 * r
    else
        f_centrifugal = 0.0
    end

    # Instantaneous Forces (Unit Height) #Based on this, radial is inward positive and tangential is in direction of rotation positive
    Ab = chord * 1.0 # planform area Assuming unit section height
    q_loc = 0.5 * rho * Ab * Vloc^2 # From Eq. 11
    aero_Np = finite_span * cn * q_loc
    aero_Tp = finite_span * ct * q_loc
    M25 = finite_span * cm_af * q_loc * chord

    Np = aero_Np + -F_addedmass_Np
    Rp = Np + F_buoy[2] + f_centrifugal# Ning Eq. 27 # Negate to match cactus frame of reference, note that delta cancels out
    Zp = Np * tan(delta) + F_buoy[3] # Ning Eq. 27 # Negate to match cactus frame of reference
    Tp = -rotation * (aero_Tp + -F_addedmass_Tp) / cos(delta) + F_buoy[1] # TODO: verify direction Ning Eq. 27 # Negate to match cactus frame of reference

    Th = finite_span * q_loc * (ct * cos(theta) + cn * sin(theta) / cos(delta)) # Eq. 11 but with delta correction

    Q = finite_span * q_loc * r * -ct # Eq. 12 but with Local radius for local torque, Negate the force for reaction torque, in the power frame of reference?

    Ast = 1.0 * r * dtheta * abs(sin(theta)) # Section 2.0, local radius for local area, however do not allow negative area, so use absolute value

    q_inf = 0.5 * rho * Ast * (V_xt^2 + V_yt^2) # this is in the radial frame of reference V_x.^2 , also (sqrt(V_x^2+V_y^2)).^2 is reduced

    CT = (k * B / (2 * pi) * dtheta * Th) ./ q_inf # Eq. 13

    if output_all
        return Th,
        Q,
        Rp,
        Tp,
        Zp,
        Vloc,
        CD,
        CT,
        alpha,
        cl,
        cd_af,
        Re,
        M_addedmass_Np,
        M_addedmass_Tp,
        F_addedmass_Np,
        F_addedmass_Tp,
        F_buoy,
        cm_af,
        M25
    else
        return CD - CT # Residual, section 2.4
    end
end

"""
DMS(turbine, env; w=0, idx_RPI=1:turbine.ntheta, solve=true)

see ?steady for detailed i/o description

Double multiple streamtube model.

`finite_span_factor` may be a nonnegative scalar or length-`ntheta` vector. It
defaults to `1.0` and scales aerodynamic blade loads, torque, thrust, induction
source terms, and aerodynamic moment without scaling added mass, buoyancy, or
centrifugal terms.
"""
function DMS(
    turbine,
    env;
    w = 0,
    idx_RPI = 1:turbine.ntheta,
    solve = true,
    finite_span_factor = 1.0,
)
    #TODO: Inputs documentation
    a_in = w
    ntheta = turbine.ntheta
    #Convert global F.O.R. winds to turbine frame
    windangle = env.windangle

    V_xtemp = env.V_x .* cos(windangle) + env.V_y .* sin(windangle) # Vinf is V_x, t is for turbine direction f.o.r.
    V_ytemp = -env.V_x .* sin(windangle) + env.V_y .* cos(windangle)
    env_turbine = Environment(
        env.rho,
        env.mu,
        V_xtemp,
        V_ytemp,
        env.V_z,
        env.V_twist,
        env.windangle,
        env.DynamicStallModel,
        env.AeroModel,
        env.Aero_AddedMass_Active,
        env.Aero_Buoyancy_Active,
        env.Aero_RotAccel_Active,
        env.AddedMass_Coeff_Ca,
        env.centrifugal_force_flag,
        env.aw_warm,
        env.steplast,
        env.idx_RPI,
        env.V_wake_old,
        env.BV_DynamicFlagL,
        env.BV_DynamicFlagD,
        env.alpha_last,
        env.suction,
        env.accel_flap,
        env.accel_edge,
        env.gravity,
    )

    # Unpack the rest

    dtheta = 2 * pi / (ntheta)
    finite_span_factor = _validated_finite_span_factor(finite_span_factor, ntheta)
    # thetavec = collect(dtheta/2:dtheta:2*pi-dtheta/2)
    thetavec = collect((dtheta/2):dtheta:(2*pi))

    astar = zeros(Real, ntheta * 2)#env.aw_warm[1:ntheta] #zeros(ntheta)

    if a_in == 0
        idx_RPI = 1:ntheta
    else
        astar[1:ntheta] = a_in[1:ntheta]
    end
    Th = zeros(Real, ntheta)
    Q = zeros(Real, ntheta)
    Rp = zeros(Real, ntheta)
    Tp = zeros(Real, ntheta)
    Zp = zeros(Real, ntheta)
    Vloc = zeros(Real, ntheta)
    CD = zeros(Real, ntheta)
    CT = zeros(Real, ntheta)
    alpha = zeros(Real, ntheta)
    cl = zeros(Real, ntheta)
    cd_af = zeros(Real, ntheta)
    cm_af = zeros(Real, ntheta)
    M25 = zeros(Real, ntheta)
    Re = zeros(Real, ntheta)
    M_addedmass_Np = zeros(Real, ntheta)
    M_addedmass_Tp = zeros(Real, ntheta)
    F_addedmass_Np = zeros(Real, ntheta)
    F_addedmass_Tp = zeros(Real, ntheta)
    F_buoy = zeros(Real, ntheta, 3)

    # For All Upper
    iter_RPI = 1
    for i = 1:Int(ntheta/2)
        # Solve Residual
        if idx_RPI[iter_RPI] == i
            iter_RPI += 1
            if solve
                resid(a) = streamtube(
                    a,
                    thetavec[i],
                    turbine,
                    env_turbine;
                    solvestep = true,
                    finite_span_factor,
                ) #solvestep makes this independent solve not mess up the dynamic stall variables during the solve
                astar[i], _ = FLOWMath.brent(resid, 0, 0.999)
            end

            Th[i],
            Q[i],
            Rp[i],
            Tp[i],
            Zp[i],
            Vloc[i],
            CD[i],
            CT[i],
            alpha[i],
            cl[i],
            cd_af[i],
            Re[i],
            M_addedmass_Np[i],
            M_addedmass_Tp[i],
            F_addedmass_Np[i],
            F_addedmass_Tp[i],
            F_buoy[i, :],
            cm_af[i],
            M25[i] = streamtube(
                astar[i],
                thetavec[i],
                turbine,
                env_turbine;
                output_all = true,
                finite_span_factor,
            )
        end
    end

    Vxsave = mean(env_turbine.V_x)

    # For All Lower
    for i = Int(ntheta/2+1):ntheta
        # Update Downstream Inflow based on upstream solution
        if idx_RPI[iter_RPI] == i
            iter_RPI += 1
            a_used = astar[Int(ntheta/2)-(i-Int(ntheta/2))+1] # We index in a circle, but the streamtubes go straight through, swapping at 180 deg
            if env_turbine.suction
                Vxwake = (1.0 - a_used) / (1.0 + a_used) * env_turbine.V_x[i] # Eq. 6
            else
                Vxwake = (1.0 - 2.0 * a_used) * env_turbine.V_x[i] # Eq. 3
            end

            # Solve Residual
            if solve
                resid(a) = streamtube(
                    a,
                    thetavec[i],
                    turbine,
                    env_turbine;
                    Vxwake,
                    solvestep = true,
                    finite_span_factor,
                ) #Solve wrt negative theta on the back side?
                astar[i], _ = FLOWMath.brent(resid, 0, 0.999)
            end

            Th[i],
            Q[i],
            Rp[i],
            Tp[i],
            Zp[i],
            Vloc[i],
            CD[i],
            CT[i],
            alpha[i],
            cl[i],
            cd_af[i],
            Re[i],
            M_addedmass_Np[i],
            M_addedmass_Tp[i],
            F_addedmass_Np[i],
            F_addedmass_Tp[i],
            F_buoy[i, :],
            cm_af[i],
            M25[i] = streamtube(
                astar[i],
                thetavec[i],
                turbine,
                env_turbine;
                output_all = true,
                Vxwake,
                finite_span_factor,
            )
        end
    end

    # Aggregate Turbine Performance
    k = 1.0
    CP = sum(
        (k * turbine.B / (2 * pi) * dtheta * Q * mean(abs.(turbine.omega))) /
        (0.5 * env_turbine.rho * 1.0 * 2 * turbine.R * Vxsave^3),
    ) # Eq. 14, normalized by nominal radius R

    return CP,
    Th,
    Q,
    Rp,
    Tp,
    Zp,
    Vloc,
    CD,
    CT,
    mean(astar[1:ntheta]),
    astar,
    alpha,
    cl,
    cd_af,
    thetavec .- windangle,
    Re,
    M_addedmass_Np,
    M_addedmass_Tp,
    F_addedmass_Np,
    F_addedmass_Tp,
    F_buoy,
    cm_af,
    M25
end
