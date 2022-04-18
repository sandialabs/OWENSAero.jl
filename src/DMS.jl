# This module is based on https://doi.org/10.5194/wes-2019-44 A Double Multiple
# Streamtube model for Vertical Axis Wind Turbines of arbitrary rotor loading
# The equation references align to this report
# Some equations also taken from (are explicitely referenced) doi:10.5194/wes-1-327-2016 Actuator cylinder theory for multiple vertical axis wind turbines

import FLOWMath
import Statistics
import Statistics:mean

"""
INTERNAL streamtube(a,theta,turbine,env;output_all=false,Vxwake=nothing,solvestep=false)

Double multiple streamtube individual streamtube calculation

Output:

if output_all
    return Th, Q, Rp, Tp, Zp, Vloc, CD, CT, alpha, cl, cd_af, Re
else
    return CD-CT # Residual, section 2.4
end

"""
function streamtube(a,theta,turbine,env;output_all=false,Vxwake=nothing,solvestep=false)

    # Unpack Vars
    B = turbine.B
    k = 1.0
    af = turbine.af
    chord = turbine.chord[1]
    ntheta = turbine.ntheta
    suction = false
    rho = env.rho
    mu = env.mu
    # winddir = env.winddir

    dtheta = 2*pi/(ntheta) #Assuming discretization is fixed equidistant (but omega can change between each point)

    if length(turbine.r)>1
        idx = round(Int, (theta+dtheta/2)/dtheta)
        r = turbine.r[idx]
        twist = turbine.twist[idx]
        delta = turbine.delta[idx]
        omega = turbine.omega[idx]
        rotation = sign(turbine.omega[idx])
        V_xt = env.V_x[idx] # Vinf is V_xt, t is for turbine direction f.o.r.
        V_yt = env.V_y[idx]
        V_zt = env.V_z[idx]
        V_twist = env.V_twist[idx]
        # V_delta = env.V_delta[idx] # Does not apply since the model calculation is centered around the point of rotation
        # V_sweep = env.V_sweep[idx] # Does not apply since the model calculation is centered around the point of rotation

    else
        idx = 1
        r = turbine.r
        twist = turbine.twist
        delta = turbine.delta
        omega = turbine.omega
        rotation = sign(mean(turbine.omega))
        V_x = env.V_x[1]
        V_y = env.V_y[1]
        V_z = env.V_z[1]
        V_twist = env.V_twist[1]
        # V_delta = env.V_delta[1] # Does not apply since the model calculation is centered around the point of rotation
        # V_sweep = env.V_sweep[1] # Does not apply since the model calculation is centered around the point of rotation
    end

    if Vxwake!=nothing
        V_xt = Vxwake
    end

    if suction
        # Drag Coefficient: Eq. 5
        if a<=0.7
            CD = 4/3 * a * (3-a)/(1+a)
        elseif a > 0.7 && a <=1
            CD = 0.889 - (0.0203-(a-0.143)^2)/0.6427
        else
            error("a > 1 is undefined")
        end
    else
        # Drag Coefficient: Eq. 2
        if a<=0.4
            CD = 4*a*(1-a)
        else
            CD = 0.889 - (0.0203-(a-0.143)^2)/0.6427
        end
    end

    # Cycle-average Thrust Coefficient
    # TSR = abs(omega)*r/sqrt(V_x^2+V_y^2) # Tip Speed Ratio, section 2.3, however since this is being used for relative local velocities, we use the local radius r, and not the nominal radius R
    # Vt = V_x*(1-a)*cos(theta)+V_y*sin(theta)+V_y*(a)*cos(theta)+TSR

    Vn = (V_xt*(1-a)*sin(theta)-V_yt*cos(theta)+V_yt*(a)*sin(theta))*cos.(delta) + V_zt*sin(delta)
    Vt = V_xt*((1-a)*cos(theta))+V_yt*sin(theta)+V_yt*(1-a)*cos(theta)+abs(omega)*r
    Vloc = sqrt(Vn^2 + Vt^2)
    phi = atan(Vn, Vt)
    alpha = phi - twist

    # Original Equations: I broke them out above to be more readable and to know where to put the cos(delta)
    # phi = atan(((1-a)*sin(theta)).*cos.(delta) , ((1-a)*cos(theta)+TSR)) # Eq. 7 but correct for slope from the 3D to 2D frame of reverence, also must use atan2 to determine which quadrant and get the signs correct
    # Vloc = Vinf .* ((1-a)^2 + 2*(1-a)*TSR*cos(theta).*cos.(delta) + TSR^2).^0.5 # Eq. 8

    Re = rho*Vloc*chord/mu
    dt = dtheta/abs.(omega)
    v_sound = 343.0 #m/s #TODO: calculate this using Atmosphere.jl
    mach = Vloc/v_sound
    if env.DSModel == "BV"
        cl, cd_af = af(alpha,Re,mach,env,V_twist,chord,dt,Vloc;solvestep,idx)
    elseif env.DSModel == "LB"
        error("LB Dynamic Stall Model Not Implemented Yet")
    else
        cl, cd_af = af(alpha,Re,mach)
    end
    ct = cd_af*cos(phi) - cl*sin(phi) # Eq. 9
    cn = cd_af*sin(phi) + cl*cos(phi) # Eq. 10

    Ab = chord * 1.0 # Assuming unit section height

    # Instantaneous Forces (Unit Height) #Based on this, radial is inward and tangential is in direction of rotation
    q_loc = 0.5 * rho * Ab * Vloc^2 # From Eq. 11
    Rp = cn.*q_loc # Ning Eq. 27 # Negate to match cactus frame of reference
    Tp = -rotation*ct.*q_loc/cos(delta) # Ning Eq. 27 # Negate to match cactus frame of reference
    Zp = cn.*q_loc*tan(delta) # Ning Eq. 27 # Negate to match cactus frame of reference

    Th = q_loc * (ct*cos(theta) + cn*sin(theta)/cos(delta)) # Eq. 11 but with delta correction

    Q = q_loc * r * -ct # Eq. 12 but with Local radius for local torque, Negate the force for reaction torque, in the power frame of reference?

    Ast = 1.0 * r * dtheta *abs(sin(theta)) # Section 2.0, local radius for local area, however do not allow negative area, so use absolute value

    q_inf = 0.5 * rho * Ast * (V_xt^2+V_yt^2) # this is in the radial frame of reference V_x.^2 , also (sqrt(V_x^2+V_y^2)).^2 is reduced

    CT = (k * B/(2*pi) * dtheta*Th) ./ q_inf # Eq. 13

    if output_all
        return Th, Q, Rp, Tp, Zp, Vloc, CD, CT, alpha, cl, cd_af, Re
    else
        return CD-CT # Residual, section 2.4
    end
end

"""
DMS(turbine, env; w=0, idx_RPI=1:turbine.ntheta, solve=true)

see ?steady for detailed i/o description

Double multiple streamtube model
"""
function DMS(turbine, env; w=0, idx_RPI=1:turbine.ntheta, solve=true)
    #TODO: Inputs documentation
    a_in = w
    ntheta = turbine.ntheta
    #Convert global F.O.R. winds to turbine frame
    windangle = env.windangle

    V_xtemp = env.V_x.*cos(windangle)+env.V_y.*sin(windangle) # Vinf is V_x, t is for turbine direction f.o.r.
    V_ytemp = -env.V_x.*sin(windangle)+env.V_y.*cos(windangle)

    env.V_x[:] = V_xtemp #TODO: ensure this doesn't mess up unsteady methods with persistent backend memory
    env.V_y[:] = V_ytemp

    # Unpack the rest

    dtheta = 2*pi/(ntheta)
    # thetavec = collect(dtheta/2:dtheta:2*pi-dtheta/2)
    thetavec = collect(dtheta/2:dtheta:2*pi)

    astar = zeros(Real,ntheta*2)#env.aw_warm[1:ntheta] #zeros(ntheta)

    if a_in == 0
        idx_RPI = 1:ntheta
    else
        astar[1:ntheta] = a_in[1:ntheta]
    end
    Th = zeros(Real,ntheta)
    Q = zeros(Real,ntheta)
    Rp = zeros(Real,ntheta)
    Tp = zeros(Real,ntheta)
    Zp = zeros(Real,ntheta)
    Vloc = zeros(Real,ntheta)
    CD = zeros(Real,ntheta)
    CT = zeros(Real,ntheta)
    alpha = zeros(Real,ntheta)
    cl = zeros(Real,ntheta)
    cd_af = zeros(Real,ntheta)
    Re = zeros(Real,ntheta)

    # For All Upper
    iter_RPI = 1
    for i = 1:Int(ntheta/2)
        # Solve Residual
        if idx_RPI[iter_RPI] == i
            iter_RPI+=1
            if solve
                resid(a) = streamtube(a,thetavec[i],turbine,env;solvestep=true) #solvestep makes this independent solve not mess up the dynamic stall variables during the solve
                astar[i], _ = FLOWMath.brent(resid, 0, .999)
            end

            Th[i], Q[i], Rp[i], Tp[i], Zp[i], Vloc[i], CD[i], CT[i], alpha[i], cl[i], cd_af[i], Re[i]  = streamtube(astar[i],thetavec[i],turbine,env;output_all=true)
        end
    end

    Vxsave = mean(env.V_x)

    # For All Lower
    for i = Int(ntheta/2+1):ntheta
        # Update Downstream Inflow based on upstream solution
        if idx_RPI[iter_RPI] == i
            iter_RPI+=1
            a_used = astar[Int(ntheta/2)-(i-Int(ntheta/2))+1] # We index in a circle, but the streamtubes go straight through, swapping at 180 deg
            if env.suction
                Vxwake = (1.0-a_used)/(1.0+a_used)*env.V_x[i] # Eq. 6
            else
                Vxwake = (1.0-2.0*a_used)*env.V_x[i] # Eq. 3
            end

            # Solve Residual
            if solve
                resid(a) = streamtube(a,thetavec[i],turbine,env;Vxwake,solvestep=true) #Solve wrt negative theta on the back side?
                astar[i], _ = FLOWMath.brent(resid, 0, .999)
            end

            Th[i], Q[i], Rp[i], Tp[i], Zp[i], Vloc[i], CD[i], CT[i], alpha[i], cl[i], cd_af[i], Re[i]  = streamtube(astar[i],thetavec[i],turbine,env;output_all=true,Vxwake)
        end
    end

    # Aggregate Turbine Performance
    k = 1.0
    CP = sum((k * turbine.B/(2*pi) * dtheta*abs.(Q)*mean(abs.(turbine.omega))) / (0.5*env.rho*1.0*2*turbine.R*Vxsave^3)) # Eq. 14, normalized by nominal radius R

    return CP, Th, Q, Rp, Tp, Zp, Vloc, CD, CT, mean(astar[1:ntheta]), astar, alpha, cl, cd_af, thetavec.-windangle, Re
end
