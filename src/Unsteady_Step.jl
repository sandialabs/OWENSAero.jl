# Common RPI/Filtered Unsteady Method
"""
    Unsteady_Step(turbine,env,us_param,mystep)

calls inflow wind init

# Inputs
* `turbine::Turbine`: turbine input for slice see ?Turbine
* `env::Env`: environment input for slice see ?Env
* `us_param::UnsteadyParams`: unsteady inputs for slice see ?UnsteadyParams
* `mystep::int`: continuous index cooresponding to the azimuthal discretation - i.e. for ntheta of 30 step 1 is the first step of rev 1, sep 31 is the first step of rev 2, etc.  Keeps track of temporal locaion

# Outputs:
* `CP`: This slice's coefficient of performance at this step
* `Th`: This slice's thrust coefficient at this step
* `Q`: Torque (N0m) at this step
* `Rp`: Radial force per height (N) at this step
* `Tp`: Tangential force per height (N) at this step
* `Zp`: Vertical force per height (N) at this step
* `Vloc`: Local velocity array for each azimuthal position (includes induction) (m/s) at this step
* `CD`: This slice's drag coefficient at this step
* `CT`: This slice's thrust coefficient (should equal drag, but may no depending on usage or solver status) at this step
* `amean`: Mean turbine induction in the streamwise direction at this step
* `astar`: Solved induction factors for each azimuthal location. First half are streamwise (u), second are cross-steam (v) at this step
* `alpha`: Local angle of attack array for each azimuthal position (includes induction) (rad) at this step
* `cl`: Local lift coefficient used for each azimuthal position at this step
* `cd_af`: Local drag coefficient used for each azimuthal position at this step
* `thetavec`: Azimuthal location of each discretization (rad)
* `Re`: Reynolds number for each azimuthal position at this step

"""
function Unsteady_Step(turbine,env,us_param,mystep)
    #TODO: RPI run with negative omega has a hard time converging, fix this
    # Unpack
    RPI = us_param.RPI
    ifw = us_param.ifw
    ntheta = turbine.ntheta
    dtheta = 2*pi/(ntheta)
    thetavec = collect(dtheta/2:dtheta:2*pi-dtheta/2)
    AeroModel = env.AeroModel
    rotation = sign(mean(turbine.omega))

    # Check that ntheta is divisible by nblades and 2
    if ntheta%turbine.B!=0 || ntheta%2!=0
        error("ntheta must be wholly divisible by the number of blades and by 2, e.g for 3 blades, the options would be 3*even numbers, or 6,12,18,24,30...")
    end

    tau = us_param.tau
    R = turbine.R
    Vinf = env.V_x
    aw_warm = zeros(Real,ntheta*2)
    aw_warm[:] = env.aw_warm[:] #Mutate, do not link
    V_wake_old = zeros(Real,1,1)
    V_wake_old[1] = max(0,env.V_wake_old[1]) # TODO: smooth max? Prevent numerical issues associated with negative wake

    # setup
    awnew = zeros(Real,size(aw_warm))

    # Calculate the RPI indices
    circular_mystep = mystep-floor((mystep-1)/ntheta*turbine.B)*ntheta/turbine.B

    if RPI
        idx_RPI = Int(circular_mystep):Int(ntheta./turbine.B):Int(ntheta*2-ntheta/turbine.B+1+circular_mystep)
    else
        idx_RPI = 1:ntheta*2
    end
    idx_RPI_half = idx_RPI[1:Int(length(idx_RPI)/2)]

    # Non-steady wind
    dt = 1.0 ./(abs.(turbine.omega) * float(ntheta) / (2*pi))

    if us_param.IECgust || ifw
        ele_x = rotation*sin.((1:float(ntheta))/float(ntheta)*2*pi).*turbine.r./R
        ele_y = rotation*cos.((1:float(ntheta))/float(ntheta)*2*pi).*turbine.r./R
        for i_theta = 1:ntheta
            #TODO: other simpler turbulence models, that may be made continuous

            if us_param.IECgust
                dt_norm = dt[i_theta]*us_param.nominalVinf/R
                gustT = us_param.gustT * us_param.nominalVinf / R
                tr = mystep*dt_norm .- ele_x .- us_param.gustX0
                if (tr[i_theta] >= 0) && (tr[i_theta]<=gustT)
                    IECGustFactor = 1.0 - 0.37 * us_param.G_amp/us_param.nominalVinf * sin(3*pi*tr[i_theta]/gustT)  * (1.0 - cos(2*pi*tr[i_theta]/gustT))
                    env.V_x[i_theta] = us_param.nominalVinf*IECGustFactor
                else
                    env.V_x[i_theta] = us_param.nominalVinf
                end
            end

            if ifw
                time = mystep*dt[i_theta] #TODO: make the input time instead of step and solve up to a given time considering a specified timestep
                velocity = OWENSOpenFASTWrappers.ifwcalcoutput([ele_x[i_theta],ele_y[i_theta],turbine.z],time)
                env.V_x[i_theta] = velocity[1]
                env.V_y[i_theta] = velocity[2]
                env.V_z[i_theta] = velocity[3]
            end
        end
    end

    # Solve for the instantaneous induction factors and calculate new unfiltered induction factors and wake
    awnew[:] = aw_warm[:] #Mutate, do not link

    _,_,_,_,_,_,_,_,_, a, w = steady(turbine, env; w=awnew, idx_RPI,ifw)
    awnew[:] = w[:]

    tau_near = tau[1]*R / V_wake_old[1]
    tau_far = tau[2]*R / V_wake_old[1]

    V_wake_filt = V_wake_old[1] * exp(-mean(dt) / tau_far) + (mean(env.V_x) * (1.0 - 2.0 * a)) * (1.0 - exp(-mean(dt) / tau_far))

    aw_filtered = aw_warm[:] .* exp.(-[dt;dt] ./ tau_near) .+ awnew[:] .* (1.0 .- exp.(-[dt;dt] ./ tau_near))

    # Get the actual performance using the filtered induction factors
    CP, Th, Q_temp, Rp_temp, Tp_temp, Zp_temp, Vloc, CD, CT, a, awstar, alpha, cl, cd, thetavec, Re, M_addedmass_Np, M_addedmass_Tp, F_addedmass_Np, F_addedmass_Tp = steady(turbine, env; w=aw_filtered,idx_RPI,solve=false,ifw)

    if env.steplast[1] != mystep
        env.aw_warm[:] = aw_filtered[:] #Mutate, do not link
        env.V_wake_old[1] = Main.ForwardDiff.value(V_wake_filt[1])
        env.steplast[1] = mystep
    end
    # env.idx_RPI[:] = idx_RPI_half #Mutate, do not link, TODO: fix this, is used for indexing outputs
    theta_temp = thetavec[idx_RPI_half[1]]

    return CP,Th,Q_temp, Rp_temp, Tp_temp, Zp_temp, Vloc, CD, CT, a, awstar, alpha, cl, cd, thetavec, Re, M_addedmass_Np, M_addedmass_Tp, F_addedmass_Np, F_addedmass_Tp
end
