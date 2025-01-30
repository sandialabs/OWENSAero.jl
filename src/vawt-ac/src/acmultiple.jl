import NLsolve
import LsqFit
import HDF5
import QuadGK
import Statistics

# --- Influence Coefficients ---

"""
applies for both Ay and Rx depending on which function ifunc(x, y, phi)
is passed in
"""
function panelIntegration(deltavec,rvec,xvec, yvec, thetavec, ifunc)

    # initialize
    nx = length(xvec)
    ntheta = length(thetavec)
    dtheta = thetavec[2] - thetavec[1]  # assumes equally spaced
    A = zeros(nx, ntheta)

    for i in eachindex(xvec)
        # redefine function so it has one parameter for use in quadgk
        integrand(phi) = ifunc(deltavec[i],rvec[i],xvec[i], yvec[i], phi)

        for j in eachindex(thetavec)
            # an Adaptive Gauss-Kronrod quadrature integration.  Tried trapz but this was faster.
            A[i, j], error = QuadGK.quadgk(integrand, thetavec[j]-dtheta/2.0, thetavec[j]+dtheta/2.0, atol=1e-10)
        end

    end

    return A
end


"""
integrand used for computing Dx
"""
function Dxintegrand(deltan,rn, x, y, phi)
    v1 = x + rn*sin(phi)
    v2 = y - rn*cos(phi)
    # v1 and v2 should never both be zero b.c. we never integrate self.  RxII handles that case.
    return (v1*sin(phi)*cos(deltan) - v2*cos(phi)*cos(deltan))*rn/(2*pi*(v1*v1 + v2*v2))
end


"""
integrand used for computing Ay
"""
function Ayintegrand(deltan,rn, x, y, phi)
    v1 = x + rn*sin(phi)
    v2 = y - rn*cos(phi)
    if abs(v1) < 1e-12 && abs(v2) < 1e-12  # occurs when integrating self, function symmetric around singularity, should integrate to zero
        return 0.0
    end
    return (v1*cos(phi)*cos(deltan) + v2*sin(phi)*cos(deltan))*rn/(2*pi*(v1*v1 + v2*v2))
end

"""
integrand used for computing AIJ
"""
function AyIJ(deltavec,rvec,xvec, yvec, thetavec)
    return panelIntegration(deltavec,rvec,xvec, yvec, thetavec, Ayintegrand)
end

"""
integrand used for computing DxIJ
"""
function DxIJ(deltavec,rvec,xvec, yvec, thetavec)
    return panelIntegration(deltavec,rvec,xvec, yvec, thetavec, Dxintegrand)
end

"""
integrand used for computing WxIJ
"""
function WxIJ(xvec, yvec, thetavec)

    # initialize
    nx = length(xvec)
    ntheta = length(thetavec)
    dtheta = thetavec[2] - thetavec[1]  # assumes equally spaced
    Wx = zeros(Real,nx, ntheta)

    for i in eachindex(xvec)
        if yvec[i] >= -1.0 && yvec[i] <= 1.0 && xvec[i] >= 0.0 && xvec[i]^2 + yvec[i]^2 >= 1.0
            # if yvec[i] >= -1.0 && yvec[i] <= 1.0 && (xvec[i] >= 0.0 || (xvec[i] >= -1 && xvec[i]^2 + yvec[i]^2 <= 1.0))
            thetak = acos(yvec[i])
            k = findfirst(thetavec + dtheta/2 .> thetak)  # index of intersection
            Wx[i, k] = -1.0
            Wx[i, ntheta-k+1] = 1.0
        end
    end

    return Wx
end

"""
integrand used for computing DxII
"""
function DxII(thetavec)

    # initialize
    ntheta = length(thetavec)
    dtheta = thetavec[2] - thetavec[1]  # assumes equally spaced
    Rx = dtheta/(4*pi)*ones(ntheta, ntheta)

    for i in eachindex(thetavec)
        if i <= ntheta/2
            Rx[i, i] = (-1 + 1.0/ntheta)/2.0
        else
            Rx[i, i] = (1 + 1.0/ntheta)/2.0
        end
    end

    return Rx
end

"""
integrand used for computing WxII
"""
function WxII(thetavec)

    # initialize
    ntheta = length(thetavec)
    Wx = zeros(ntheta, ntheta)

    for i = div(ntheta,2)+1:ntheta
        Wx[i, ntheta+1-i] = -1
    end

    return Wx
end

"""
Internal, precomputes influence coefficient matricies and saves them as HDF5 files
"""
function precomputeMatrices(deltavec,rvec,ntheta,file)

    # precompute self influence matrices

    # setup discretization (all the same, and uniformly spaced in theta)
    dtheta = 2*pi/ntheta
    theta = collect(dtheta/2:dtheta:2*pi)

    Dxself = DxII(theta)
    Wxself = WxII(theta)
    Ayself = AyIJ(deltavec,rvec,-sin.(theta), cos.(theta), theta)

    # write to file
    HDF5.h5open(file, "w") do file
        HDF5.write(file, "theta", theta)
        HDF5.write(file, "Dx", Dxself)
        HDF5.write(file, "Wx", Wxself)
        HDF5.write(file, "Ay", Ayself)
    end

end

"""
Internal, assembles the matrices of multiple turbine systems into a combined system
centerX, centerY: array of x,y coordinates for centers of the VAWTs in the farm
radii: corresponding array of their radii
"""
function matrixAssemble(turbine, centerX, centerY, radii, ntheta)

    if length(turbine.delta) == 1
        deltavec = ones(Real,ntheta)*turbine.delta
        rvec = ones(Real,ntheta)*turbine.r./mean(turbine.r)
    else
        deltavec = turbine.delta
        rvec = turbine.r./mean(turbine.r)
    end

    if turbine.r_delta_influence
        file = "theta-r-delta-$ntheta.h5"
    else
        file = "theta-$ntheta.h5"
        deltavec = zeros(Real,length(deltavec))
        rvec = ones(Real,length(rvec))
    end

    if !isfile(file) || turbine.r_delta_influence
        precomputeMatrices(deltavec,rvec,ntheta,file)
    end

    theta = HDF5.h5read(file, "theta")
    Dxself = HDF5.h5read(file, "Dx")
    Wxself = HDF5.h5read(file, "Wx")
    Ayself = HDF5.h5read(file, "Ay")
    # for row = 1:length(Ayself[:,1])
    #     println((Ayself[row,:]))
    # end
    # initialize global matrices
    nturbines = length(radii)
    Dx = zeros(Real,nturbines*ntheta, nturbines*ntheta)
    Wx = zeros(Real,nturbines*ntheta, nturbines*ntheta)
    Ay = zeros(Real,nturbines*ntheta, nturbines*ntheta)

    # iterate through turbines
    for I in eachindex(radii)
        for J in eachindex(radii)

            # find normalized i locations relative to center of turbine J
            x = (centerX[I].-radii[I]*sin.(theta) .- centerX[J])/radii[J]
            y = (centerY[I].+radii[I]*cos.(theta) .- centerY[J])/radii[J]

            # self-influence is precomputed
            if I == J
                Dxsub = Dxself
                Wxsub = Wxself
                Aysub = Ayself

                # pairs can be mapped for same radius
            elseif J < I && radii[I] == radii[J]

                # grab cross-diagonal I,J -> J,I matrix
                Dxsub = Dx[(J-1)*ntheta+1:J*ntheta, (I-1)*ntheta+1:I*ntheta]
                Aysub = Ay[(J-1)*ntheta+1:J*ntheta, (I-1)*ntheta+1:I*ntheta]

                # mapping index for coefficients that are the same
                idx = [div(ntheta,2)+1:ntheta; 1:div(ntheta,2)]

                # directly map over
                Dxsub = Dxsub[idx, idx]
                Aysub = Aysub[idx, idx]

                # wake term must be recomptued
                Wxsub = WxIJ(x, y, theta)

                # # if VAWTs are very far apart we can approximate some of the influence coefficients
                # elseif approxfar && sqrt((centerX[I]-centerX[J])^2 + (centerY[I]-centerY[J])^2) > 10*radii[I]
                #     println("far apart")
                #     xc = (centerX[I] - centerX[J])/radii[J]
                #     yc = (centerY[I] - centerY[J])/radii[J]

                #     Rxsub = RxIJFar(xc, yc, theta)
                #     Wxsub = zeros(ntheta, ntheta)  # should have negligible wake impact
                #     Aysub = AyIJFar(xc, yc, theta)

            else
                Dxsub = DxIJ(deltavec,rvec, x, y, theta)
                Wxsub = WxIJ(x, y, theta)
                Aysub = AyIJ(deltavec,rvec, x, y, theta)
            end

            # assemble into global matrix
            Dx[(I-1)*ntheta+1:I*ntheta, (J-1)*ntheta+1:J*ntheta] = Dxsub
            Wx[(I-1)*ntheta+1:I*ntheta, (J-1)*ntheta+1:J*ntheta] = Wxsub
            Ay[(I-1)*ntheta+1:I*ntheta, (J-1)*ntheta+1:J*ntheta] = Aysub

        end
    end
    Ax = Dx + Wx

    return Ax, Ay, theta
end
# ---------------------------

# ------- Force Coefficients ---------
"""
Internal, calculates the radial force used in the residual function as well as the turbine performance when converged
"""
function radialforce(uvec, vvec, thetavec, turbine, env)

    # unpack
    ntheta = turbine.ntheta
    R = turbine.R #reference radius
    r = turbine.r #deformed (potentially) radius at each azimuthal location
    chord = turbine.chord#[1]
    thickness = turbine.thick#[idx]
    twist = turbine.twist
    delta = turbine.delta
    B = turbine.B
    Omega = turbine.omega

    rho = env.rho
    mu = env.mu
    V_x = env.V_x
    V_y = env.V_y
    V_wind = sqrt.(V_x.^2 + V_y.^2)
    V_z = env.V_z
    V_twist = env.V_twist # in rad/s
    # V_delta = env.V_delta # Does not apply since the model calculation is centered around the point of rotation
    # V_sweep = env.V_sweep # Same as V_delta

    # set the rotation direction
    rotation = sign(mean(Omega))

    # velocity components and angles
    Vn = (V_x.*(1.0 .+ uvec).*sin.(thetavec) .- V_x.*vvec.*cos.(thetavec) .-
        V_y.*(1.0 .+ vvec).*cos.(thetavec) .+ V_y.*(vvec).*sin.(thetavec)
        ).*cos.(delta) .+ V_z.*sin.(delta) #multiply by slope to get 3D to 2D transformation normal to radial
    Vt = (V_x.*(1.0 .+ uvec).*cos.(thetavec) .+ V_x.*vvec.*sin.(thetavec) .+
        V_y.*(1.0 .+ vvec).*sin.(thetavec) .+ V_y.*uvec.*cos.(thetavec)) .+ abs.(Omega).*r
    W = sqrt.(Vn.^2 + Vt.^2)
    phi = atan.(Vn, Vt)
    alpha = phi .- twist
    Re = rho*W*chord/mu  # currently no Re dependence

    # airfoil
    dtheta = 2*pi/ntheta
    dt = dtheta./abs.(Omega)
    v_sound = 343.0 #m/s #TODO: calculate this using Atmosphere.jl
    mach = W/v_sound
    if env.DynamicStallModel == "BV"
        cl = zeros(Real,length(alpha))
        cd = zeros(Real,length(alpha))
        for ii = 1:length(alpha)
            cl[ii], cd[ii] = turbine.af(alpha[ii],Re[ii],mach[ii],env,V_twist[ii],chord,dt[ii],W[ii])
        end
    elseif env.DynamicStallModel == "LB"
        error("LB Dynamic Stall Model Not Coupled Yet")
    else
        cl, cd = turbine.af(alpha,Re,mach)
    end

    # rotate force coefficients
    cn = cl.*cos.(phi) + cd.*sin.(phi)
    ct = cl.*sin.(phi) - cd.*cos.(phi)

    # radial force
    sigma = B*chord./r
    q = sigma./(4*pi).*cn./cos.(delta).*(W./V_wind).^2 #divide by slope to get 3D to 2D transformation normal to radial

    # Added Mass
    if env.Aero_AddedMass_Active
        if length(turbine.r) > 1  # TODO CM: theta/thetavec
            # dtheta = 2 * pi / (ntheta) #Assuming discretization is fixed equidistant (but omega can change between each point)
            # println("theta: $theta")
            # println("dtheta: $dtheta")
            # println("ntheta: $ntheta")
            # idx = round(Int, (theta + dtheta / 2) / dtheta)

            twist = turbine.twist#[idx]
            omega = turbine.omega#[idx]
            r = turbine.r#[idx]
            accel_flap = env.accel_flap#[idx]
            accel_edge = env.accel_edge#[idx]
        else # TODO CM: is this ever used?
            # twist = turbine.twist
            # omega = turbine.omega
            # r = turbine.r
            # accel_flap = env.accel_flap[1]
            # accel_edge = env.accel_edge[1]
        end
        chord = turbine.chord#[1]
        thickness = turbine.thick#[idx]
        rho = env.rho

        Vol_flap = @. pi * (chord / 2)^2 * 1.0
        Vol_edge = @. pi * ((thickness / 10) / 2)^2 * 1.0

        if env.Aero_RotAccel_Active
            accel_rot = @. omega^2 * r
        else
            accel_rot = 0.0
        end

        M_addedmass_flap = @. rho * env.AddedMass_Coeff_Ca * Vol_flap
        M_addedmass_edge = @. rho * env.AddedMass_Coeff_Ca * Vol_edge

        F_addedmass_flap = @. M_addedmass_flap * (accel_flap + accel_rot)
        F_addedmass_edge = @. M_addedmass_edge * (accel_edge + accel_rot)

        M_addedmass_Np = @. M_addedmass_flap * cos(twist) + M_addedmass_edge * sin(twist) # Go from the beam frame of reference to the normal and tangential direction #TODO: verify the directions
        M_addedmass_Tp = @. M_addedmass_edge * cos(twist) - M_addedmass_flap * sin(twist)

        F_addedmass_Np = @. F_addedmass_flap * cos(twist) + F_addedmass_edge * sin(twist) # Go from the beam frame of reference to the normal and tangential direction #TODO: verify the directions
        F_addedmass_Tp = @. F_addedmass_edge * cos(twist) - F_addedmass_flap * sin(twist)
    else
        M_addedmass_Np = zero(alpha)
        M_addedmass_Tp = zero(alpha)
        F_addedmass_Np = zero(alpha)
        F_addedmass_Tp = zero(alpha)
    end

    # Buoyancy
    F_buoy = zeros(3, ntheta)
    section_area = chord * thickness / 2 * 1.0 # per unit length TODO: input volume
    mass = -env.gravity .* (rho * section_area .- turbine.rhoA) # buoyancy mass minus structural mass since added mass requires moving the gravity here
    if env.Aero_Buoyancy_Active
        for (itheta, theta) in enumerate(thetavec)
            dcm = [
                cos(theta) -sin(theta) 0
                sin(theta) cos(theta) 0
                0 0 1
            ]
            F_buoy[:, itheta] = dcm * mass
        end
    end

    # instantaneous forces #Based on this, radial is inward and tangential is in direction of rotation
    qdyn = 0.5*rho*W.^2
    Rp = cn.*qdyn*chord - F_addedmass_Np .+ F_buoy[2, :] # TODO CM: correct?
    Tp = (rotation*ct.*qdyn.*chord + F_addedmass_Tp)./cos.(delta)  .+ F_buoy[1, :] # TODO CM: correct?
    Zp = (cn.*qdyn.*chord - F_addedmass_Np).*tan.(delta)  .+ F_buoy[3, :] # TODO CM: correct?

    # nonlinear correction factor
    integrand = (W./V_wind).^2 .* (cn.*sin.(thetavec) - rotation*ct.*cos.(thetavec)./cos.(delta))
    CT = mean(sigma./(4*pi) .* pInt(thetavec, integrand))
    if CT > 2.0
        a = 0.5*(1.0 + sqrt(1.0 + CT))
        ka = 1.0 / (a-1)

    elseif CT > 0.96
        a = 1.0/7*(1 + 3.0*sqrt(7.0/2*CT - 3))
        ka = 18.0*a / (7*a^2 - 2*a + 4)

    else
        a = 0.5*(1 - sqrt(1.0 - CT))
        ka = 1.0 / (1-a)
    end

    # power coefficient
    H = 1.0  # per unit height
    Sref = 2*R*H
    Q = r.*Tp*rotation
    P = abs(mean(Omega))*B/(2*pi)*pInt(thetavec, Q)
    CP = P / (0.5*rho*mean(V_wind)^3 * Sref)

    return q, ka, CT, CP, Rp, Tp, Zp, a, alpha, cl, cd, Vn, Vt, Re, Q, M_addedmass_Np, M_addedmass_Tp, F_addedmass_Np, F_addedmass_Tp
end

# -----------------------------------------

# ------ Solve System --------------
"""
Internal, sets up the residual function
"""
function residual(w, A, theta, k, turbines, env;w_RPI=zeros(Real,1),idx_RPI=idx_RPI=1:2*turbines[1].ntheta)
    if length(w_RPI)>1
        w[idx_RPI] = w_RPI
    end
    # setup
    ntheta = length(theta)
    nturbines = Int(length(w)/2/ntheta)
    q = zeros(Real,ntheta*nturbines)
    ka = 0.0

    for i = 1:nturbines
        idx = collect((i-1)*ntheta+1:i*ntheta)

        u = w[idx]
        v = w[ntheta*nturbines .+ idx]

        q[idx], ka, _, _, _, _, _, _, _, cl, cd, Vn, Vt, Re = radialforce(u, v, theta, turbines[i], env)
    end

    if nturbines == 1  # if only one turbine use the k from the analysis
        k = [ka]
    end  # otherwise, use k that was input to this function

    # reformat to multiply in correct locations
    kmult = repeat(k, inner=[ntheta])
    kmult = [kmult; kmult]

    output = (A*q).*kmult - w

    return output[:]
end

"""
AC(turbines, env; w=zeros(Real,2*turbines[1].ntheta), idx_RPI=1:2*turbine.ntheta, solve=true, ifw=false)

see ?steady for detailed i/o description

Double multiple streamtube model
"""
function AC(turbines, env; w=zeros(Real,2*turbines[1].ntheta), idx_RPI=1:2*turbines[1].ntheta,solve=true,ifw=false)
    #TODO: make these modifications work for more than one turbine!
    if w==zeros(Real,2*turbines[1].ntheta)
        w[:].=0.0
    end
    if length(turbines) != 1
        error("Multiple AC Turbines with RPI/Unsteady not yet implemented!")
    end

    # list comprehensions
    centerX = [turbine.centerX for turbine in turbines]
    centerY = [turbine.centerY for turbine in turbines]
    radii = [turbine.R for turbine in turbines]

    ntheta = turbines[1].ntheta

    # assemble global matrices
    Ax, Ay, theta = matrixAssemble(turbines[1], centerX, centerY, radii, ntheta)

    # setup
    # ntheta = length(theta)
    nturbines = length(turbines)
    tol = 1e-6
    CT = zeros(Real,nturbines)
    CP = zeros(Real,nturbines)
    Rp = zeros(Real,ntheta, nturbines)
    Tp = zeros(Real,ntheta, nturbines)
    Zp = zeros(Real,ntheta, nturbines)

    q = zeros(Real,ntheta)

    # compute nonlinear correction factors (each turbine individaully)
    k = zeros(Real,nturbines)

    # for i = 1:length(turbines)
    i = 1
    if solve
        w0 = w[idx_RPI]

        idx = collect((i-1)*ntheta+1:i*ntheta)

        if length(idx_RPI)<ntheta*2
            resid_RPI(x) = residual(w,[Ax[idx, idx]; Ay[idx, idx]], theta, [1.0], turbines, env;w_RPI=x,idx_RPI)
            if ifw
                autodiff=:central
            else
                autodiff=:central #TODO: this runs most of the time, but sometimes not
            end
            result = LsqFit.lmfit(resid_RPI, Float64.(w0), Float64[];min_step_quality=1e-5,autodiff=autodiff) #TODO: figure out real vectors/how to do this solve in the midst of automatic gradient calcs
            w_RPI = result.param
            if !result.converged
                println("LMFIT terminated prematurely. Max Residual = ", maximum(abs.(result.resid)))
            end
        else
            resid_ALL(x) = residual(x,[Ax[idx, idx]; Ay[idx, idx]], theta, [1.0], turbines, env)
            if minimum(abs.(turbines[1].delta))<pi/4 || ifw
                autodiff=:central
            else #Weird error, haven't been able to figure out the type issue; when I debug it goes away... and... if I run it through the lmfit solver it isn't an issue
                autodiff=:central
            end
            result = NLsolve.nlsolve(resid_ALL, Float64.(w0), ftol=1e-3,xtol=1e-3,iterations=10,autodiff=autodiff)
            w_RPI = result.zero
            if !NLsolve.converged(result)
                println("NLsolve terminated prematurely. Residual Norm = ", result.residual_norm)
            end
        end

        w[idx_RPI] = w_RPI
    end #solve

    idx = collect(1:ntheta)
    u = w[idx]
    v = w[ntheta .+ idx]
    q, k, CT, CP, Rp, Tp, Zp, a, alpha, cl, cd, Vn, Vt, Re, Q, M_addedmass_Np, M_addedmass_Tp, F_addedmass_Np, F_addedmass_Tp = radialforce(u, v, theta, turbines[i], env)

    return CP, q ,Q, Rp, Tp, Zp, sqrt.(Vn.^2 .+ Vt.^2), CT, CT, a, w, alpha, cl, cd, theta, Re, M_addedmass_Np, M_addedmass_Tp, F_addedmass_Np, F_addedmass_Tp

end
# TODO: keep this as it shows how to handle multiple turbines, might also look at original code
# function AC_steady(turbines, env)
#     ntheta = turbines[1].ntheta
#
#     #Convert global F.O.R. winds to turbine frame
#     windangle = env.windangle
#
#     V_xtemp = env.V_x.*cos(windangle)+env.V_y.*sin(windangle) # Vinf is V_x, t is for turbine direction f.o.r.
#     V_ytemp = -env.V_x.*sin(windangle)+env.V_y.*cos(windangle)
#
#     env.V_x[:] = V_xtemp #TODO: ensure this doesn't mess up unsteady methods with persistent backend memory
#     env.V_y[:] = V_ytemp
#
#     # list comprehensions
#     centerX = [turbine.centerX for turbine in turbines]
#     centerY = [turbine.centerY for turbine in turbines]
#     radii = [turbine.R for turbine in turbines]
#
#     # assemble global matrices
#     Ax, Ay, theta = matrixAssemble(turbines[1],centerX, centerY, radii, ntheta)
#
#     # setup
#     ntheta = length(theta)
#     nturbines = length(turbines)
#     tol = 1e-6
#     CT = zeros(Real,nturbines)
#     CP = zeros(Real,nturbines)
#     Rp = zeros(Real,ntheta, nturbines)
#     Tp = zeros(Real,ntheta, nturbines)
#     Zp = zeros(Real,ntheta, nturbines)
#     a = zeros(Real, nturbines)
#     alpha = zeros(Real,ntheta, nturbines)
#     cl = zeros(Real,ntheta,nturbines)
#     cd = zeros(Real,ntheta,nturbines)
#     Vn = zeros(Real,ntheta,nturbines)
#     Vt = zeros(Real,ntheta,nturbines)
#     Re = zeros(Real,ntheta,nturbines)
#     wsave = zeros(Real,ntheta*2,nturbines)
#     q = zeros(Real,ntheta)
#
#     # compute nonlinear correction factors (each turbine individaully)
#     k = zeros(Real,nturbines)
#
#     for i = 1:length(turbines)
#         w0 = env.aw_warm # zeros(ntheta*2)
#
#         idx = collect((i-1)*ntheta+1:i*ntheta)
#
#         resid_single(x) = residual(x,[Ax[idx, idx]; Ay[idx, idx]], theta, [1.0], turbines, env)
#         result = NLsolve.nlsolve(resid_single, w0, ftol=tol,autodiff=:forwarddiff)
#         w = result.zero
#         if !NLsolve.converged(result)
#             println("NLsolve terminated prematurely.") # info = ", result)
#         end
#
#         idx = collect(1:ntheta)
#         u = w[idx]
#         v = w[ntheta .+ idx]
#         q, k[i], CT[i], CP[i], Rp[:, i], Tp[:, i], Zp[:, i], a[i], alpha[:, i], cl[:, i], cd[:, i], Vn[:, i], Vt[:, i], Re[:, i] = radialforce(u, v, theta, turbines[i], env)
#         wsave[:,i] = w
#     end
#
#
#     if nturbines == 1
#              # CP,      Th,   Q,      Rp, Tp, Zp, Vloc,                  CD,      CT,astar,                     _,                       alpha, cl, cd_af, thetavec, Re = OWENSAero.steady(turbine, env)
#         return CP, nothing ,nothing, Rp, Tp, Zp, sqrt.(Vn.^2 .+ Vt.^2), nothing, CT, a, wsave, alpha, cl, cd, theta.-windangle, Re
#     end
#
#     # Solve coupled system
#     w0 = zeros(Real,nturbines*ntheta*2) #TODO: warmstart from wsave flattened
#
#     resid_multiple(x) = residual(x,[Ax; Ay], theta, k, turbines, env)
#     result = NLsolve.nlsolve(resid_multiple, w0, ftol=tol)
#     w = result.zero
#     if !NLsolve.converged(result)
#         println("hybrd terminated prematurely. info = ", info)
#     end
#
#     for i in eachindex(turbines)
#         idx = collect((i-1)*ntheta+1:i*ntheta)
#
#         u = w[idx]
#         v = w[ntheta*nturbines .+ idx]
#         _, _, CT[i], CP[i], Rp[:, i], Tp[:, i], Zp[:, i], _, alpha[:, i], cl[:, i], cd[:, i], Vn[:, i], Vt[:, i], Re[:, i] = radialforce(u, v, theta, turbines[i], env)
#
#     end
#     @warn "interface for multiple turbines not unified yet"
#     return CT, CP, Rp, Tp, Zp, theta, alpha, cl, cd, Vn, Vt, w, Re
# end

# -----------------------------------------

# ---------- helper methods --------------

"""
Internal, trapezoidal integration of y w.r.t. x
"""
function trapz(x, y)  # integrate y w.r.t. x

    integral = 0.0
    for i = 1:length(x)-1
        integral += (x[i+1]-x[i])*0.5*(y[i] + y[i+1])
    end
    return integral
end

"""
Internal, integration for a periodic function where end points don't reach ends (uses trapezoidal method)
"""
function pInt(theta, f)

    integral = trapz(theta, f)

    # add end points
    dtheta = 2*theta[1]  # assumes equally spaced, starts at 0
    integral += dtheta * 0.5*(f[1] + f[end])

    return integral
end
