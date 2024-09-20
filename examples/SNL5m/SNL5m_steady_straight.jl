
import PyPlot
PyPlot.close("all")
import Statistics:mean
import DelimitedFiles
import Dierckx
import QuadGK
import FLOWMath
import VAWTAero

close("all")
path,_ = splitdir(@__FILE__)
# include("$(path)/../../../src/VAWTAero.jl")

PyPlot.rc("figure", figsize=(4, 3))
PyPlot.rc("font", size=10.0)
PyPlot.rc("lines", linewidth=1.5)
PyPlot.rc("lines", markersize=3.0)
PyPlot.rc("legend", frameon=false)
PyPlot.rc("axes.spines", right=false, top=false)
PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
PyPlot.rc("figure",max_open_warning=500)
# rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]
PyPlot.ion()
function runme(A_model,iscurved)
    # Main Function
    DS_model = "BV"
    r_delta_infl = false
    ntheta = 30
    R = 5.0/2 #m
    H = 1.02*R*2; #m
    chord = 0.1524 #m
    plots = true
    eta = 0.4 #blade mount point, or distance from leading edge that the rotation plane passes through.

    #                      #image rotation correction
    shapeX_raw_temp = 1 / cos(12.38*pi/180) .* [ -6.29E-03, 1.33E-01, 2.93E-01, 4.45E-01, 5.71E-01, 6.56E-01, 7.22E-01, 7.81E-01, 8.29E-01, 8.77E-01, 9.16E-01, 9.46E-01, 9.65E-01, 9.77E-01, 9.83E-01, 9.83E-01, 9.83E-01, 9.79E-01, 9.68E-01, 9.54E-01, 9.34E-01, 9.11E-01, 8.93E-01, 8.68E-01, 8.45E-01, 8.13E-01, 7.76E-01, 7.30E-01, 6.73E-01, 6.27E-01, 5.34E-01, 4.13E-01, 2.92E-01, 1.76E-01, 5.88E-02, 1.74E-03]
    shapeX_raw = R * (shapeX_raw_temp./maximum(shapeX_raw_temp))
    shapeY_raw_temp = [-4.12E-03,3.25E-02,8.62E-02,1.36E-01,1.79E-01,2.08E-01,2.35E-01,2.63E-01,2.90E-01,3.19E-01,3.56E-01,3.92E-01,4.24E-01,4.57E-01,4.84E-01,5.01E-01,5.20E-01,5.47E-01,5.84E-01,6.09E-01,6.41E-01,6.65E-01,6.84E-01,7.05E-01,7.20E-01,7.40E-01,7.58E-01,7.83E-01,8.02E-01,8.17E-01,8.44E-01,8.82E-01,9.20E-01,9.56E-01,9.93E-01,1.01E+00]
    shapeY_raw = H * shapeY_raw_temp./maximum(shapeY_raw_temp)
    shapeX_spline = FLOWMath.Akima(shapeY_raw, shapeX_raw)
    n_slices = 30#length(shapeX)-1;
    shapeY = collect(LinRange(0,H,n_slices+1))
    shapeX = R.*(1.0.-4.0.*(shapeY/H.-.5).^2)#shapeX_spline(shapeY)
    h_frac = (shapeY[2:end] - shapeY[1:end-1])./shapeY[end];
    h_elem = (shapeY[2:end] - shapeY[1:end-1])
    h = (shapeY[2:end] + shapeY[1:end-1])/2.0;

    RefArea_half, error = QuadGK.quadgk(shapeX_spline, 0, H, atol=1e-10)
    RefArea = RefArea_half*2
    delta_xs = shapeX[2:end] - shapeX[1:end-1]
    delta_zs = shapeY[2:end] - shapeY[1:end-1]

    delta3D = atan.(delta_xs./delta_zs)

    if iscurved=="straight"
        delta3D = delta3D.*0.0
    end

    element_planf_A = sqrt.(delta_xs.^2+delta_zs.^2)*chord
    element_planf_L = sqrt.(delta_xs.^2+delta_zs.^2)

    r3D = (shapeX[2:end,1]+shapeX[1:end-1,1])/2.0
    aerocenter_dist = (eta-.25)*chord

    twist3D = -atan.(aerocenter_dist./r3D)#ones(n_slices)*-0.4*pi/180

    if DS_model=="BV"
        af3D = VAWTAero.readaerodyn_BV("$(path)/airfoils/NACA_0015_RE3E5.dat")
    else
        af3D = VAWTAero.readaerodyn("$(path)/airfoils/NACA_0015_RE3E5.dat")
    end

    # Unchanging Parameters
    B = 3
    omega = -ones(ntheta) * 150.0 / 60.0 * 2*pi # RPM -> radians/sec
    rho = 1.225*0.8 #80percent at albuquerque elevation
    mu = 1.7894e-5


    # tsrvec = [2.1,3.1,4.2,4.7,5.2,6.0,7.1] #FOR PLOTTING
    tsrvec = [5.2] #FOR TESTING
    n_tsr = length(tsrvec)
    CPvec = zeros(n_tsr)
    CTvec = zeros(n_tsr)
    Qvec = zeros(n_tsr)
    Vinf_used = zeros(n_tsr)
    CT = zeros(n_slices,n_tsr)
    CP = zeros(n_slices,n_tsr)
    Qave = zeros(n_slices,n_tsr)
    Rp = zeros(ntheta,n_slices,n_tsr)
    Tp = zeros(ntheta,n_slices,n_tsr)
    Rp_out = zeros(n_slices,n_tsr)
    Tp_out = zeros(n_slices,n_tsr)
    Zp_out = zeros(n_slices,n_tsr)
    Zp = zeros(ntheta,n_slices,n_tsr)
    alpha = zeros(ntheta,n_slices,n_tsr)
    cl = zeros(ntheta,n_slices,n_tsr)
    cd = zeros(ntheta,n_slices,n_tsr)
    W = zeros(ntheta,n_slices,n_tsr)
    Re = zeros(ntheta,n_slices,n_tsr)
    w = zeros(ntheta*2,n_slices,n_tsr)
    thetavec = zeros(ntheta)

    Rp_dist_cactus = zeros(n_slices,n_tsr)
    Tp_dist_cactus = zeros(n_slices,n_tsr)
    Zp_dist_cactus = zeros(n_slices,n_tsr)
    CP_dist_cactus = zeros(n_slices,n_tsr)
    te_cactus = zeros(n_slices,n_tsr)
    cactus_CP = zeros(n_tsr)

    println("Aero_model: $(A_model), DS_model: $(DS_model), r_delta_infl: $(r_delta_infl)")
    dataout = Array{Any}(undef,n_tsr)
    for tsr_idx = 1:n_tsr

        tsr_digit1 = Int(floor(tsrvec[tsr_idx]))
        tsr_digit2 = Int(round((tsrvec[tsr_idx]-tsr_digit1)*10.0))
        Vinf = abs.(omega)/tsrvec[tsr_idx]*R

        cactus_Rev_data = DelimitedFiles.readdlm("$(path)/data/steady/TestVAWT$(tsr_digit1)_$(tsr_digit2)_troposkein_RevData.csv",',',Float64,skipstart = 1)
        cactus_CP[tsr_idx] = cactus_Rev_data[end,2]

        data = DelimitedFiles.readdlm("$(path)/data/steady/TestVAWT$(tsr_digit1)_$(tsr_digit2)_troposkein_ElementData.csv", ',',skipstart = 1)
        #azimuthal discretization, data columns, slices
        data_full_rev = zeros(ntheta,length(data[1,:]),Int(maximum(data[:,4])))
        elapsed = zeros(n_slices)
        for slice = 1:n_slices
            # start = time()
            # println(slice)
            # Potentially Aerostructural Parameters that my change
            r = ones(ntheta)*r3D[slice] #m
            twist = ones(ntheta)*twist3D[slice] #rad
            delta = ones(ntheta)*delta3D[slice] #rad
            af = af3D

            turbine = VAWTAero.Turbine(R,r,chord,twist,delta,omega,B,af,ntheta,r_delta_infl)
            if slice == 1 && tsr_idx != 1
                aw_warm = w[:,slice,tsr_idx-1]
            elseif tsr_idx == 1
                aw_warm = zeros(ntheta*2)
            else
                aw_warm = w[:,slice-1,tsr_idx]
            end
            env = VAWTAero.Environment(rho,mu,Vinf,DS_model,A_model,aw_warm)

            # Option to Plot the model
            # VAWTAero.plot_domain(turbine)

            start = time()

            # Juno.@enter VAWTAero.steady(turbine, env)
            CP_in, _, _, Rp_in, Tp_in, Zp_in,W[:,slice,tsr_idx], _, _, _, w_in, alpha_in,
            cl_in, cd_in, thetavec[:], Re_in = VAWTAero.steady(turbine, env)

            CP[slice,tsr_idx] = CP_in[1]
            Tp[:,slice,tsr_idx] = Tp_in
            Rp[:,slice,tsr_idx] = Rp_in
            Zp[:,slice,tsr_idx] = Zp_in
            alpha[:,slice,tsr_idx] = alpha_in
            cl[:,slice,tsr_idx] = cl_in
            cd[:,slice,tsr_idx] = cd_in
            Re[:,slice,tsr_idx] = Re_in
            w[1:ntheta,slice,tsr_idx] = w_in[1:ntheta]
            Rp_out[slice,tsr_idx] = B/(2*pi)*VAWTAero.pInt(thetavec, abs.(Rp[:,slice,tsr_idx])) #average absolute force over the revolution
            Tp_out[slice,tsr_idx] = B/(2*pi)*VAWTAero.pInt(thetavec, Tp_in) #average NON-absolute force over the revolution
            Zp_out[slice,tsr_idx] = B/(2*pi)*VAWTAero.pInt(thetavec, abs.(Zp[:,slice,tsr_idx])) #average absolute force over the revolution

            elapsed[slice] = time() - start

            # println("$(A_model)_$(DS_model)_TSR: $(tsrvec[tsr_idx]) Slice $(slice) Time: $(time()-start)")
            # Blade 1
            #                                  Next to Last revolution           Blade 1                Slice #
            used_data_logic = ((data[:,5].==maximum(data[:,5])-1) .& (data[:,3].==1.0)) .& (data[:,4].==slice)
            data_full_rev_b1 = data[used_data_logic,:]
            data_full_rev[:,:,slice] = data_full_rev_b1

            q_loc_b1 = 0.5*rho*element_planf_A[slice].*(data_full_rev_b1[:,15].*Vinf).^2
            Qp_cactus_b1 = q_loc_b1.*data_full_rev_b1[:,25]./element_planf_L[slice].*r3D[slice]

            # Blade 2
            used_data_logic = ((data[:,5].==maximum(data[:,5])-1.0) .& (data[:,3].==2.0)) .& (data[:,4].==slice)
            data_full_rev_b2 = data[used_data_logic,:]
            q_loc_b2 = 0.5*rho*element_planf_A[slice].*(data_full_rev_b2[:,15].*Vinf).^2
            Qp_cactus_b2 = q_loc_b2.*data_full_rev_b2[:,25]./element_planf_L[slice].*r3D[slice]

            # Blade 3
            used_data_logic = ((data[:,5].==maximum(data[:,5])-1.0) .& (data[:,3].==3.0)) .& (data[:,4].==slice)
            data_full_rev_b3 = data[used_data_logic,:]
            q_loc_b3 = 0.5*rho*element_planf_A[slice].*(data_full_rev_b3[:,15].*Vinf).^2
            Qp_cactus_b3 = q_loc_b3.*data_full_rev_b3[:,25]./element_planf_L[slice].*r3D[slice]

            P_cactus = abs(mean(abs.(omega)))/(2*pi)*VAWTAero.pInt(thetavec, Qp_cactus_b1+Qp_cactus_b2+Qp_cactus_b3)

            H_cactus = 1.0  # per unit height
            Sref_cactus = 2*R*H_cactus

            CP_dist_cactus[slice,tsr_idx] = -P_cactus / (0.5*rho*mean(Vinf)^3 * RefArea)
            # convert local tangential coefficient to turbine tangential coefficient to force averaged over the revolution (ntheta)
            Rp_dist_cactus[slice,tsr_idx] = sum(0.5*rho*chord*((data_full_rev_b1[:,15].*Vinf).^2 .*abs.(data_full_rev_b1[:,24])+(data_full_rev_b2[:,15].*Vinf).^2 .*abs.(data_full_rev_b2[:,24])+(data_full_rev_b3[:,15].*Vinf).^2 .*abs.(data_full_rev_b3[:,24])))/ntheta * cos(delta3D[slice])
            Tp_dist_cactus[slice,tsr_idx] = sum(0.5*rho*chord*((data_full_rev_b1[:,15].*Vinf).^2 .*(data_full_rev_b1[:,25])+(data_full_rev_b2[:,15].*Vinf).^2 .*(data_full_rev_b2[:,25])+(data_full_rev_b3[:,15].*Vinf).^2 .*(data_full_rev_b3[:,25])))/ntheta # Average NON-absolute
            Zp_dist_cactus[slice,tsr_idx] = sum(0.5*rho*chord*((data_full_rev_b1[:,15].*Vinf).^2 .*abs.(data_full_rev_b1[:,24])+(data_full_rev_b2[:,15].*Vinf).^2 .*abs.(data_full_rev_b2[:,24])+(data_full_rev_b3[:,15].*Vinf).^2 .*abs.(data_full_rev_b3[:,24])))/ntheta * sin(delta3D[slice])
            te_cactus[slice,tsr_idx] = sum(data_full_rev_b1[:,29]+data_full_rev_b2[:,29]+data_full_rev_b3[:,29])/ntheta #average over revolution by dividing by ntheta

            q_loc = 0.5*rho.*(data_full_rev[:,15,slice].*Vinf).^2

            # if plots
            #     # Angle of Attack
            #     PyPlot.figure()
            #     PyPlot.plot((data_full_rev[:,2,slice].-minimum(data_full_rev[:,2,slice]))*180/pi, data_full_rev[:,9,slice],"k")
            #     PyPlot.plot((data_full_rev[:,2,slice].-minimum(data_full_rev[:,2,slice]))*180/pi, data_full_rev[:,10,slice],"k--")
            #     PyPlot.plot((data_full_rev[:,2,slice].-minimum(data_full_rev[:,2,slice]))*180/pi, data_full_rev[:,11,slice],"k-.")
            #     PyPlot.plot(thetavec*180/pi, alpha[:,slice,tsr_idx]*180/pi,color=plot_cycle[1])
            #     PyPlot.xlabel("Theta (deg)")
            #     PyPlot.ylabel("AOA (deg)")
            #     PyPlot.legend(["Cactus 25","Cactus 50","Cactus 75","$A_model"])
            #     # PyPlot.legend(["Cactus $(sum(q_loc.*data_full_rev[:,25,slice]))","$A_model $(sum(h_elem[slice].*Tp[:,slice,tsr_idx]))"])
            #     PyPlot.savefig("$(path)/figs/troposkein/SNL5m_AOA_$(A_model)_$(DS_model)_TSR_$(tsr_digit1)_$(tsr_digit2)_slice_$slice.pdf",transparent = true)
            #
            #     # cl
            #     PyPlot.figure()
            #     PyPlot.plot((data_full_rev[:,2,slice].-minimum(data_full_rev[:,2,slice]))*180/pi, data_full_rev[:,20,slice],"k")
            #     PyPlot.plot(thetavec*180/pi, cl[:,slice,tsr_idx],color=plot_cycle[1])
            #     PyPlot.xlabel("Theta (deg)")
            #     PyPlot.ylabel("cl")
            #     PyPlot.legend(["Cactus ","$A_model"])
            #     # PyPlot.legend(["Cactus $(sum(q_loc.*data_full_rev[:,25,slice]))","$A_model $(sum(h_elem[slice].*Tp[:,slice,tsr_idx]))"])
            #     PyPlot.savefig("$(path)/figs/troposkein/SNL5m_cl_$(A_model)_$(DS_model)_TSR_$(tsr_digit1)_$(tsr_digit2)_slice_$slice.pdf",transparent = true)
            #
            #     # Re
            #     PyPlot.figure()
            #     PyPlot.plot((data_full_rev[:,2,slice].-minimum(data_full_rev[:,2,slice]))*180/pi, data_full_rev[:,13,slice],"k")
            #     PyPlot.plot(thetavec*180/pi, Re[:,slice,tsr_idx],color=plot_cycle[1])
            #     PyPlot.xlabel("Theta (deg)")
            #     PyPlot.ylabel("Re")
            #     PyPlot.legend(["Cactus ","$A_model"])
            #     # PyPlot.legend(["Cactus $(sum(q_loc.*data_full_rev[:,25,slice]))","$A_model $(sum(h_elem[slice].*Tp[:,slice,tsr_idx]))"])
            #     PyPlot.savefig("$(path)/figs/troposkein/SNL5m_Re_$(A_model)_$(DS_model)_TSR_$(tsr_digit1)_$(tsr_digit2)_slice_$slice.pdf",transparent = true)
            #
            #     # cd
            #     PyPlot.figure()
            #     PyPlot.plot((data_full_rev[:,2,slice].-minimum(data_full_rev[:,2,slice]))*180/pi, data_full_rev[:,21,slice],"k")
            #     PyPlot.plot(thetavec*180/pi, cd[:,slice,tsr_idx],color=plot_cycle[1])
            #     PyPlot.xlabel("Theta (deg)")
            #     PyPlot.ylabel("cd")
            #     PyPlot.legend(["Cactus ","$A_model"])
            #     # PyPlot.legend(["Cactus $(sum(q_loc.*data_full_rev[:,25,slice]))","$A_model $(sum(h_elem[slice].*Tp[:,slice,tsr_idx]))"])
            #     PyPlot.savefig("$(path)/figs/troposkein/SNL5m_cd_$(A_model)_$(DS_model)_TSR_$(tsr_digit1)_$(tsr_digit2)_slice_$slice.pdf",transparent = true)
            #
            #     # W
            #     PyPlot.figure()
            #     PyPlot.plot((data_full_rev[:,2,slice].-minimum(data_full_rev[:,2,slice]))*180/pi, data_full_rev[:,15,slice].*Vinf,"k")
            #     PyPlot.plot(thetavec*180/pi, W[:,slice,tsr_idx],color=plot_cycle[1])
            #     PyPlot.xlabel("Theta (deg)")
            #     PyPlot.ylabel("Local Velocity")
            #     PyPlot.legend(["Cactus ","$A_model"])
            #     # PyPlot.legend(["Cactus $(sum(q_loc.*data_full_rev[:,25,slice]))","$A_model $(sum(h_elem[slice].*Tp[:,slice,tsr_idx]))"])
            #     PyPlot.savefig("$(path)/figs/troposkein/SNL5m_W_$(A_model)_$(DS_model)_TSR_$(tsr_digit1)_$(tsr_digit2)_slice_$slice.pdf",transparent = true)
            #
            #     # W - omega*r - Vinf*cos(thetavec)
            #     PyPlot.figure()
            #     PyPlot.plot((data_full_rev[:,2,slice].-minimum(data_full_rev[:,2,slice]))*180/pi, data_full_rev[:,15,slice].*Vinf.- omega[1]*r3D[slice] .- Vinf.*cos.(thetavec),"k")
            #     PyPlot.plot(thetavec*180/pi, W[:,slice,tsr_idx] .- omega[1]*r3D[slice] .- Vinf.*cos.(thetavec),color=plot_cycle[1])
            #     PyPlot.xlabel("Theta (deg)")
            #     PyPlot.ylabel("Local Velocity minus rotation and wind speed")
            #     PyPlot.legend(["Cactus ","$A_model"])
            #     # PyPlot.legend(["Cactus $(sum(q_loc.*data_full_rev[:,25,slice]))","$A_model $(sum(h_elem[slice].*Tp[:,slice,tsr_idx]))"])
            #     PyPlot.savefig("$(path)/figs/troposkein/SNL5m_W2_$(A_model)_$(DS_model)_TSR_$(tsr_digit1)_$(tsr_digit2)_slice_$slice.pdf",transparent = true)
            #
            #     # u
            #     PyPlot.figure()
            #     PyPlot.plot((data_full_rev[:,2,slice].-minimum(data_full_rev[:,2,slice]))*180/pi, data_full_rev[:,16,slice].*Vinf,"k")
            #     PyPlot.plot(thetavec*180/pi, w[1:ntheta,slice,tsr_idx]*180/pi,color=plot_cycle[1])
            #     PyPlot.xlabel("Theta (deg)")
            #     PyPlot.ylabel("u induced velocity")
            #     PyPlot.legend(["Cactus ","$A_model"])
            #     # PyPlot.legend(["Cactus $(sum(q_loc.*data_full_rev[:,25,slice]))","$A_model $(sum(h_elem[slice].*Tp[:,slice,tsr_idx]))"])
            #     PyPlot.savefig("$(path)/figs/troposkein/SNL5m_u_$(A_model)_$(DS_model)_TSR_$(tsr_digit1)_$(tsr_digit2)_slice_$slice.pdf",transparent = true)
            #
            #     # v
            #     PyPlot.figure()
            #     PyPlot.plot((data_full_rev[:,2,slice].-minimum(data_full_rev[:,2,slice]))*180/pi, data_full_rev[:,17,slice].*Vinf,"k")
            #     PyPlot.plot(thetavec*180/pi, w[ntheta+1:ntheta*2,slice,tsr_idx]*180/pi,color=plot_cycle[1])
            #     PyPlot.xlabel("Theta (deg)")
            #     PyPlot.ylabel("v induced velocity")
            #     PyPlot.legend(["Cactus ","$A_model"])
            #     # PyPlot.legend(["Cactus $(sum(q_loc.*data_full_rev[:,25,slice]))","$A_model $(sum(h_elem[slice].*Tp[:,slice,tsr_idx]))"])
            #     PyPlot.savefig("$(path)/figs/troposkein/SNL5m_v_$(A_model)_$(DS_model)_TSR_$(tsr_digit1)_$(tsr_digit2)_slice_$slice.pdf",transparent = true)
            #
            #     # Tangential Force
            #     PyPlot.figure()
            #     PyPlot.plot((data_full_rev[:,2,slice].-minimum(data_full_rev[:,2,slice]))*180/pi, q_loc.*data_full_rev[:,25,slice]*chord,"k")
            #     PyPlot.plot(thetavec*180/pi, Tp[:,slice,tsr_idx].*cos.(delta3D[slice]),color=plot_cycle[1])
            #     PyPlot.xlabel("Theta (deg)")
            #     PyPlot.ylabel("Tangential Force per span (N/m)")
            #     PyPlot.legend(["Cactus","$A_model"])
            #     # PyPlot.legend(["Cactus $(sum(q_loc.*data_full_rev[:,25,slice]))","$A_model $(sum(h_elem[slice].*Tp[:,slice,tsr_idx]))"])
            #     PyPlot.savefig("$(path)/figs/troposkein/SNL5m_Tangential_Force_$(A_model)_$(DS_model)_TSR_$(tsr_digit1)_$(tsr_digit2)_slice_$slice.pdf",transparent = true)
            #
            #
            #     # Radial Force
            #     PyPlot.figure()
            #     PyPlot.plot((data_full_rev[:,2,slice].-minimum(data_full_rev[:,2,slice]))*180/pi, q_loc.*data_full_rev[:,24,slice].*chord.*cos(delta3D[slice]),"k")
            #     PyPlot.plot(thetavec*180/pi, Rp[:,slice,tsr_idx].*cos.(delta3D[slice]),color=plot_cycle[1])
            #     PyPlot.xlabel("Theta (deg)")
            #     PyPlot.ylabel("Radial Force per span (N/m)")
            #     PyPlot.legend(["Cactus","$A_model"])
            #     PyPlot.savefig("$(path)/figs/troposkein/SNL5m_Radial_Force_$(A_model)_$(DS_model)_TSR_$(tsr_digit1)_$(tsr_digit2)_slice_$slice.pdf",transparent = true)
            #
            #     # # Torque Force
            #     # PyPlot.figure()
            #     # PyPlot.plot((data_full_rev[:,2,slice].-minimum(data_full_rev[:,2,slice]))*180/pi, q_loc.*data_full_rev[:,24,slice],"k")
            #     # PyPlot.plot(thetavec*180/pi, h_elem[slice].*-Rp[:,slice,tsr_idx],color=plot_cycle[1])
            #     # PyPlot.xlabel("Theta (deg)")
            #     # PyPlot.ylabel("Torque")
            #     # PyPlot.legend(["Cactus",A_model])
            #     # PyPlot.savefig("$(path)/figs/troposkein/SNL5m_Torque_$(A_model)_$(DS_model)_TSR_5_2_slice_$slice.pdf",transparent = true)
            #
            #     # Vertical Force
            #     PyPlot.figure()
            #     PyPlot.plot((data_full_rev[:,2,slice].-minimum(data_full_rev[:,2,slice]))*180/pi, q_loc.*data_full_rev[:,24,slice].*chord.*sin(delta3D[slice]),"k") #Normal to vertical by sin(delta)
            #     PyPlot.plot(thetavec*180/pi, Zp[:,slice,tsr_idx].*cos.(delta3D[slice]),color=plot_cycle[1]) #remove the blade span length by multiplying by cos(delta)
            #     PyPlot.xlabel("Theta (deg)")
            #     PyPlot.ylabel("Vertical Force per span (N/m)")
            #     PyPlot.legend(["Cactus","$A_model"])
            #     PyPlot.savefig("$(path)/figs/troposkein/SNL5m_Vertical_Force_$(A_model)_$(DS_model)_TSR_$(tsr_digit1)_$(tsr_digit2)_slice_$slice.pdf",transparent = true)
            #
            #     PyPlot.close("all")
            # end

        end
        mean_elapsed = mean(elapsed)
        println("$(mean_elapsed) seconds")

        # CPvec[tsr_idx] = sum(CP[:,tsr_idx].*h_frac)
        CPvec[tsr_idx] = VAWTAero.trapz(h,CP[:,tsr_idx].*2*R./RefArea) #Undo local normalization and normalize by full turbine

        dataout[tsr_idx] = data_full_rev
    end

    # Validation against Experimental, Cactus for CP
    experimental = DelimitedFiles.readdlm("$(path)/data/steady/SNL5m_CP_150RPM_0015.txt")

    # if plots
    #     PyPlot.figure()
    #     PyPlot.plot(experimental[1:end-1,1],experimental[1:end-1,2],"k.")
    #     PyPlot.plot(tsrvec,cactus_CP,"-",color = plot_cycle[1])
    #     PyPlot.plot(tsrvec,CPvec,"-",color = plot_cycle[2])
    #     PyPlot.xlabel("Tip Speed Ratio")
    #     PyPlot.ylabel("Coefficient of Performance")
    #     PyPlot.legend(["Experimental", "Cactus", "$(A_model)_$DS_model"])
    #     PyPlot.savefig("$(path)/figs/troposkein/SNL5m_CP3D_$(A_model)_$DS_model.pdf",transparent = true)
    # end


    return delta3D,experimental,cactus_CP,tsrvec,CPvec,shapeX,shapeY,element_planf_L,Rp_out,Tp_out,Zp_out,CP,h,h_frac,Rp_dist_cactus,Tp_dist_cactus,Zp_dist_cactus,CP_dist_cactus,te_cactus,r3D,R,omega,B,thetavec,RefArea, Tp, Rp, Zp, dataout,tsrvec

end

# Juno.@enter runme("DMS","noDS")
delta3D,experimental,cactus_CP,tsrvec,CPvec_DMS_straight,shapeX,shapeY,element_planf_L,Rp_DMS_straight,Tp_DMS_straight,Zp_DMS_straight,Cp_1DMS_straight,h,h_frac,Rp_dist_cactus,Tp_dist_cactus,Zp_dist_cactus,CP_dist_cactus,te_cactus,r3D,R,omega,B,thetavec,RefArea,Tp_DMS_straight_theta,Rp_DMS_straight_theta,Zp_DMS_straight_theta,dataout,tsrvec = runme("DMS","straight")
delta3D,experimental,cactus_CP,tsrvec,CPvec_DMS_curved,shapeX,shapeY,element_planf_L,Rp_DMS_curved,Tp_DMS_curved,Zp_DMS_curved,Cp_1DMS_curved,h,h_frac,Rp_dist_cactus,Tp_dist_cactus,Zp_dist_cactus,CP_dist_cactus,te_cactus,r3D,R,omega,B,thetavec,RefArea,Tp_DMS_curved_theta,Rp_DMS_curved_theta,Zp_DMS_curved_theta,dataout,tsrvec = runme("DMS","curved")
delta3D,experimental,cactus_CP,tsrvec,CPvec_AC_straight,shapeX,shapeY,element_planf_L,Rp_AC_straight,Tp_AC_straight,Zp_AC_straight,Cp_1AC_straight,h,h_frac,Rp_dist_cactus,Tp_dist_cactus,Zp_dist_cactus,CP_dist_cactus,te_cactus,r3D,R,omega,B,thetavec,RefArea,Tp_AC_straight_theta,Rp_AC_straight_theta,Zp_AC_straight_theta,dataout,tsrvec = runme("AC","straight")
delta3D,experimental,cactus_CP,tsrvec,CPvec_AC_curved,shapeX,shapeY,element_planf_L,Rp_AC_curved,Tp_AC_curved,Zp_AC_curved,Cp_1AC_curved,h,h_frac,Rp_dist_cactus,Tp_dist_cactus,Zp_dist_cactus,CP_dist_cactus,te_cactus,r3D,R,omega,B,thetavec,RefArea,Tp_AC_curved_theta,Rp_AC_curved_theta,Zp_AC_curved_theta,dataout,tsrvec = runme("AC","curved")

h_elem = (shapeY[2:end] - shapeY[1:end-1])

CP_total_dist_cactus2 = zeros(length(tsrvec))
CP_total_AC_straight = zeros(length(tsrvec))
CP_total_AC_curved = zeros(length(tsrvec))
CP_total_DMS_straight = zeros(length(tsrvec))
CP_total_DMS_curved = zeros(length(tsrvec))
CP_total_dist_cactus3 = zeros(length(tsrvec))
chord = 0.1524
rho = 1.225*0.8
# for tsr_idx = 1:length(tsrvec)
tsr_idx = 1
tsr_digit1 = Int(floor(tsrvec[tsr_idx]))
tsr_digit2 = Int(round((tsrvec[tsr_idx]-tsr_digit1)*10.0))
data_full_rev = dataout[tsr_idx]
R = 5.0/2
Vinf = abs.(omega)/tsrvec[tsr_idx]*R


############################################
# Force normalizaion
# Tangential Force

# Calculate the integral torque

Qint_AC_curved = VAWTAero.trapz(h,-Tp_AC_curved[:,tsr_idx].*cos.(delta3D))
Qint_DMS_curved = VAWTAero.trapz(h,-Tp_DMS_curved[:,tsr_idx].*cos.(delta3D))
Qint_AC_straight = VAWTAero.trapz(h,-Tp_AC_straight[:,tsr_idx].*cos.(delta3D))
Qint_DMS_straight = VAWTAero.trapz(h,-Tp_DMS_straight[:,tsr_idx].*cos.(delta3D))
Qint_cactus = VAWTAero.trapz(h,te_cactus[:,tsr_idx].*(0.5*1.225*0.8*mean(Vinf)^2 * RefArea * R)./r3D./h_elem.*cos.(delta3D))

err_Qint_AC_curved = (Qint_AC_curved-Qint_cactus)/Qint_cactus*100
err_Qint_DMS_curved = (Qint_DMS_curved-Qint_cactus)/Qint_cactus*100
err_Qint_AC_straight = (Qint_AC_straight-Qint_cactus)/Qint_cactus*100
err_Qint_DMS_straight = (Qint_DMS_straight-Qint_cactus)/Qint_cactus*100

println("err_Qint_AC_curved: $err_Qint_AC_curved %")
println("err_Qint_DMS_curved: $err_Qint_DMS_curved %")
println("err_Qint_AC_straight: $err_Qint_AC_straight %")
println("err_Qint_DMS_straight: $err_Qint_DMS_straight %")


# Calculate peaks
n_slices = 30
for islice = 1:Int(n_slices/2)
islice = 7
q_loc = 0.5*rho.*(data_full_rev[:,15,islice].*Vinf).^2
peak_cactus = maximum(-q_loc.*data_full_rev[:,25,islice]*chord)
peak_AC_curved = maximum(-Tp_AC_curved_theta[:,islice,tsr_idx].*cos.(delta3D[islice]))
peak_DMS_curved = maximum(-Tp_DMS_curved_theta[:,islice,tsr_idx].*cos.(delta3D[islice]))
peak_AC_straight = maximum(-Tp_AC_straight_theta[:,islice,tsr_idx].*cos.(delta3D[islice]))
peak_DMS_straight = maximum(-Tp_DMS_straight_theta[:,islice,tsr_idx].*cos.(delta3D[islice]))

minpeak_cactus = abs(minimum(-q_loc.*data_full_rev[:,25,islice]*chord))
minpeak_AC_curved = abs(minimum(-Tp_AC_curved_theta[:,islice,tsr_idx].*cos.(delta3D[islice])))
minpeak_DMS_curved = abs(minimum(-Tp_DMS_curved_theta[:,islice,tsr_idx].*cos.(delta3D[islice])))
minpeak_AC_straight = abs(minimum(-Tp_AC_straight_theta[:,islice,tsr_idx].*cos.(delta3D[islice])))
minpeak_DMS_straight = abs(minimum(-Tp_DMS_straight_theta[:,islice,tsr_idx].*cos.(delta3D[islice])))


err_peak_AC_curved = abs(peak_AC_curved-peak_cactus)/peak_cactus*100
err_peak_DMS_curved = abs(peak_DMS_curved-peak_cactus)/peak_cactus*100
err_peak_AC_straight = abs(peak_AC_straight-peak_cactus)/peak_cactus*100
err_peak_DMS_straight = abs(peak_DMS_straight-peak_cactus)/peak_cactus*100

err_minpeak_AC_curved = abs(minpeak_AC_curved-minpeak_cactus)/minpeak_cactus*100
err_minpeak_DMS_curved = abs(minpeak_DMS_curved-minpeak_cactus)/minpeak_cactus*100
err_minpeak_AC_straight = abs(minpeak_AC_straight-minpeak_cactus)/minpeak_cactus*100
err_minpeak_DMS_straight = abs(minpeak_DMS_straight-minpeak_cactus)/minpeak_cactus*100

println(islice)
println("err_peak_AC_curved $(err_peak_AC_curved) %")
println("err_peak_DMS_curved $(err_peak_DMS_curved) %")
println("err_peak_AC_straight $(err_peak_AC_straight) %")
println("err_peak_DMS_straight $(err_peak_DMS_straight) %")

println("err_minpeak_AC_curved $(err_minpeak_AC_curved) %")
println("err_minpeak_DMS_curved $(err_minpeak_DMS_curved) %")
println("err_minpeak_AC_straight $(err_minpeak_AC_straight) %")
println("err_minpeak_DMS_straight $(err_minpeak_DMS_straight) %")

Qint_AC_curved = VAWTAero.trapz(h,-Tp_AC_curved[:,tsr_idx].*cos.(delta3D).*r3D)
Qint_DMS_curved = VAWTAero.trapz(h,-Tp_DMS_curved[:,tsr_idx].*cos.(delta3D).*r3D)
Qint_AC_straight = VAWTAero.trapz(h,-Tp_AC_straight[:,tsr_idx].*cos.(delta3D).*r3D)
Qint_DMS_straight = VAWTAero.trapz(h,-Tp_DMS_straight[:,tsr_idx].*cos.(delta3D).*r3D)
Qint_cactus = VAWTAero.trapz(h,te_cactus[:,tsr_idx].*(0.5*1.225*0.8*mean(Vinf)^2 * RefArea * R)./h_elem.*cos.(delta3D))

cactmean = te_cactus[islice,tsr_idx].*(0.5*1.225*0.8*mean(Vinf)^2 * RefArea * R)./r3D[islice]./h_elem[islice].*cos.(delta3D[islice])

err_meanTpQ_AC_curved = (-Tp_AC_curved[islice,tsr_idx].*cos.(delta3D[islice])-cactmean)/cactmean*100
err_meanTpQ_DMS_curved = (-Tp_DMS_curved[islice,tsr_idx].*cos.(delta3D[islice])-cactmean)/cactmean*100
err_meanTpQ_AC_straight = (-Tp_AC_straight[islice,tsr_idx].*cos.(delta3D[islice])-cactmean)/cactmean*100
err_meanTpQ_DMS_straight = (-Tp_DMS_straight[islice,tsr_idx].*cos.(delta3D[islice])-cactmean)/cactmean*100

println("err_meanTpQ_AC_curved: $err_meanTpQ_AC_curved %")
println("err_meanTpQ_DMS_curved: $err_meanTpQ_DMS_curved %")
println("err_meanTpQ_AC_straight: $err_meanTpQ_AC_straight %")
println("err_meanTpQ_DMS_straight: $err_meanTpQ_DMS_straight %")

end
PyPlot.rc("figure.subplot", left=.12, bottom=.2, top=0.95, right=.95)
PyPlot.rc("figure", figsize=(5, 3.5))
PyPlot.figure()
R = 5.0/2
PyPlot.plot([0.0,0.0],[0.0,0.0],"-",color = plot_cycle[1])
# PyPlot.plot(te_cactus[:,tsr_idx].*(0.5*1.225*0.8*mean(Vinf)^2 * RefArea * R)./r3D./h_elem.*cos.(delta3D),h,color = plot_cycle[1])
# PyPlot.plot(-Tp_dist_cactus[:,tsr_idx],h,"-",color = plot_cycle[1])  #cactus gives torque over the element, not the unitized
PyPlot.plot(-Tp_AC_curved[:,tsr_idx].*cos.(delta3D),h,color = plot_cycle[2])
PyPlot.plot(-Tp_DMS_curved[:,tsr_idx].*cos.(delta3D),h,color = plot_cycle[3])
PyPlot.plot(-Tp_AC_straight[:,tsr_idx].*cos.(delta3D),h,"--",color = plot_cycle[2]) #multiply by cos(delta) to get,rid of element length from angle
PyPlot.plot(-Tp_DMS_straight[:,tsr_idx].*cos.(delta3D),h,"--",color = plot_cycle[3])
PyPlot.plot(te_cactus[:,tsr_idx].*(0.5*1.225*0.8*mean(Vinf)^2 * RefArea * R)./h_elem./r3D.*cos.(delta3D),h,color = plot_cycle[1])
# PyPlot.plot(shapeX,shapeY,"k")
# PyPlot.plot(-shapeX,shapeY,"k")

# Continuing, Plot slices cylinders
n_slices = 30
R = 1.0#5.0/2 #m
H = 1.02*R*2; #m
hr = H/R
n_slices = 30
shapeY1 = collect(LinRange(0,H,n_slices+1))
shapeX1 = R.*(1.0.-4.0.*(shapeY1/H.-.5).^2)

shapeY2 = collect(LinRange(0,H,Int(n_slices/2)+1))
shapeX2 = R.*(1.0.-4.0.*(shapeY2/H.-.5).^2)

shapeYslice = zeros(Int(n_slices/2)*2-1)
shapeXslice = zeros(Int(n_slices/2)*2-1)

shapeXslice[1:2:end] = shapeX2[1:end-1]
shapeXslice[2:2:end] = shapeX2[2:end-1]
shapeXslice = [shapeXslice;0.0]
shapeXslice[2:3] = shapeXslice[2:3].-0.02*2.6
shapeXslice[end-2:end-1] = shapeXslice[end-2:end-1].-0.02*2.6

shapeYslice[1:2:end] = shapeY2[1:end-1]
shapeYslice[2:2:end] = shapeY2[1:end-2]
shapeYslice = [shapeYslice;shapeY2[end-1]]
shapeYslice[1:2] = shapeYslice[1:2].-0.025*2.6
shapeYslice[3:4] = shapeYslice[3:4].-0.025/2*2.6
shapeYslice[end-1:end] = shapeYslice[end-1:end].+0.025*2.6
shapeYslice[end-3:end-2] = shapeYslice[end-3:end-2].+0.025/2*2.6
slice = 22
bottomslice = 1
lenarrow = 0.4

# Blade Shape
PyPlot.plot(-shapeXslice*2.5,(shapeYslice.+0.035*2).*2.5,"-",color ="0.5")
PyPlot.plot(-shapeX1[bottomslice:end]*2.5,shapeY1[bottomslice:end]*2.5,"-k")
PyPlot.plot(shapeXslice*2.5,(shapeYslice.+0.035*2).*2.5,"-",color ="0.5")
PyPlot.plot(shapeX1[bottomslice:end]*2.5,shapeY1[bottomslice:end]*2.5,"-k")

# for ii = 1:2:length(shapeXslice)
#     PyPlot.plot([-shapeXslice[ii], shapeXslice[ii]].*2.5,[shapeYslice[ii].+0.035*2, shapeYslice[ii].+0.035*2].*2.5,"--",color ="0.5")
# end

PyPlot.plot(0,7,"w.")
PyPlot.xlim([-5,25])
PyPlot.ylim([-2,8])
PyPlot.xlabel("Azimuthally Averaged\nTangential Force per span (N/m)")#"\n Lateral Location (m)")
PyPlot.ylabel("Vertical Location (m)")
PyPlot.legend(["CACTUS","AC Curved","DMS Curved","AC Cylinders","DMS Cylinders"])#, "VAWT Shape"])
PyPlot.axis("equal")
PyPlot.savefig("$(path)/figs/straight_comparison/SNL5m_shape_Tp_$(tsr_digit1)_$(tsr_digit2).pdf",transparent = true)


###############################
# Torque
PyPlot.rc("figure.subplot", left=.12, bottom=.2, top=0.95, right=.95)
PyPlot.rc("figure", figsize=(5, 3.5))
PyPlot.figure()
R = 5.0/2
PyPlot.plot(te_cactus[:,tsr_idx].*(0.5*1.225*0.8*mean(Vinf)^2 * RefArea * R)./h_elem.*cos.(delta3D),h,color = plot_cycle[1])
# PyPlot.plot(-Tp_dist_cactus[:,tsr_idx],h,"-",color = plot_cycle[1])  #cactus gives torque over the element, not the unitized
PyPlot.plot(-Tp_AC_curved[:,tsr_idx].*cos.(delta3D).*r3D,h,color = plot_cycle[2])
PyPlot.plot(-Tp_DMS_curved[:,tsr_idx].*cos.(delta3D).*r3D,h,color = plot_cycle[3])
PyPlot.plot(-Tp_AC_straight[:,tsr_idx].*cos.(delta3D).*r3D,h,"--",color = plot_cycle[2]) #multiply by cos(delta) to get,rid of element length from angle
PyPlot.plot(-Tp_DMS_straight[:,tsr_idx].*cos.(delta3D).*r3D,h,"--",color = plot_cycle[3])

# Blade Shape
PyPlot.plot(-shapeXslice*2.5,(shapeYslice.+0.035*2).*2.5,"-",color ="0.5")
PyPlot.plot(-shapeX1[bottomslice:end]*2.5,shapeY1[bottomslice:end]*2.5,"-k")
PyPlot.plot(shapeXslice*2.5,(shapeYslice.+0.035*2).*2.5,"-",color ="0.5")
PyPlot.plot(shapeX1[bottomslice:end]*2.5,shapeY1[bottomslice:end]*2.5,"-k")

# for ii = 1:2:length(shapeXslice)
#     PyPlot.plot([-shapeXslice[ii], shapeXslice[ii]].*2.5,[shapeYslice[ii].+0.035*2, shapeYslice[ii].+0.035*2].*2.5,"--",color ="0.5")
# end

PyPlot.plot(0,7,"w.")
# PyPlot.xlim([-5,25])
# PyPlot.ylim([-2,8])
PyPlot.xlabel("Mean Torque per Span (N/m)")#"\n Lateral Location (m)")
PyPlot.ylabel("Vertical Location (m)")
PyPlot.legend(["CACTUS","AC Curved","DMS Curved","AC Cylinders","DMS Cylinders"])#, "VAWT Shape"])
# PyPlot.axis("equal")
PyPlot.savefig("$(path)/figs/straight_comparison/SNL5m_shape_Qp_$(tsr_digit1)_$(tsr_digit2).pdf",transparent = true)


# Force normalizaion
# Torque Force
PyPlot.rc("figure.subplot", left=.12, bottom=.2, top=0.95, right=.95)
PyPlot.rc("figure", figsize=(5, 3.5))
PyPlot.figure()
R = 5.0/2
cact = te_cactus[:,tsr_idx].*(0.5*1.225*0.8*mean(Vinf)^2 * RefArea * R)./h_elem.*cos.(delta3D)
ac_cyl = -Tp_AC_straight[:,tsr_idx].*cos.(delta3D).*r3D
ac_cur = -Tp_AC_curved[:,tsr_idx].*cos.(delta3D).*r3D
dms_cur = -Tp_DMS_curved[:,tsr_idx].*cos.(delta3D).*r3D
dms_cyl = -Tp_DMS_straight[:,tsr_idx].*cos.(delta3D).*r3D
ac_cyl_error = abs.(cact.-ac_cyl)./abs.(cact).*100#.*cact./maximum(cact).*100
ac_cur_error = abs.(cact.-ac_cur)./abs.(cact).*100#.*cact./maximum(cact).*100
dms_cyl_error = abs.(cact.-dms_cyl)./abs.(cact).*100#.*cact./maximum(cact).*100
dms_cur_error = abs.(cact.-dms_cur)./abs.(cact).*100#.*cact./maximum(cact).*100

PyPlot.plot(abs.(ac_cur_error),h./h[end],color = plot_cycle[2],label="AC Curved")
PyPlot.plot(abs.(dms_cur_error),h./h[end],color = plot_cycle[3],label="DMS Curved")
PyPlot.plot(abs.(ac_cyl_error),h./h[end],"--",color = plot_cycle[2],label="AC Cylinder") #multiply by cos(delta) to get,rid of element length from angle
PyPlot.plot(abs.(dms_cyl_error),h./h[end],"--",color = plot_cycle[3],label="DMS Cylinder")

# Blade Shape
# PyPlot.plot(-shapeXslice*2.5,(shapeYslice.+0.035*2).*2.5,"-",color ="0.5")
# PyPlot.plot(-shapeX1[bottomslice:end]*2.5,shapeY1[bottomslice:end]*2.5,"-k")
# PyPlot.plot(shapeXslice*2.5,(shapeYslice.+0.035*2).*2.5,"-",color ="0.5")
# PyPlot.plot(shapeX1[bottomslice:end]*2.5,shapeY1[bottomslice:end]*2.5,"-k")

# for ii = 1:2:length(shapeXslice)
#     PyPlot.plot([-shapeXslice[ii], shapeXslice[ii]].*2.5,[shapeYslice[ii].+0.035*2, shapeYslice[ii].+0.035*2].*2.5,"--",color ="0.5")
# end

# PyPlot.plot(0,7,"w.")
PyPlot.xlim([0,50])
# PyPlot.ylim([-2,8])
PyPlot.xlabel("Average Torque Error (%)")#"\n Lateral Location (m)")
PyPlot.ylabel("Vertical Location (m)")
PyPlot.legend()
# PyPlot.axis("equal")
PyPlot.savefig("$(path)/figs/straight_comparison/SNL5m_shape_Tp_error$(tsr_digit1)_$(tsr_digit2).pdf",transparent = true)


# Tangential Force Surface Plot
meshgrid(x,y) = (repeat(x',length(y),1),repeat(y,1,length(x)))

#Cactus
q_loc2 = 0.5*rho.*(dataout[1][:,15,:].*Vinf).^2
PyPlot.figure()
X,Y=meshgrid(thetavec*180/pi,h/h[end])

CactTpplot = q_loc2.*data_full_rev[:,25,7]*chord.*r3D
levels = LinRange(minimum(-CactTpplot),maximum(-CactTpplot),20)
cp = PyPlot.contourf(X, Y, -CactTpplot', linestyles="solid",cmap="viridis",levels)
cb = PyPlot.colorbar(cp)
cb.set_label("Torque per height (N-m/m)")
PyPlot.xlabel("Azimuth Angle (deg)")
PyPlot.ylabel("Normalized Height")
# PyPlot.title("using PMSG_arms from GeneratorSE \n blank spots didn't converge")
PyPlot.savefig("$path/figs/straight_comparison/SNL5m_surface_CACTUS_Qp_$(tsr_digit1)_$(tsr_digit2).png",transparent = true)


#AC Curved
PyPlot.figure()
X,Y=meshgrid(thetavec*180/pi,h/h[end])
ACcurvTpplot = Tp_AC_curved_theta[:,:,tsr_idx].*cos.(delta3D[:]).*r3D
# levels = LinRange(minimum(-ACcurvTpplot),maximum(-ACcurvTpplot),20)
cp = PyPlot.contourf(X, Y, -ACcurvTpplot', linestyles="solid",cmap="viridis",levels)
cb = PyPlot.colorbar(cp)
cb.set_label("Torque per height (N-m/m)")
PyPlot.xlabel("Azimuth Angle (deg)")
PyPlot.ylabel("Normalized Height")
# PyPlot.title("using PMSG_arms from GeneratorSE \n blank spots didn't converge")
PyPlot.savefig("$path/figs/straight_comparison/SNL5m_surface_AC_Qp_$(tsr_digit1)_$(tsr_digit2).png",transparent = true)

#AC Straight
PyPlot.figure()
X,Y=meshgrid(thetavec*180/pi,h/h[end])
ACstraTpplot = Tp_AC_straight_theta[:,:,tsr_idx].*cos.(delta3D[:]).*r3D
# levels = LinRange(minimum(-ACstraTpplot),maximum(-ACstraTpplot),20)
cp = PyPlot.contourf(X, Y, -ACstraTpplot', linestyles="solid",cmap="viridis",levels)
cb = PyPlot.colorbar(cp)
cb.set_label("Torque Force per height (N/m)")
PyPlot.xlabel("Azimuth Angle (deg)")
PyPlot.ylabel("Normalized Height")
# PyPlot.title("using PMSG_arms from GeneratorSE \n blank spots didn't converge")
PyPlot.savefig("$path/figs/straight_comparison/SNL5m_surface_AC_straight_Qp_$(tsr_digit1)_$(tsr_digit2).png",transparent = true)

#AC difference
PyPlot.figure()
X,Y=meshgrid(thetavec*180/pi,h/h[end])
ACdiffplot = (Tp_AC_straight_theta[:,:,tsr_idx].*cos.(delta3D[:]).*r3D.-Tp_AC_curved_theta[:,:,tsr_idx].*cos.(delta3D[:]).*r3D)
levels = LinRange(minimum(ACdiffplot),maximum(ACdiffplot),20)
cp = PyPlot.contourf(X, Y, ACdiffplot', linestyles="solid",cmap="viridis",levels)
cb = PyPlot.colorbar(cp)
cb.set_label("Difference in Torque per height (N/m)")
PyPlot.xlabel("Azimuth Angle (deg)")
PyPlot.ylabel("Normalized Height")
# PyPlot.title("using PMSG_arms from GeneratorSE \n blank spots didn't converge")
PyPlot.savefig("$path/figs/straight_comparison/SNL5m_surface_AC_diff_Tp_$(tsr_digit1)_$(tsr_digit2).png",transparent = true)

#AC difference
PyPlot.figure()
X,Y=meshgrid(thetavec*180/pi,h/h[end])
TpCurvPlot = (CactTpplot.-ACcurvTpplot)
levels = LinRange(minimum(TpCurvPlot),maximum(TpCurvPlot),20)
cp = PyPlot.contourf(X, Y, TpCurvPlot', linestyles="solid",cmap="viridis",levels)
cb = PyPlot.colorbar(cp)
cb.set_label("Tangential Force per height (N/m)")
PyPlot.xlabel("Azimuth Angle (deg)")
PyPlot.ylabel("Normalized Height")
# PyPlot.title("using PMSG_arms from GeneratorSE \n blank spots didn't converge")
PyPlot.savefig("$path/figs/straight_comparison/SNL5m_surface_AC_Cact_diff_Qp_$(tsr_digit1)_$(tsr_digit2).png",transparent = true)


# Radial Force
PyPlot.rc("figure.subplot", left=.15, bottom=.25, top=0.95, right=.95)
PyPlot.figure()
PyPlot.plot(Rp_dist_cactus[:,tsr_idx]./10,h,color = plot_cycle[1])  #cactus gives torque over the element, not the unitized
PyPlot.plot(Rp_AC_straight[:,tsr_idx].*cos.(delta3D)./10,h,"--",color = plot_cycle[2])
PyPlot.plot(Rp_AC_curved[:,tsr_idx].*cos.(delta3D)./10,h,color = plot_cycle[2])
PyPlot.plot(Rp_DMS_straight[:,tsr_idx].*cos.(delta3D)./10,h,"--",color = plot_cycle[3])
PyPlot.plot(Rp_DMS_curved[:,tsr_idx].*cos.(delta3D)./10,h,color = plot_cycle[3])
PyPlot.plot(shapeX,shapeY,"k")
PyPlot.plot(-shapeX,shapeY,"k")
PyPlot.xlabel("Azimuthally Averaged Absolute\nRadial Force per span (N/m) /10\n Lateral Location (m)")
PyPlot.ylabel("Vertical Location (m)")
PyPlot.legend(["CACTUS","AC Straight","AC Curved","DMS Straight","DMS Curved", "VAWT Shape"])
PyPlot.axis("equal")
PyPlot.savefig("$(path)/figs/straight_comparison/SNL5m_shape_Rp_$(tsr_digit1)_$(tsr_digit2).pdf",transparent = true)

# Vertical Force
PyPlot.rc("figure.subplot", left=.15, bottom=.25, top=0.95, right=.95)
PyPlot.figure()
PyPlot.plot(abs.(Zp_dist_cactus[:,tsr_idx])/10,h,color = plot_cycle[1])  #cactus gives torque over the element, not the unitized
PyPlot.plot(Zp_AC_straight[:,tsr_idx].*cos.(delta3D)/10,h,"--",color = plot_cycle[2])
PyPlot.plot(Zp_AC_curved[:,tsr_idx].*cos.(delta3D)/10,h,color = plot_cycle[2])
PyPlot.plot(Zp_DMS_straight[:,tsr_idx].*cos.(delta3D)/10,h,"--",color = plot_cycle[3])
PyPlot.plot(Zp_DMS_curved[:,tsr_idx].*cos.(delta3D)/10,h,color = plot_cycle[3])
PyPlot.plot(shapeX,shapeY,"k")
PyPlot.plot(-shapeX,shapeY,"k")
PyPlot.xlabel("Azimuthally Averaged Absolute\nVertical Force per span (N/m) /10\n Lateral Location (m)")
PyPlot.ylabel("Vertical Location (m)")
PyPlot.legend(["CACTUS","AC Straight","AC Curved","DMS Straight","DMS Curved", "VAWT Shape"])
PyPlot.axis("equal")
PyPlot.savefig("$(path)/figs/straight_comparison/SNL5m_shape_Zp_$(tsr_digit1)_$(tsr_digit2).pdf",transparent = true)

for islice = 1:n_slices
    # Tangential Force ALL
    PyPlot.figure()
    # q_loc = 0.5*rho.*(data_full_rev[:,15,15].*Vinf).^2
    # PyPlot.plot((data_full_rev[:,2,15].-minimum(data_full_rev[:,2,15]))*180/pi, q_loc.*data_full_rev[:,25,15]*chord,"k")
    q_loc = 0.5*rho.*(data_full_rev[:,15,islice].*Vinf).^2
    PyPlot.plot((data_full_rev[:,2,islice].-minimum(data_full_rev[:,2,islice]))*180/pi, -q_loc.*data_full_rev[:,25,islice]*chord,color=plot_cycle[1],label="Cactus")
    PyPlot.plot(thetavec*180/pi, -Tp_AC_curved_theta[:,islice,tsr_idx].*cos.(delta3D[islice]),color=plot_cycle[2],label="AC Curved") #remove relative span len2th by multiplying by cos(delta)
    # PyPlot.plot(thetavec*180/pi, -Tp_DMS_curved_theta[:,islice,tsr_idx].*cos.(delta3D[islice]),color=plot_cycle[3],label="DMS Curved")
    PyPlot.plot(thetavec*180/pi, -Tp_AC_straight_theta[:,islice,tsr_idx].*cos.(delta3D[islice]),"--",color=plot_cycle[2],label="AC Cylinders")
    # PyPlot.plot(thetavec*180/pi, -Tp_DMS_straight_theta[:,islice,tsr_idx].*cos.(delta3D[islice]),"--",color=plot_cycle[3],label="DMS Cylinders")
    PyPlot.xlabel("Theta (deg)")
    # Blade Shape
    scalex = 365/12.0
    scaley = 1.3
    offsetx = 340.0
    offsety = 6.5
    if islice==7
        PyPlot.plot(-shapeXslice*scalex.+offsetx,(shapeYslice.+0.035*2).*scaley.+offsety,"-",color ="0.5")
        PyPlot.plot(-shapeX1[bottomslice:end]*scalex.+offsetx,shapeY1[bottomslice:end]*scaley.+offsety,"-k")
        PyPlot.plot(shapeXslice*scalex.+offsetx,(shapeYslice.+0.035*2).*scaley.+offsety,"-",color ="0.5")
        PyPlot.plot(shapeX1[bottomslice:end]*scalex.+offsetx,shapeY1[bottomslice:end]*scaley.+offsety,"-k")

        for ii = [6,8]#1:2:length(shapeXslice)
            PyPlot.plot(([-shapeXslice[ii], shapeXslice[ii]])*scalex.+offsetx,([shapeYslice[ii], shapeYslice[ii]])*scaley.+offsety,"-",color=color=plot_cycle[2])
        end
    end
    PyPlot.ylabel("Tangential Force per span (N/m)")
    PyPlot.legend(loc=1,bbox_to_anchor=(0.84,0.97))
    PyPlot.savefig("$(path)/figs/straight_comparison/SNL5m_Tangential_Force_slice_$(islice)_ALL_TSR_$(tsr_digit1)_$(tsr_digit2).pdf",transparent = true)

    # Radial Force ALL
    PyPlot.figure()
    # q_loc = 0.5*rho.*(data_full_rev[:,15,15].*Vinf).^2
    # PyPlot.plot((data_full_rev[:,2,15].-minimum(data_full_rev[:,2,15]))*180/pi, q_loc.*data_full_rev[:,24,15]*chord.*cos(delta3D[15]),"k") #Normal to vertical by cos(delta)
    q_loc = 0.5*rho.*(data_full_rev[:,15,islice].*Vinf).^2
    PyPlot.plot((data_full_rev[:,2,islice].-minimum(data_full_rev[:,2,islice]))*180/pi, q_loc.*data_full_rev[:,24,islice]*chord.*cos(delta3D[islice]),color=plot_cycle[1],label="Cactus")
    PyPlot.plot(thetavec*180/pi, Rp_AC_curved_theta[:,islice,tsr_idx].*cos.(delta3D[islice]),color=plot_cycle[2],label="AC Curved")
    # PyPlot.plot(thetavec*180/pi, Rp_DMS_curved_theta[:,islice,tsr_idx].*cos.(delta3D[islice]),color=plot_cycle[3],label="DMS Curved")
    PyPlot.plot(thetavec*180/pi, Rp_AC_straight_theta[:,islice,tsr_idx].*cos.(delta3D[islice]),"--",color=plot_cycle[2],label="AC Cylinders")
    # PyPlot.plot(thetavec*180/pi, Rp_DMS_straight_theta[:,islice,tsr_idx].*cos.(delta3D[islice]),"--",color=plot_cycle[3],label="DMS Cylinders")
    PyPlot.xlabel("Theta (deg)")
    PyPlot.ylabel("Radial Force per span (N/m)")
    scalex = 365/12.0
    scaley = 9.0
    offsetx = 345.0
    offsety = 30.0
    if islice==7
        PyPlot.plot(-shapeXslice*scalex.+offsetx,(shapeYslice.+0.035*2).*scaley.+offsety,"-",color ="0.5")
        PyPlot.plot(-shapeX1[bottomslice:end]*scalex.+offsetx,shapeY1[bottomslice:end]*scaley.+offsety,"-k")
        PyPlot.plot(shapeXslice*scalex.+offsetx,(shapeYslice.+0.035*2).*scaley.+offsety,"-",color ="0.5")
        PyPlot.plot(shapeX1[bottomslice:end]*scalex.+offsetx,shapeY1[bottomslice:end]*scaley.+offsety,"-k")

        for ii = [6,8]#1:2:length(shapeXslice)
            PyPlot.plot(([-shapeXslice[ii], shapeXslice[ii]])*scalex.+offsetx,([shapeYslice[ii], shapeYslice[ii]])*scaley.+offsety,"-",color=color=plot_cycle[2])
        end
    end
    PyPlot.legend(loc=1,bbox_to_anchor=(0.83,0.96))
    PyPlot.savefig("$(path)/figs/straight_comparison/SNL5m_Radial_Force_slice_$(islice)_ALL_TSR_$(tsr_digit1)_$(tsr_digit2).pdf",transparent = true)

    # Vertical Force ALL
    PyPlot.figure()
    # q_loc = 0.5*rho.*(data_full_rev[:,15,15].*Vinf).^2
    # PyPlot.plot((data_full_rev[:,2,15].-minimum(data_full_rev[:,2,15]))*180/pi, q_loc.*data_full_rev[:,24,15]*chord.*sin(delta3D[15]),"k") #Normal to vertical by sin(delta)
    q_loc = 0.5*rho.*(data_full_rev[:,15,islice].*Vinf).^2
    PyPlot.plot((data_full_rev[:,2,islice].-minimum(data_full_rev[:,2,islice]))*180/pi, q_loc.*data_full_rev[:,24,islice]*chord.*sin(delta3D[islice]),color=plot_cycle[1],label="Cactus")
    PyPlot.plot(thetavec*180/pi, Zp_AC_curved_theta[:,islice,tsr_idx].*cos.(delta3D[islice]),color=plot_cycle[2],label="AC Curved")
    # PyPlot.plot(thetavec*180/pi, Zp_DMS_curved_theta[:,islice,tsr_idx].*cos.(delta3D[islice]),color=plot_cycle[3],label="DMS Curved")
    PyPlot.plot(thetavec*180/pi, Zp_AC_straight_theta[:,islice,tsr_idx].*cos.(delta3D[islice]),"--",color=plot_cycle[2],label="AC Cylinders")
    # PyPlot.plot(thetavec*180/pi, Zp_DMS_straight_theta[:,islice,tsr_idx].*cos.(delta3D[islice]),"--",color=plot_cycle[3],label="DMS Cylinders")
    PyPlot.xlabel("Theta (deg)")
    PyPlot.ylabel("Vertical Force per span (N/m)")
    scalex = 365/12.0
    scaley = 8.0
    offsetx = 340.0
    offsety = 23.0
    if islice==7
        PyPlot.plot(-shapeXslice*scalex.+offsetx,(shapeYslice.+0.035*2).*scaley.+offsety,"-",color ="0.5")
        PyPlot.plot(-shapeX1[bottomslice:end]*scalex.+offsetx,shapeY1[bottomslice:end]*scaley.+offsety,"-k")
        PyPlot.plot(shapeXslice*scalex.+offsetx,(shapeYslice.+0.035*2).*scaley.+offsety,"-",color ="0.5")
        PyPlot.plot(shapeX1[bottomslice:end]*scalex.+offsetx,shapeY1[bottomslice:end]*scaley.+offsety,"-k")

        for ii = [6,8]#1:2:length(shapeXslice)
            PyPlot.plot(([-shapeXslice[ii], shapeXslice[ii]])*scalex.+offsetx,([shapeYslice[ii], shapeYslice[ii]])*scaley.+offsety,"-",color=color=plot_cycle[2])
        end
    end
    PyPlot.legend(loc=1,bbox_to_anchor=(0.83,.95))
    PyPlot.savefig("$(path)/figs/straight_comparison/SNL5m_Vertical_Force_slice_$(islice)_ALL_TSR_$(tsr_digit1)_$(tsr_digit2).pdf",transparent = true)
end

# end

# PyPlot.figure()
# PyPlot.plot(experimental[1:end-1,1],experimental[1:end-1,2],"k.")
# PyPlot.plot(tsrvec,cactus_CP,color = plot_cycle[1])
# PyPlot.plot(tsrvec,CPvec_DMS_straight,"--",color = plot_cycle[2])
# PyPlot.plot(tsrvec,CPvec_DMS_curved,color = plot_cycle[2])
# PyPlot.plot(tsrvec,CPvec_AC_straight,"--",color = plot_cycle[3])
# PyPlot.plot(tsrvec,CPvec_AC_curved,color = plot_cycle[3])
# # PyPlot.plot(tsrvec,CP_total_dist_cactus2,color = plot_cycle[1])
# # PyPlot.plot(tsrvec,CP_total_dist_cactus3,".-",color = plot_cycle[1])
# # PyPlot.plot(tsrvec,CP_total_AC_straight,color = plot_cycle[4])
# # PyPlot.plot(tsrvec,CP_total_AC_curved,color = plot_cycle[4])
# # PyPlot.plot(tsrvec,CP_total_DMS_straight,color = plot_cycle[5])
# # PyPlot.plot(tsrvec,CP_total_DMS_curved,color = plot_cycle[5])
# PyPlot.xlabel("Tip Speed Ratio")
# PyPlot.ylabel("Coefficient of Performance")
# PyPlot.legend(["Experimental", "Cactus","DMS Straight","DMS Curved","AC Straight","AC Curved","Cactus CT","Cactus te","AC Straight","AC Curved","DMS Straight","DMS Curved"])
# PyPlot.savefig("$(path)/figs/straight_comparison/SNL5m_CP3D.pdf",transparent = true)
