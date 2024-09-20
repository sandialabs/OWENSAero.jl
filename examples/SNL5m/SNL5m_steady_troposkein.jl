
import PyPlot
PyPlot.close("all")
PyPlot.ion()
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

using Test
import Statistics:mean
import DelimitedFiles
import Dierckx
import QuadGK
import FLOWMath
import HDF5
import OWENSAero

path,_ = splitdir(@__FILE__)
# include("$(path)/../../../src/OWENSAero.jl")


function runme(A_model,DS_model,r_delta_infl)
    # Main Function
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

    element_planf_A = sqrt.(delta_xs.^2+delta_zs.^2)*chord
    element_planf_L = sqrt.(delta_xs.^2+delta_zs.^2)

    r3D = (shapeX[2:end,1]+shapeX[1:end-1,1])/2.0
    aerocenter_dist = (eta-.25)*chord

    twist3D = -atan.(aerocenter_dist./r3D)#ones(n_slices)*-0.4*pi/180

    if DS_model=="BV"
        af3D = OWENSAero.readaerodyn_BV("$(path)/airfoils/NACA_0015_RE3E5.dat")
    else
        af3D = OWENSAero.readaerodyn("$(path)/airfoils/NACA_0015_RE3E5.dat")
    end

    # Unchanging Parameters
    B = 3
    omega = -ones(ntheta) * 150.0 / 60.0 * 2*pi # RPM -> radians/sec
    rho = 1.225*0.8 #80percent at albuquerque elevation
    mu = 1.7894e-5


    tsrvec = [2.1,3.1,4.2,4.7,5.2,6.0,7.1] #FOR PLOTTING
    # tsrvec = [5.2] #FOR TESTING
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

        tsr_digit1 = 5#Int(floor(tsrvec[tsr_idx]))
        tsr_digit2 = 2#Int(round((tsrvec[tsr_idx]-tsr_digit1)*10.0))
        Vinf = abs.(omega)/tsrvec[tsr_idx]*R

        cactus_Rev_data = DelimitedFiles.readdlm("$(path)/data/steady/TestVAWT$(tsr_digit1)_$(tsr_digit2)_troposkein_RevData.csv",',',Float64,skipstart = 1)
        cactus_CP[tsr_idx] = cactus_Rev_data[end,2]

        data = DelimitedFiles.readdlm("$(path)/data/steady/TestVAWT$(tsr_digit1)_$(tsr_digit2)_troposkein_ElementData.csv", ',',skipstart = 1)
        #azimuthal discretization, data columns, slices
        data_full_rev = zeros(ntheta,length(data[1,:]),Int(maximum(data[:,4])))
        elapsed = zeros(n_slices)
        for slice = 1:n_slices
            # start = time()
            # println("$i : $slice")
            # Potentially Aerostructural Parameters that my change
            r = ones(ntheta)*r3D[slice] #m
            twist = ones(ntheta)*twist3D[slice] #rad
            delta = ones(ntheta)*delta3D[slice] #rad
            af = af3D

            turbine = OWENSAero.Turbine(R,r,chord,twist,delta,omega,B,af,ntheta,r_delta_infl)
            if slice == 1 && tsr_idx != 1
                aw_warm = w[:,slice,tsr_idx-1]
            elseif tsr_idx == 1
                aw_warm = zeros(ntheta*2)
            else
                aw_warm = w[:,slice-1,tsr_idx]
            end
            env = OWENSAero.Environment(rho,mu,Vinf,DS_model,A_model,aw_warm)

            # Option to Plot the model
            # OWENSAero.plot_domain(turbine)

            # Juno.@enter OWENSAero.DMS(turbine, env, Vinf)
            # CP[slice,i], _, _, Rp[:,slice,i], Tp[:,slice,i], Zp[:,slice,i], _, _, CT[slice,i], _, alpha[:,slice,i], Vinf_used[:,slice,i], _, _, thetavec
            start = time()

            CP_in, _, _, Rp_in, Tp_in, Zp_in,W[:,slice,tsr_idx], _, _, _, w_in, alpha_in,
            cl_in, cd_in, thetavec[:], Re_in = OWENSAero.steady(turbine, env)

            CP[slice,tsr_idx] = CP_in[1]
            Tp[:,slice,tsr_idx] = Tp_in
            Rp[:,slice,tsr_idx] = Rp_in
            Zp[:,slice,tsr_idx] = Zp_in
            alpha[:,slice,tsr_idx] = alpha_in
            cl[:,slice,tsr_idx] = cl_in
            cd[:,slice,tsr_idx] = cd_in
            Re[:,slice,tsr_idx] = Re_in
            w[1:ntheta,slice,tsr_idx] = w_in[1:ntheta]
            Rp_out[slice,tsr_idx] = B/(2*pi)*OWENSAero.pInt(thetavec, abs.(Rp[:,slice,tsr_idx])) #average absolute force over the revolution
            Tp_out[slice,tsr_idx] = B/(2*pi)*OWENSAero.pInt(thetavec, Tp_in) #average NON-absolute force over the revolution
            Zp_out[slice,tsr_idx] = B/(2*pi)*OWENSAero.pInt(thetavec, abs.(Zp[:,slice,tsr_idx])) #average absolute force over the revolution

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

            P_cactus = abs(mean(abs.(omega)))/(2*pi)*OWENSAero.pInt(thetavec, Qp_cactus_b1+Qp_cactus_b2+Qp_cactus_b3)

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
        CPvec[tsr_idx] = OWENSAero.trapz(h,CP[:,tsr_idx].*2*R./RefArea) #Undo local normalization and normalize by full turbine

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
println("Starting")
r_delta_infl = false
start = time()
delta3D,experimental,cactus_CP,tsrvec,CPvec_DMS_noDS,shapeX,shapeY,element_planf_L,Rp_DMS_noDS,Tp_DMS_noDS,Zp_DMS_noDS,Cp_1DMS_noDS,h,h_frac,Rp_dist_cactus,Tp_dist_cactus,Zp_dist_cactus,CP_dist_cactus,te_cactus,r3D,R,omega,B,thetavec,RefArea,Tp_DMS_noDS_theta,Rp_DMS_noDS_theta,Zp_DMS_noDS_theta,dataout,tsrvec = runme("DMS","noDS",r_delta_infl)
elapsed_DMS_noDS = time() - start
println("elapsed_DMS_noDS $elapsed_DMS_noDS")
start = time()
delta3D,experimental,cactus_CP,tsrvec,CPvec_DMS_BV,shapeX,shapeY,element_planf_L,Rp_DMS_BV,Tp_DMS_BV,Zp_DMS_BV,Cp_1DMS_BV,h,h_frac,Rp_dist_cactus,Tp_dist_cactus,Zp_dist_cactus,CP_dist_cactus,te_cactus,r3D,R,omega,B,thetavec,RefArea,Tp_DMS_BV_theta,Rp_DMS_BV_theta,Zp_DMS_BV_theta,dataout,tsrvec = runme("DMS","BV",r_delta_infl)
elapsed_DMS_BV = time() - start
println("elapsed_DMS_BV $elapsed_DMS_BV")
start = time()
delta3D,experimental,cactus_CP,tsrvec,CPvec_AC_noDS,shapeX,shapeY,element_planf_L,Rp_AC_noDS,Tp_AC_noDS,Zp_AC_noDS,Cp_1AC_noDS,h,h_frac,Rp_dist_cactus,Tp_dist_cactus,Zp_dist_cactus,CP_dist_cactus,te_cactus,r3D,R,omega,B,thetavec,RefArea,Tp_AC_noDS_theta,Rp_AC_noDS_theta,Zp_AC_noDS_theta,dataout,tsrvec = runme("AC","noDS",r_delta_infl)
elapsed_AC_noDS = time() - start
println("elapsed_AC_noDS $elapsed_AC_noDS")
start = time()
delta3D,experimental,cactus_CP,tsrvec,CPvec_AC_BV,shapeX,shapeY,element_planf_L,Rp_AC_BV,Tp_AC_BV,Zp_AC_BV,Cp_1AC_BV,h,h_frac,Rp_dist_cactus,Tp_dist_cactus,Zp_dist_cactus,CP_dist_cactus,te_cactus,r3D,R,omega,B,thetavec,RefArea,Tp_AC_BV_theta,Rp_AC_BV_theta,Zp_AC_BV_theta,dataout,tsrvec = runme("AC","BV",r_delta_infl)
elapsed_AC_BV = time() - start
println("elapsed_AC_BV $elapsed_AC_BV")

QCx = [-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00,-0.00000e+00 ]
QCy = [ 0.00000e+00, 6.80000e-02, 1.36000e-01, 2.04000e-01, 2.72000e-01, 3.40000e-01, 4.08000e-01, 4.76000e-01, 5.44000e-01, 6.12000e-01, 6.80000e-01, 7.48000e-01, 8.16000e-01, 8.84000e-01, 9.52000e-01, 1.02000e+00, 1.08800e+00, 1.15600e+00, 1.22400e+00, 1.29200e+00, 1.36000e+00, 1.42800e+00, 1.49600e+00, 1.56400e+00, 1.63200e+00, 1.70000e+00, 1.76800e+00, 1.83600e+00, 1.90400e+00, 1.97200e+00, 2.04000e+00 ]
QCz = [-1.08950e-02,-1.39191e-01,-2.41913e-01,-3.43988e-01,-4.48635e-01,-5.49196e-01,-6.50442e-01,-7.36009e-01,-8.06550e-01,-8.66658e-01,-9.12919e-01,-9.44895e-01,-9.70493e-01,-9.87413e-01,-9.98089e-01,-1.00000e+00,-9.97602e-01,-9.89086e-01,-9.72320e-01,-9.51171e-01,-9.18891e-01,-8.80189e-01,-8.25827e-01,-7.60379e-01,-6.66251e-01,-5.51174e-01,-4.41618e-01,-3.32788e-01,-2.22668e-01,-1.13861e-01,-1.77009e-03 ]
ECtoR = [6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02,6.09570e-02 ]
EAreaR = [8.88308e-03   8.40765e-03   7.94073e-03   7.48392e-03   7.03918e-03   6.60895e-03   6.19626e-03   5.80484e-03   5.43929e-03   5.10517e-03   4.80904e-03   4.55831e-03   4.36080e-03   4.22401e-03   4.15392e-03   4.15392e-03   4.22401e-03   4.36080e-03   4.55831e-03   4.80904e-03   5.10517e-03   5.43929e-03   5.80484e-03   6.19626e-03   6.60895e-03   7.03918e-03   7.48392e-03   7.94073e-03   8.40765e-03   8.88308e-03 ]
myEAreaR = (element_planf_L*0.1524)./(R^2)

PyPlot.figure()
PyPlot.plot(EAreaR',h)
PyPlot.plot(myEAreaR,h)
PyPlot.plot(diff(shapeY)./cos.(delta3D)*0.1524./(R^2),h)
PyPlot.legend(["CACTUS EAreaR", "My EAreaR", "test"])
PyPlot.savefig("$(path)/figs/troposkein/EAreaR.pdf",transparent = true)

PyPlot.figure()
PyPlot.plot(QCz,QCy,".-",color=plot_cycle[1])
PyPlot.plot(-shapeX./R,shapeY./R,".-",color=plot_cycle[2])
PyPlot.legend(["SNL 5m Image", "Parabolic"])
PyPlot.savefig("$(path)/figs/troposkein/QCzy.pdf",transparent = true)

h_elem = (shapeY[2:end] - shapeY[1:end-1])

CP_total_dist_cactus2 = zeros(length(tsrvec))
CP_total_AC_noDS = zeros(length(tsrvec))
CP_total_AC_BV = zeros(length(tsrvec))
CP_total_DMS_noDS = zeros(length(tsrvec))
CP_total_DMS_BV = zeros(length(tsrvec))
CP_total_dist_cactus3 = zeros(length(tsrvec))
chord = 0.1524
rho = 1.225*0.8
println("Aggregate Plots")
for tsr_idx = 1:length(tsrvec)
    println("TSR: $(tsrvec[tsr_idx])")
    local tsr_digit1 = Int(floor(tsrvec[tsr_idx]))
    local tsr_digit2 = Int(round((tsrvec[tsr_idx]-tsr_digit1)*10.0))
    data_full_rev = dataout[tsr_idx]
    local Vinf = abs.(omega)/tsrvec[tsr_idx]*R

    # Tangential Force ALL
    PyPlot.figure()

    q_loc = 0.5*rho.*(data_full_rev[:,15,15].*Vinf).^2
    cactus_15_T = -q_loc.*data_full_rev[:,25,15]*chord

    PyPlot.plot((data_full_rev[:,2,15].-minimum(data_full_rev[:,2,15]))*180/pi,cactus_15_T ,"-",color = plot_cycle[1])

    q_loc = 0.5*rho.*(data_full_rev[:,15,7].*Vinf).^2
    cactus_7_T = -q_loc.*data_full_rev[:,25,7]*chord

    PyPlot.plot((data_full_rev[:,2,7].-minimum(data_full_rev[:,2,7]))*180/pi, cactus_7_T,"--",color = plot_cycle[1])
    PyPlot.plot(thetavec*180/pi, -Tp_AC_BV_theta[:,15,tsr_idx].*cos.(delta3D[15]),color=plot_cycle[2]) #remove relative span length by multiplying by cos(delta)
    PyPlot.plot(thetavec*180/pi, -Tp_AC_BV_theta[:,7,tsr_idx].*cos.(delta3D[7]),"--",color=plot_cycle[2])
    PyPlot.plot(thetavec*180/pi, -Tp_DMS_BV_theta[:,15,tsr_idx].*cos.(delta3D[15]),color=plot_cycle[3])
    PyPlot.plot(thetavec*180/pi, -Tp_DMS_BV_theta[:,7,tsr_idx].*cos.(delta3D[7]),"--",color=plot_cycle[3])
    PyPlot.xlabel("Theta (deg)")
    PyPlot.ylabel("Tangential Force per span (N/m)")
    PyPlot.legend(["Cactus 49% height","Cactus 22% height","AC 49% height","AC 22% height","DMS 49% height","DMS 22% height"])
    PyPlot.savefig("$(path)/figs/troposkein/SNL5m_Tangential_Force_ALL_TSR_$(tsr_digit1)_$(tsr_digit2)_rdelta_$(r_delta_infl).pdf",transparent = true)

    println("Tangential")
    # Error 15
    twod_15 = -Tp_DMS_BV_theta[:,15,tsr_idx].*cos.(delta3D[15])
    maxerr = maximum(abs.((cactus_15_T .- twod_15)./cactus_15_T))
    meanerr = mean(abs.((cactus_15_T .- twod_15)./cactus_15_T))
    println("Max at 15: $(maxerr)")
    println("Mean at 15: $(meanerr)")

    # Error 7
    twod_7 = -Tp_DMS_BV_theta[:,7,tsr_idx].*cos.(delta3D[7])
    maxerr = maximum(abs.((cactus_7_T .- twod_7)./cactus_7_T))
    meanerr = mean(abs.((cactus_7_T .- twod_7)./cactus_7_T))
    println("Max at 7: $(maxerr)")
    println("Mean at 7: $(meanerr)")

    # Radial Force ALL
    PyPlot.figure()
    q_loc = 0.5*rho.*(data_full_rev[:,15,15].*Vinf).^2
    cactus_15_R = q_loc.*data_full_rev[:,24,15]*chord.*cos(delta3D[15])

    PyPlot.plot((data_full_rev[:,2,15].-minimum(data_full_rev[:,2,15]))*180/pi, cactus_15_R,"-",color = plot_cycle[1]) #Normal to vertical by cos(delta)
    q_loc = 0.5*rho.*(data_full_rev[:,15,7].*Vinf).^2
    cactus_7_R = q_loc.*data_full_rev[:,24,7]*chord.*cos(delta3D[7])

    PyPlot.plot((data_full_rev[:,2,7].-minimum(data_full_rev[:,2,7]))*180/pi, cactus_7_R,"--",color = plot_cycle[1])
    PyPlot.plot(thetavec*180/pi, Rp_AC_BV_theta[:,15,tsr_idx].*cos.(delta3D[15]),color=plot_cycle[2])
    PyPlot.plot(thetavec*180/pi, Rp_AC_BV_theta[:,7,tsr_idx].*cos.(delta3D[7]),"--",color=plot_cycle[2])
    PyPlot.plot(thetavec*180/pi, Rp_DMS_BV_theta[:,15,tsr_idx].*cos.(delta3D[15]),color=plot_cycle[3])
    PyPlot.plot(thetavec*180/pi, Rp_DMS_BV_theta[:,7,tsr_idx].*cos.(delta3D[7]),"--",color=plot_cycle[3])
    PyPlot.xlabel("Theta (deg)")
    PyPlot.ylabel("Radial Force per span (N/m)")
    PyPlot.legend(["Cactus 49% height","Cactus 22% height","AC 49% height","AC 22% height","DMS 49% height","DMS 22% height"])
    PyPlot.savefig("$(path)/figs/troposkein/SNL5m_Radial_Force_ALL_TSR_$(tsr_digit1)_$(tsr_digit2)_rdelta_$(r_delta_infl).pdf",transparent = true)

    println("Radial")
    # Error 15
    twod_15 = Rp_DMS_BV_theta[:,15,tsr_idx].*cos.(delta3D[15])
    maxerr = maximum(abs.((cactus_15_R .- twod_15)./cactus_15_R))
    meanerr = mean(abs.((cactus_15_R .- twod_15)./cactus_15_R))
    println("Max at 15: $(maxerr)")
    println("Mean at 15: $(meanerr)")

    # Error 7
    twod_7 = Rp_DMS_BV_theta[:,7,tsr_idx].*cos.(delta3D[7])
    maxerr = maximum(abs.((cactus_7_R .- twod_7)./cactus_7_R))
    meanerr = mean(abs.((cactus_7_R .- twod_7)./cactus_7_R))
    println("Max at 7: $(maxerr)")
    println("Mean at 7: $(meanerr)")

    # Vertical Force ALL
    PyPlot.figure()
    q_loc = 0.5*rho.*(data_full_rev[:,15,15].*Vinf).^2
    cactus_15_V = q_loc.*data_full_rev[:,24,15]*chord.*sin(delta3D[15])

    PyPlot.plot((data_full_rev[:,2,15].-minimum(data_full_rev[:,2,15]))*180/pi, cactus_15_V,"-",color = plot_cycle[1]) #Normal to vertical by sin(delta)
    q_loc = 0.5*rho.*(data_full_rev[:,15,7].*Vinf).^2
    cactus_7_V = q_loc.*data_full_rev[:,24,7]*chord.*sin(delta3D[7])

    PyPlot.plot((data_full_rev[:,2,7].-minimum(data_full_rev[:,2,7]))*180/pi, cactus_7_V,"--",color = plot_cycle[1])
    PyPlot.plot(thetavec*180/pi, Zp_AC_BV_theta[:,15,tsr_idx].*cos.(delta3D[15]),color=plot_cycle[2])
    PyPlot.plot(thetavec*180/pi, Zp_AC_BV_theta[:,7,tsr_idx].*cos.(delta3D[7]),"--",color=plot_cycle[2])
    PyPlot.plot(thetavec*180/pi, Zp_DMS_BV_theta[:,15,tsr_idx].*cos.(delta3D[15]),color=plot_cycle[3])
    PyPlot.plot(thetavec*180/pi, Zp_DMS_BV_theta[:,7,tsr_idx].*cos.(delta3D[7]),"--",color=plot_cycle[3])
    PyPlot.xlabel("Theta (deg)")
    PyPlot.ylabel("Vertical Force per span (N/m)")
    PyPlot.legend(["Cactus 49% height","Cactus 22% height","AC 49% height","AC 22% height","DMS 49% height","DMS 22% height"])
    PyPlot.savefig("$(path)/figs/troposkein/SNL5m_Vertical_Force_ALL_TSR_$(tsr_digit1)_$(tsr_digit2)_rdelta_$(r_delta_infl).pdf",transparent = true)

    println("Vertical")
    # Error 15
    twod_15 = Zp_DMS_BV_theta[:,15,tsr_idx].*cos.(delta3D[15])
    maxerr = maximum(abs.((cactus_15_V .- twod_15)./cactus_15_V))
    meanerr = mean(abs.((cactus_15_V .- twod_15)./cactus_15_V))
    println("Max at 15: $(maxerr)")
    println("Mean at 15: $(meanerr)")

    # Error 7
    twod_7 = Zp_DMS_BV_theta[:,7,tsr_idx].*cos.(delta3D[7])
    maxerr = maximum(abs.((cactus_7_V .- twod_7)./cactus_7_V))
    meanerr = mean(abs.((cactus_7_V .- twod_7)./cactus_7_V))
    println("Max at 7: $(maxerr)")
    println("Mean at 7: $(meanerr)")


    # Tangential Force
    PyPlot.rc("figure.subplot", left=.12, bottom=.22, top=0.95, right=.95)
    PyPlot.rc("figure", figsize=(5, 3.5))
    PyPlot.figure()
    PyPlot.plot(te_cactus[:,tsr_idx].*(0.5*1.225*0.8*mean(Vinf)^2 * RefArea * R)./r3D./h_elem.*cos.(delta3D),h,color = plot_cycle[1])
    # PyPlot.plot(-Tp_dist_cactus[:,tsr_idx],h,"-",color = plot_cycle[1])  #cactus gives torque over the element, not the unitized
    PyPlot.plot(-Tp_AC_noDS[:,tsr_idx].*cos.(delta3D),h,"--",color = plot_cycle[2]) #multiply by cos(delta) to get rid of element length from angle
    PyPlot.plot(-Tp_AC_BV[:,tsr_idx].*cos.(delta3D),h,color = plot_cycle[2])
    PyPlot.plot(-Tp_DMS_noDS[:,tsr_idx].*cos.(delta3D),h,"--",color = plot_cycle[3])
    PyPlot.plot(-Tp_DMS_BV[:,tsr_idx].*cos.(delta3D),h,color = plot_cycle[3])
    PyPlot.plot(shapeX,shapeY,"k")
    PyPlot.plot(-shapeX,shapeY,"k")
    PyPlot.plot(0,7,"w.")
    PyPlot.xlabel("Azimuthally Averaged\nTangential Force per span (N/m)")#"\n Lateral Location (m)")
    PyPlot.ylabel("Vertical Location (m)")
    PyPlot.legend(["CACTUS","AC no DS","AC BV","DMS no DS","DMS BV"])#, "VAWT Shape"])
    PyPlot.axis("equal")
    PyPlot.savefig("$(path)/figs/troposkein/SNL5m_shape_Tp_$(tsr_digit1)_$(tsr_digit2)_rdelta_$(r_delta_infl).pdf",transparent = true)

    # Radial Force
    PyPlot.rc("figure.subplot", left=.15, bottom=.25, top=0.95, right=.95)
    PyPlot.figure()
    PyPlot.plot(Rp_dist_cactus[:,tsr_idx]./10,h,color = plot_cycle[1])  #cactus gives torque over the element, not the unitized
    PyPlot.plot(Rp_AC_noDS[:,tsr_idx].*cos.(delta3D)./10,h,color = plot_cycle[2])
    PyPlot.plot(Rp_AC_BV[:,tsr_idx].*cos.(delta3D)./10,h,"--",color = plot_cycle[2])
    PyPlot.plot(Rp_DMS_noDS[:,tsr_idx].*cos.(delta3D)./10,h,color = plot_cycle[3])
    PyPlot.plot(Rp_DMS_BV[:,tsr_idx].*cos.(delta3D)./10,h,"--",color = plot_cycle[3])
    PyPlot.plot(shapeX,shapeY,"k")
    PyPlot.plot(-shapeX,shapeY,"k")
    PyPlot.xlabel("Azimuthally Averaged Absolute\nRadial Force per span (N/m) /10\n Lateral Location (m)")
    PyPlot.ylabel("Vertical Location (m)")
    PyPlot.legend(["CACTUS","AC no DS","AC BV","DMS no DS","DMS BV", "VAWT Shape"])
    PyPlot.axis("equal")
    PyPlot.savefig("$(path)/figs/troposkein/SNL5m_shape_Rp_$(tsr_digit1)_$(tsr_digit2)_rdelta_$(r_delta_infl).pdf",transparent = true)

    # Vertical Force
    PyPlot.rc("figure.subplot", left=.15, bottom=.25, top=0.95, right=.95)
    PyPlot.figure()
    PyPlot.plot(abs.(Zp_dist_cactus[:,tsr_idx])./10.0,h,color = plot_cycle[1])  #cactus gives torque over the element, not the unitized
    PyPlot.plot(Zp_AC_noDS[:,tsr_idx].*cos.(delta3D)./10.0,h,color = plot_cycle[2])
    PyPlot.plot(Zp_AC_BV[:,tsr_idx].*cos.(delta3D)./10.0,h,"--",color = plot_cycle[2])
    PyPlot.plot(Zp_DMS_noDS[:,tsr_idx].*cos.(delta3D)./10.0,h,color = plot_cycle[3])
    PyPlot.plot(Zp_DMS_BV[:,tsr_idx].*cos.(delta3D)./10.0,h,"--",color = plot_cycle[3])
    PyPlot.plot(shapeX,shapeY,"k")
    PyPlot.plot(-shapeX,shapeY,"k")
    PyPlot.xlabel("Azimuthally Averaged Absolute\nVertical Force per span (N/m) /10\n Lateral Location (m)")
    PyPlot.ylabel("Vertical Location (m)")
    PyPlot.legend(["CACTUS","AC no DS","AC BV","DMS no DS","DMS BV", "VAWT Shape"])
    PyPlot.axis("equal")
    PyPlot.savefig("$(path)/figs/troposkein/SNL5m_shape_Zp_$(tsr_digit1)_$(tsr_digit2)_$(r_delta_infl).pdf",transparent = true)

    # Calculate own integral CP
    Q_dist_cactus2 = r3D.*Tp_dist_cactus
    P_dist_cactus2 = abs(mean(omega)) * Q_dist_cactus2
    P_dist_cactus3 = abs(mean(omega)).* te_cactus.*(0.5*1.225*0.8*mean(Vinf)^2 * RefArea * R)./h_elem.*cos.(delta3D)
    P_total_dist_cactus2 = abs(mean(omega)) * - OWENSAero.trapz(h,Q_dist_cactus2[:,tsr_idx]./cos.(delta3D))
    P_total_dist_cactus3 = OWENSAero.trapz(h,P_dist_cactus3[:,tsr_idx]./cos.(delta3D))
    CP_dist_cactus2 = P_dist_cactus2 ./ (0.5*1.225*0.8*mean(Vinf)^3 .* RefArea)
    CP_dist_cactus3 = P_dist_cactus3 ./ (0.5*1.225*0.8*mean(Vinf)^3 .* RefArea)
    CP_total_dist_cactus2[tsr_idx] = P_total_dist_cactus2 ./ (0.5*1.225*0.8*mean(Vinf)^3 .* RefArea)
    CP_total_dist_cactus3[tsr_idx] = P_total_dist_cactus3 ./ (0.5*1.225*0.8*mean(Vinf)^3 .* RefArea)

    Q_DMS_noDS = r3D.*Tp_DMS_noDS
    P_DMS_noDS = abs(mean(omega)) * Q_DMS_noDS
    P_total_DMS_noDS = abs(mean(omega)) * - OWENSAero.trapz(h,Q_DMS_noDS[:,tsr_idx])
    CP_DMS_noDS = P_DMS_noDS ./ (0.5*1.225*0.8*mean(Vinf)^3 .* RefArea)
    CP_total_DMS_noDS[tsr_idx] = P_total_DMS_noDS ./ (0.5*1.225*0.8*mean(Vinf)^3 .* RefArea)

    Q_DMS_BV = r3D.*Tp_DMS_BV
    P_DMS_BV = abs(mean(omega)) * Q_DMS_BV
    P_total_DMS_BV = abs(mean(omega)) * - OWENSAero.trapz(h,Q_DMS_BV[:,tsr_idx])
    CP_DMS_BV = P_DMS_BV ./ (0.5*1.225*0.8*mean(Vinf)^3 .* RefArea)
    CP_total_DMS_BV[tsr_idx] = P_total_DMS_BV ./ (0.5*1.225*0.8*mean(Vinf)^3 .* RefArea)

    Q_AC_noDS = r3D.*Tp_AC_noDS
    P_AC_noDS = abs(mean(omega)) * Q_AC_noDS
    P_total_AC_noDS = abs(mean(omega)) * - OWENSAero.trapz(h,Q_AC_noDS[:,tsr_idx])
    CP_AC_noDS = P_AC_noDS ./ (0.5*1.225*0.8*mean(Vinf)^3 .* RefArea)
    CP_total_AC_noDS[tsr_idx] = P_total_AC_noDS ./ (0.5*1.225*0.8*mean(Vinf)^3 .* RefArea)

    Q_AC_BV = r3D.*Tp_AC_BV
    P_AC_BV = abs(mean(omega)) * Q_AC_BV
    P_total_AC_BV = abs(mean(omega)) * - OWENSAero.trapz(h,Q_AC_BV[:,tsr_idx])
    CP_AC_BV = P_AC_BV ./ (0.5*1.225*0.8*mean(Vinf)^3 .* RefArea)
    CP_total_AC_BV[tsr_idx] = P_total_AC_BV ./ (0.5*1.225*0.8*mean(Vinf)^3 .* RefArea)

    # Distributed Torque
    PyPlot.figure()
    PyPlot.plot(te_cactus[:,tsr_idx].*(0.5*1.225*0.8*mean(Vinf)^2 * RefArea * R)./h_elem.*cos.(delta3D)./10,h,"--",color = plot_cycle[1])
    PyPlot.plot(Q_dist_cactus2[:,tsr_idx]./10,h,"-",color = plot_cycle[1])  #cactus gives torque over the element, not the unitized
    PyPlot.plot(Q_AC_noDS[:,tsr_idx].*cos.(delta3D)./10,h,color = plot_cycle[2])
    PyPlot.plot(Q_AC_BV[:,tsr_idx].*cos.(delta3D)./10,h,"--",color = plot_cycle[2])
    PyPlot.plot(Q_DMS_noDS[:,tsr_idx].*cos.(delta3D)./10,h,color = plot_cycle[3])
    PyPlot.plot(Q_DMS_BV[:,tsr_idx].*cos.(delta3D)./10,h,"--",color = plot_cycle[3])
    PyPlot.plot(shapeX,shapeY,"k")
    PyPlot.plot(-shapeX,shapeY,"k")
    PyPlot.xlabel("Azimuthally Averaged\nTorque per span (N-m/m) /10\n Lateral Location (m)")
    PyPlot.ylabel("Vertical Location (m)")
    PyPlot.legend(["CACTUS te","CACTUS CT","AC no DS","AC BV","DMS no DS","DMS BV", "VAWT Shape"])
    # PyPlot.legend(["DMS no DS","DMS BV","AC no DS","AC BV"])
    PyPlot.axis("equal")
    PyPlot.savefig("$(path)/figs/troposkein/SNL5m_shape_Qp_$(tsr_digit1)_$(tsr_digit2).pdf",transparent = true)

    # Distributed Power
    PyPlot.figure()
    PyPlot.plot(P_dist_cactus3[:,tsr_idx]./100,h,"--",color = plot_cycle[1])
    PyPlot.plot(P_dist_cactus2[:,tsr_idx]./100,h,"-",color = plot_cycle[1])  #cactus gives torque over the element, not the unitized
    PyPlot.plot(P_AC_noDS[:,tsr_idx].*cos.(delta3D)./100,h,color = plot_cycle[2])
    PyPlot.plot(P_AC_BV[:,tsr_idx].*cos.(delta3D)./100,h,"--",color = plot_cycle[2])
    PyPlot.plot(P_DMS_noDS[:,tsr_idx].*cos.(delta3D)./100,h,color = plot_cycle[3])
    PyPlot.plot(P_DMS_BV[:,tsr_idx].*cos.(delta3D)./100,h,"--",color = plot_cycle[3])
    PyPlot.plot(shapeX,shapeY,"k")
    PyPlot.plot(-shapeX,shapeY,"k")
    PyPlot.xlabel("Azimuthally Averaged\nPower per span (W/m) /100\n Lateral Location (m)")
    PyPlot.ylabel("Vertical Location (m)")
    PyPlot.legend(["CACTUS te: $(round(P_total_dist_cactus3,digits=3)) W","CACTUS CT: $(round(P_total_dist_cactus2,digits=3)) W","AC no DS: $(round(P_total_AC_noDS,digits=3)) W","AC BV: $(round(P_total_AC_BV,digits=3)) W","DMS no DS: $(round(P_total_DMS_noDS,digits=3)) W","DMS BV: $(round(P_total_DMS_BV,digits=3)) W", "VAWT Shape"])
    # PyPlot.legend(["DMS no DS","DMS BV","AC no DS","AC BV"])
    PyPlot.axis("equal")
    PyPlot.savefig("$(path)/figs/troposkein/SNL5m_shape_Power_$(tsr_digit1)_$(tsr_digit2).pdf",transparent = true)

    # Distributed CP as recalculated
    PyPlot.figure()
    PyPlot.plot(CP_dist_cactus3[:,tsr_idx]*100,h,"--",color = plot_cycle[1]) #CACTUS te
    PyPlot.plot(CP_dist_cactus2[:,tsr_idx]*100,h,"-",color = plot_cycle[1]) #CACTUS CT
    PyPlot.plot(CP_AC_noDS[:,tsr_idx].*cos.(delta3D)*100,h,color = plot_cycle[2])
    PyPlot.plot(CP_AC_BV[:,tsr_idx].*cos.(delta3D)*100,h,color = plot_cycle[2],"--")
    PyPlot.plot(CP_DMS_noDS[:,tsr_idx].*cos.(delta3D)*100,h,color = plot_cycle[3])
    PyPlot.plot(CP_DMS_BV[:,tsr_idx].*cos.(delta3D)*100,h,color = plot_cycle[3],"--")
    # PyPlot.plot(-r3D,h,"b.")
    # PyPlot.plot(r3D,h,"b.")
    PyPlot.plot(shapeX,shapeY,"k")
    PyPlot.plot(-shapeX,shapeY,"k")
    PyPlot.xlabel("X")
    PyPlot.ylabel("Y")
    PyPlot.title("CP distribution per span x100\n CACTUS Nominal CP: $(round(cactus_CP[tsr_idx],digits=3))")
    PyPlot.legend(["CACTUS te: $(round(CP_total_dist_cactus3[tsr_idx],digits=3))","CACTUS CT: $(round(CP_total_dist_cactus2[tsr_idx],digits=3))","AC no DS: $(round(CP_total_AC_noDS[tsr_idx],digits=3))","DMS BV: $(round(CP_total_DMS_BV[tsr_idx],digits=3))","DMS no DS: $(round(CP_total_DMS_noDS[tsr_idx],digits=3))","AC BV: $(round(CP_total_AC_BV[tsr_idx],digits=3))"])
    PyPlot.axis("equal")
    PyPlot.savefig("$(path)/figs/troposkein/SNL5m_shape_CP_$(tsr_digit1)_$(tsr_digit2).pdf",transparent = true)

    # Distributed CP as calculated within the models
    PyPlot.figure()
    PyPlot.plot(CP_dist_cactus3[:,tsr_idx]*100,h,"--",color = plot_cycle[1]) #CACTUS te
    PyPlot.plot(CP_dist_cactus2[:,tsr_idx]*100,h,"-",color = plot_cycle[1]) #CACTUS CT
    PyPlot.plot(Cp_1AC_noDS[:,tsr_idx].*2*R./RefArea*100,h,color = plot_cycle[2])
    PyPlot.plot(Cp_1AC_BV[:,tsr_idx].*2*R./RefArea*100,h,color = plot_cycle[2],"--")
    PyPlot.plot(Cp_1DMS_noDS[:,tsr_idx].*2*R./RefArea*100,h,color = plot_cycle[3])
    PyPlot.plot(Cp_1DMS_BV[:,tsr_idx].*2*R./RefArea*100,h,color = plot_cycle[3],"--")
    # PyPlot.plot(-r3D,h,"b.")
    # PyPlot.plot(r3D,h,"b.")
    PyPlot.plot(shapeX,shapeY,"k")
    PyPlot.plot(-shapeX,shapeY,"k")
    PyPlot.xlabel("X")
    PyPlot.ylabel("Y")
    PyPlot.title("CP distribution per height x100")
    # PyPlot.legend(["CACTUS: $(round(cactus_CP[tsr_idx],digits=3))","CACTUS CT: $(round(CP_total_dist_cactus2[tsr_idx],digits=3))","AC no DS: $(round(CPvec_AC_noDS[tsr_idx]),digits=3))","AC BV: $(round(CPvec_AC_BV[tsr_idx]),digits=3))","DMS no DS: $(round(CPvec_DMS_noDS[tsr_idx]),digits=3))","DMS BV: $(round(CPvec_DMS_BV[tsr_idx]),digits=3))"])
    PyPlot.legend(["CACTUS te: $(round(CP_total_dist_cactus3[tsr_idx],digits=3))","CACTUS CT: $(round(CP_total_dist_cactus2[tsr_idx],digits=3))","AC no DS: $(round(CPvec_AC_noDS[tsr_idx],digits=3))","AC BV: $(round(CPvec_AC_BV[tsr_idx],digits=3))","DMS no DS: $(round(CPvec_DMS_noDS[tsr_idx],digits=3))","DMS BV: $(round(CPvec_DMS_BV[tsr_idx],digits=3))"])
    PyPlot.axis("equal")
    PyPlot.savefig("$(path)/figs/troposkein/SNL5m_shape_as_calculated_CP_$(tsr_digit1)_$(tsr_digit2).pdf",transparent = true)
end


cactus_Time_data = DelimitedFiles.readdlm("$(path)/data/steady/TestVAWT5_2_TimeData.csv",',',Float64,skipstart = 1)

test = cactus_Time_data[cactus_Time_data[:,3].==maximum(cactus_Time_data[:,3])-1,:] #All data that aligns with revolution 15
println("Average CP: $(mean(test[:,5]))")

PyPlot.rc("figure", figsize=(4, 2.5))
PyPlot.pygui(true)
PyPlot.figure()
PyPlot.plot(experimental[1:end-1,1],experimental[1:end-1,2],"k.",label="Exp.")
PyPlot.plot(tsrvec,cactus_CP,".-",markersize = 6,linewidth = 1.4,color = plot_cycle[1],label="Cactus")
# PyPlot.plot(tsrvec,CPvec_AC_noDS,"-",color = plot_cycle[2],label="CPvec_AC_noDS")
PyPlot.plot(tsrvec,CPvec_AC_BV,".-",markersize = 6,linewidth = 1.4,color = plot_cycle[2],label="AC")
# PyPlot.plot(tsrvec,CP_total_DMS_noDS,"-",color = plot_cycle[3],label="CP_total_DMS_noDS")
PyPlot.plot(tsrvec,CP_total_DMS_BV,".-",markersize = 6,linewidth = 1.4,color = plot_cycle[3],label="DMS")
# PyPlot.plot(tsrvec,CPvec_DMS_noDS,"-",color = plot_cycle[2],label="CPvec_DMS_noDS")
PyPlot.plot(tsrvec,CPvec_DMS_BV,"--",color = plot_cycle[2],label="CPvec_DMS_BV")
# PyPlot.plot(tsrvec,CP_total_AC_noDS,"-",color = plot_cycle[4],label="CP_total_AC_noDS")
PyPlot.plot(tsrvec,CP_total_AC_BV,"--",color = plot_cycle[4],label="CP_total_AC_BV")
PyPlot.xlabel("Tip Speed Ratio")
PyPlot.ylabel("Coefficient of Performance")
PyPlot.legend()#["Experimental", "Cactus","AC","DMS"],loc = "lower center")#"AC noDS","AC BV","DMS noDS","DMS BV","Cactus CT","Cactus te","AC noDS","AC BV","DMS noDS","DMS BV"])
# PyPlot.savefig("$(path)/figs/troposkein/SNL5m_CP3D.pdf",transparent = true)

#
#
# Influence of R-Delta
#
#

r_delta_infl = true
start = time()
delta3D,experimental,cactus_CP,tsrvec,CPvec_AC_BV_rdel,shapeX,shapeY,element_planf_L,Rp_AC_BV_rdel,Tp_AC_BV_rdel,Zp_AC_BV_rdel,Cp_1AC_BV_rdel,h,h_frac,Rp_dist_cactus,Tp_dist_cactus,Zp_dist_cactus,CP_dist_cactus,te_cactus,r3D,R,omega,B,thetavec,RefArea,Tp_AC_BV_rdel_theta,Rp_AC_BV_rdel_theta,Zp_AC_BV_rdel_theta,dataout,tsrvec = runme("AC","BV",r_delta_infl)
elapsed_AC_BV_rdel = time() - start
println("elapsed_AC_BV $elapsed_AC_BV_rdel")

h_elem = (shapeY[2:end] - shapeY[1:end-1])

# Tangential Force
tsr_idx = 1#5
tsr_digit1 = Int(floor(tsrvec[tsr_idx]))
tsr_digit2 = Int(round((tsrvec[tsr_idx]-tsr_digit1)*10.0))
Vinf = abs.(omega)/tsrvec[tsr_idx]*R

PyPlot.rc("figure.subplot", left=.12, bottom=.22, top=0.95, right=.95)
PyPlot.rc("figure", figsize=(5, 3.5))
PyPlot.figure()
PyPlot.plot(te_cactus[:,tsr_idx].*(0.5*1.225*0.8*mean(Vinf)^2 * RefArea * R)./r3D./h_elem.*cos.(delta3D),h,color = plot_cycle[1])
# PyPlot.plot(-Tp_dist_cactus[:,tsr_idx],h,"-",color = plot_cycle[1])  #cactus gives torque over the element, not the unitized
PyPlot.plot(-Tp_AC_BV_rdel[:,tsr_idx].*cos.(delta3D),h,"--",color = plot_cycle[2]) #multiply by cos(delta) to get rid of element length from angle
PyPlot.plot(-Tp_AC_BV[:,tsr_idx].*cos.(delta3D),h,color = plot_cycle[2])
PyPlot.plot(shapeX,shapeY,"k")
PyPlot.plot(-shapeX,shapeY,"k")
# PyPlot.plot(0,5,"w.")
PyPlot.xlabel("Azimuthally Averaged\nTangential Force per span (N/m)")#"\n Lateral Location (m)")
PyPlot.ylabel("Vertical Location (m)")
PyPlot.legend(["CACTUS","AC","AC Influence Coefficient Slope"])#, "VAWT Shape"])
PyPlot.axis("equal")
PyPlot.savefig("$(path)/figs/troposkein/SNL5m_shape_Tp_$(tsr_digit1)_$(tsr_digit2)_rdelta_BOTH.pdf",transparent = true)

# Difference
PyPlot.figure()
normalized_force = abs.(Tp_AC_BV[:,tsr_idx]).*cos.(delta3D) ./ (maximum(abs.(Tp_AC_BV[:,tsr_idx]).*cos.(delta3D)))
myerror = (Tp_AC_BV_rdel[:,tsr_idx].*cos.(delta3D).*r3D - Tp_AC_BV[:,tsr_idx].*cos.(delta3D).*r3D)./Tp_AC_BV[:,tsr_idx].*cos.(delta3D).*r3D
PyPlot.plot(myerror*100,h,"-",color = plot_cycle[2]) #multiply by cos(delta) to get rid of element length from angle
# PyPlot.plot(myerror.*normalized_force*100,h,"--",color = plot_cycle[2])
PyPlot.plot(shapeX,shapeY,"k")
PyPlot.plot(-shapeX,shapeY,"k")
# PyPlot.plot(-35,0,"w.")
PyPlot.xlabel("Torque Difference (%)")#"\n Lateral Location (m)")
PyPlot.ylabel("Vertical Location (m)")
PyPlot.legend(["Difference"])#, "VAWT Shape"])
PyPlot.axis("equal")
PyPlot.savefig("$(path)/figs/troposkein/SNL5m_shape_Tp_$(tsr_digit1)_$(tsr_digit2)_rdelta_BOTH_ERROR.pdf",transparent = true)
