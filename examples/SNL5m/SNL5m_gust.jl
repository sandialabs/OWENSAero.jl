
import PyPlot
PyPlot.ion()
PyPlot.rc("figure", figsize=(4, 3))
PyPlot.rc("font", size=10.0)
PyPlot.rc("lines", linewidth=1.5)
PyPlot.rc("lines", markersize=3.0)
PyPlot.rc("legend", frameon=false)
PyPlot.rc("axes.spines", right=false, top=false)
PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
# rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

import Statistics:mean
import DelimitedFiles
import HDF5
import FLOWMath
import QuadGK
using Test
import OWENSAero

path,_ = splitdir(@__FILE__)
# include("$(path)/../../../OWENSAero.jl/src/OWENSAero.jl")

function runme(AeroModel,slice,bnum)
    # AeroModel = "DMS"
    # slice = 7
    bnum = 1
    ntheta = 30
    # AeroModel = "AC"
    if AeroModel == "AC"
        N_aw = ntheta*2 #Number of either "a" (DMS induction factor) or w (u and v induction factors)
    elseif AeroModel == "DMS"
        N_aw = ntheta
    else
        error("Aeromodel not recognized, choose AC or DMS")
    end
    rotation = -1.0
    R = 5.0/2 #m
    H = 1.02*R*2 #m
    chord = 0.1524 #m
    B = 3
    # bnum=1 #for plotting
    k = 1.0
    rho = 1.225*0.8
    mu = 1.7894e-5
    suction = false
    eta = 0.4 #blade mount point, or distance from leading edge that the rotation plane passes through.
    # Unsteady Params
    RPI = false
    tau = [0.7,0.2]
    N_Rev = 40#20
    G_amp = 5.0
    gustX0 = 26.0#13.0
    steplast = zeros(Int,1,1)
    idx_RPI = zeros(Int,ntheta*2,1)
    awsave = zeros(Real,N_aw)
    gustT = 0.8

    function affun(alpha, Re, M;V_twist=nothing,chord=nothing,dt=nothing,Vloc=nothing)

        cl = 6.2*alpha
        cd = 0.008 .- 0.003.*cl + 0.01.*cl.*cl

        return cl, cd
    end

    DS_model = "BV" #TODO: resolve the convergence issues when RPI is used with BV
    # if DS_model=="BV"
        # af = OWENSAero.readaerodyn_BV("$(path)/airfoils/NACA_0015_RE3E5.dat")
    # else
    #     af = OWENSAero.readaerodyn("$(path)/airfoils/NACA_0015_RE3E5.dat")
    # end
    af = OWENSAero.readaerodyn_BV_NEW("$(path)/airfoils/NACA_0015.dat",DynamicStallModel="BV")
    awwarm = [0.0734080061048192, 0.12251242922733005, 0.1741623096669431, 0.21463048262337026, 0.24095130886069838, 0.2525947743402928, 0.2506267437330584, 0.23777593110718653, 0.2183280835908269, 0.1966037927062711, 0.17431172191618466, 0.14941630864787642, 0.1179208332087899, 0.07703049955652619, 0.025895573316087175, 0.006295173150749581, 0.05855375428983589, 0.11585393975035356, 0.17692022645228292, 0.23679806306693904, 0.29188854267978814, 0.3388985559576654, 0.3716153405175147, 0.38102986264893646, 0.3602220274361156, 0.3094479766496759, 0.23763375800321238, 0.1587359710414027, 0.08511863905660223, 0.031421675399937114, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    awwarm=zero(awwarm)
    shapeX_raw_temp = 1 / cos(12.38*pi/180) .* [ -6.29E-03, 1.33E-01, 2.93E-01, 4.45E-01, 5.71E-01, 6.56E-01, 7.22E-01, 7.81E-01, 8.29E-01, 8.77E-01, 9.16E-01, 9.46E-01, 9.65E-01, 9.77E-01, 9.83E-01, 9.83E-01, 9.83E-01, 9.79E-01, 9.68E-01, 9.54E-01, 9.34E-01, 9.11E-01, 8.93E-01, 8.68E-01, 8.45E-01, 8.13E-01, 7.76E-01, 7.30E-01, 6.73E-01, 6.27E-01, 5.34E-01, 4.13E-01, 2.92E-01, 1.76E-01, 5.88E-02, 1.74E-03]
    shapeX_raw = R * (shapeX_raw_temp./maximum(shapeX_raw_temp))
    shapeY_raw_temp = [-4.12E-03,3.25E-02,8.62E-02,1.36E-01,1.79E-01,2.08E-01,2.35E-01,2.63E-01,2.90E-01,3.19E-01,3.56E-01,3.92E-01,4.24E-01,4.57E-01,4.84E-01,5.01E-01,5.20E-01,5.47E-01,5.84E-01,6.09E-01,6.41E-01,6.65E-01,6.84E-01,7.05E-01,7.20E-01,7.40E-01,7.58E-01,7.83E-01,8.02E-01,8.17E-01,8.44E-01,8.82E-01,9.20E-01,9.56E-01,9.93E-01,1.01E+00]
    shapeY_raw = H * shapeY_raw_temp./maximum(shapeY_raw_temp)
    shapeX_spline = FLOWMath.Akima(shapeY_raw, shapeX_raw)
    n_slices = 30#length(shapeX)-1;
    shapeY = LinRange(0,H,n_slices+1)
    shapeX = R.*(1.0.-4.0.*(shapeY/H.-.5).^2)#shapeX_spline(shapeY)
    n_slices = length(shapeX)-1;
    h_frac = (shapeY[2:end] - shapeY[1:end-1])./shapeY[end];
    h = (shapeY[2:end] + shapeY[1:end-1])/2.0;

    RefArea_half, _ = QuadGK.quadgk(shapeX_spline, 0, H, atol=1e-10)
    RefArea = RefArea_half*2

    delta_xs = shapeX[2:end] - shapeX[1:end-1]
    delta_zs = shapeY[2:end] - shapeY[1:end-1]

    delta3D = atan.(delta_xs./delta_zs)

    element_planf_A = sqrt.(delta_xs.^2+delta_zs.^2)*chord
    element_planf_L = sqrt.(delta_xs.^2+delta_zs.^2)

    r3D = (shapeX[2:end,1]+shapeX[1:end-1,1])/2.0
    aerocenter_dist = (eta-.25)*chord

    twist3D = -atan.(aerocenter_dist./r3D)#ones(Real,n_slices)*-0.4*pi/180


    # Evaluate the slice
    r = ones(Real,ntheta)*r3D[Int(slice)] #m
    chord = 0.1524 #m
    twist = ones(Real,ntheta)*twist3D[Int(slice)] #rad
    delta = ones(Real,ntheta)*delta3D[Int(slice)] #rad
    omega = rotation .* ones(Real,ntheta) * 150.0 / 60.0 * 2*pi # RPM -> radians/sec
    tsrvec = [2.1,3.1,4.2,5.2,6.0,7.0]
    Vinf = abs.(omega)/4.0*maximum(r3D)
    V_wake_old = ones(Real,1).*mean(Vinf)
    nominalVinf = Vinf[1]

    IECgust = true
    ifw = false
    substep_level = 1

        ###
    ### Solve 1st order filter
    ###
    env = OWENSAero.Environment(rho,mu,Vinf,DS_model,AeroModel,awwarm)
    turbine2D = OWENSAero.Turbine(R,r,chord,twist,delta,omega,B,af,ntheta,false)
    us_param = OWENSAero.UnsteadyParams(RPI,tau,ifw,IECgust,nominalVinf,G_amp,gustX0,gustT)

    CP_filt = zeros(Real,ntheta*N_Rev*substep_level)
    Qsave = zeros(Real,ntheta)
    Vinf_used_filt = zeros(Real,ntheta)
    Vinfhist_filt = zeros(Real,ntheta*N_Rev*substep_level)
    Rp_filt = zeros(Real,ntheta*N_Rev*substep_level)
    Rp_filt_blades = zeros(Real,ntheta*N_Rev*substep_level,B)
    Tp_filt = zeros(Real,ntheta*N_Rev*substep_level)
    Tp_filt_blades = zeros(Real,ntheta*N_Rev*substep_level,B)
    Zp_filt = zeros(Real,ntheta*N_Rev*substep_level)
    alpha_filt = zeros(Real,ntheta*N_Rev*substep_level)
    cl_filt = zeros(Real,ntheta*N_Rev*substep_level)
    cd_filt = zeros(Real,ntheta*N_Rev*substep_level)
    Vloc_filt = zeros(Real,ntheta*N_Rev*substep_level)
    Re_filt = zeros(Real,ntheta*N_Rev*substep_level)

    println("Running First Order Filter Method")
    start = time()

    for t = 1:ntheta*(N_Rev)*substep_level
        # println(t/(ntheta*(N_Rev)*substep_level))
        step1 = Int(ceil(t/substep_level))
        rev_step1 = Int(step1-floor(Int,(step1-1)/ntheta)*ntheta + (1-1)*ntheta/B) #TODO: include other blades simultaneously for aerostructural deformation (i.e. updating r for each blade as it passes)
        rev_step2 = Int(step1-floor(Int,(step1-1)/ntheta)*ntheta + (2-1)*ntheta/B) #TODO: include other blades simultaneously for aerostructural deformation (i.e. updating r for each blade as it passes)
        rev_step3 = Int(step1-floor(Int,(step1-1)/ntheta)*ntheta + (3-1)*ntheta/B) #TODO: include other blades simultaneously for aerostructural deformation (i.e. updating r for each blade as it passes)
        if rev_step1>ntheta
            rev_step1 -= ntheta
        end

        if rev_step2>ntheta
            rev_step2 -= ntheta
        end

        if rev_step3>ntheta
            rev_step3 -= ntheta
        end
        # Juno.@enter OWENSAero.Unsteady_Step(turbine2D,env,us_param,step1)
        _, _, Q_temp, Rp_temp, Tp_temp, Zp_temp, Vloc_temp,_,_,_,_, alpha_temp, cl_temp, cd_temp, thetavec, Re_temp = OWENSAero.Unsteady_Step(turbine2D,env,us_param,step1)

        Rp_filt[step1] = Rp_temp[rev_step1]#.+Rp_temp[rev_step2].+Rp_temp[rev_step3]
        Tp_filt[step1] = Tp_temp[rev_step1]#.+Tp_temp[rev_step2].+Tp_temp[rev_step3]
        Tp_filt_blades[step1,1] = Tp_temp[rev_step1]
        Tp_filt_blades[step1,2] = Tp_temp[rev_step2]
        Tp_filt_blades[step1,3] = Tp_temp[rev_step3]
        Rp_filt_blades[step1,1] = Rp_temp[rev_step1]
        Rp_filt_blades[step1,2] = Rp_temp[rev_step2]
        Rp_filt_blades[step1,3] = Rp_temp[rev_step3]
        Zp_filt[step1] = Zp_temp[rev_step1]#.+Zp_temp[rev_step2].+Zp_temp[rev_step3]
        alpha_filt[step1] = alpha_temp[rev_step1]
        cl_filt[step1] = cl_temp[rev_step1]#.+cl_temp[rev_step2].+cl_temp[rev_step3]
        cd_filt[step1] = cd_temp[rev_step1]#.+cd_temp[rev_step2].+cd_temp[rev_step3]
        Vloc_filt[step1] = Vloc_temp[rev_step1]#.+Vloc_temp[rev_step2].+Vloc_temp[rev_step3]
        Re_filt[step1] = Re_temp[rev_step1]#.+Re_temp[rev_step2].+Re_temp[rev_step3]
        Qsave[:] = Q_temp[:]
        Vinf_used_filt[:] = env.V_x[:]
        # println(env.V_x[1])
        # println(us_param.nominalVinf)
        Vinfhist_filt[step1] = env.V_x[rev_step1]

        CP_filt[step1] = sum((turbine2D.B/(2*pi) * 2*pi/(turbine2D.ntheta)*Qsave*mean(abs.(turbine2D.omega))) / (0.5*env.rho*1.0*2*turbine2D.R*mean(Vinf_used_filt)^3)) # Eq. 14, normalized by nominal radius R
    end


    elapsed = time() - start
    println("Time 1st Order Filter:")
    println(elapsed)

    CP_RPI = zeros(Real,ntheta*N_Rev*substep_level)
    Qsave = zeros(Real,ntheta)
    Vinf_used_RPI = zeros(Real,ntheta)
    Vinfhist_RPI = zeros(Real,ntheta*N_Rev*substep_level)
    Rp_RPI = zeros(Real,ntheta*N_Rev*substep_level)
    Rp_RPI_blades = zeros(Real,ntheta*N_Rev*substep_level,B)
    Tp_RPI = zeros(Real,ntheta*N_Rev*substep_level)
    Tp_RPI_blades = zeros(Real,ntheta*N_Rev*substep_level,B)
    Tp_RPI2 = zeros(Real,ntheta*N_Rev*substep_level)
    Zp_RPI = zeros(Real,ntheta*N_Rev*substep_level)
    alpha_RPI = zeros(Real,ntheta*N_Rev*substep_level)
    cl_RPI = zeros(Real,ntheta*N_Rev*substep_level)
    cd_RPI = zeros(Real,ntheta*N_Rev*substep_level)
    Vloc_RPI = zeros(Real,ntheta*N_Rev*substep_level)
    Re_RPI = zeros(Real,ntheta*N_Rev*substep_level)

    idx_RPI = zeros(Int,B*2,1)
    awsave = zeros(Real,N_aw)
    V_wake_old = ones(Real,1).*Vinf
    nominalVinf = Vinf[1]
    tau = [0.3,3.0]
    RPI = true
    env = OWENSAero.Environment(rho,mu,Vinf,DS_model,AeroModel,awwarm)
    turbine2D = OWENSAero.Turbine(R,r,chord,twist,delta,omega,B,af,ntheta,false)
    us_param = OWENSAero.UnsteadyParams(RPI,tau,ifw,IECgust,nominalVinf,G_amp,gustX0,gustT)
    start = time()
    turbines = Array{OWENSAero.Turbine}(undef,1)
    turbines[1] = turbine2D

    println("Running RPI Method")

    for t = 1:ntheta*(N_Rev)*substep_level
        # t = 1
        # println(t/(ntheta*(N_Rev)*substep_level))
        step1 = Int(ceil(t/substep_level))
        rev_step1 = Int(step1-floor(Int,(step1-1)/ntheta)*ntheta + (1-1)*ntheta/B) #TODO: include other blades simultaneously for aerostructural deformation (i.e. updating r for each blade as it passes)
        rev_step2 = Int(step1-floor(Int,(step1-1)/ntheta)*ntheta + (2-1)*ntheta/B) #TODO: include other blades simultaneously for aerostructural deformation (i.e. updating r for each blade as it passes)
        rev_step3 = Int(step1-floor(Int,(step1-1)/ntheta)*ntheta + (3-1)*ntheta/B) #TODO: include other blades simultaneously for aerostructural deformation (i.e. updating r for each blade as it passes)
        if rev_step1>ntheta
            rev_step1 -= ntheta
        end

        if rev_step2>ntheta
            rev_step2 -= ntheta
        end

        if rev_step3>ntheta
            rev_step3 -= ntheta
        end

        # Juno.@enter OWENSAero.Unsteady_Step(turbine2D,env,us_param,step1)
        _, _, Q_temp, Rp_temp, Tp_temp, Zp_temp, Vloc_temp,_,_,_,_, alpha_temp, cl_temp, cd_temp, thetavec, Re_temp = OWENSAero.Unsteady_Step(turbine2D,env,us_param,step1)

        Rp_RPI[step1] = Rp_temp[rev_step1]#.+Rp_temp[rev_step2].+Rp_temp[rev_step3]
        Tp_RPI[step1] = Tp_temp[rev_step1]#.+Tp_temp[rev_step2].+Tp_temp[rev_step3]
        Tp_RPI_blades[step1,1] = Tp_temp[rev_step1]
        Tp_RPI_blades[step1,2] = Tp_temp[rev_step2]
        Tp_RPI_blades[step1,3] = Tp_temp[rev_step3]
        Rp_RPI_blades[step1,1] = Rp_temp[rev_step1]
        Rp_RPI_blades[step1,2] = Rp_temp[rev_step2]
        Rp_RPI_blades[step1,3] = Rp_temp[rev_step3]
        Zp_RPI[step1] = Zp_temp[rev_step1]#.+Zp_temp[rev_step2].+Zp_temp[rev_step3]
        alpha_RPI[step1] = alpha_temp[rev_step1]
        cl_RPI[step1] = cl_temp[rev_step1]#.+cl_temp[rev_step2].+cl_temp[rev_step3]
        cd_RPI[step1] = cd_temp[rev_step1]#.+cd_temp[rev_step2].+cd_temp[rev_step3]
        Vloc_RPI[step1] = Vloc_temp[rev_step1].+Vloc_temp[rev_step2].+Vloc_temp[rev_step3]
        Re_RPI[step1] = Re_temp[rev_step1]#.+Re_temp[rev_step2].+Re_temp[rev_step3]

        Qsave[:] = Q_temp[:]
        Vinf_used_RPI[:] = env.V_x[:]
        Vinfhist_RPI[step1] = env.V_x[Int(ntheta/2)]

        CP_RPI[step1] = sum((turbine2D.B/(2*pi) * 2*pi/(turbine2D.ntheta)*Qsave*mean(abs.(turbine2D.omega))) / (0.5*env.rho*1.0*2*turbine2D.R*mean(Vinf_used_RPI)^3)) # Eq. 14, normalized by nominal radius R
    end
    elapsed = time() - start
    println("Time RPI:")
    println(elapsed)


    # PyPlot.figure()
    # rev_array = LinRange(2*pi/(turbine2D.ntheta)/2,N_Rev*2.0*pi,length(Tp_filt))/(2*pi)
    # PyPlot.plot(rev_array,Vinfhist_RPI)
    # PyPlot.plot(1:ntheta,Vinf_used_RPI)
    # PyPlot.savefig("$(path)/figs/gust/SNL5m_gustVel_$(slice)_gtime$(gustT).pdf",transparent = true)

    if false#true


        PyPlot.close("all")

        t = collect(1:ntheta*(N_Rev)*substep_level)

        data = DelimitedFiles.readdlm("$(path)/data/gust/RerunTestVAWT5_2_Gust2_Troposkein_ElementData.csv", ',',skipstart = 1)
        revdata = DelimitedFiles.readdlm("$(path)/data/gust/RerunTestVAWT5_2_Gust2_Troposkein_RevData.csv", ',',skipstart = 1)
        # Extract Data for Blade bnum
        used_data_logic = (data[:,3].==bnum) .& (data[:,4].==slice)
        data_full_rev1 = data[used_data_logic,:]

        rev_array = LinRange(2*pi/(turbine2D.ntheta)/2,N_Rev*2.0*pi,length(Tp_filt))/(2*pi)
        q_loc = 0.5*env.rho.*(data_full_rev1[:,15]*mean(Vinf)).^2

        r_start = 20.0;
        r_end = 25.0;

        # CP
        PyPlot.figure()
        PyPlot.plot(revdata[:,1],revdata[:,2],"-",color=plot_cycle[1],linewidth = 1.0)
        PyPlot.plot(rev_array,CP_filt,color = plot_cycle[2],linewidth = 1.0)
        PyPlot.plot(rev_array,CP_RPI,color = plot_cycle[1],linewidth = 1.0)
        PyPlot.ylim([0.0,0.6])
        PyPlot.xlabel("Revolution")
        PyPlot.ylabel("CP")
        # PyPlot.xlim([r_start,r_end])
        PyPlot.legend(["CACTUS", "$AeroModel 1st Order Filter","$AeroModel RPI"])
        # PyPlot.savefig("$(path)/figs/gust/SNL5m_$(AeroModel)Transient_CP_$(slice)_gtime$(gustT).pdf",transparent = true)

        # Angle of Attack
        PyPlot.figure()
        PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), data_full_rev1[:,9],"",color=plot_cycle[1],linewidth = 0.75)
        PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), data_full_rev1[:,10],"--",color=plot_cycle[1],linewidth = 0.75)
        PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), data_full_rev1[:,11],"-.",color=plot_cycle[1],linewidth = 0.75)
        PyPlot.plot(rev_array, alpha_filt*180/pi,color=plot_cycle[2],linewidth = 0.75)
        PyPlot.plot(rev_array, alpha_RPI*180/pi,color=plot_cycle[3],linewidth = 0.75)
        PyPlot.xlabel("Theta (deg)")
        PyPlot.xlim([r_start,r_end])
        PyPlot.ylabel("AOA (deg)")
        PyPlot.legend(["CACTUS 25","CACTUS 50","CACTUS 75", "$AeroModel 1st Order Filter","$AeroModel RPI"])
        # PyPlot.savefig("$(path)/figs/gust/SNL5m_$(AeroModel)Transient_alpha$(bnum)_$(slice)_gtime$(gustT).pdf",transparent = true)

        # Rp
        PyPlot.figure()
        PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), -q_loc.*(data_full_rev1[:,24]).*chord.*cos.(turbine2D.delta[1]),"-",color=plot_cycle[1],linewidth = 0.75)
        PyPlot.plot(rev_array,-Rp_filt.*cos.(turbine2D.delta[1]),color = plot_cycle[2],linewidth = 0.75)
        PyPlot.plot(rev_array,-Rp_RPI.*cos.(turbine2D.delta[1]),color = plot_cycle[3],linewidth = 0.75)
        PyPlot.xlabel("Revolution")
        PyPlot.ylabel("Radial Force per Span (N/m)")
        PyPlot.xlim([r_start,r_end])
        # PyPlot.legend(["CACTUS", "$AeroModel 1st Order Filter","$AeroModel RPI"])
        # PyPlot.savefig("$(path)/figs/gust/SNL5m_$(AeroModel)Transient_Rp$(bnum)_$(slice)_gtime$(gustT).pdf",transparent = true)

        # Tp
        PyPlot.figure()
        PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), -q_loc.*(data_full_rev1[:,25]).*chord,"-",color=plot_cycle[1],linewidth = 0.75)
        PyPlot.plot(rev_array, -Tp_filt.*cos.(turbine2D.delta[1]),color = plot_cycle[2],linewidth = 0.75)
        PyPlot.plot(rev_array, -Tp_RPI.*cos.(turbine2D.delta[1]),color = plot_cycle[3],linewidth = 0.75)
        PyPlot.xlabel("Revolution")
        PyPlot.ylabel("Tangential Force per Span (N/m)")
        PyPlot.xlim([r_start,r_end])
        PyPlot.legend(["CACTUS", "$AeroModel 1st Order Filter","$AeroModel RPI"])
        # PyPlot.savefig("$(path)/figs/gust/SNL5m_$(AeroModel)Transient_Tp$(bnum)_$(slice)_gtime$(gustT).pdf",transparent = true)

        # Zp
        PyPlot.figure()
        PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), -q_loc.*(data_full_rev1[:,24]).*chord.*sin.(turbine2D.delta[1]),"-",color=plot_cycle[1],linewidth = 0.75)
        PyPlot.plot(rev_array,-Zp_filt.*cos.(turbine2D.delta[1]),color = plot_cycle[2],linewidth = 0.75)
        PyPlot.plot(rev_array,-Zp_RPI.*cos.(turbine2D.delta[1]),color = plot_cycle[3],linewidth = 0.75)
        PyPlot.xlabel("Revolution")
        PyPlot.ylabel("Vertical Force per Span (N/m)")
        PyPlot.xlim([r_start,r_end])
        # PyPlot.legend(["CACTUS", "$AeroModel 1st Order Filter","$AeroModel RPI"])
        # PyPlot.savefig("$(path)/figs/gust/SNL5m_$(AeroModel)Transient_Zp$(bnum)_$(slice)_gtime$(gustT).pdf",transparent = true)
        # end
        # #
        # function myplots()
        # Plot gust velocity
        mystep = collect(1:ntheta*N_Rev)
        dt = 1.0./(abs(turbine2D.omega[1]) .* ntheta ./ (2.0.*pi))
        dt_norm = dt.*us_param.nominalVinf./turbine2D.R
        gustT = us_param.gustT .* us_param.nominalVinf ./ turbine2D.R
        tr = mystep.*dt_norm .- us_param.gustX0
        IECGustFactor = 1.0 .- 0.37 .* us_param.G_amp/us_param.nominalVinf .* sin.(3*pi*tr./gustT)  .* (1.0 .- cos.(2*pi*tr./gustT))

        IECGustFactor[tr .<= 0] .= 1.0
        IECGustFactor[tr .>= gustT] .= 1.0

        Vinf_used = zeros(Real,turbine2D.B)
        Vinf_used_long = zeros(Real,length(mystep),turbine2D.B)
        for ii = 1:length(mystep)
            ele_x = rotation*sin.(float((mystep[ii]:ntheta/turbine2D.B:mystep[ii]+ntheta-ntheta/turbine2D.B+1)./ntheta*2*pi)).*turbine2D.r[1]./turbine2D.R
            tr = mystep[ii]*dt_norm .- ele_x .- gustX0
            for bld_i = 1:turbine2D.B
                if (tr[bld_i] >= 0) && (tr[bld_i]<=gustT)

                    IECGustFactor_blade = 1.0 - 0.37 * G_amp/us_param.nominalVinf * sin(3*pi*tr[bld_i]/gustT)  * (1.0 - cos(2*pi*tr[bld_i]/gustT))
                    Vinf_used[bld_i] = us_param.nominalVinf*IECGustFactor_blade

                else
                    Vinf_used[bld_i] = us_param.nominalVinf
                end
            end
            Vinf_used_long[ii,:] = Vinf_used[:]
        end

        # PyPlot.figure()
        # PyPlot.plot(rev_array, us_param.nominalVinf*IECGustFactor,"k-")
        # PyPlot.plot(rev_array, Vinf_used_long[:,1],"-",color = plot_cycle[1])
        # # PyPlot.plot(rev_array, Vinf_used_long[:,2],"-",color = plot_cycle[2])
        # # PyPlot.plot(rev_array, Vinf_used_long[:,3],"-",color = plot_cycle[3])
        # PyPlot.xlim([r_start,r_end])
        # PyPlot.ylim([5,12])
        # PyPlot.xlabel("Revolution")
        # PyPlot.ylabel("Freestream Velocity (m/s)")
        # PyPlot.legend(["Hub","Blade 1","Blade 2","Blade 3"],loc = "upper right")
        # # PyPlot.savefig("$(path)/figs/gust/SNL5m_blade_slice_$(slice)_gust_gtime$(us_param.gustT).pdf",transparent = true)

    end

    return RefArea,substep_level,ntheta,us_param,env,turbine2D,chord,Vinf,CP_filt,CP_RPI,alpha_filt,alpha_RPI,Rp_filt,Rp_RPI,Tp_filt,Tp_RPI,Zp_filt,Zp_RPI,Vinf_used_filt,Vinf_used_RPI,Tp_filt_blades,Tp_RPI_blades,Rp_filt_blades,Rp_RPI_blades

end

slice = 15
# function runme2(slice)
bnum = 1 # For plots
println("Running AC")
RefArea,substep_level,ntheta,us_param,env,turbine2D,chord,Vinf,CP_filt_AC,CP_RPI_AC,alpha_filt_AC,alpha_RPI_AC,Rp_filt_AC,Rp_RPI_AC,Tp_filt_AC,Tp_RPI_AC,Zp_filt_AC,Zp_RPI_AC,Vinf_used_filt_AC,Vinf_used_RPI_AC,Tp_filt_blades_AC,Tp_RPI_blades_AC,Rp_filt_blades_AC,Rp_RPI_blades_AC = runme("AC",slice,bnum)
println("Running DMS")
RefArea,substep_level,ntheta,us_param,env,turbine2D,chord,Vinf,CP_filt_DMS,CP_RPI_DMS,alpha_filt_DMS,alpha_RPI_DMS,Rp_filt_DMS,Rp_RPI_DMS,Tp_filt_DMS,Tp_RPI_DMS,Zp_filt_DMS,Zp_RPI_DMS,Vinf_used_filt_DMS,Vinf_used_RPI_DMS,Tp_filt_blades_DMS,Tp_RPI_blades_DMS,Rp_filt_blades_DMS,Rp_RPI_blades_DMS = runme("DMS",slice,bnum)

N_Rev = 40
# if true
println("Slice: $slice")

PyPlot.close("all")

t = collect(1:ntheta*(N_Rev)*substep_level)

data = DelimitedFiles.readdlm("$(path)/data/gust/RerunTestVAWT5_2_Gust2_Troposkein_ElementData.csv", ',',skipstart = 1)
revdata = DelimitedFiles.readdlm("$(path)/data/gust/RerunTestVAWT5_2_Gust2_Troposkein_RevData.csv", ',',skipstart = 1)
# Extract Data for Blade bnum
used_data_logic1 = (data[:,3].==1) .& (data[:,4].==slice)
used_data_logic2 = (data[:,3].==2) .& (data[:,4].==slice)
used_data_logic3 = (data[:,3].==3) .& (data[:,4].==slice)
data_full_rev1 = data[used_data_logic1,:]
data_full_rev2 = data[used_data_logic2,:]
data_full_rev3 = data[used_data_logic3,:]

rev_array = LinRange(2*pi/(turbine2D.ntheta)/2,N_Rev*2.0*pi,length(Tp_filt_AC))/(2*pi)
q_loc1 = 0.5*env.rho.*(data_full_rev1[:,15]*mean(Vinf)).^2
q_loc2 = 0.5*env.rho.*(data_full_rev2[:,15]*mean(Vinf)).^2
q_loc3 = 0.5*env.rho.*(data_full_rev3[:,15]*mean(Vinf)).^2

r_start = 20.0;
r_end = 25.0;

# # CP
# PyPlot.figure()
# PyPlot.plot(revdata[:,1],revdata[:,2],"-",color=plot_cycle[1],linewidth = 1.0)
# PyPlot.plot(rev_array,CP_filt,color = plot_cycle[2],linewidth = 1.0)
# PyPlot.plot(rev_array,CP_RPI,color = plot_cycle[1],linewidth = 1.0)
# PyPlot.ylim([0.0,0.6])
# PyPlot.xlabel("Revolution")
# PyPlot.ylabel("CP")
# # PyPlot.xlim([r_start,r_end])
# PyPlot.legend(["CACTUS", "$AeroModel 1st Order Filter","$AeroModel RPI"])
# PyPlot.savefig("$(path)/figs/gust/SNL5m_ALL_Transient_CP_$(slice)_gtime$(gustT).pdf",transparent = true)

# Angle of Attack
# PyPlot.figure()
# PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), data_full_rev1[:,9],"",color=plot_cycle[1],linewidth = 0.75)
# PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), data_full_rev1[:,10],"--",color=plot_cycle[1],linewidth = 0.75)
# PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), data_full_rev1[:,11],"-.",color=plot_cycle[1],linewidth = 0.75)
# PyPlot.plot(rev_array, alpha_filt_AC*180/pi,color=plot_cycle[2],linewidth = 0.75)
# PyPlot.plot(rev_array, alpha_RPI_AC*180/pi,"--",color=plot_cycle[2],linewidth = 0.75)
# PyPlot.plot(rev_array, alpha_filt_DMS*180/pi,color=plot_cycle[3],linewidth = 0.75)
# PyPlot.plot(rev_array, alpha_RPI_DMS*180/pi,"--",color=plot_cycle[3],linewidth = 0.75)
# PyPlot.xlabel("Theta (deg)")
# PyPlot.xlim([r_start,r_end])
# PyPlot.ylabel("AOA (deg)")
# if bnum == 1
#     PyPlot.legend(["CACTUS 25","CACTUS 50","CACTUS 75", "AC 1st Order Filter","AC RPI", "DMS 1st Order Filter", "DMS RPI"])
# end
# PyPlot.savefig("$(path)/figs/gust/SNL5m_ALL_Transient_alpha$(bnum)_$(slice)_gtime$(us_param.gustT).pdf",transparent = true)

# Rp
# PyPlot.figure()
# rpcactus = -q_loc1.*(data_full_rev1[:,24])#.-q_loc2.*(data_full_rev2[:,24]).-q_loc3.*(data_full_rev3[:,24])
# PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), rpcactus.*chord.*cos.(turbine2D.delta[1]),"-",color=plot_cycle[1],linewidth = 0.75) #Normal to radial by cos(delta)
# PyPlot.plot(rev_array,-Rp_filt_AC.*cos.(turbine2D.delta[1]),color = plot_cycle[2],linewidth = 0.75) #Remove span length by multiplying by cos(delta)
# PyPlot.plot(rev_array,-Rp_RPI_AC.*cos.(turbine2D.delta[1]),"--",color = plot_cycle[2],linewidth = 0.75) #Remove span length by multiplying by cos(delta)
# PyPlot.plot(rev_array,-Rp_filt_DMS.*cos.(turbine2D.delta[1]),color = plot_cycle[3],linewidth = 0.75) #Remove span length by multiplying by cos(delta)
# PyPlot.plot(rev_array,-Rp_RPI_DMS.*cos.(turbine2D.delta[1]),"--",color = plot_cycle[3],linewidth = 0.75) #Remove span length by multiplying by cos(delta)
# PyPlot.xlabel("Revolution")
# PyPlot.ylabel("Radial Force per Span (N/m)")
# PyPlot.xlim([r_start,r_end])
# PyPlot.legend(["CACTUS", "AC 1st Order Filter","AC RPI", "DMS 1st Order Filter", "DMS RPI"])
# PyPlot.savefig("$(path)/figs/gust/SNL5m_ALL_Transient_Rp$(bnum)_$(slice)_gtime$(us_param.gustT).pdf",transparent = true)

println("Radial")
idx_start = round(Int, length(Rp_filt_AC)/2)
twod_filt = -Rp_filt_AC[idx_start:end].*cos.(turbine2D.delta[1])
twod_RPI = -Rp_RPI_AC[idx_start:end].*cos.(turbine2D.delta[1])
maxerr = maximum(abs.((twod_filt .- twod_RPI)./twod_filt))
meanerr = mean(abs.((twod_filt .- twod_RPI)./twod_filt))
println("Max: $(maxerr)")
println("Mean: $(meanerr)")

# Rp Blades
PyPlot.figure()
PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), -q_loc1.*(data_full_rev1[:,24])*chord.*cos.(turbine2D.delta[1]),"-",color=plot_cycle[1],linewidth = 0.75,label="CACTUS1")
PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), -q_loc2.*(data_full_rev2[:,24])*chord.*cos.(turbine2D.delta[1]),"-",color=plot_cycle[2],linewidth = 0.75,label="CACTUS2")
PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), -q_loc3.*(data_full_rev3[:,24])*chord.*cos.(turbine2D.delta[1]),"-",color=plot_cycle[3],linewidth = 0.75,label="CACTUS3")
PyPlot.plot(rev_array, -Rp_filt_blades_AC[:,1].*cos.(turbine2D.delta[1]),".-",color = plot_cycle[1],linewidth = 0.75,label="AC 1st Order Filter1") #Remove span length by multiplying by cos(delta)
PyPlot.plot(rev_array, -Rp_filt_blades_AC[:,2].*cos.(turbine2D.delta[1]),".-",color = plot_cycle[2],linewidth = 0.75,label="AC 1st Order Filter2") #Remove span length by multiplying by cos(delta)
PyPlot.plot(rev_array, -Rp_filt_blades_AC[:,3].*cos.(turbine2D.delta[1]),".-",color = plot_cycle[3],linewidth = 0.75,label="AC 1st Order Filter3") #Remove span length by multiplying by cos(delta)
# PyPlot.plot(rev_array, -Rp_RPI_AC.*cos.(turbine2D.delta[1]),color = plot_cycle[2],linewidth = 0.75,label="AC RPI") #Remove span length by multiplying by cos(delta)
PyPlot.plot(rev_array, -Rp_filt_blades_DMS[:,1].*cos.(turbine2D.delta[1]),"--",color = plot_cycle[1],linewidth = 0.75,label="DMS 1st Order Filter1") #Remove span length by multiplying by cos(delta)
PyPlot.plot(rev_array, -Rp_filt_blades_DMS[:,2].*cos.(turbine2D.delta[1]),"--",color = plot_cycle[2],linewidth = 0.75,label="DMS 1st Order Filter2") #Remove span length by multiplying by cos(delta)
PyPlot.plot(rev_array, -Rp_filt_blades_DMS[:,3].*cos.(turbine2D.delta[1]),"--",color = plot_cycle[3],linewidth = 0.75,label="DMS 1st Order Filter3") #Remove span length by multiplying by cos(delta)
# PyPlot.plot(rev_array, -Rp_RPI_DMS.*cos.(turbine2D.delta[1]),color = plot_cycle[3],linewidth = 0.75,label="DMS RPI") #Remove span length by multiplying by cos(delta)
PyPlot.xlabel("Revolution")
PyPlot.ylabel("Radial Force per Span (N/m)")
PyPlot.xlim([r_start,r_end])
if slice == 15
    PyPlot.ylim([-2,40])
elseif slice == 7
    # PyPlot.ylim([-2,20])
end
if bnum == 1
    PyPlot.legend(loc=1,bbox_to_anchor=(0.7,1.16))
end
# PyPlot.savefig("$(path)/figs/gust/SNL5m_ALL_Transient_Rpblades_$(slice)_gtime$(us_param.gustT).pdf",transparent = true)


# Tp
PyPlot.rc("figure", figsize=(4, 3.1))
PyPlot.figure()
tpcactus = -q_loc1.*(data_full_rev1[:,25])#.-q_loc2.*(data_full_rev2[:,25]).-q_loc3.*(data_full_rev3[:,25])
# PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), tpcactus*chord,"-",color=plot_cycle[1],linewidth = 0.75,label="CACTUS")
PyPlot.plot(rev_array, -Tp_filt_AC.*cos.(turbine2D.delta[1]),"--",color = plot_cycle[2],linewidth = 0.75,label="AC 1st") #Remove span length by multiplying by cos(delta)
PyPlot.plot(rev_array, -Tp_RPI_AC.*cos.(turbine2D.delta[1]),color = plot_cycle[2],linewidth = 0.75,label="AC RPI") #Remove span length by multiplying by cos(delta)
PyPlot.plot(rev_array, -Tp_filt_DMS.*cos.(turbine2D.delta[1]),"--",color = plot_cycle[3],linewidth = 0.75,label="DMS 1st") #Remove span length by multiplying by cos(delta)
PyPlot.plot(rev_array, -Tp_RPI_DMS.*cos.(turbine2D.delta[1]),color = plot_cycle[3],linewidth = 0.75,label="DMS RPI") #Remove span length by multiplying by cos(delta)
PyPlot.xlabel("Revolution")
PyPlot.ylabel("Tangential Force per Span (N/m)")
PyPlot.xlim([r_start,r_end])
if slice == 15
    PyPlot.ylim([-2,40])
elseif slice == 7
    PyPlot.ylim([-1,12])
end
if bnum == 1
    PyPlot.legend(loc=1,bbox_to_anchor=(.42,1.05))
end
# # PyPlot.savefig("$(path)/figs/gust/SNL5m_ALL_Transient_Tp$(bnum)_$(slice)_gtime$(us_param.gustT).pdf",transparent = true)
# PyPlot.rc("figure", figsize=(4, 3))
# # Tp Blades
# PyPlot.figure()
# tpcactus = -q_loc1.*(data_full_rev1[:,25])#.-q_loc2.*(data_full_rev2[:,25]).-q_loc3.*(data_full_rev3[:,25])
# PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), -q_loc1.*(data_full_rev1[:,25])*chord,"-",color=plot_cycle[1],linewidth = 0.75,label="CACTUS1")
# PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), -q_loc2.*(data_full_rev2[:,25])*chord,"-",color=plot_cycle[2],linewidth = 0.75,label="CACTUS2")
# PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), -q_loc3.*(data_full_rev3[:,25])*chord,"-",color=plot_cycle[3],linewidth = 0.75,label="CACTUS3")
# PyPlot.plot(rev_array, -Tp_filt_blades_AC[:,1].*cos.(turbine2D.delta[1]),".-",color = plot_cycle[1],linewidth = 0.75,label="AC 1st Order Filter1") #Remove span length by multiplying by cos(delta)
# PyPlot.plot(rev_array, -Tp_filt_blades_AC[:,2].*cos.(turbine2D.delta[1]),".-",color = plot_cycle[2],linewidth = 0.75,label="AC 1st Order Filter2") #Remove span length by multiplying by cos(delta)
# PyPlot.plot(rev_array, -Tp_filt_blades_AC[:,3].*cos.(turbine2D.delta[1]),".-",color = plot_cycle[3],linewidth = 0.75,label="AC 1st Order Filter3") #Remove span length by multiplying by cos(delta)
# # PyPlot.plot(rev_array, -Tp_RPI_AC.*cos.(turbine2D.delta[1]),color = plot_cycle[2],linewidth = 0.75,label="AC RPI") #Remove span length by multiplying by cos(delta)
# PyPlot.plot(rev_array, -Tp_filt_blades_DMS[:,1].*cos.(turbine2D.delta[1]),"--",color = plot_cycle[1],linewidth = 0.75,label="DMS 1st Order Filter1") #Remove span length by multiplying by cos(delta)
# PyPlot.plot(rev_array, -Tp_filt_blades_DMS[:,2].*cos.(turbine2D.delta[1]),"--",color = plot_cycle[2],linewidth = 0.75,label="DMS 1st Order Filter2") #Remove span length by multiplying by cos(delta)
# PyPlot.plot(rev_array, -Tp_filt_blades_DMS[:,3].*cos.(turbine2D.delta[1]),"--",color = plot_cycle[3],linewidth = 0.75,label="DMS 1st Order Filter3") #Remove span length by multiplying by cos(delta)
# # PyPlot.plot(rev_array, -Tp_RPI_DMS.*cos.(turbine2D.delta[1]),color = plot_cycle[3],linewidth = 0.75,label="DMS RPI") #Remove span length by multiplying by cos(delta)
# PyPlot.xlabel("Revolution")
# PyPlot.ylabel("Tangential Force per Span (N/m)")
# PyPlot.xlim([r_start,r_end])
# if slice == 15
#     PyPlot.ylim([-2,40])
# elseif slice == 7
#     # PyPlot.ylim([-2,20])
# end
# if bnum == 1
#     PyPlot.legend(loc=1,bbox_to_anchor=(0.7,1.16))
# end
# # PyPlot.savefig("$(path)/figs/gust/SNL5m_ALL_Transient_Tpblades_$(slice)_gtime$(us_param.gustT).pdf",transparent = true)
#
#
# println("Tangential")
# twod_filt = -Tp_filt_AC[idx_start:end].*cos.(turbine2D.delta[1])
# twod_RPI = -Tp_RPI_AC[idx_start:end].*cos.(turbine2D.delta[1])
# maxerr = maximum(abs.((twod_filt .- twod_RPI)./twod_filt))
# meanerr = mean(abs.((twod_filt .- twod_RPI)./twod_filt))
# println("Max: $(maxerr)")
# println("Mean: $(meanerr)")
#
# # Zp
# PyPlot.figure()
# zpcactus = -q_loc1.*(data_full_rev1[:,24])#.-q_loc2.*(data_full_rev2[:,24]).-q_loc3.*(data_full_rev3[:,24])
# PyPlot.plot((data_full_rev1[:,2] .- minimum(data_full_rev1[:,2]))/(2*pi), zpcactus.*chord.*sin.(turbine2D.delta[1]),"-",color=plot_cycle[1],linewidth = 0.75) #Normal to vertical by sin(delta)
# PyPlot.plot(rev_array,-Zp_filt_AC.*cos.(turbine2D.delta[1]),color = plot_cycle[2],linewidth = 0.75) #Remove span length by multiplying by cos(delta)
# PyPlot.plot(rev_array,-Zp_RPI_AC.*cos.(turbine2D.delta[1]),"--",color = plot_cycle[2],linewidth = 0.75) #Remove span length by multiplying by cos(delta)
# PyPlot.plot(rev_array,-Zp_filt_DMS.*cos.(turbine2D.delta[1]),color = plot_cycle[3],linewidth = 0.75) #Remove span length by multiplying by cos(delta)
# PyPlot.plot(rev_array,-Zp_RPI_DMS.*cos.(turbine2D.delta[1]),"--",color = plot_cycle[3],linewidth = 0.75) #Remove span length by multiplying by cos(delta)
# PyPlot.xlabel("Revolution")
# PyPlot.ylabel("Vertical Force per Span (N/m)")
# PyPlot.xlim([r_start,r_end])
# # PyPlot.legend(["CACTUS", "AC 1st Order Filter","AC RPI", "DMS 1st Order Filter", "DMS RPI"])
# # PyPlot.savefig("$(path)/figs/gust/SNL5m_ALL_Transient_Zp$(bnum)_$(slice)_gtime$(us_param.gustT).pdf",transparent = true)
#
# println("Vertical")
# twod_filt = -Zp_filt_AC[idx_start:end].*cos.(turbine2D.delta[1])
# twod_RPI = -Zp_RPI_AC[idx_start:end].*cos.(turbine2D.delta[1])
# maxerr = maximum(abs.((twod_filt .- twod_RPI)./twod_filt))
# meanerr = mean(abs.((twod_filt .- twod_RPI)./twod_filt))
# println("Max: $(maxerr)")
# println("Mean: $(meanerr)")
#
# # end
# # end
# # println("starting...")
# #
# # println("running slice 7")
# # runme2(7)
# #
# # println("running slice 15")
# # runme2(15)
# # # Juno.@enter runme2(7)
