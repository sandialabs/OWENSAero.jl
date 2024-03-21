using Test
import DelimitedFiles
import Dierckx: Spline1D, evaluate
import HDF5
# import PyPlot
# PyPlot.close("all")

import OWENSAero
# include("../src/OWENSAero.jl")
path,_ = splitdir(@__FILE__)
tol = 1e-4

@testset "NACA 0012 Boeing-Vertol" begin

    # Data taken from https://doi.org/10.1007/s40430-019-1907-4 figure 10
    function runme()
        k = 0.048 # k = omega*c/(2*U)
        Re = 3.8e6 # Re = rho * l * V / mu
        rho = 1.225
        mu = 1.7894e-5
        M = 0.3
        v_sound = 343.0 #m/s
        U = M * v_sound
        c = Re * mu / (rho * U)
        w = k*2*U/c
        a_amp = 10.0*pi/180
        a_mean = 10.0*pi/180
        dt = 0.002
        tmax = 3.0
        aoaStallPos = 10*pi/180
        aoaStallNeg = -10*pi/180
        AOA0 = 0.0
        tc = 0.12
        BV_DynamicFlagL = false
        BV_DynamicFlagD = false
        n_samples = round(Int,tmax/dt)+1
        alpha_BV = ones(n_samples+1).*a_mean
        CL_BV = zeros(n_samples)
        CD_BV = zeros(n_samples)
        CM_BV = zeros(n_samples)

        # Load Airfoil Data
        af_re3_6e6 = DelimitedFiles.readdlm("$(path)/data/dynstall/NACA_0012.dat", ' ',Float64,skipstart = 13)
        # Data taken from https://doi.org/10.1007/s40430-019-1907-4 figure 1: The airfoil data used significantly changes the stall characteristics, must use same airfoil data as article for 1:1 comparison.
        cl_sadr = DelimitedFiles.readdlm("$(path)/data/dynstall/CL_NACA0012_Sadr.txt", ',',Float64,skipstart = 0)

        full_alpha_cl = [-reverse(cl_sadr[2:end,1]);cl_sadr[:,1]]
        full_cl = [-reverse(cl_sadr[2:end,2]);cl_sadr[:,2]]

        cm_sadr = DelimitedFiles.readdlm("$(path)/data/dynstall/CM_NACA0012_Sadr.txt", ',',Float64,skipstart = 0)

        full_alpha_cm = [-reverse(cm_sadr[2:end,1]);cm_sadr[:,1]]
        full_cm = [-reverse(cm_sadr[2:end,2]);cm_sadr[:,2]]

        # PyPlot.figure()
        # PyPlot.plot(full_alpha_cl,full_cl)
        # PyPlot.plot(full_alpha_cm,full_cm)

        afcl = Spline1D(full_alpha_cl*pi/180, full_cl, s=0.1)
        afcm = Spline1D(full_alpha_cm*pi/180, full_cm, s=0.1)
        afcd = Spline1D(af_re3_6e6[:,1]*pi/180, af_re3_6e6[:,3], s=0.01)

        function af(alpha,Re,mach,family_factor)

            cl = evaluate(afcl, alpha)
            cm = evaluate(afcm, alpha)
            cd = evaluate(afcd, alpha)

            return cl, cd, cm
        end

        # Run the BV Model
        i = 1
        for t = 0:dt:tmax
            alpha_BV[i+1] = a_mean+a_amp*sin(w*t)
            dalpha=alpha_BV[i+1]-alpha_BV[i]
            adotnorm=dalpha/dt*c/(2.0*max(U,0.001))
            CL_BV[i], CD_BV[i], CM_BV[i], BV_DynamicFlagL, BV_DynamicFlagD = OWENSAero.Boeing_Vertol(af,alpha_BV[i+1],adotnorm,M,Re,aoaStallPos,aoaStallNeg,AOA0,tc,BV_DynamicFlagL,BV_DynamicFlagD, family_factor = 0.0)
            i += 1
        end
        return af_re3_6e6,full_alpha_cl,full_cl,full_alpha_cm,full_cm, alpha_BV[2:end], CL_BV, CD_BV, CM_BV
    end

    # Juno.@enter runme()
    af_re3_6e6,full_alpha_cl,full_cl,full_alpha_cm,full_cm, alpha_BV, CL_BV, CD_BV, CM_BV = runme()

    # Load the Experimental and Comparative BV Results
    CL_exp = DelimitedFiles.readdlm("$(path)/data/dynstall/Fig10_Exp_CL.txt", ',',skipstart = 1)
    CL_BV_paper = DelimitedFiles.readdlm("$(path)/data/dynstall/Fig10_BV_CL.txt", ',',skipstart = 1)

    CM_exp = DelimitedFiles.readdlm("$(path)/data/dynstall/CM_Fig10_EXP_Sadr.txt", ',',skipstart = 0)
    CM_BV_paper = DelimitedFiles.readdlm("$(path)/data/dynstall/CM_Fig10_BV_Sadr.txt", ',',skipstart = 0)

    #Unit Data
    filename = "$path/data/dynstall/BV_unit_data.h5"
    # HDF5.h5open(filename, "w") do file
    #     HDF5.write(file,"alpha_BV",alpha_BV)
    #     HDF5.write(file,"CL_BV",CL_BV)
    #     HDF5.write(file,"CD_BV",CD_BV)
    #     HDF5.write(file,"CM_BV",CM_BV)
    # end

    alpha_BV_old = HDF5.h5read(filename,"alpha_BV")
    CL_BV_old = HDF5.h5read(filename,"CL_BV")
    CD_BV_old = HDF5.h5read(filename,"CD_BV")
    CM_BV_old = HDF5.h5read(filename,"CM_BV")

    # Perform Testing Comparisons
    # PyPlot.figure()
    # PyPlot.plot(LinRange(0,1,length(alpha_BV)),alpha_BV,"r-")
    # PyPlot.plot(LinRange(0,1,length(alpha_BV_old)),alpha_BV_old,"b-")

    # PyPlot.figure()
    # PyPlot.plot(LinRange(0,1,length(CL_BV)),CL_BV,"r-")
    # PyPlot.plot(LinRange(0,1,length(CL_BV_old)),CL_BV_old,"b-")

    # PyPlot.figure()
    # PyPlot.plot(LinRange(0,1,length(CD_BV)),CD_BV,"r-")
    # PyPlot.plot(LinRange(0,1,length(CD_BV_old)),CD_BV_old,"b-")

    # PyPlot.figure()
    # PyPlot.plot(LinRange(0,1,length(CM_BV)),CM_BV,"r-")
    # PyPlot.plot(LinRange(0,1,length(CM_BV_old)),CM_BV_old,"b-")

    for ii = 1:length(alpha_BV_old)
        @test isapprox(alpha_BV_old[ii],alpha_BV[ii],atol=tol)
        @test isapprox(CL_BV_old[ii],CL_BV[ii],atol=tol)
        @test isapprox(CD_BV_old[ii],CD_BV[ii],atol=tol)
        @test isapprox(CM_BV_old[ii],CM_BV[ii],atol=tol)
    end

#     # Plot
#     plots = true
#     if plots
#         # import PyPlot
#         PyPlot.rc("figure", figsize=(4, 3))
#         PyPlot.rc("font", size=10.0)
#         PyPlot.rc("lines", linewidth=1.5)
#         PyPlot.rc("lines", markersize=3.0)
#         PyPlot.rc("legend", frameon=false)
#         PyPlot.rc("axes.spines", right=false, top=false)
#         PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
#         # rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
#         plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]
#
#         PyPlot.figure()
#         PyPlot.plot(full_alpha_cl,full_cl,"k.-")
#         PyPlot.plot(CL_exp[:,1],CL_exp[:,2],".-",color = plot_cycle[1])
#         PyPlot.plot(CL_BV_paper[:,1],CL_BV_paper[:,2],".-",color = plot_cycle[2])
#         PyPlot.plot(alpha_BV*180/pi, CL_BV,".-",color = plot_cycle[3])
#         PyPlot.arrow(7.0,1.1,3.5,0.4,head_width=0.05,width=0.015,head_length=0.7,overhang=0.5,head_starts_at_zero="true",facecolor="k")
#         PyPlot.arrow(18.0,0.4,-5.0,0.0,head_width=0.05,width=0.015,head_length=0.7,overhang=0.5,head_starts_at_zero="true",facecolor="k")
#         PyPlot.xlabel("Angle of Attack (deg)")
#         PyPlot.ylabel("Lift Coefficient")
#         PyPlot.xlim([-22,22])
#         PyPlot.legend(["Static Sadr","Exp Dynamic","Sadr BV","CACTUS BV"])
#         PyPlot.savefig("$(path)/../doc/paper/analyses/figs/dynstall/CL_Fig10_BV.pdf",transparent = true)
#
#         PyPlot.figure()
#         PyPlot.plot(full_alpha_cm,full_cm,"k.-")
#         PyPlot.plot(CM_exp[:,1],CM_exp[:,2],".-",color = plot_cycle[1])
#         PyPlot.plot(CM_BV_paper[:,1],CM_BV_paper[:,2],".-",color = plot_cycle[2])
#         PyPlot.plot(alpha_BV*180/pi, CM_BV,".-",color = plot_cycle[3])
#         PyPlot.xlabel("Angle of Attack (deg)")
#         PyPlot.ylabel("Moment Coefficient")
#         PyPlot.xlim([0,22])
#         PyPlot.legend(["Static Sadr","Exp Dynamic","Sadr BV","CACTUS BV"])
#         PyPlot.savefig("$(path)/../doc/paper/analyses/figs/dynstall/CM_Fig10_BV.pdf",transparent = true)
#
#         PyPlot.figure()
#         PyPlot.plot(af_re3_6e6[:,1],af_re3_6e6[:,3],"k.-")
#         # PyPlot.plot(CD_exp[:,1],CD_exp[:,2],".-",color = plot_cycle[1])
#         # PyPlot.plot(CD_BV_paper[:,1],CD_BV_paper[:,2],".-",color = plot_cycle[2])
#         PyPlot.plot(alpha_BV*180/pi, CD_BV,".-",color = plot_cycle[3])
#         PyPlot.arrow(17.0,0.1,2.0,0.1,head_width=0.02,width=0.005,head_length=0.3,overhang=0.3,head_starts_at_zero="true",facecolor="k")
#         PyPlot.arrow(14.0,0.28,-2.0,-0.1,head_width=0.02,width=0.005,head_length=0.3,overhang=0.3,head_starts_at_zero="true",facecolor="k")
#         PyPlot.xlabel("Angle of Attack (deg)")
#         PyPlot.ylabel("Drag Coefficient")
#         PyPlot.xlim([0,22])
#         PyPlot.ylim([0,0.4])
#         PyPlot.legend(["Static Xfoil","CACTUS BV"])
#         PyPlot.savefig("$(path)/../doc/paper/analyses/figs/dynstall/CD_Fig10_BV.pdf",transparent = true)
#     end
#
end

# @testset "NACA 0012 Leishman-Beddoes" begin
#
#     # Data taken from https://doi.org/10.1007/s40430-019-1907-4 figure 10
#     function runme2()
#         k = 0.048 # k = omega*c/(2*U)
#         Re = 3.8e6 # Re = rho * l * V / mu
#         rho = 1.225
#         mu = 1.7894e-5
#         M = 0.3
#         v_sound = 343.0 #m/s
#         U = M * v_sound
#         c = Re * mu / (rho * U)
#         w = k*2*U/c
#         a_amp = 10.0*pi/180
#         a_mean = 10.0*pi/180
#         dt = 0.002
#         tmax = 3.0
#         CLCritP = 13.25 * pi/180 # stall angle
#         CLCritN = -13.25 * pi/180 # stall angle
#         cv_Last = 0.0
#         CLa = 2.0*pi
#         AOA0 = 0.0
#         tc = 0.12
#         ds = 2.0*U*dt/c
#         CLRefLE_Last = 0.0
#         CLRef_Last = 0.0
#         Fstat_Last = 0.0
#         cv_Last = 0.0
#         dp = 0.0
#         dF = 0.0
#         dCNv = 0.0
#         LESepState = 0
#         sLEv = 0.0
#         n_samples = round(Int,tmax/dt)+1
#         alpha_LB = zeros(n_samples)
#         CL_LB = zeros(n_samples)
#         CD_LB = zeros(n_samples)
#         CM_LB = zeros(n_samples)
#
#         # Load Airfoil Data
#         af_re3_6e6 = DelimitedFiles.readdlm("$(path)/data/dynstall/NACA_0012.dat", ' ',Float64,skipstart = 13)
#         # Data taken from https://doi.org/10.1007/s40430-019-1907-4 figure 1: The airfoil data used significantly changes the stall characteristics, must use same airfoil data as article for 1:1 comparison.
#         cl_sadr = DelimitedFiles.readdlm("$(path)/data/dynstall/CL_NACA0012_Sadr.txt", ',',Float64,skipstart = 0)
#         full_alpha_cl = [-reverse(cl_sadr[2:end,1]);cl_sadr[:,1]]
#         full_cl = [-reverse(cl_sadr[2:end,2]);cl_sadr[:,2]]
#
#         cm_sadr = DelimitedFiles.readdlm("$(path)/data/dynstall/CM_NACA0012_Sadr.txt", ',',Float64,skipstart = 0)
#
#         full_alpha_cm = [-reverse(cm_sadr[2:end,1]);cm_sadr[:,1]]
#         full_cm = [-reverse(cm_sadr[2:end,2]);cm_sadr[:,2]]
#
#         afcl = Spline1D(full_alpha_cl*pi/180, full_cl, s=0.1)
#         afcm = Spline1D(full_alpha_cm*pi/180, full_cm, s=0.1)
#         afcd = Spline1D(af_re3_6e6[:,1]*pi/180, af_re3_6e6[:,3], s=0.01)
#
#         function af(alpha,Re,mach,family_factor)
#
#             cl = evaluate(afcl, alpha)
#             cm = evaluate(afcm, alpha)
#             cd = evaluate(afcd, alpha)
#
#             return cl, cd, cm
#         end
#
#         # Run the LB Model
#         i = 1
#         for t = 0:dt:tmax
#             alpha_LB[i] = a_mean+a_amp*sin(w*t)
#             # Juno.@enter OWENSAero.Leishman_Beddoes!(af,alpha_LB[i],alpha_LB[i],ds,M,Re,CLa,AOA0,CLCritP,CLCritN,CLRefLE_Last, CLRef_Last,Fstat_Last, cv_Last, dp, dF, dCNv, LESepState, Tp=1.7, family_factor = 0.0)
#             CL_LB[i], CD_LB[i], CM_LB[i], CLRefLE_Last,CLRef_Last,Fstat_Last,cv_Last, dp, dF, dCNv, LESepState, sLEv = OWENSAero.Leishman_Beddoes(af,alpha_LB[i],alpha_LB[i],ds,M,Re,CLa,AOA0,CLCritP,CLCritN,CLRefLE_Last, CLRef_Last,Fstat_Last, cv_Last, dp, dF, dCNv, LESepState, sLEv, Tp=1.7, family_factor = 0.0)
#             i += 1
#         end
#         return af_re3_6e6,full_alpha_cl,full_cl,full_alpha_cm,full_cm, alpha_LB, CL_LB, CD_LB, CM_LB
#     end
#
#     af_re3_6e6,full_alpha_cl,full_cl,full_alpha_cm,full_cm, alpha_LB, CL_LB, CD_LB, CM_LB = runme2()
#
#     # Load the Experimental and Comparative BV Results
#     CL_exp = DelimitedFiles.readdlm("$(path)/data/dynstall/Fig10_Exp_CL.txt", ',',skipstart = 1)
#     CL_LB_paper = DelimitedFiles.readdlm("$(path)/data/dynstall/Fig10_LB_CL.txt", ',',skipstart = 1)
#
#     CM_exp = DelimitedFiles.readdlm("$(path)/data/dynstall/CM_Fig10_EXP_Sadr.txt", ',',skipstart = 0)
#     CM_LB_paper = DelimitedFiles.readdlm("$(path)/data/dynstall/CM_Fig10_LB_Sadr.txt", ',',skipstart = 0)
#
#     #Unit Data
#     filename = "$path/data/dynstall/LB_unit_data.h5"
#     # HDF5.h5open(filename, "w") do file
#     #     HDF5.write(file,"alpha_LB",alpha_LB)
#     #     HDF5.write(file,"CL_LB",CL_LB)
#     #     HDF5.write(file,"CD_LB",CD_LB)
#     #     HDF5.write(file,"CM_LB",CM_LB)
#     # end
#
#     alpha_LB_old = HDF5.h5read(filename,"alpha_LB")
#     CL_LB_old = HDF5.h5read(filename,"CL_LB")
#     CD_LB_old = HDF5.h5read(filename,"CD_LB")
#     CM_LB_old = HDF5.h5read(filename,"CM_LB")
#
#     # Perform Testing Comparisons
#     @test isapprox(alpha_LB_old,alpha_LB,atol=tol)
#     @test isapprox(CL_LB_old,CL_LB,atol=tol)
#     @test isapprox(CD_LB_old,CD_LB,atol=tol)
#     @test isapprox(CM_LB_old,CM_LB,atol=tol)
#
#     # Plot
#     plots = true
#     if plots
#         # import PyPlot
#         PyPlot.rc("figure", figsize=(4, 3))
#         PyPlot.rc("font", size=10.0)
#         PyPlot.rc("lines", linewidth=1.5)
#         PyPlot.rc("lines", markersize=3.0)
#         PyPlot.rc("legend", frameon=false)
#         PyPlot.rc("axes.spines", right=false, top=false)
#         PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
#         # rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
#         plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]
#
#         PyPlot.figure()
#         PyPlot.plot(full_alpha_cl,full_cl,"k.-")
#         PyPlot.plot(CL_exp[:,1],CL_exp[:,2],".-",color = plot_cycle[1])
#         PyPlot.plot(CL_LB_paper[:,1],CL_LB_paper[:,2],".-",color = plot_cycle[2])
#         PyPlot.plot(alpha_LB[100:end]*180/pi, CL_LB[100:end],".-",color = plot_cycle[3])
#         PyPlot.xlabel("Angle of Attack (deg)")
#         PyPlot.ylabel("Lift Coefficient")
#         PyPlot.xlim([-22,22])
#         PyPlot.legend(["Static","Exp Dynamic","Sadr LB","CACTUS LB"])
#         PyPlot.savefig("$(path)/../doc/paper/analyses/figs/dynstall/CL_Fig10_LB.pdf",transparent = true)
#
#         PyPlot.figure()
#         PyPlot.plot(full_alpha_cm,full_cm,"k.-")
#         PyPlot.plot(CM_exp[:,1],CM_exp[:,2],".-",color = plot_cycle[1])
#         PyPlot.plot(CM_LB_paper[:,1],CM_LB_paper[:,2],".-",color = plot_cycle[2])
#         PyPlot.plot(alpha_LB[100:end]*180/pi, CM_LB[100:end],".-",color = plot_cycle[3])
#         PyPlot.xlabel("Angle of Attack (deg)")
#         PyPlot.ylabel("Moment Coefficient")
#         PyPlot.xlim([0,22])
#         PyPlot.legend(["Static","Exp Dynamic","Sadr LB","CACTUS LB"])
#         PyPlot.savefig("$(path)/../doc/paper/analyses/figs/dynstall/CM_Fig10_LB.pdf",transparent = true)
#
#         PyPlot.figure()
#         PyPlot.plot(af_re3_6e6[:,1],af_re3_6e6[:,3],"k.-")
#         # PyPlot.plot(CD_exp[:,1],CD_exp[:,2],".-",color = plot_cycle[1])
#         # PyPlot.plot(CD_LB_paper[:,1],CD_LB_paper[:,2],".-",color = plot_cycle[2])
#         PyPlot.plot(alpha_LB[100:end]*180/pi, CD_LB[100:end],".-",color = plot_cycle[3])
#         PyPlot.xlabel("Angle of Attack (deg)")
#         PyPlot.ylabel("Drag Coefficient")
#         PyPlot.xlim([0,22])
#         PyPlot.ylim([0,0.4])
#         PyPlot.legend(["Static Xfoil","CACTUS LB"])
#         PyPlot.savefig("$(path)/../doc/paper/analyses/figs/dynstall/CD_Fig10_LB.pdf",transparent = true)
#     end
#
# end
