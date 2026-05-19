using Test
import DelimitedFiles
import HDF5
import FLOWMath
import Statistics: mean

import OWENSAero
# include("../src/OWENSAero.jl")
path,_ = splitdir(@__FILE__)
tol = 1e-4

function akima_without_duplicate_abscissae(x, y)
    xvec = Float64.(collect(x))
    yvec = Float64.(collect(y))
    order = sortperm(xvec)
    xsorted = xvec[order]
    ysorted = yvec[order]
    unique_x = Float64[]
    unique_y = Float64[]
    i = firstindex(xsorted)
    while i <= lastindex(xsorted)
        j = i
        while j < lastindex(xsorted) && xsorted[j+1] == xsorted[i]
            j += 1
        end
        push!(unique_x, xsorted[i])
        push!(unique_y, sum(ysorted[i:j]) / (j - i + 1))
        i = j + 1
    end
    return FLOWMath.Akima(unique_x, unique_y)
end

function interpolate_monotonic(x, y, x_query)
    idx = searchsortedlast(x, x_query)
    idx == length(x) && return y[end]
    idx == 0 && return y[1]
    fraction = (x_query - x[idx]) / (x[idx+1] - x[idx])
    return y[idx] + fraction * (y[idx+1] - y[idx])
end

function branch_validation_metrics(alpha_model_deg, coefficient_model, reference_data; branch)
    max_alpha_idx = argmax(alpha_model_deg)
    min_alpha_idx = argmin(alpha_model_deg)
    reference_max_alpha_idx = argmax(reference_data[:, 1])

    if branch === :upstroke
        model_alpha = alpha_model_deg[1:max_alpha_idx]
        model_coefficient = coefficient_model[1:max_alpha_idx]
        reference_branch = reference_data[1:reference_max_alpha_idx, :]
    elseif branch === :downstroke
        model_alpha = reverse(alpha_model_deg[max_alpha_idx:min_alpha_idx])
        model_coefficient = reverse(coefficient_model[max_alpha_idx:min_alpha_idx])
        reference_branch = reference_data[reference_max_alpha_idx:end, :]
    else
        throw(ArgumentError("branch must be :upstroke or :downstroke"))
    end

    overlap = findall(
        alpha -> first(model_alpha) <= alpha <= last(model_alpha),
        reference_branch[:, 1],
    )
    isempty(overlap) && throw(ArgumentError("reference branch does not overlap model branch"))

    reference_alpha = reference_branch[overlap, 1]
    reference_coefficient = reference_branch[overlap, 2]
    model_on_reference = [
        interpolate_monotonic(model_alpha, model_coefficient, alpha)
        for alpha in reference_alpha
    ]
    error = model_on_reference .- reference_coefficient
    abs_error = abs.(error)
    return (
        n = length(reference_alpha),
        rmse = sqrt(mean(error .^ 2)),
        mean_bias = mean(error),
        max_abs_error = maximum(abs_error),
    )
end

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

        afcl = akima_without_duplicate_abscissae(full_alpha_cl*pi/180, full_cl)
        afcm = akima_without_duplicate_abscissae(full_alpha_cm*pi/180, full_cm)
        afcd = akima_without_duplicate_abscissae(af_re3_6e6[:,1]*pi/180, af_re3_6e6[:,3])

        function af(alpha,Re,mach,family_factor)

            cl = afcl(alpha)
            cm = afcm(alpha)
            cd = afcd(alpha)

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

    CL_BV_paper = DelimitedFiles.readdlm(
        "$(path)/data/dynstall/Fig10_BV_CL.txt",
        ',',
        Float64,
    )
    CM_BV_paper = DelimitedFiles.readdlm(
        "$(path)/data/dynstall/CM_Fig10_BV_Sadr.txt",
        ',',
        Float64,
    )

    # Perform pinned regression comparisons against the checked-in dynamic-stall
    # fixture. These are intentionally direct array checks rather than loose
    # "is real" assertions, because small changes in the dynamic flags or prior
    # alpha state can otherwise pass unnoticed.
    @test alpha_BV isa Vector{Float64}
    @test CL_BV isa Vector{Float64}
    @test CD_BV isa Vector{Float64}
    @test CM_BV isa Vector{Float64}
    @test length(alpha_BV) == 1501
    @test length(CL_BV) == 1501
    @test length(CD_BV) == 1501
    @test length(CM_BV) == 1501
    @test all(isfinite, CL_BV)
    @test all(isfinite, CD_BV)
    @test all(isfinite, CM_BV)
    @test isapprox(alpha_BV_old, alpha_BV; atol=tol, rtol=0)
    @test isapprox(CL_BV_old, CL_BV; atol=tol, rtol=0)
    @test isapprox(CD_BV_old, CD_BV; atol=tol, rtol=0)
    @test isapprox(CM_BV_old, CM_BV; atol=tol, rtol=0)
    @test maximum(CL_BV) ≈ 1.6453597607493362 atol=tol
    @test minimum(CL_BV) ≈ 0.003243925119865854 atol=tol
    @test maximum(CD_BV) ≈ 0.3136503672765179 atol=tol
    @test minimum(CM_BV) ≈ -0.10459421783230137 atol=tol
    @test CL_BV[1] ≈ 1.1116093542418823 atol=tol
    @test CD_BV[1] ≈ 0.0184 atol=tol
    @test CM_BV[end] ≈ -0.001917210444783475 atol=tol

    alpha_BV_deg = alpha_BV .* 180 ./ pi
    cl_upstroke_metrics =
        branch_validation_metrics(alpha_BV_deg, CL_BV, CL_BV_paper; branch = :upstroke)
    cl_downstroke_metrics =
        branch_validation_metrics(alpha_BV_deg, CL_BV, CL_BV_paper; branch = :downstroke)
    cm_upstroke_metrics =
        branch_validation_metrics(alpha_BV_deg, CM_BV, CM_BV_paper; branch = :upstroke)
    cm_downstroke_metrics =
        branch_validation_metrics(alpha_BV_deg, CM_BV, CM_BV_paper; branch = :downstroke)

    @test cl_upstroke_metrics.n == 15
    @test cl_upstroke_metrics.rmse ≈ 0.04276727373228908 atol=tol
    @test cl_upstroke_metrics.mean_bias ≈ -0.018794652588474445 atol=tol
    @test cl_upstroke_metrics.max_abs_error ≈ 0.09946657595093722 atol=tol

    @test cl_downstroke_metrics.n == 14
    @test cl_downstroke_metrics.rmse ≈ 0.019404808461011173 atol=tol
    @test cl_downstroke_metrics.mean_bias ≈ 0.012941709005760067 atol=tol
    @test cl_downstroke_metrics.max_abs_error ≈ 0.039470874818617885 atol=tol

    @test cm_upstroke_metrics.n == 18
    @test cm_upstroke_metrics.rmse ≈ 0.070528150967298 atol=tol
    @test cm_upstroke_metrics.mean_bias ≈ 0.048736483519612094 atol=tol
    @test cm_upstroke_metrics.max_abs_error ≈ 0.14983040863490035 atol=tol

    @test cm_downstroke_metrics.n == 25
    @test cm_downstroke_metrics.rmse ≈ 0.0545389720609115 atol=tol
    @test cm_downstroke_metrics.mean_bias ≈ -0.026397059261107584 atol=tol
    @test cm_downstroke_metrics.max_abs_error ≈ 0.12274142264021343 atol=tol

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
