using Test
import HDF5
import ForwardDiff
import FiniteDiff
import Statistics:mean
import OWENSOpenFASTWrappers
# import RollingFunctions
#
# ENV["MPLBACKEND"]="qt5agg"
# import PyPlot
# PyPlot.close("all")
# PyPlot.rc("figure", figsize=(4, 3))
# PyPlot.rc("font", size=10.0)
# PyPlot.rc("lines", linewidth=1.5)
# PyPlot.rc("lines", markersize=3.0)
# PyPlot.rc("legend", frameon=false)
# PyPlot.rc("axes.spines", right=false, top=false)
# PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
# PyPlot.rc("figure",max_open_warning=500)
# # rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

const UNSTEADY_BASELINES = Dict(
    ("AC", true) => (Rp1=14.61798502718863, Tp1=-1.3578451054187242, Rpmean=5.478111754622732, Tpmean=6.599005743637965, Vloc1=42.65570346549178, Re1=445031.47870641906),
    ("AC", false) => (Rp1=12.529927964477698, Tp1=-1.408546292257496, Rpmean=4.660093508900539, Tpmean=4.321887489076866, Vloc1=42.66169039146707, Re1=445093.9409401469),
    ("DMS", true) => (Rp1=17.8724152206306, Tp1=-1.2526644527903148, Rpmean=6.733379018517041, Tpmean=8.114827894681484, Vloc1=42.5286710408061, Re1=443706.1359454618),
    ("DMS", false) => (Rp1=17.8724152206306, Tp1=-1.2526644527903148, Rpmean=18.419804953792056, Tpmean=4.136206532176461, Vloc1=42.5286710408061, Re1=443706.1359454618),
)

const IFW_UNSTEADY_BASELINES = Dict(
    "AC" => (Rp1=-2.6331031705280084, Tp1=-1.2559507588560983, Rpmean=2.3980398649584904, Tpmean=0.30031142795009114, Vloc1=38.56407465412821, Re1=402343.0813221861),
    "DMS" => (Rp1=-3.958783356071784, Tp1=-1.272566646170957, Rpmean=4.65480460237821, Tpmean=-0.017477189310114483, Vloc1=38.95222326170589, Re1=406392.6769156071),
)

function test_unsteady_outputs(AeroModel, RPI, Rp, Tp, Zp, Vloc, Re, baseline)
    expected_steps = 2 * ntheta

    @test length(Rp) == expected_steps
    @test length(Tp) == expected_steps
    @test length(Zp) == expected_steps
    @test length(Vloc) == expected_steps
    @test length(Re) == expected_steps
    @test all(isfinite, Float64.(Rp))
    @test all(isfinite, Float64.(Tp))
    @test all(isfinite, Float64.(Zp))
    @test all(isfinite, Float64.(Vloc))
    @test all(isfinite, Float64.(Re))
    @test minimum(Float64.(Vloc)) > 0
    @test minimum(Float64.(Re)) > 0
    @test maximum(abs.(Float64.(Rp))) > 1
    @test maximum(abs.(Float64.(Tp))) > 1

    @test Rp[1] ≈ baseline.Rp1 rtol=1e-7 atol=1e-7
    @test Tp[1] ≈ baseline.Tp1 rtol=1e-7 atol=1e-7
    @test mean(Rp) ≈ baseline.Rpmean rtol=1e-7 atol=1e-7
    @test mean(Tp) ≈ baseline.Tpmean rtol=1e-7 atol=1e-7
    @test Vloc[1] ≈ baseline.Vloc1 rtol=1e-7 atol=1e-7
    @test Re[1] ≈ baseline.Re1 rtol=1e-7 atol=1e-4
end

path,_ = splitdir(@__FILE__)
import OWENSAero
# include("$(path)/../src/OWENSAero.jl")

function aerowrapper(x;RPI=true,returnall=false,windangle_D=0.0,AeroModel="DMS",steady=true,ifw=false,DynamicStallModel="BV")

    ntheta = 30
    RPM = x[1]
    R = 4.55/2 #m
    r = ones(Real,ntheta)*R #m
    chord = 0.1524 #m
    twist = ones(Real,ntheta) * 0.0 #rad
    delta = ones(Real,ntheta) * 0.0*pi/180 #rad
    omega = ones(Real,ntheta) * RPM / 60.0 * 2*pi # RPM -> radians/sec
    B = 3

    function affun(alpha, Re, M;V_twist=nothing,chord=nothing,dt=nothing,Vloc=nothing)
        cl = 6.2*alpha
        cd = 0.008 .- 0.003.*cl + 0.01.*cl.*cl
        return cl, cd
    end

    if DynamicStallModel=="BV" #Both are gradient safe
        af = OWENSAero.readaerodyn_BV("$(path)/airfoils/NACA_0015_RE3E5.dat")
    elseif DynamicStallModel=="none"
        af = OWENSAero.readaerodyn("$(path)/airfoils/NACA_0015_RE3E5.dat")
    else
        af = affun
    end

    rho = 1.225
    mu = 1.7894e-5
    tsrvec = [2.1,3.1,4.2,x[2],6.0,7.0]
    V_xtemp = omega/tsrvec[4]*R
    V_ytemp = zeros(Real,size(V_xtemp))
    V_z = zeros(Real,size(V_xtemp))
    V_twist = zeros(Real,size(V_xtemp))
    windangle = windangle_D * pi/180
    V_x = V_xtemp*cos(windangle)-V_ytemp*sin(windangle)
    V_y = V_xtemp*sin(windangle)+V_ytemp*cos(windangle)

    z = 1.0
    turbine = OWENSAero.Turbine(R,r,z,chord,twist,delta,omega,B,af,ntheta,false)
    env = OWENSAero.Environment(rho,mu,V_x,V_y,V_z,V_twist,windangle,DynamicStallModel,AeroModel,zeros(Real,ntheta*2))

    if steady
        CP, Th, Q, Rp, Tp, Zp, Vloc, CD, CT, a, awstar, alpha, cl, cd_af, thetavec, Re = OWENSAero.steady(turbine, env)
    else
        N_Rev = 2
        tau = [0.3,3.0]

        CP = zeros(Real,ntheta*N_Rev)
        Rp = zeros(Real,ntheta*N_Rev)
        Tp = zeros(Real,ntheta*N_Rev)
        Zp = zeros(Real,ntheta*N_Rev)
        alpha = zeros(Real,ntheta*N_Rev)
        cl = zeros(Real,ntheta*N_Rev)
        cd_af = zeros(Real,ntheta*N_Rev)
        Vloc = zeros(Real,ntheta*N_Rev)
        Re = zeros(Real,ntheta*N_Rev)
        thetavec = zeros(Real,ntheta*N_Rev)
        Q = 0.0
        Th = 0.0
        CD = 0.0
        CT = 0.0
        a = 0.0
        awstar = 0.0
        bnum = 1
        us_param = OWENSAero.UnsteadyParams(RPI,tau,ifw)

        if ifw
            turbsim_filename = "$path/data/ifw/test.bts"
            OWENSOpenFASTWrappers.ifwinit(;inflowlib_filename=nothing,turbsim_filename)
        end

        for step1 = 1:ntheta*N_Rev
            rev_step = Int(step1-floor(Int,(step1-1)/ntheta)*ntheta + (bnum-1)*ntheta/B)
            if rev_step>ntheta
                rev_step -= ntheta
            end

            CP_temp, Th_temp, Q_temp, Rp_temp, Tp_temp, Zp_temp, Vloc_temp,
            CD_temp, CT_temp, a_temp, awstar_temp, alpha_temp, cl_temp, cd_temp,
            thetavec_temp, Re_temp = OWENSAero.Unsteady_Step(turbine,env,us_param,step1)

            Rp[step1] = Rp_temp[rev_step]
            Tp[step1] = Tp_temp[rev_step]
            Zp[step1] = Zp_temp[rev_step]
            alpha[step1] = alpha_temp[rev_step]
            if returnall
                cl[step1] = cl_temp[rev_step]
                cd_af[step1] = cd_temp[rev_step]
                Vloc[step1] = Vloc_temp[rev_step]
                Re[step1] = Re_temp[rev_step]
            end
            thetavec[step1] =
                thetavec_temp[rev_step] + 2 * pi * floor(Int, (step1 - 1) / ntheta)

        end
        if ifw
            OWENSOpenFASTWrappers.ifwend()
        end
    end

    if returnall
        return CP, Th, Q, Rp, Tp, Zp, Vloc, CD, CT, a, awstar, alpha, cl, cd_af, thetavec, Re
    else
        return [CP;Tp;Rp;Zp]
    end

end

ntheta = 30
Rp_us = zeros(ntheta) #initialize scope
Tp_us = zeros(ntheta) #initialize scope
for AeroModel in ["AC","DMS"]
    # AeroModel = "AC"
    println(AeroModel)
    RPM = 150.0
    TSR = 5.2
    xin = [RPM,TSR]
    filterwindow = 1*ntheta
    ##########################################
    @testset "TEST Unsteady Method with ifw $AeroModel" begin
        ##########################################

        CP, Th, Q, Rp, Tp, Zp, Vloc, CD, CT, a, awstar, alpha, cl_af, cd_af, thetavec, Re = aerowrapper(xin;returnall=true,AeroModel,steady=false,ifw=true)
        test_unsteady_outputs(AeroModel, true, Rp, Tp, Zp, Vloc, Re, IFW_UNSTEADY_BASELINES[AeroModel])
        # PyPlot.figure()
        # PyPlot.plot(thetavec./ntheta,Tp,label="Turbulent",color=plot_cycle[1])
        # Tp_rollingave = RollingFunctions.runmean(Tp,filterwindow)
        # PyPlot.plot(thetavec./ntheta,Tp_rollingave,"--",label="Turbulent Averaged",color=plot_cycle[1])
    end
    ################################
    @testset "TEST Unsteady Method $AeroModel" begin
        #################################
        # PyPlot.figure()
        i=1
        for RPI in [true,false]
            i+=1
            global Rp_us
            global Tp_us
            CP, Th, Q, Rp_us, Tp_us, Zp, Vloc, CD, CT, a, awstar, alpha, cl_af, cd_af, thetavec, Re = aerowrapper(xin;RPI,returnall=true,AeroModel,steady=false)
            test_unsteady_outputs(AeroModel, RPI, Rp_us, Tp_us, Zp, Vloc, Re, UNSTEADY_BASELINES[(AeroModel, RPI)])
            @test thetavec[1] ≈ pi / ntheta atol=1e-14
            @test thetavec[ntheta] ≈ 2 * pi - pi / ntheta atol=1e-14
            @test thetavec[ntheta+1] ≈ 2 * pi + pi / ntheta atol=1e-14
            @test thetavec[end] ≈ 4 * pi - pi / ntheta atol=1e-14

            # PyPlot.plot(thetavec./ntheta,Tp_us,label="Steady Wind RPI: $RPI",color=plot_cycle[i])
            # Tp_us_rollingave = RollingFunctions.runmean(Tp_us,filterwindow)
            # PyPlot.plot(thetavec./ntheta,Tp_us_rollingave,"--",label="Steady Wind Averaged RPI: $RPI",color=plot_cycle[i])
        end
        # PyPlot.xlabel("Revolution")
        # PyPlot.ylabel("Tangential Force per Height")
        # PyPlot.legend()
    end
    ###############################


    ################################
    @testset "TEST Steady GRADIENTS $AeroModel" begin
        ################################

        if AeroModel=="DMS"
            aerowrapper0(x) = aerowrapper(x;AeroModel)
            J = ForwardDiff.jacobian(aerowrapper0, xin)

            function aerowrapper2(x)
                return Float64.(aerowrapper(x))
            end

            J2 = FiniteDiff.finite_difference_jacobian(aerowrapper2,xin)

            for i1 = 1:length(J[:,1])
                for i2 = 1:length(J[1,:])
                    @test isapprox(J[i1,i2],J2[i1,i2];atol=1e-5)
                end
            end
        else
            @test_skip "AC analytic gradients need to make it past NLsolve to work"
        end
    end
    ################################
    @testset "TEST Wind Direction $AeroModel" begin
        #################################
        if AeroModel=="DMS"
            windangle_D = collect(0.0:30.0:90.0)
            CP_ref, Rp_ref, Tp_ref, Vloc_ref, alpha_ref, thetavec_ref = nothing, nothing, nothing, nothing, nothing, nothing
            # PyPlot.figure()
            for iwind = 1:length(windangle_D)

                local CP, Th, Q, Rp, Tp, Zp, Vloc, CD, CT, a, awstar, alpha, cl_af, cd_af, thetavec, Re = aerowrapper(xin;AeroModel,returnall=true,windangle_D=windangle_D[iwind])
                @test isfinite(CP)
                @test all(isfinite, Float64.(Rp))
                @test all(isfinite, Float64.(Tp))
                @test all(isfinite, Float64.(Vloc))
                @test all(isfinite, Float64.(alpha))
                if iwind == 1
                    CP_ref = CP
                    Rp_ref = copy(Rp)
                    Tp_ref = copy(Tp)
                    Vloc_ref = copy(Vloc)
                    alpha_ref = copy(alpha)
                    thetavec_ref = copy(thetavec)
                else
                    @test CP ≈ CP_ref atol=1e-14
                    @test Rp ≈ Rp_ref atol=1e-12
                    @test Tp ≈ Tp_ref atol=1e-12
                    @test Vloc ≈ Vloc_ref atol=1e-12
                    @test alpha ≈ alpha_ref atol=1e-14
                    @test thetavec ≈ thetavec_ref .- windangle_D[iwind] * pi / 180 atol=1e-14
                end
                rel_ang_off = 0.0#windangle_D[iwind]*pi/180 #turn on to align loads
                # PyPlot.plot(thetavec[1:ntheta].+rel_ang_off,Tp,label="Angle $(windangle_D[iwind])")
            end
            # PyPlot.legend()
        else
            for windangle_D in (0.0, 30.0, 60.0)
                local CP, Th, Q, Rp, Tp, Zp, Vloc, CD, CT, a, awstar, alpha, cl_af, cd_af, thetavec, Re = aerowrapper(xin;AeroModel,returnall=true,windangle_D)
                @test isfinite(CP)
                @test CP > 0
                @test all(isfinite, Float64.(Rp))
                @test all(isfinite, Float64.(Tp))
                @test all(isfinite, Float64.(Vloc))
                @test minimum(Float64.(Vloc)) > 0
                @test thetavec[1] ≈ pi / ntheta atol=1e-14
                @test thetavec[end] ≈ 2 * pi - pi / ntheta atol=1e-14
            end
            @test_skip "AC wind-angle frame invariance is not implemented yet"
        end
    end
    ################################
    @testset "TEST Accuracy $AeroModel" begin
        #################################

        CP, Th, Q, Rp, Tp, Zp, Vloc, CD, CT, a, awstar, alpha, cl_af, cd_af, thetavec, Re = aerowrapper(xin;returnall=true,AeroModel)

        filename = "$path/data/$(AeroModel)_simple_unit_data.h5"
        # HDF5.h5open(filename, "w") do file
        #     HDF5.write(file,"Rp_us_old",Float64.(Rp_us))
        #     HDF5.write(file,"Tp_us_old",Float64.(Tp_us))
        #     HDF5.write(file,"CP_old",Float64.(CP))
        #     HDF5.write(file,"Rp_old",Float64.(Rp))
        #     HDF5.write(file,"Tp_old",Float64.(Tp))
        #     HDF5.write(file,"Zp_old",Float64.(Zp))
        #     HDF5.write(file,"Vloc_old",Float64.(Vloc))
        #     HDF5.write(file,"CT_old",Float64.(CT))
        #     HDF5.write(file,"a_old",Float64.(a))
        #     HDF5.write(file,"alpha_old",Float64.(alpha))
        #     HDF5.write(file,"cl_af_old",Float64.(cl_af))
        #     HDF5.write(file,"cd_af_old",Float64.(cd_af))
        #     HDF5.write(file,"thetavec_old",Float64.(thetavec))
        # end

        CP_old = HDF5.h5read(filename,"CP_old")
        Rp_old = HDF5.h5read(filename,"Rp_old")
        Tp_old = HDF5.h5read(filename,"Tp_old")
        Rp_us_old = HDF5.h5read(filename,"Rp_us_old")
        Tp_us_old = HDF5.h5read(filename,"Tp_us_old")
        Zp_old = HDF5.h5read(filename,"Zp_old")
        CT_old = HDF5.h5read(filename,"CT_old")
        Vloc_old = HDF5.h5read(filename,"Vloc_old")
        a_old = HDF5.h5read(filename,"a_old")
        alpha_old = HDF5.h5read(filename,"alpha_old")
        cl_af_old = HDF5.h5read(filename,"cl_af_old")
        cd_af_old = HDF5.h5read(filename,"cd_af_old")
        thetavec_old = HDF5.h5read(filename,"thetavec_old")

        fixture_rtol = AeroModel == "AC" ? 1e-3 : 1e-10
        @test size(Rp) == size(Rp_old)
        @test size(Tp) == size(Tp_old)
        @test size(Zp) == size(Zp_old)
        @test size(Vloc) == size(Vloc_old)
        @test size(alpha) == size(alpha_old)
        @test size(cl_af) == size(cl_af_old)
        @test size(cd_af) == size(cd_af_old)
        @test size(thetavec) == size(thetavec_old)
        @test isfinite(CP)
        @test CP ≈ CP_old rtol=fixture_rtol atol=1e-12
        @test Rp[1] ≈ Rp_old[1] rtol=fixture_rtol atol=1e-12
        @test Tp[1] ≈ Tp_old[1] rtol=fixture_rtol atol=1e-12
        @test Vloc[1] ≈ Vloc_old[1] rtol=fixture_rtol atol=1e-12
        @test thetavec[1] ≈ thetavec_old[1] atol=1e-14
        @test thetavec[end] ≈ thetavec_old[end] atol=1e-14
        @test all(isfinite, Float64.(Rp))
        @test all(isfinite, Float64.(Tp))
        @test all(isfinite, Float64.(Vloc))
        @test all(isfinite, Float64.(alpha))
        @test all(isfinite, Float64.(Re))
        @test minimum(Float64.(Vloc)) > 0
        @test minimum(Float64.(Re)) > 0

        tol = 1e-2

        atol = maximum(abs.(Rp_us_old))*0.02
        for ii = 1:length(Rp_us)
            @test isapprox(Rp_us[ii],Rp_us_old[ii];atol)
        end

        atol = maximum(abs.(Tp_us_old))*0.04
        for ii = 1:length(Tp_us)
            @test isapprox(Tp_us[ii],Tp_us_old[ii];atol)
        end
        for ii = 1:length(CP)
            @test isapprox(CP[ii],CP_old[ii][1];atol = tol)
        end
        # PyPlot.figure()
        # PyPlot.plot(LinRange(0,1,length(Rp)),Rp,"r-")
        # PyPlot.plot(LinRange(0,1,length(Rp)),Rp_old,"k-")
        atol = maximum(abs.(Rp_old))*0.02
        for ii = 1:length(Rp)
            @test isapprox(Rp[ii],Rp_old[ii];atol)
        end
        # PyPlot.figure()
        # PyPlot.plot(LinRange(0,1,length(Tp)),Tp,"r-")
        # PyPlot.plot(LinRange(0,1,length(Tp)),Tp_old,"k-")
        atol = maximum(abs.(Tp_old))*0.02
        for ii = 1:length(Tp)
            @test isapprox(Tp[ii],Tp_old[ii];atol)
        end
        atol = maximum(abs.(Zp_old))*0.02
        for ii = 1:length(Zp)
            @test isapprox(Zp[ii],Zp_old[ii];atol)
        end
        atol = maximum(abs.(CT_old))*0.02
        for ii = 1:length(CT)
            @test isapprox(CT[ii],CT_old[ii];atol)
        end
        atol = maximum(abs.(Vloc_old))*0.02
        for ii = 1:length(Vloc)
            @test isapprox(Vloc[ii],Vloc_old[ii];atol)
        end
        atol = maximum(abs.(a_old))*0.02
        for ii = 1:length(a)
            @test isapprox(a[ii],a_old[ii];atol)
        end
        atol = maximum(abs.(alpha_old))*0.02
        for ii = 1:length(alpha)
            @test isapprox(alpha[ii],alpha_old[ii];atol)
        end
        atol = maximum(abs.(cl_af_old))*0.02
        for ii = 1:length(cl_af)
            @test isapprox(cl_af[ii],cl_af_old[ii];atol)
        end
        atol = maximum(abs.(cd_af_old))*0.02
        for ii = 1:length(cd_af)
            @test isapprox(cd_af[ii],cd_af_old[ii];atol)
        end
        atol = maximum(abs.(thetavec_old))*0.02
        for ii = 1:length(thetavec)
            @test isapprox(thetavec[ii],thetavec_old[ii];atol)
        end
    end
end
