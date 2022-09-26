using Test
import HDF5
import ForwardDiff
import FiniteDiff
import Statistics:mean
import OpenFASTWrappers
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

path,_ = splitdir(@__FILE__)
import VAWTAero
# include("$(path)/../src/VAWTAero.jl")

function aerowrapper(x;RPI=true,returnall=false,windangle_D=0.0,AModel="DMS",steady=true,ifw=false,DSModel="BV")

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

    if DSModel=="BV" #Both are gradient safe
        af = VAWTAero.readaerodyn_BV("$(path)/airfoils/NACA_0015_RE3E5.dat")
    elseif DSModel=="none"
        af = VAWTAero.readaerodyn("$(path)/airfoils/NACA_0015_RE3E5.dat")
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
    turbine = VAWTAero.Turbine(R,r,z,chord,twist,delta,omega,B,af,ntheta,false)
    env = VAWTAero.Environment(rho,mu,V_x,V_y,V_z,V_twist,windangle,DSModel,AModel,zeros(Real,ntheta*2))

    if steady
        CP, Th, Q, Rp, Tp, Zp, Vloc, CD, CT, a, awstar, alpha, cl, cd_af, thetavec, Re = VAWTAero.steady(turbine, env)
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
        us_param = VAWTAero.UnsteadyParams(RPI,tau,ifw)

        if ifw
            turbsim_filename = "$path/data/ifw/test.bts"
            OpenFASTWrappers.ifwinit(;turbsim_filename)
        end

        for step1 = 1:ntheta*N_Rev
            rev_step = Int(step1-floor(Int,(step1-1)/ntheta)*ntheta + (bnum-1)*ntheta/B)
            if rev_step>ntheta
                rev_step -= ntheta
            end

            CP_temp, Th_temp, Q_temp, Rp_temp, Tp_temp, Zp_temp, Vloc_temp,
            CD_temp, CT_temp, a_temp, awstar_temp, alpha_temp, cl_temp, cd_temp,
            thetavec_temp, Re_temp = VAWTAero.Unsteady_Step(turbine,env,us_param,step1)

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
            thetavec[step1] = step1 #TODO: fix this

        end
        if ifw
            OpenFASTWrappers.ifwend()
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
for AModel in ["AC","DMS"]
    # AModel = "AC"
    println(AModel)
    RPM = 150.0
    TSR = 5.2
    xin = [RPM,TSR]
    filterwindow = 1*ntheta
    ##########################################
    @testset "TEST Unsteady Method with ifw $AModel" begin
        ##########################################

        CP, Th, Q, Rp, Tp, Zp, Vloc, CD, CT, a, awstar, alpha, cl_af, cd_af, thetavec, Re = aerowrapper(xin;returnall=true,AModel,steady=false,ifw=true)
        @test true #check that it runs
        # PyPlot.figure()
        # PyPlot.plot(thetavec./ntheta,Tp,label="Turbulent",color=plot_cycle[1])
        # Tp_rollingave = RollingFunctions.runmean(Tp,filterwindow)
        # PyPlot.plot(thetavec./ntheta,Tp_rollingave,"--",label="Turbulent Averaged",color=plot_cycle[1])
    end
    ################################
    @testset "TEST Unsteady Method $AModel" begin
        #################################
        # PyPlot.figure()
        i=1
        for RPI in [true,false]
            i+=1
            global Rp_us
            global Tp_us
            CP, Th, Q, Rp_us, Tp_us, Zp, Vloc, CD, CT, a, awstar, alpha, cl_af, cd_af, thetavec, Re = aerowrapper(xin;RPI,returnall=true,AModel,steady=false)
            @test true #check that it runs

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
    @testset "TEST Steady GRADIENTS $AModel" begin
        ################################

        if AModel=="DMS"
            aerowrapper0(x) = aerowrapper(x;AModel)
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
            @warn "AC analytic gradients need to make it past NLsolve to work, still to do"
        end
    end
    ################################
    @testset "TEST Wind Direction $AModel" begin
        #################################
        if AModel=="DMS"
            windangle_D = collect(0.0:30.0:90.0)
            # PyPlot.figure()
            for iwind = 1:length(windangle_D)

                local CP, Th, Q, Rp, Tp, Zp, Vloc, CD, CT, a, awstar, alpha, cl_af, cd_af, thetavec, Re = aerowrapper(xin;AModel,returnall=true,windangle_D=windangle_D[iwind])
                rel_ang_off = 0.0#windangle_D[iwind]*pi/180 #turn on to align loads
                # PyPlot.plot(thetavec[1:ntheta].+rel_ang_off,Tp,label="Angle $(windangle_D[iwind])")
            end
            # PyPlot.legend()
        else
            @warn "AC nonzero wind angle still in work"
        end
    end
    ################################
    @testset "TEST Accuracy $AModel" begin
        #################################

        CP, Th, Q, Rp, Tp, Zp, Vloc, CD, CT, a, awstar, alpha, cl_af, cd_af, thetavec, Re = aerowrapper(xin;returnall=true,AModel)

        filename = "$path/data/$(AModel)_simple_unit_data.h5"
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
