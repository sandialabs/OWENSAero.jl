import PyPlot
PyPlot.rc("figure", figsize=(4, 3))
PyPlot.rc("font", size=10.0)
PyPlot.rc("lines", linewidth=1.5)
PyPlot.rc("lines", markersize=3.0)
PyPlot.rc("legend", frameon=false)
PyPlot.rc("axes.spines", right=false, top=false)
PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
PyPlot.rc("figure",max_open_warning=500)
color_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]
PyPlot.pygui(true)
# PyPlot.close("all")
import OWENSAero
using Test
import HDF5
import FLOWMath
using Statistics:mean
import DelimitedFiles

path,_ = splitdir(@__FILE__)

radius = 0.538
height = 0.807
chordmid = 0.066667
chordtip = 0.04
B = 3

RE_d = [0.9,1.3]
# for (ire,Vinf) in enumerate([0.8,1.2])
    ire = 2
    Vinf = 1.2
    rho=1000.0
    mu=1.792E-3
    eta=0.5

    ntheta = 30
    Nslices = 30
    ifw=false
    AeroModel = "AC"

    Aero_AddedMass_Active = true
    Aero_RotAccel_Active = true

    chord = FLOWMath.akima([0.0,0.5,1.0],[chordtip,chordmid,chordtip],LinRange(0,1,Nslices))

    # PyPlot.figure()
    # PyPlot.plot(collect(LinRange(0,1,Nslices)),chord)

    shapeX = ones(Nslices).*radius
    shapeZ = LinRange(0,height,Nslices)


    OWENSAero.setupTurb(shapeX,shapeZ,B,chord,10.0,Vinf;rho,mu,eta,afname = "$(path)/airfoils/NACA_0021.dat",DynamicStallModel="BV",Nslices,Aero_AddedMass_Active,Aero_RotAccel_Active)

    TSRrange = LinRange(1.0,5.0,15)
    CP = zeros(length(TSRrange))
    RpSteady = zeros(B,Nslices,ntheta,length(TSRrange))
    TpSteady = zeros(B,Nslices,ntheta,length(TSRrange))
    alphaSteady = zeros(B,Nslices,ntheta,length(TSRrange))
    for (iTSR,TSR) in enumerate(collect(TSRrange))
        # iTSR = 1
        # TSR = 2.2
        #Steady State Test
        omega = Vinf/radius*TSR
        
        RPM = omega * 60 / (2*pi)

        CPSteady,
        RpSteady[:,:,:,iTSR],
        TpSteady[:,:,:,iTSR],
        ZpSteady,
        alphaSteady[:,:,:,iTSR],
        cl_afSteady,
        cd_afSteady,
        VlocSteady,
        ReSteady,
        thetavecSteady,
        nstepSteady,
        Fx_baseSteady,
        Fy_baseSteady,
        Fz_baseSteady,
        Mx_baseSteady,
        My_baseSteady,
        Mz_baseSteady,
        powerSteady,
        power2Steady,torque,z3Dnorm,delta,Mz_base2,M_addedmass_Np,
        M_addedmass_Tp,F_addedmass_Np,F_addedmass_Tp = OWENSAero.steadyTurb(;omega,Vinf)

        println("RPM: $RPM")
        # println(mean(ReSteady))
        area = height * radius*2

        CP[iTSR] = powerSteady/(0.5*rho*Vinf^3*area)
    end

    PyPlot.figure("CP")
    PyPlot.plot(TSRrange,CP,"-",color=color_cycle[2],label="OWENS Aero, $(RE_d[ire]) RE_d (No Added Mass)") #,color=color_cycle[2]
# end
RM2_0_538D_RE_D_0_9E6 = DelimitedFiles.readdlm("$(path)/RM2_0.538D_RE_D_0.9E6.csv", ',',Float64)
RM2_0_538D_RE_D_1_3E6 = DelimitedFiles.readdlm("$(path)/RM2_0.538D_RE_D_1.3E6.csv", ',',Float64)
# PyPlot.plot(RM2_0_538D_RE_D_0_9E6[:,1],RM2_0_538D_RE_D_0_9E6[:,2],"k--",label="Exp. 0.9e6 RE_d")
PyPlot.plot(RM2_0_538D_RE_D_1_3E6[:,1],RM2_0_538D_RE_D_1_3E6[:,2],"k-",label="Exp. 1.2e6 RE_d")
# PyPlot.plot(right_TSR,right_CP,"ko",label="Right Only Exp.")
PyPlot.legend()
PyPlot.xlabel("TSR")
PyPlot.ylabel("Cp")

PyPlot.figure("Rp")
PyPlot.plot((1:length(RpSteady[1,15,:,7]))./length(RpSteady[1,15,:,7]).*360,RpSteady[1,15,:,7],".-",label="No Added Mass")
PyPlot.legend()
PyPlot.ylabel("Rp")
PyPlot.xlabel("Azimuth")

PyPlot.figure("Tp")
PyPlot.plot((1:length(TpSteady[1,15,:,7]))./length(TpSteady[1,15,:,7]).*360,TpSteady[1,15,:,7],".-",label="No Added Mass")
PyPlot.legend()
PyPlot.ylabel("Tp")
PyPlot.xlabel("Azimuth")

PyPlot.figure("Alpha")
PyPlot.plot((1:length(alphaSteady[1,15,:,7]))./length(alphaSteady[1,15,:,7]).*360,alphaSteady[1,15,:,7]./2.0./pi.*360,".-",label="Added Mass")
PyPlot.legend()
PyPlot.ylabel("Alpha")
PyPlot.xlabel("Azimuth")


