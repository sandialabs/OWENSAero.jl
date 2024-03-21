# import PyPlot
# PyPlot.close("all")
import OWENSAero
using Test
import HDF5

path,_ = splitdir(@__FILE__)
# include("$(path)/../src/OWENSAero.jl")

global Radius

function aerowrapper(chord;TSR=6.0,ifw=false,bladelength=110.0)
    ntheta = 30#450
    Nslices = 30
    Vinf = 13.0
    println("Running")
    xrotor = [0.0,10,10.0017885627159,11.1553560041493,12.2884687527063,13.4014193201263,14.4944950133051,15.5679780084628,16.6221454239871,17.6572693919714,18.6736171284654,19.6714510024561,20.6510286035987,21.6126028087127,22.5564218470618,23.4827293644340,24.3917644860386,25.2837618782360,26.1589518091164,27.0175602079435,27.8598087234779,28.6859147811952,29.4960916394141,30.2905484443493,31.0694902841018,31.8331182416024,32.5816294465209,33.3152171261553,34.0340706553126,34.7383756051961,35.4283137913100,36.1040633203956,36.7657986364086,37.4136905655517,38.0479063603738,38.6686097429453,39.2759609471231,39.8701167599149,40.4512305619533,41.0194523670911,41.5749288611271,42.1178034396727,42.6482162451696,43.1663042030671,43.6722010571695,44.1660374041621,44.6479407273238,45.1180354294375,45.5764428649045,46.0232813710714,46.4586662987795,46.8827100421419,47.2955220675583,47.6972089419733,48.0878743603867,48.4676191726225,48.8365414093626,49.1947363074538,49.5422963344927,49.8793112126966,50.2058679420645,50.5220508228369,50.8279414772568,51.1236188706413,51.4091593317660,51.6846365725691,51.9501217071802,52.2056832702783,52.4513872347840,52.6872970288901,52.9134735524361,53.1299751926284,53.3368578391141,53.5341748984079,53.7219773076794,53.9003135479022,54.0692296563697,54.2287692385790,54.3789734794877,54.5198811541461,54.6515286377068,54.7739499148146,54.8871765883801,54.9912378877372,55.0861606761897,55.1719694579451,55.2486863844405,55.3163312600615,55.3749215472542,55.4244723710330,55.4649965228855,55.4965044640745,55.5190043283384,55.5325019239909,55.5370007354208,55.5325019239909,55.5190043283384,55.4965044640745,55.4649965228855,55.4244723710330,55.3749215472542,55.3163312600615,55.2486863844405,55.1719694579451,55.0861606761897,54.9912378877372,54.8871765883801,54.7739499148146,54.6515286377068,54.5198811541461,54.3789734794877,54.2287692385790,54.0692296563697,53.9003135479022,53.7219773076794,53.5341748984079,53.3368578391141,53.1299751926284,52.9134735524361,52.6872970288901,52.4513872347840,52.2056832702783,51.9501217071802,51.6846365725691,51.4091593317660,51.1236188706413,50.8279414772568,50.5220508228369,50.2058679420645,49.8793112126966,49.5422963344927,49.1947363074538,48.8365414093626,48.4676191726225,48.0878743603867,47.6972089419733,47.2955220675583,46.8827100421419,46.4586662987795,46.0232813710714,45.5764428649045,45.1180354294375,44.6479407273238,44.1660374041621,43.6722010571695,43.1663042030671,42.6482162451696,42.1178034396727,41.5749288611271,41.0194523670911,40.4512305619533,39.8701167599149,39.2759609471231,38.6686097429453,38.0479063603738,37.4136905655517,36.7657986364086,36.1040633203956,35.4283137913100,34.7383756051961,34.0340706553126,33.3152171261553,32.5816294465209,31.8331182416024,31.0694902841018,30.2905484443493,29.4960916394141,28.6859147811952,27.8598087234779,27.0175602079435,26.1589518091164,25.2837618782360,24.3917644860386,23.4827293644340,22.5564218470618,21.6126028087127,20.6510286035987,19.6714510024561,18.6736171284654,17.6572693919714,16.6221454239871,15.5679780084628,14.4944950133051,13.4014193201263,12.2884687527063,11.1553560041493,10.0017885627159,8.82746863631162,7.63209307561193,6.41535329580485,5.17693519693007,5,0.0]

    zrotor = [1.0,4.47914708496519,4.48000000000000,5.04000000000000,5.60000000000000,6.16000000000000,6.72000000000000,7.28000000000000,7.84000000000000,8.40000000000000,8.96000000000000,9.52000000000000,10.0800000000000,10.6400000000000,11.2000000000000,11.7600000000000,12.3200000000000,12.8800000000000,13.4400000000000,14,14.5600000000000,15.1200000000000,15.6800000000000,16.2400000000000,16.8000000000000,17.3600000000000,17.9200000000000,18.4800000000000,19.0400000000000,19.6000000000000,20.1600000000000,20.7200000000000,21.2800000000000,21.8400000000000,22.4000000000000,22.9600000000000,23.5200000000000,24.0800000000000,24.6400000000000,25.2000000000000,25.7600000000000,26.3200000000000,26.8800000000000,27.4400000000000,28,28.5600000000000,29.1200000000000,29.6800000000000,30.2400000000000,30.8000000000000,31.3600000000000,31.9200000000000,32.4800000000000,33.0400000000000,33.6000000000000,34.1600000000000,34.7200000000000,35.2800000000000,35.8400000000000,36.4000000000000,36.9600000000000,37.5200000000000,38.0800000000000,38.6400000000000,39.2000000000000,39.7600000000000,40.3200000000000,40.8800000000000,41.4400000000000,42,42.5600000000000,43.1200000000000,43.6800000000000,44.2400000000000,44.8000000000000,45.3600000000000,45.9200000000000,46.4800000000000,47.0400000000000,47.6000000000000,48.1600000000000,48.7200000000000,49.2800000000000,49.8400000000000,50.4000000000000,50.9600000000000,51.5200000000000,52.0800000000000,52.6400000000000,53.2000000000000,53.7600000000000,54.3200000000000,54.8800000000000,55.4400000000000,56,56.5600000000000,57.1200000000000,57.6800000000000,58.2400000000000,58.8000000000000,59.3600000000000,59.9200000000000,60.4800000000000,61.0400000000000,61.6000000000000,62.1600000000000,62.7200000000000,63.2800000000000,63.8400000000000,64.4000000000000,64.9600000000000,65.5200000000000,66.0800000000000,66.6400000000000,67.2000000000000,67.7600000000000,68.3200000000000,68.8800000000000,69.4400000000000,70,70.5600000000000,71.1200000000000,71.6800000000000,72.2400000000000,72.8000000000000,73.3600000000000,73.9200000000000,74.4800000000000,75.0400000000000,75.6000000000000,76.1600000000000,76.7200000000000,77.2800000000000,77.8400000000000,78.4000000000000,78.9600000000000,79.5200000000000,80.0800000000000,80.6400000000000,81.2000000000000,81.7600000000000,82.3200000000000,82.8800000000000,83.4400000000000,84,84.5600000000000,85.1200000000000,85.6800000000000,86.2400000000000,86.8000000000000,87.3600000000000,87.9200000000000,88.4800000000000,89.0400000000000,89.6000000000000,90.1600000000000,90.7200000000000,91.2800000000000,91.8400000000000,92.4000000000000,92.9600000000000,93.5200000000000,94.0800000000000,94.6400000000000,95.2000000000000,95.7600000000000,96.3200000000000,96.8800000000000,97.4400000000000,98,98.5600000000000,99.1200000000000,99.6800000000000,100.240000000000,100.800000000000,101.360000000000,101.920000000000,102.480000000000,103.040000000000,103.600000000000,104.160000000000,104.720000000000,105.280000000000,105.840000000000,106.400000000000,106.960000000000,107.520000000000,108.080000000000,108.640000000000,109.200000000000,109.760000000000,109.838611903775,109.9]

    shapeX_raw = xrotor./150*bladelength
    shapeY_raw = zrotor./150*bladelength
    B = 3

    OWENSAero.setupTurb(shapeX_raw,shapeY_raw,B,chord,TSR,Vinf;afname = "$(path)/airfoils/NACA_0015.dat",DSModel="BV")

    global Radius = maximum(shapeX_raw)
    omega = Vinf/Radius*TSR
    RPM = omega/2/pi*60
    N_Rev = 15
    mytime = N_Rev/RPM*60
    n_steps = ntheta*N_Rev#120
    CP = zeros(Nslices,n_steps)
    Rp = zeros(B,Nslices,n_steps)
    Tp = zeros(B,Nslices,n_steps)
    Zp = zeros(B,Nslices,n_steps)
    alpha = zeros(B,Nslices,n_steps)
    cl = zeros(Nslices,n_steps)
    cd_af = zeros(Nslices,n_steps)
    Vloc = zeros(Nslices,n_steps)
    Re = zeros(Nslices,n_steps)
    thetavec = zeros(n_steps)

    # Base Loads
    Fx_base = zeros(n_steps)
    Fy_base = zeros(n_steps)
    Fz_base = zeros(n_steps)
    Mx_base = zeros(n_steps)
    My_base = zeros(n_steps)
    Mz_base = zeros(n_steps)
    power = zeros(n_steps)
    power2 = zeros(n_steps)

    for (ii,tnew) in enumerate([mytime/2,mytime])

        CP_temp,
        Rp_temp,
        Tp_temp,
        Zp_temp,
        alpha_temp,
        cl_temp,
        cd_af_temp,
        Vloc_temp,
        Re_temp,
        thetavec_temp,
        nstep_temp,
        Fx_base_temp,
        Fy_base_temp,
        Fz_base_temp,
        Mx_base_temp,
        My_base_temp,
        Mz_base_temp,
        power_temp,
        power2_temp = OWENSAero.advanceTurb(tnew)

        start = Int((ii-1)*n_steps/2+1)
        stop = Int((ii)*n_steps/2)

        CP[:,start:stop] = CP_temp[:,1:stop-start+1]
        Rp[:,:,start:stop] = Rp_temp[:,:,1:stop-start+1]
        Tp[:,:,start:stop] = Tp_temp[:,:,1:stop-start+1]
        Zp[:,:,start:stop] = Zp_temp[:,:,1:stop-start+1]
        alpha[:,:,start:stop] = alpha_temp[:,:,1:stop-start+1]
        cl[:,start:stop] = cl_temp[:,1:stop-start+1]
        cd_af[:,start:stop] = cd_af_temp[:,1:stop-start+1]
        Vloc[:,start:stop] = Vloc_temp[1,:,1:stop-start+1]
        Re[:,start:stop] = Re_temp[:,1:stop-start+1]
        thetavec[start:stop] = thetavec_temp[1,1:stop-start+1]
        Fx_base[start:stop] = Fx_base_temp[1:stop-start+1]
        Fy_base[start:stop] = Fy_base_temp[1:stop-start+1]
        Fz_base[start:stop] = Fz_base_temp[1:stop-start+1]
        Mx_base[start:stop] = Mx_base_temp[1:stop-start+1]
        My_base[start:stop] = My_base_temp[1:stop-start+1]
        Mz_base[start:stop] = Mz_base_temp[1:stop-start+1]
        power[start:stop] = power_temp[1:stop-start+1]
        power2[start:stop] = power2_temp[1:stop-start+1]
    end

    if ifw
        OWENSAero.ifwend()
    end

    Q = 0.0
    Th = 0.0
    CD = 0.0
    CT = 0.0
    a = 0.0
    awstar = 0.0

    return CP, Th, Q, Rp, Tp, Zp, Vloc, CD, CT, a, awstar, alpha, cl, cd_af,
    thetavec, Re, Fx_base,Fy_base,Fz_base,Mx_base,My_base,Mz_base,power,
    power2,omega[1]

end



# for AModel in ["DMS"]#,"AC"]
AModel = "DMS"
# for bladelength in [330.0]#,220.0,150.0]
bladelength = 330.0
TSR = 6.0

if bladelength==150.0
    chord = 2.0
elseif bladelength==220.0
    chord = 3.5
elseif bladelength==330.0
    chord = 5.0
end

ntheta = 30#450
filterwindow = 3*ntheta

##########################################
######## Unsteady Method with ifw ########
##########################################
ifw=false
# Juno.@enter aerowrapper(chord;TSR,ifw,bladelength)
CP, Th, Q, Rp, Tp, Zp, Vloc, CD, CT, a, awstar, alpha, cl_af, cd_af, thetavec, Re,
Fx_base,Fy_base,Fz_base,Mx_base,My_base,Mz_base,power,power2,omega = aerowrapper(chord;TSR,
ifw,bladelength)

RPS = omega/(2*pi)
timetemp = thetavec[:]/ntheta/RPS #1:round(Int,ntheta*N_Rev)
idx_start = 1#round(Int,length(Fx_base)/5)

mytime = timetemp[idx_start:end] .- timetemp[idx_start]

# filename = "$path/data/unsteadyFullTurb_ORIGINAL.h5"
#
# HDF5.h5open(filename, "w") do file
#     HDF5.write(file,"CP",Float64.(CP[:,idx_start:end]))
#     HDF5.write(file,"Rp",Float64.(Rp[:,:,idx_start:end]))
#     HDF5.write(file,"Tp",Float64.(Tp[:,:,idx_start:end]))
#     HDF5.write(file,"Zp",Float64.(Zp[:,:,idx_start:end]))
#     HDF5.write(file,"Vloc",Float64.(Vloc[:,idx_start:end]))
#     HDF5.write(file,"alpha",Float64.(alpha[:,:,idx_start:end]))
#     HDF5.write(file,"cl_af",Float64.(cl_af[:,idx_start:end]))
#     HDF5.write(file,"cd_af",Float64.(cd_af[:,idx_start:end]))
#     HDF5.write(file,"thetavec",Float64.(thetavec[idx_start:end]))
#     HDF5.write(file,"Re",Float64.(Re[:,idx_start:end]))
#     HDF5.write(file,"Fx_base",Float64.(-Fx_base[idx_start:end]))
#     HDF5.write(file,"Fy_base",Float64.(-Fy_base[idx_start:end]))
#     HDF5.write(file,"Fz_base",Float64.(-Fz_base[idx_start:end]))
#     HDF5.write(file,"Mx_base",Float64.(-Mx_base[idx_start:end]))
#     HDF5.write(file,"My_base",Float64.(-My_base[idx_start:end]))
#     HDF5.write(file,"Mz_base",Float64.(-Mz_base[idx_start:end]))
#     HDF5.write(file,"power",Float64.(power[idx_start:end]))
#     HDF5.write(file,"power2",Float64.(power2[idx_start:end]))
#     HDF5.write(file,"time",Float64.(mytime))
# end

file2 = "$path/data/unsteadyFullTurb_ORIGINAL.h5"
CP2 = HDF5.h5read(file2,"CP")
Rp2 = HDF5.h5read(file2,"Rp")
Tp2 = HDF5.h5read(file2,"Tp")
Zp2 = HDF5.h5read(file2,"Zp")
Vloc2 = HDF5.h5read(file2,"Vloc")
alpha2 = HDF5.h5read(file2,"alpha")
cl_af2 = HDF5.h5read(file2,"cl_af")
cd_af2 = HDF5.h5read(file2,"cd_af")
thetavec2 = HDF5.h5read(file2,"thetavec")
Re2 = HDF5.h5read(file2,"Re")
Fx_base2 = HDF5.h5read(file2,"Fx_base")
Fy_base2 = HDF5.h5read(file2,"Fy_base")
Fz_base2 = HDF5.h5read(file2,"Fz_base")
Mx_base2 = HDF5.h5read(file2,"Mx_base")
My_base2 = HDF5.h5read(file2,"My_base")
Mz_base2 = HDF5.h5read(file2,"Mz_base")
power2 = HDF5.h5read(file2,"power")
power22 = HDF5.h5read(file2,"power2")
mytime2 = HDF5.h5read(file2,"time")

# # Blade Loads
#
# PyPlot.figure()
# PyPlot.plot(1:length(Rp2[1,1,:]),Rp2[1,1,:],label="old")
# PyPlot.plot(1:length(Rp[1,1,:]),Rp[1,1,:],".-",label="new")
# PyPlot.legend()
# PyPlot.ylabel("Rp")
#
# PyPlot.figure()
# PyPlot.plot(1:length(Tp2[1,1,:]),Tp2[1,1,:],label="old")
# PyPlot.plot(1:length(Tp[1,1,:]),Tp[1,1,:],".-",label="new")
# PyPlot.legend()
# PyPlot.ylabel("Tp")
#
# PyPlot.figure()
# PyPlot.plot(1:length(Zp2[1,1,:]),Zp2[1,1,:],label="old")
# PyPlot.plot(1:length(Zp[1,1,:]),Zp[1,1,:],".-",label="new")
# PyPlot.legend()
# PyPlot.ylabel("Zp")
#
# # Base Loads
# PyPlot.figure()
# PyPlot.plot(1:length(Fx_base2),Fx_base2,label="old")
# PyPlot.plot(1:length(Fx_base),-Fx_base,".-",label="new")
# PyPlot.legend()
# PyPlot.ylabel("Fx_base")
#
# PyPlot.figure()
# PyPlot.plot(1:length(Fy_base2),Fy_base2,label="old")
# PyPlot.plot(1:length(Fy_base),-Fy_base,".-",label="new")
# PyPlot.legend()
# PyPlot.ylabel("Fy_base")
#
# PyPlot.figure()
# PyPlot.plot(1:length(Fz_base2),Fz_base2,label="old")
# PyPlot.plot(1:length(Fz_base),-Fz_base,".-",label="new")
# PyPlot.legend()
# PyPlot.ylabel("Fz_base")
#
# PyPlot.figure()
# PyPlot.plot(1:length(Mx_base2),Mx_base2,label="old")
# PyPlot.plot(1:length(Mx_base),-Mx_base,".-",label="new")
# PyPlot.legend()
# PyPlot.ylabel("Mx_base")
#
# PyPlot.figure()
# PyPlot.plot(1:length(My_base2),My_base2,label="old")
# PyPlot.plot(1:length(My_base),-My_base,".-",label="new")
# PyPlot.legend()
# PyPlot.ylabel("My_base")
#
# PyPlot.figure()
# PyPlot.plot(1:length(Mz_base2),Mz_base2,label="old")
# PyPlot.plot(1:length(Mz_base),-Mz_base,".-",label="new")
# PyPlot.legend()
# PyPlot.ylabel("Mz_base")

atol = 1e-6
for ii = 1:length(CP2)
    @test isapprox(float.(CP[ii]),CP2[ii];atol)
end
for ii = 1:length(Rp2)
    @test isapprox(float.(Rp[ii]),Rp2[ii];atol)
end
for ii = 1:length(Tp2)
    @test isapprox(float.(Tp[ii]),Tp2[ii];atol)
end
for ii = 1:length(Zp2)
    @test isapprox(float.(Zp[ii]),Zp2[ii];atol)
end
for ii = 1:length(Vloc2)
    @test isapprox(float.(Vloc[ii]),Vloc2[ii];atol)
end
for ii = 1:length(alpha2)
    @test isapprox(float.(alpha[ii]),alpha2[ii];atol)
end
for ii = 1:length(cl_af2)
    @test isapprox(float.(cl_af[ii]),cl_af2[ii];atol)
end
for ii = 1:length(cd_af2)
    @test isapprox(float.(cd_af[ii]),cd_af2[ii];atol)
end
for ii = 1:length(thetavec2)
    @test isapprox(float.(thetavec[ii]),thetavec2[ii];atol)
end
for ii = 1:length(Re2)
    @test isapprox(float.(Re[ii]),Re2[ii];atol)
end
for ii = 1:length(Fx_base2)
    @test isapprox(float.(-Fx_base[ii]),Fx_base2[ii];atol)
end
for ii = 1:length(Fy_base2)
    @test isapprox(float.(-Fy_base[ii]),Fy_base2[ii];atol)
end
for ii = 1:length(Fz_base2)
    @test isapprox(float.(-Fz_base[ii]),Fz_base2[ii];atol)
end
for ii = 1:length(Mx_base2)
    @test isapprox(float.(-Mx_base[ii]),Mx_base2[ii];atol)
end
for ii = 1:length(My_base2)
    @test isapprox(float.(-My_base[ii]),My_base2[ii];atol)
end
for ii = 1:length(Mz_base2)
    @test isapprox(float.(-Mz_base[ii]),Mz_base2[ii];atol)
end
for ii = 1:length(power2)
    @test isapprox(float.(power[ii]),power2[ii];atol)
end
for ii = 1:length(mytime2)
    @test isapprox(float.(mytime[ii]),mytime2[ii];atol)
end

#Steady State Test
global Radius
Vinf = 13.0
omega = Vinf/Radius*TSR

CPSteady,
RpSteady,
TpSteady,
ZpSteady,
alphaSteady,
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
power2Steady = OWENSAero.steadyTurb(;omega,Vinf)

idx_start = 1
filename = "$path/data/steadyFullTurb.h5"

# HDF5.h5open(filename, "w") do file
#     HDF5.write(file,"CP",Float64.(CPSteady[:,idx_start:end]))
#     HDF5.write(file,"Rp",Float64.(RpSteady[:,:,idx_start:end]))
#     HDF5.write(file,"Tp",Float64.(TpSteady[:,:,idx_start:end]))
#     HDF5.write(file,"Zp",Float64.(ZpSteady[:,:,idx_start:end]))
#     HDF5.write(file,"Vloc",Float64.(VlocSteady[:,idx_start:end]))
#     HDF5.write(file,"alpha",Float64.(alphaSteady[:,:,idx_start:end]))
#     HDF5.write(file,"cl_af",Float64.(cl_afSteady[:,idx_start:end]))
#     HDF5.write(file,"cd_af",Float64.(cd_afSteady[:,idx_start:end]))
#     HDF5.write(file,"thetavec",Float64.(thetavecSteady[idx_start:end]))
#     HDF5.write(file,"Re",Float64.(ReSteady[:,idx_start:end]))
#     HDF5.write(file,"Fx_base",Float64.(Fx_baseSteady[idx_start:end]))
#     HDF5.write(file,"Fy_base",Float64.(Fy_baseSteady[idx_start:end]))
#     HDF5.write(file,"Fz_base",Float64.(Fz_baseSteady[idx_start:end]))
#     HDF5.write(file,"Mx_base",Float64.(Mx_baseSteady[idx_start:end]))
#     HDF5.write(file,"My_base",Float64.(My_baseSteady[idx_start:end]))
#     HDF5.write(file,"Mz_base",Float64.(Mz_baseSteady[idx_start:end]))
#     HDF5.write(file,"power",Float64.(powerSteady))
#     HDF5.write(file,"power2",Float64.(power2Steady))
# end

fileSteadyOld = filename
CPSteadyOld = HDF5.h5read(fileSteadyOld,"CP")
RpSteadyOld = HDF5.h5read(fileSteadyOld,"Rp")
TpSteadyOld = HDF5.h5read(fileSteadyOld,"Tp")
ZpSteadyOld = HDF5.h5read(fileSteadyOld,"Zp")
VlocSteadyOld = HDF5.h5read(fileSteadyOld,"Vloc")
alphaSteadyOld = HDF5.h5read(fileSteadyOld,"alpha")
cl_afSteadyOld = HDF5.h5read(fileSteadyOld,"cl_af")
cd_afSteadyOld = HDF5.h5read(fileSteadyOld,"cd_af")
thetavecSteadyOld = HDF5.h5read(fileSteadyOld,"thetavec")
ReSteadyOld = HDF5.h5read(fileSteadyOld,"Re")
Fx_baseSteadyOld = HDF5.h5read(fileSteadyOld,"Fx_base")
Fy_baseSteadyOld = HDF5.h5read(fileSteadyOld,"Fy_base")
Fz_baseSteadyOld = HDF5.h5read(fileSteadyOld,"Fz_base")
Mx_baseSteadyOld = HDF5.h5read(fileSteadyOld,"Mx_base")
My_baseSteadyOld = HDF5.h5read(fileSteadyOld,"My_base")
Mz_baseSteadyOld = HDF5.h5read(fileSteadyOld,"Mz_base")
powerSteadyOld = HDF5.h5read(fileSteadyOld,"power")
power2SteadyOld = HDF5.h5read(fileSteadyOld,"power2")

# PyPlot.figure()
# PyPlot.plot(1:length(RpSteadyOld[1,15,:]),RpSteadyOld[1,15,:],label="old")
# PyPlot.plot(1:length(RpSteady[1,15,:]),RpSteady[1,15,:],".-",label="new")
# PyPlot.legend()
# PyPlot.ylabel("Rp")
#
# PyPlot.figure()
# PyPlot.plot(1:length(TpSteadyOld[1,15,:]),TpSteadyOld[1,15,:],label="old")
# PyPlot.plot(1:length(TpSteady[1,15,:]),TpSteady[1,15,:],".-",label="new")
# PyPlot.legend()
# PyPlot.ylabel("Tp")
#
# PyPlot.figure()
# PyPlot.plot(1:length(ZpSteadyOld[1,15,:]),ZpSteadyOld[1,15,:],label="old")
# PyPlot.plot(1:length(ZpSteady[1,15,:]),ZpSteady[1,15,:],".-",label="new")
# PyPlot.legend()
# PyPlot.ylabel("Zp")
#
# PyPlot.figure()
# PyPlot.plot(1:length(alphaSteadyOld[1,15,:]),alphaSteadyOld[1,15,:].*180/pi,label="old")
# PyPlot.plot(1:length(alphaSteady[1,15,:]),alphaSteady[1,15,:].*180/pi,".-",label="new")
# PyPlot.legend()
# PyPlot.ylabel("alpha")

atol = 1e-4
for ii = 1:length(CPSteadyOld)
    @test isapprox(float.(CPSteady[ii]),CPSteadyOld[ii];atol=atol*abs(CPSteadyOld[ii]))
end
for ii = 1:length(RpSteadyOld)
    @test isapprox(float.(RpSteady[ii]),RpSteadyOld[ii];atol=atol*abs(RpSteadyOld[ii]))
end
for ii = 1:length(TpSteadyOld)
    @test isapprox(float.(TpSteady[ii]),TpSteadyOld[ii];atol=atol*abs(TpSteadyOld[ii]))
end
for ii = 1:length(ZpSteadyOld)
    @test isapprox(float.(ZpSteady[ii]),ZpSteadyOld[ii];atol=atol*abs(ZpSteadyOld[ii]))
end
for ii = 1:length(VlocSteadyOld)
    @test isapprox(float.(VlocSteady[ii]),VlocSteadyOld[ii];atol=atol*abs(VlocSteadyOld[ii]))
end
for ii = 1:length(alphaSteadyOld)
    @test isapprox(float.(alphaSteady[ii]),alphaSteadyOld[ii];atol=atol*abs(alphaSteadyOld[ii]))
end
for ii = 1:length(cl_afSteadyOld)
    @test isapprox(float.(cl_afSteady[ii]),cl_afSteadyOld[ii];atol=atol*abs(cl_afSteadyOld[ii]))
end
for ii = 1:length(cd_afSteadyOld)
    @test isapprox(float.(cd_afSteady[ii]),cd_afSteadyOld[ii];atol=atol*abs(cd_afSteadyOld[ii]))
end
for ii = 1:length(thetavecSteadyOld)
    @test isapprox(float.(thetavecSteady[ii]),thetavecSteadyOld[ii];atol=atol*abs(thetavecSteadyOld[ii]))
end
for ii = 1:length(ReSteadyOld)
    @test isapprox(float.(ReSteady[ii]),ReSteadyOld[ii];atol=atol*abs(ReSteadyOld[ii]))
end
for ii = 1:length(Fx_baseSteadyOld)
    @test isapprox(float.(Fx_baseSteady[ii]),Fx_baseSteadyOld[ii];atol=atol*abs(Fx_baseSteadyOld[ii]))
end
for ii = 1:length(Fy_baseSteadyOld)
    @test isapprox(float.(Fy_baseSteady[ii]),Fy_baseSteadyOld[ii];atol=atol*abs(Fy_baseSteadyOld[ii]))
end
for ii = 1:length(Fz_baseSteadyOld)
    @test isapprox(float.(Fz_baseSteady[ii]),Fz_baseSteadyOld[ii];atol=atol*abs(Fz_baseSteadyOld[ii]))
end
for ii = 1:length(Mx_baseSteadyOld)
    @test isapprox(float.(Mx_baseSteady[ii]),Mx_baseSteadyOld[ii];atol=atol*abs(Mx_baseSteadyOld[ii]))
end
for ii = 1:length(My_baseSteadyOld)
    @test isapprox(float.(My_baseSteady[ii]),My_baseSteadyOld[ii];atol=atol*abs(My_baseSteadyOld[ii]))
end
for ii = 1:length(Mz_baseSteadyOld)
    @test isapprox(float.(Mz_baseSteady[ii]),Mz_baseSteadyOld[ii];atol=atol*abs(Mz_baseSteadyOld[ii]))
end
for ii = 1:length(powerSteadyOld)
    @test isapprox(float.(powerSteady[ii]),powerSteadyOld[ii];atol=atol*abs(powerSteadyOld[ii]))
end

# PyPlot.figure()
# PyPlot.plot(1:length(Rp[1,15,:]),Rp[1,15,:])
# PyPlot.plot(length(Rp[1,15,:])-length(RpSteady[1,15,:])+1:length(Rp[1,15,:]),RpSteady[1,15,:])

# end
# end
