
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


A_model = "DMS"
DS_model = "BV"
r_delta_infl = false

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
Nslices = 30#length(shapeX)-1;
shapeY = collect(LinRange(0,H,Nslices+1))
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

twist3D = -atan.(aerocenter_dist./r3D)#ones(Nslices)*-0.4*pi/180

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
CT = zeros(Nslices,n_tsr)
CP = zeros(Nslices,n_tsr)
Qave = zeros(Nslices,n_tsr)
Rp = zeros(ntheta,Nslices,n_tsr)
Tp = zeros(ntheta,Nslices,n_tsr)
Rp_out = zeros(Nslices,n_tsr)
Tp_out = zeros(Nslices,n_tsr)
Zp_out = zeros(Nslices,n_tsr)
Zp = zeros(ntheta,Nslices,n_tsr)
alpha = zeros(ntheta,Nslices,n_tsr)
cl = zeros(ntheta,Nslices,n_tsr)
cd = zeros(ntheta,Nslices,n_tsr)
W = zeros(ntheta,Nslices,n_tsr)
Re = zeros(ntheta,Nslices,n_tsr)
w = zeros(ntheta*2,Nslices,n_tsr)
thetavec = zeros(ntheta)

println("Aero_model: $(A_model), DS_model: $(DS_model), r_delta_infl: $(r_delta_infl)")

for tsr_idx = 1:n_tsr

    tsr_digit1 = 5#Int(floor(tsrvec[tsr_idx]))
    tsr_digit2 = 2#Int(round((tsrvec[tsr_idx]-tsr_digit1)*10.0))
    Vinf = abs.(omega)/tsrvec[tsr_idx]*R

    #azimuthal discretization, data columns, slices
    elapsed = zeros(Nslices)
    for slice = 1:Nslices
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

    end
    mean_elapsed = mean(elapsed)
    println("$(mean_elapsed) seconds")

    # CPvec[tsr_idx] = sum(CP[:,tsr_idx].*h_frac)
    CPvec[tsr_idx] = OWENSAero.trapz(h,CP[:,tsr_idx].*2*R./RefArea) #Undo local normalization and normalize by full turbine

end

experimental = DelimitedFiles.readdlm("$(path)/data/steady/SNL5m_CP_150RPM_0015.txt")






# New Interface
Vinf = 10.0
omega = 150.0 / 60.0 * 2*pi 
OWENSAero.setupTurb(shapeX,shapeY,B,chord,omega,Vinf;rho,mu,eta,afname = "$(path)/../airfoils/NACA_0015.dat",DSModel="BV",Nslices)

TSRrange = LinRange(1.0,8.0,20)
CP = zeros(length(TSRrange))

for (iTSR,TSR) in enumerate(collect(TSRrange))
    # iTSR = 1
    # TSR = 2.2
    #Steady State Test
    # omega = Vinf/radius*TSR
    radius = R
    Vinf = omega*radius/TSR
    
    RPM = omega * 60 / (2*pi)

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
    power2Steady,torque,z3Dnorm,delta,Mz_base2,M_addedmass_Np,M_addedmass_Tp,F_addedmass_Np,F_addedmass_Tp = OWENSAero.steadyTurb(;omega,Vinf)

    # println("Re: $(rho/mu*omega*radius*chord)")
    # println(mean(ReSteady))
    height = H
    area = height * radius*2
    
    windpower = 0.5*rho*Vinf.^3*RefArea

    CP[iTSR] = powerSteady/windpower
end

PyPlot.rc("figure", figsize=(4, 2.5))
PyPlot.pygui(true)
PyPlot.figure()
PyPlot.plot(experimental[1:end-1,1],experimental[1:end-1,2],"k.",label="Exp.")
PyPlot.plot(tsrvec,CPvec,".-",markersize = 6,linewidth = 1.4,color = plot_cycle[1],label="DMS")
PyPlot.plot(TSRrange,CP,".-",markersize = 6,linewidth = 1.4,color = plot_cycle[2],label="DMS_new")
PyPlot.xlabel("Tip Speed Ratio")
PyPlot.ylabel("Coefficient of Performance")
PyPlot.legend()