import Statistics:mean
import OWENSAero
using Test
import HDF5

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

close("all")

path,_ = splitdir(@__FILE__)
# include("$(path)/../src/OWENSAero.jl")
# Main Function

function runme()
    AModel = "AC"
    plots = false
    ntheta = 30
    R = 5.0/2#1.73#slice 7#5.0/2 #m
    # r = ones(ntheta)*R #m
    chord = 0.1524 #m
    twist = ones(ntheta)*0.0 #rad
    delta = ones(ntheta)*0.0#23.0*pi/180 #slice 7 #rad
    omega = ones(ntheta) * 150.0 / 60.0 * 2*pi # RPM -> radians/sec
    B = 3
    k = 1.0
    af = OWENSAero.readaerodyn("$(path)/airfoils/NACA_0015_RE3E5.dat")

    rho = 1.225
    mu = 1.7894e-5
    suction = false
    tsrvec = [2.1,3.1,4.2,5.2,6.0,7.0]
    Vinf = (omega)/tsrvec[4]*R

    # Pairs where the mean radius of the ellipse is 1
    a = R.*[1.0,0.913,0.793,0.7135,0.5965]#[1.0,0.99,0.955,0.913,0.845,0.793,0.7135,0.5965]
    b = R.*[1.0,1.1,1.3,1.5,2.0]#[1.0,1.01,0.05,1.1,1.2,1.3,1.5,2.0]
    dtheta = 2*pi/ntheta
    theta = collect(dtheta/2:dtheta:2*pi)
    thetad = theta*180/pi
    Tp = zeros(ntheta,length(a))
    Tp_1 = zeros(ntheta,length(a))
    error = zeros(ntheta,length(a))
    max_error = zeros(length(a))
    mean_error = zeros(length(a))
    mymean = zeros(length(a))
    max1 = zeros(length(a))
    elapsed_reffect = zeros(length(a))
    elapsed_no_reffect = zeros(length(a))

for shift = 11
    PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)

    PyPlot.rc("figure", figsize=(4, 3))
    PyPlot.figure(1)#Rp
    PyPlot.figure(2)#Tp
    PyPlot.figure(3)#Zp


    # WITHOUT Radius Effect in Influence Coefficient Matrix
    for ii = 1:length(a)
        r1=a[ii]*b[ii] ./ ((b[ii]*cos.(theta)).^2+(a[ii]*sin.(theta)).^2).^0.5
        r = zero(r1)
        circshift!(r,r1,shift)
        env = OWENSAero.Environment(rho,mu,Vinf,"NONE",AModel,zeros(ntheta*2))
        turbine = OWENSAero.Turbine(R,r.*R,chord,twist,delta,omega,B,af,ntheta,false)

        start = time()
        _, _, _, Rp, Tp_1[:,ii], Zp = OWENSAero.steady(turbine, env)
        elapsed_no_reffect[ii] = time() - start

        PyPlot.figure(1)#Rp
        PyPlot.plot(thetad,Rp,"-",color = plot_cycle[ii])
        PyPlot.figure(2)#Tp
        PyPlot.plot(thetad,-Tp_1[:,ii],"-",color = plot_cycle[ii])
        PyPlot.figure(3)#Zp
        PyPlot.plot(thetad,Zp,"-",color = plot_cycle[ii])

    end

    # With Radius Effect in Influence Coefficient Matrix

    for ii = 1:length(a)
        r1=a[ii]*b[ii] ./ ((b[ii]*cos.(theta)).^2+(a[ii]*sin.(theta)).^2).^0.5
        r = zero(r1)
        circshift!(r,r1,shift)
        turbine = OWENSAero.Turbine(R,r.*R,chord,twist,delta,omega,B,af,ntheta,true)
        env = OWENSAero.Environment(rho,mu,Vinf,"NONE",AModel,zeros(ntheta*2))

        start = time()

        _, _, _, Rp, Tp[:,ii], Zp = OWENSAero.steady(turbine, env)
        elapsed_reffect[ii] = time() - start

        PyPlot.figure(1)#Rp
        PyPlot.plot(thetad,Rp,"--",color = plot_cycle[ii])
        PyPlot.figure(2)#Tp
        PyPlot.plot(thetad,-Tp[:,ii],"--",color = plot_cycle[ii])
        PyPlot.figure(3)#Zp
        PyPlot.plot(thetad,Zp,"--",color = plot_cycle[ii])

        # println("Time Fraction with/without: $(round(elapsed_reffect[ii]/elapsed_no_reffect[ii],digits=3))")
        max_error[ii] = maximum(abs.((Tp_1[:,ii].-Tp[:,ii])) ./ mean(abs.(Tp_1[:,ii]))) .* 100
        mymean[ii] = mean(abs.(Tp_1[:,ii]))
        max1[ii] = minimum(Tp_1[:,ii])
        # mean_error[ii] = mean(abs.((Tp_1[:,ii].-Tp[:,ii])) ./ mean(abs.(Tp_1[:,ii]))) .* 100
        mean_error[ii] = sqrt(sum((Tp_1[:,ii].-Tp[:,ii]).^2)/length(Tp_1[:,ii]))
        error[:,ii] = Tp_1[:,ii].-Tp[:,ii]
        # PyPlot.figure(4)#Rp
        # PyPlot.plot(thetad,Rp,"--",color = plot_cycle[ii])
        PyPlot.figure(5)#Tp
        PyPlot.plot(thetad,(Tp_1[:,ii].-Tp[:,ii]),"--",color = plot_cycle[ii])
        # PyPlot.figure(6)#Zp
        # PyPlot.plot(thetad,Zp,"--",color = plot_cycle[ii])

    end

    PyPlot.figure(1)#Rp
    PyPlot.legend(["0%","10%","30%","50%","100%"])
    PyPlot.xlabel("Azimuth Angle (deg)")
    PyPlot.ylabel("Radial Force Per Unit Height (N)")
    PyPlot.savefig("$(path)/figs/radius_effect/Rp_$shift.pdf",transparent = true)
    PyPlot.figure(2)#Tp
    PyPlot.legend(["0%","10%","30%","50%","100%"])
    PyPlot.xlabel("Azimuth Angle (deg)")
    PyPlot.ylabel("Tangential Force Per Unit Height (N)")
    PyPlot.savefig("$(path)/figs/radius_effect/Tp_$shift.pdf",transparent = true)
    PyPlot.figure(3)#Zp
    PyPlot.legend(["0%","10%","30%","50%","100%"])
    PyPlot.xlabel("Azimuth Angle (deg)")
    PyPlot.ylabel("Vertical Force Per Unit Height (N)")
    PyPlot.savefig("$(path)/figs/radius_effect/Zp_$shift.pdf",transparent = true)

    # PyPlot.figure(4)#Rp
    # PyPlot.legend(["0%","10%","30%","50%","100%"])
    # PyPlot.xlabel("Azimuth Angle (deg)")
    # PyPlot.ylabel("Difference in Radial Force Per Unit Height (N)")
    # PyPlot.savefig("$(path)/figs/radius_effect/Rp_$shift.pdf",transparent = true)
    PyPlot.figure(5)#Tp
    PyPlot.legend(["0%","10%","30%","50%","100%"])
    PyPlot.xlabel("Azimuth Angle (deg)")
    PyPlot.ylabel("Difference in Tangential Force Per Unit Height (N)")
    PyPlot.savefig("$(path)/figs/radius_effect/Tp_diff_$shift.pdf",transparent = true)
    # PyPlot.figure(6)#Zp
    # PyPlot.legend(["0%","10%","30%","50%","100%"])
    # PyPlot.xlabel("Azimuth Angle (deg)")
    # PyPlot.ylabel("Difference in Vertical Force Per Unit Height (N)")
    # PyPlot.savefig("$(path)/figs/radius_effect/Zp_$shift.pdf",transparent = true)

    PyPlot.rc("figure", figsize=(4, 3.1))
    PyPlot.figure()
    # PyPlot.plot((b/R.-1.0).*100,max_error,"k-")
    PyPlot.plot((b/R.-1.0).*100,mean_error,"k--")#,color = plot_cycle[ii])
    PyPlot.xlabel("Percent Ovality")# (b=R*(1+%/100) & mean(r)==R)")
    PyPlot.ylabel("Tangential Force Difference (%)")
    PyPlot.legend(["Max","Mean"])
    PyPlot.savefig("$(path)/figs/radius_effect/error_$shift.pdf",transparent = true)

    println(shift)

    PyPlot.rc("figure", figsize=(4, 3))
    PyPlot.rc("figure.subplot", left=.15, bottom=.21, top=0.9, right=.99)
    PyPlot.figure()
    for ii = 1:length(a)
        r1=a[ii]*b[ii] ./ ((b[ii]*cos.(theta)).^2+(a[ii]*sin.(theta)).^2).^0.5
        r = zero(r1)
        circshift!(r,r1,shift)
        x = -r.*sin.(theta)
        y = r.*cos.(theta)
        angle = atan((maximum(y)-minimum(y))/(maximum(x)-minimum(x)))*180/pi
        println("angle: $angle")
        PyPlot.plot([x;x[1]],[y;y[1]],color = plot_cycle[ii],"-")
    end
    PyPlot.plot(5.0,0,"w.")
    PyPlot.axis("equal")
    PyPlot.legend(["0%","10%","30%","50%","100%"])
    PyPlot.xlabel("x")
    PyPlot.ylabel("y")
    PyPlot.savefig("$(path)/figs/radius_effect/path_$shift.pdf",transparent = true)
    # println(b./a)


    println(maximum(max_error))
    println(maximum(mean_error))
    # PyPlot.close("all")
    return mymean,error,max1
end
end
# Juno.@enter runme()
 mymean,myerror,max1 = runme()
