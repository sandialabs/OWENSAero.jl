import Statistics:mean
import OWENSAero
import FLOWMath
import ForwardDiff
using Test
import HDF5
using SNOW

import PyPlot
PyPlot.pygui(true)
PyPlot.rc("figure", figsize=(4, 3))
PyPlot.rc("font", size=10.0)
PyPlot.rc("lines", linewidth=1.5)
PyPlot.rc("lines", markersize=3.0)
PyPlot.rc("legend", frameon=false)
PyPlot.rc("axes.spines", right=false, top=false)
PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
# rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

PyPlot.close("all")

path,_ = splitdir(@__FILE__)

global iter = 1

function objcon!(con,design_vars;runplot=false)
    global iter += 1
    AModel = "DMS"
    plots = false
    ntheta = 30
    R = 5.0/2#1.73#slice 7#5.0/2 #m
    # Oval 
    # Pairs where the mean radius of the ellipse is 1
    dtheta = 2*pi/ntheta
    theta = collect(dtheta/2:dtheta:2*pi)
    thetad = theta*180/pi
    a = R.*1.0 #0.5965#[1.0,0.99,0.955,0.913,0.845,0.793,0.7135,0.5965]
    b = R.*1.0 #2.0#[1.0,1.01,0.05,1.1,1.2,1.3,1.5,2.0]
    shift = 11 # amount that the ellipse is rotated
    height = 1.5 # NOTE: this model is just for a single slice that is unitized, so we scale up the height based on that by multiplying the single slice's value
    chord = design_vars[1]#0.1524 #m
    RPM = 150.0

    delta = ones(Real,ntheta)*0.0#23.0*pi/180 #slice 7 #rad
    omega = ones(Real,ntheta) * RPM / 60.0 * 2*pi # RPM -> radians/sec
    B = 3
    k = 1.0
    af = OWENSAero.readaerodyn("$(path)/airfoils/NACA_0015_RE3E5.dat") #You'll want to use better airfoil data

    rho = 1.225
    mu = 1.7894e-5
    N_tsr = 5
    tsrvec = LinRange(0.1,4.0,N_tsr)#[2.1,3.1,4.2,5.2,6.0,7.0]

    CPvec = zeros(Real,N_tsr)

    for itsr = 1:length(tsrvec)
        Vinf = (omega)/tsrvec[itsr]*R

        # Apply the twist control points for the given TSR
        twist_control_points = design_vars[(itsr-1)+1:itsr*10+1] .* pi/180
        dtheta = 2*pi/length(twist_control_points)
        theta_control_points = collect(dtheta/2:dtheta:2*pi)
        twist = FLOWMath.akima(theta_control_points,twist_control_points,theta)#ones(Real,ntheta)*0.0 #rad


        # Radius as a function of azimuth angle and the input a and b widths
        r1=a*b ./ ((b*cos.(theta)).^2+(a*sin.(theta)).^2).^0.5
        # now rotate it discretely
        r = zero(r1)
        circshift!(r,r1,shift)

        x = -r.*sin.(theta)
        y = r.*cos.(theta)
        area = (maximum(y) - minimum(y))*height # this is the frontal area the wind sees

        env = OWENSAero.Environment(rho,mu,Vinf,"NONE",AModel,zeros(ntheta*2))
        turbine = OWENSAero.Turbine(R,r.*R,chord,twist,delta,omega,B,af,ntheta,false)

        # start = time()

        _, Th, Q, Rp, Tp, Zp, Vloc, CD, CT, mean_astar,
        astar, alpha, cl_af, cd_af, thetavec, Re, M_addedmass_Np,
        M_addedmass_Tp, F_addedmass_Np, F_addedmass_Tp= OWENSAero.steady(turbine, env)

        # time_elapsed = time() - start
        
        CPvec[itsr] = mean(Q.*omega ./ (0.5.*rho.*Vinf.^3.0.*area))

        if runplot
            # PyPlot.figure(1)#Rp
            # PyPlot.plot(thetad,Rp,"-",color = plot_cycle[itsr],label="TSR $(tsrvec[itsr])")
            # PyPlot.xlabel("Azimuth Angle (deg)")
            # PyPlot.ylabel("Radial Force Per Unit Height (N)")
            # PyPlot.legend()
            # PyPlot.savefig("$(path)/figs/Rp_$shift.pdf",transparent = true)

            # PyPlot.figure(2)#Tp
            # PyPlot.plot(thetad,-Tp,"-",color = plot_cycle[itsr],label="TSR $(tsrvec[itsr])")
            # PyPlot.xlabel("Azimuth Angle (deg)")
            # PyPlot.ylabel("Tangential Force Per Unit Height (N)")
            # PyPlot.legend()
            # PyPlot.savefig("$(path)/figs/Tp_$shift.pdf",transparent = true)

            # PyPlot.figure(3)#Zp
            # PyPlot.plot(thetad,Zp,"-",color = plot_cycle[itsr],label="TSR $(tsrvec[itsr])")
            # PyPlot.xlabel("Azimuth Angle (deg)")
            # PyPlot.ylabel("Vertical Force Per Unit Height (N)")
            # PyPlot.legend()
            # PyPlot.savefig("$(path)/figs/Zp_$shift.pdf",transparent = true)

            # PyPlot.rc("figure", figsize=(4, 3))
            # PyPlot.rc("figure.subplot", left=.15, bottom=.21, top=0.9, right=.99)
            # PyPlot.figure()
            # x = -r.*sin.(theta)
            # y = r.*cos.(theta)
            # angle = atan((maximum(y)-minimum(y))/(maximum(x)-minimum(x)))*180/pi
            # # println("angle: $angle")
            # PyPlot.plot([x;x[1]],[y;y[1]],color = plot_cycle[1],"-")
            # PyPlot.plot(5.0,0,"w.")
            # PyPlot.axis("equal")
            # PyPlot.xlabel("x")
            # PyPlot.ylabel("y")
            # PyPlot.savefig("$(path)/figs/path_$shift.pdf",transparent = true)
            # # println(b./a)
        end
            
    end

    if runplot
        PyPlot.figure(75)
        PyPlot.plot(tsrvec,CPvec,"-")
        PyPlot.xlabel("TSR")
        PyPlot.ylabel("CP")
        sleep(0.1)
        PyPlot.savefig("$(path)/figs/CP_$shift.pdf",transparent = true)
    end
    con[1] = 0 #- FLOWMath.ksmin(CPvec)  # 0 < power

    objective = -mean(CPvec)

    return objective

end


# Since we defined 5 TSRs above, and 10 control points for the pitch for each TSR, apply accordingly
x0 = zeros(10*5+1)
lx = zeros(10*5+1).-20.0
ux = zeros(10*5+1).+20.0
x0[1] = 0.1524 
lx[1] = 0.01
ux[1] = 2.0

N_contstraints = 1
myIpoptoptions = Dict{String, Any}()
myIpoptoptions["hessian_approximation"] = "limited-memory"
myIpoptoptions["limited_memory_update_type"] = "bfgs"
myIpoptoptions["print_level"] = 5
myIpoptoptions["dual_inf_tol"] = 1e-1
myIpoptoptions["constr_viol_tol"] = 1e-1
myIpoptoptions["compl_inf_tol"] = 1e-1
myIpoptoptions["tol"] = 1e-4
myIpoptoptions["max_cpu_time"] = 200.0
optionsIPOPT = Options(solver=IPOPT(myIpoptoptions),derivatives=ForwardAD())
xopt, fopt, info, out = minimize(objcon!, x0, N_contstraints, lx, ux, -Inf.*ones(N_contstraints), 0.0.*ones(N_contstraints), optionsIPOPT)

println(fopt)
x0 = zeros(10*5+1)
x0[1] = 0.1524 
# Rerun with plotting
original_fopt = objcon!(zeros(N_contstraints),x0;runplot=true)
objcon!(zeros(N_contstraints),xopt;runplot=true)

println("Percent improvement in objective: $(-(original_fopt-fopt)/original_fopt*100)")