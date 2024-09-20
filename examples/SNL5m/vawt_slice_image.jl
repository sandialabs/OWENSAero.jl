# using PyCall
# @pyimport matplotlib.patches as patches
# pygui(:qt5)

using PyPlot
# PyPlot.ioff()
using Statistics
path,_ = splitdir(@__FILE__)

function getxy(theta,r)
    circx = r*sin.(theta)
    circy = r*cos.(theta)

    return circx, circy
end

function af_placement(theta,scale,radius)
    #theta is in radians

    NACA0015x = [1.0000,0.9500,0.9000,0.8000,0.7000,0.6000,0.5000,0.4000,0.3000,0.2500,0.2000,0.1500,0.1000,0.0750,0.0500,0.0250,0.0125,0.0000,0.0125,0.0250,0.0500,0.0750,0.1000,0.1500,0.2000,0.2500,0.3000,0.4000,0.5000,0.6000,0.7000,0.8000,0.9000,0.9500,1.0000]
    NACA0015y = [0.00158,0.01008,0.01810,0.03279,0.04580,0.05704,0.06617,0.07254,0.07502,0.07427,0.07172,0.06682,0.05853,0.05250,0.04443,0.03268,0.02367,0.00000,-0.02367,-0.03268,-0.04443,-0.05250,-0.05853,-0.06682,-0.07172,-0.07427,-0.07502,-0.07254,-0.06617,-0.05704,-0.04580,-0.03279,-0.01810,-0.01008,-0.00158]


    #Scale the airfoil
    NACA0015x = NACA0015x*scale
    NACA0015y = NACA0015y*scale

    # Center the airfoil at top dead center
    NACA0015x = NACA0015x .- maximum(NACA0015x)/2.55
    NACA0015y = NACA0015y .+ radius

    #Rotate the airfoil
    afx = NACA0015x*cos(theta) - NACA0015y*sin(theta)
    afy = NACA0015x*sin(theta) + NACA0015y*cos(theta)

    return afx,afy
end

rc("figure", figsize=(4, 4))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=.1, bottom=.1, top=0.9, right=.9)
# rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]
close("all")


u =[0.0752,-0.0176,-0.1044,-0.1787,-0.2362,-0.2725,-0.2828,-0.2793,-0.2614,-0.2325,-0.1940,-0.1466,-0.0869,-0.0240,0.0406,-0.0224,-0.1922,-0.3354,-0.4555,-0.5440,-0.6111,-0.6588,-0.6869,-0.6918,-0.6735,-0.6139,-0.5166,-0.3838,-0.2140,-0.0037]
v =[0.2657,0.2826,0.2690,0.2323,0.1782,0.1128,0.0455,-0.0178,-0.0764,-0.1283,-0.1729,-0.2095,-0.2338,-0.2426,-0.2353,-0.2080,-0.1650,-0.1191,-0.0792,-0.0477,-0.0240,-0.0072,0.0042,0.0127,0.0206,0.0327,0.0558,0.0936,0.1476,0.2136]

#Draw circle
r = 1
Ntheta = 30
arrowscale = .45

#
#
#
# AC Velocities and Blades
#
#
#

PyPlot.figure()

#Circles
circx, circy = getxy(LinRange(0,2*pi,Ntheta*10),r)
PyPlot.plot(circx,circy,"k",linewidth = .75)
circx, circy = getxy(LinRange(0,2*pi-2*pi/Ntheta,Ntheta),r*1)
PyPlot.plot(circx,circy,"k.",markersize = 5)
PyPlot.plot(-2.5,.05,"w.")
circx, circy = getxy(LinRange(0,2*pi-2*pi/Ntheta,Ntheta),r*1)

#Arrows
for i = 1:Ntheta
    PyPlot.arrow(circx[i],circy[i],u[i]*arrowscale,0,overhang = .125, head_width = .02,facecolor = "k")
    PyPlot.arrow(circx[i],circy[i],0,v[i]*arrowscale,overhang = .125, head_width = .02,facecolor = "k")
end

#Airfoils
blade_offset = 35
for thetablade = [0+blade_offset,120+blade_offset,240+blade_offset]
    NACA0015x,NACA0015y = af_placement(thetablade*pi/180,.6,.89)
    PyPlot.plot(NACA0015x,NACA0015y,color = plot_cycle[1],linewidth = .75)
    PyPlot.plot(mean(NACA0015x),mean(NACA0015y),"k.")
    local circx, circy = getxy(thetablade*pi/180,r*1)
    # PyPlot.plot(circx,circy,"r.", markersize = 5)
end

#Axes
PyPlot.plot([-1.5,1.5],[0,0],"k-.",linewidth = .5)
PyPlot.plot([0,0],[-1.5,1.5],"k-.",linewidth = .5)
PyPlot.annotate(L"x,u", xy=(-1.5, 0), xytext=(1.6, -0.0500))
PyPlot.arrow(1.5,0,.003,0.0,overhang = .125, head_width = .04,facecolor = "k")
PyPlot.annotate(L"y,v", xy=(-1.5, 0), xytext=(-0.15, 1.6))
PyPlot.arrow(0,1.45,0.0,0.001,overhang = .125, head_width = .04,facecolor = "k")

#Theta
Ndots = 30
PyPlot.plot(LinRange(0,-1.5*sin(blade_offset*pi/180),Ndots),LinRange(0,1.5*cos(blade_offset*pi/180),Ndots),"k.",markersize = .5)
PyPlot.annotate(L"\theta", xy=(-1.5, 0), xytext=(-.19, .45))
circx, circy = getxy(LinRange(0,-20*pi/180,Ntheta),r*.4)
PyPlot.plot(circx,circy,"k",linewidth = .75)
PyPlot.arrow(circx[end],circy[end],-.01,-.007,overhang = .125, head_width = .04,facecolor = "k")

#Vinf
PyPlot.annotate(L"V_{\infty}", xy=(-1.5, 0), xytext=(-2.5, -0.10))
PyPlot.arrow(-2.23,.0,.55,0,overhang = .125, width = 0.018,head_width = .06,facecolor = "k")
axis("off")
axis("equal")
PyPlot.savefig("$(path)/figs/vawt_slice.pdf",transparent = true)

#
#
#
# Inflow Wind Frame of Reference
#
#
#

PyPlot.figure()

#Circles
circx, circy = getxy(LinRange(0,2*pi,Ntheta*10),r)
PyPlot.plot(circx,circy,"k",linewidth = .75)
PyPlot.plot(-2.5,.05,"w.")
circx, circy = getxy(LinRange(0,2*pi-2*pi/Ntheta,Ntheta),r*1)

#Arrows
# for i = 1:Ntheta
#     PyPlot.arrow(circx[i],circy[i],u[i]*arrowscale,0,overhang = .125, head_width = .02,facecolor = "k")
#     PyPlot.arrow(circx[i],circy[i],0,v[i]*arrowscale,overhang = .125, head_width = .02,facecolor = "k")
# end

#Airfoils
blade_offset = 35
for thetablade = [0+blade_offset,120+blade_offset,240+blade_offset]
    NACA0015x,NACA0015y = af_placement(thetablade*pi/180,.6,.89)
    PyPlot.plot(NACA0015x,NACA0015y,color = plot_cycle[1],linewidth = .75)
end

#Axes
PyPlot.plot([-1.5,1.5],[0,0],"k-.",linewidth = .5)
PyPlot.plot([0,0],[-1.5,1.5],"k-.",linewidth = .5)
PyPlot.annotate(L"x", xy=(-1.5, 0), xytext=(1.6, -0.0500))
PyPlot.arrow(1.5,0,.003,0.0,overhang = .125, head_width = .04,facecolor = "k")
PyPlot.annotate(L"y", xy=(-1.5, 0), xytext=(-0.05, 1.6))
PyPlot.arrow(0,1.5,0.0,0.003,overhang = .125, head_width = .04,facecolor = "k")

#Theta
Ndots = 30
PyPlot.plot(LinRange(0,-2.5*sin(blade_offset*pi/180),Ndots),LinRange(0,0.5*cos(blade_offset*pi/180),Ndots),"k.",markersize = .5)
PyPlot.plot(-LinRange(0,-2.5*sin(blade_offset*pi/180),Ndots),-LinRange(0,0.5*cos(blade_offset*pi/180),Ndots),"k.",markersize = .5)
# PyPlot.annotate(L"\theta", xy=(-1.5, 0), xytext=(-1.35, .1))
circx, circy = getxy(LinRange(-pi/2,-80*pi/180,Ntheta),r*1.1)
# PyPlot.plot(circx,circy,"k",linewidth = .75)
# PyPlot.arrow(circx[end],circy[end],.003,.007,overhang = .125, head_width = .04,facecolor = "k")

# Other side of Theta
PyPlot.plot(-LinRange(0,-2.5*sin(blade_offset*pi/180),Ndots),-LinRange(0,0.5*cos(blade_offset*pi/180),Ndots),"k.",markersize = .5)
PyPlot.annotate(L"\theta", xy=(-1.5, 0), xytext=(1.2, -.25))
PyPlot.plot(-circx,-circy,"k",linewidth = .75)
PyPlot.arrow(-circx[end],-circy[end],-.003,-.007,overhang = .125, head_width = .04,facecolor = "k")


#Vinf
PyPlot.annotate(L"V_{\infty}", xy=(-1.5, 0), xytext=(-2.55, 0.65))
PyPlot.arrow(-2.23,.65,.55,-.17,overhang = .125, width = 0.018,head_width = .06,facecolor = "k")
axis("off")
axis("equal")
PyPlot.savefig("$(path)/figs/inflow_wind.pdf",transparent = true)

#
#
# PLOT VERTICAL Slope
#
#

R = 1.0#5.0/2 #m
H = 1.02*R*2; #m
hr = H/R
n_slices = 30
shapeY = collect(LinRange(0,H,n_slices+1))
shapeX = R.*(1.0.-4.0.*(shapeY/H.-.5).^2)
Y = (shapeY[2:end] + shapeY[1:end-1])/2.0;
X = (shapeX[2:end] + shapeX[1:end-1])/2.0;
delta_xs = shapeX[2:end] - shapeX[1:end-1]
delta_zs = shapeY[2:end] - shapeY[1:end-1]
delta3D = atan.(delta_xs./delta_zs)
slice = 22
bottomslice = 1
lenarrow = 0.4

rc("figure", figsize=(3, 4))
PyPlot.figure()
# Force Vectors
PyPlot.annotate(L"force", xy=(-1.5, 0), xytext=(-0.78, H*0.67))
PyPlot.arrow(-X[slice],Y[slice],-lenarrow*cos(delta3D[slice]),-lenarrow*sin(delta3D[slice]),overhang = .125, width = 0.0025,head_width = .03,color = plot_cycle[1])
PyPlot.arrow(-X[slice],Y[slice],-lenarrow*cos(delta3D[slice]),0.0,overhang = .125, width = 0.0025,head_width = .03,color = plot_cycle[1])
# Dashed Delta Line
circx, circy = getxy(LinRange(-pi/2,delta3D[slice]*1.7,20),R*.3)
PyPlot.plot(circx.-X[slice].+(R*.29)/2,circy.+Y[slice],"k--",linewidth = .75)
# Force Vector Symbols
PyPlot.annotate(L"\delta", xy=(-1.5, 0.0), xytext=(-X[slice]-.25, Y[slice]+.03))
PyPlot.annotate(L"\hat{n}", xy=(-1.5, 0.0), xytext=(-X[slice]-.44, Y[slice]+.28))
PyPlot.annotate(L"\hat{r}", xy=(-1.5, 0.0), xytext=(-X[slice]-.42, Y[slice]-.03))

# Velocity Vectors
slice = 26
lenarrow = 0.3
PyPlot.annotate(L"velocity", xy=(-1.5, 0), xytext=(-0.95, H*0.87))
PyPlot.arrow(-X[slice],Y[slice],lenarrow*cos(delta3D[slice])*cos(delta3D[slice]),lenarrow*cos(delta3D[slice])*sin(delta3D[slice]),overhang = .125, width = 0.0025,head_width = .03,color = plot_cycle[2])
PyPlot.arrow(-X[slice],Y[slice],lenarrow,0.0,overhang = .125, width = 0.0025,head_width = .03,color = plot_cycle[2])
# Dashed Delta Line
circx, circy = getxy(LinRange(pi/2-delta3D[slice]*0.5,pi/2,20),R*.3)
PyPlot.plot(circx.-X[slice].+(R*-.17),circy.+Y[slice],"k--",linewidth = .75)
# Force Vector Symbols
PyPlot.annotate(L"\delta", xy=(-1.5, 0.0), xytext=(-X[slice]+.15, Y[slice]-.13))
PyPlot.annotate(L"\hat{n}", xy=(-1.5, 0.0), xytext=(-X[slice]+.15, Y[slice]-.3))
PyPlot.annotate(L"\hat{r}", xy=(-1.5, 0.0), xytext=(-X[slice]+.4, Y[slice]-.05))

# Z axis
PyPlot.annotate(L"\hat{z}", xy=(-1.5, 0), xytext=(-0.15, H*1.1))
PyPlot.arrow(0,H*1.1,0.0,0.001,overhang = .125, head_width = .04,facecolor = "k")
PyPlot.plot([0.0,0.0],[shapeY[bottomslice],H*1.1],"k--",linewidth = .75)
# Blade Shape
PyPlot.plot(-shapeX[bottomslice:end],shapeY[bottomslice:end],"-k")
PyPlot.plot(0,H*1.2,"w.")

# # Blade Shape 2
# PyPlot.plot(shapeX[bottomslice:end],shapeY[bottomslice:end],"-k")
# PyPlot.plot(0,H*1.2,"w.")

PyPlot.axis("equal")
PyPlot.axis("off")
PyPlot.savefig("$(path)/figs/vertical_delta.pdf",transparent = true)

#
#
# PLOT Disks vs Curved
#
#

R = 1.0#5.0/2 #m
H = 1.02*R*2; #m
hr = H/R
n_slices = 30
shapeY = collect(LinRange(0,H,n_slices+1))
shapeX = R.*(1.0.-4.0.*(shapeY/H.-.5).^2)

shapeY2 = collect(LinRange(0,H,Int(n_slices/2)+1))
shapeX2 = R.*(1.0.-4.0.*(shapeY2/H.-.5).^2)

shapeYslice = zeros(Int(n_slices/2)*2-1)
shapeXslice = zeros(Int(n_slices/2)*2-1)

shapeXslice[1:2:end] = shapeX2[1:end-1]
shapeXslice[2:2:end] = shapeX2[2:end-1]
shapeXslice = [shapeXslice;0.0]
shapeXslice[2:3] = shapeXslice[2:3].-0.02*2.6
shapeXslice[end-2:end-1] = shapeXslice[end-2:end-1].-0.02*2.6

shapeYslice[1:2:end] = shapeY2[1:end-1]
shapeYslice[2:2:end] = shapeY2[1:end-2]
shapeYslice = [shapeYslice;shapeY2[end-1]]
shapeYslice[1:2] = shapeYslice[1:2].-0.025*2.6
shapeYslice[3:4] = shapeYslice[3:4].-0.025/2*2.6
shapeYslice[end-1:end] = shapeYslice[end-1:end].+0.025*2.6
shapeYslice[end-3:end-2] = shapeYslice[end-3:end-2].+0.025/2*2.6
Y = (shapeY[2:end] + shapeY[1:end-1])/2.0;
X = (shapeX[2:end] + shapeX[1:end-1])/2.0;
delta_xs = shapeX[2:end] - shapeX[1:end-1]
delta_zs = shapeY[2:end] - shapeY[1:end-1]
delta3D = atan.(delta_xs./delta_zs)
slice = 22
bottomslice = 1
lenarrow = 0.4

rc("figure", figsize=(3, 4))
PyPlot.figure()
# Blade Shape
PyPlot.plot(-shapeX[bottomslice:end],shapeY[bottomslice:end],"-k")
PyPlot.plot(-shapeXslice,shapeYslice.+0.035*2,"-",color = plot_cycle[1])
PyPlot.plot(shapeX[bottomslice:end],shapeY[bottomslice:end],"-k")
PyPlot.plot(shapeXslice,shapeYslice.+0.035*2,"-",color = plot_cycle[1])

for ii = 1:2:length(shapeXslice)
    PyPlot.plot([-shapeXslice[ii], shapeXslice[ii]],[shapeYslice[ii].+0.035*2, shapeYslice[ii].+0.035*2],"-",color = plot_cycle[1])
end

PyPlot.axis("equal")
PyPlot.axis("off")
# PyPlot.legend(["Curved","Straight"])
PyPlot.savefig("$(path)/figs/stacked_cylinders2.pdf",transparent = true)

#Slice Lines

rc("figure", figsize=(3, 4))
PyPlot.figure()
# Blade Shape
PyPlot.plot(-shapeX[bottomslice:end],shapeY[bottomslice:end],"-k")
PyPlot.plot(shapeX[bottomslice:end],shapeY[bottomslice:end],"-k")

for ii = 1:length(shapeX)
    PyPlot.plot([-shapeX[ii], shapeX[ii]],[shapeY[ii], shapeY[ii]],"-",color = plot_cycle[1])
end
PyPlot.axis("equal")
PyPlot.axis("off")
# PyPlot.legend(["Curved","Straight"])
PyPlot.savefig("$(path)/figs/stacked_cylinders_lines.pdf",transparent = true)


#
#
# DMS
#
#

R = 1.0
ntheta = 30

dtheta = 2*pi/(ntheta) # Each streamtube has two parts
thetavec = [collect(dtheta/2:dtheta:2*pi+dtheta);0.0]
theta_st = thetavec.-dtheta/2

PyPlot.figure()

# Theta discretization
PyPlot.plot(R*sin.(thetavec),R*cos.(thetavec),"k.")

# Swept Disk
circx, circy = getxy(LinRange(0,2*pi,ntheta*10),R)
PyPlot.plot(circx,circy,"k",linewidth = .75)



#Upstream
for i = 1:15
    if thetavec[i] > pi/2
        offset = dtheta
        offset2 = 0
    else
        offset = 0
        offset2 = dtheta
    end
    # Streamtube boundaries
    xbndry = R*[-sin(thetavec[i]-offset),sin(thetavec[i]-offset2)]
    ybndry = R*[cos(theta_st[i]),cos(theta_st[i])]
    PyPlot.plot(xbndry,ybndry, "-",linewidth = .6,color = (0.0,0.0,0.0))

    # Streamtube disks
    xdisk = R*[-sin(thetavec[i]),-sin(thetavec[i])]
    ydisk = R*[cos(theta_st[i]),cos(theta_st[i+1])]
    PyPlot.plot(xdisk,ydisk, "k-",linewidth = .55)
end

#Downstream
for i = 16:31
    if thetavec[i] > 3pi/2
        offset = 0
        offset2 = 0
    else
        offset = dtheta
        offset2 = dtheta
    end
    # Streamtube boundaries
    xbndry = R*[-sin(thetavec[i]-offset),R*1.5]
    ybndry = R*[cos(theta_st[i]),cos(theta_st[i])]
    PyPlot.plot(xbndry,ybndry, "-",linewidth = .35,color = (0.7,0.7,0.7))

    # Streamtube disks
    xdisk = R*[-sin(thetavec[i]),-sin(thetavec[i])]
    ydisk = R*[cos(theta_st[i]),cos(theta_st[i+1])]
    PyPlot.plot(xdisk,ydisk, "k-",linewidth = .55)
end

PyPlot.axis("equal")
PyPlot.xlim(-1.5*R,1.5*R)
PyPlot.ylim(-R,R)
# PyPlot.xlabel("X")
# PyPlot.ylabel("Y")


#Airfoils
blade_offset = 35
for thetablade = [0+blade_offset,120+blade_offset,240+blade_offset]
    NACA0015x,NACA0015y = af_placement(thetablade*pi/180,.6,.89)
    PyPlot.plot(NACA0015x,NACA0015y,color = plot_cycle[1],linewidth = .75)
    PyPlot.plot(mean(NACA0015x),mean(NACA0015y),"k.")
    local circx, circy = getxy(thetablade*pi/180,r*1)
    # PyPlot.plot(circx,circy,"r.", markersize = 5)
end

#Axes
PyPlot.plot([-1.5,1.5],[0,0],"k-.",linewidth = .5)
PyPlot.plot([0,0],[-1.5,1.5],"k-.",linewidth = .5)
PyPlot.annotate(L"x,u", xy=(-1.5, 0), xytext=(1.6, -0.0500))
PyPlot.arrow(1.5,0,.003,0.0,overhang = .125, head_width = .04,facecolor = "k")
PyPlot.annotate(L"y,v", xy=(-1.5, 0), xytext=(-0.15, 1.6))
PyPlot.arrow(0,1.45,0.0,0.001,overhang = .125, head_width = .04,facecolor = "k")

#Theta
Ndots = 30
PyPlot.plot(LinRange(0,-1.5*sin(blade_offset*pi/180),Ndots),LinRange(0,1.5*cos(blade_offset*pi/180),Ndots),"k.",markersize = .5)
PyPlot.annotate(L"\theta", xy=(-1.5, 0), xytext=(-.19, .45))
local circx, circy = getxy(LinRange(0,-20*pi/180,Ntheta),r*.4)
PyPlot.plot(circx,circy,"k",linewidth = .75)
PyPlot.arrow(circx[end],circy[end],-.01,-.007,overhang = .125, head_width = .04,facecolor = "k")

#Vinf
PyPlot.annotate(L"V_{\infty}", xy=(-1.5, 0), xytext=(-2.5, -0.10))
PyPlot.arrow(-2.23,.0,.55,0,overhang = .125, width = 0.018,head_width = .06,facecolor = "k")
PyPlot.plot(-2.0,0.0,"w")
axis("off")
axis("equal")
PyPlot.savefig("$(path)/figs/vawt_dms_slice.pdf",transparent = true)
