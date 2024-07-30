# module AirfoilRead: From the repository associated with: doi:10.5194/wes-1-327-2016 Actuator cylinder theory for multiple vertical axis wind turbines

# export AirfoilData, readaerodyn

# import Interpolations: interpolate, Gridded, Linear, GriddedInterpolation
import FLOWMath
import Dierckx: Spline1D, evaluate

"""
    readaerodyn(filename)

create airfoil lookup for a file with only one reynolds number

# Inputs
* `filename::string`: file path/name to airfoil file formatted like in the test folder

# Outputs:
* `af::function`: cl, cd = af(alpha,re,mach) with alpha in rad

"""
function readaerodyn(filename)
    """currently only reads one Reynolds number if multiple exist"""

    alpha = Float64[]
    cl = Float64[]
    cd = Float64[]

    open(filename) do f

        # skip header
        for i = 1:13
            readline(f)
        end

        # read until EOT
        while true
            line = readline(f)
            if occursin(line, "EOT")
                break
            end
            parts = split(line)
            push!(alpha, parse(Float64,parts[1]))
            push!(cl, parse(Float64,parts[2]))
            push!(cd, parse(Float64,parts[3]))
        end
    end

    # af = AirfoilData(alpha*pi/180.0, cl, cd)

    # 1D interpolations for now.  ignoring Re dependence (which is very minor)
    # afcl = interpolate((alpha*pi/180.0,), cl, Gridded(Linear()))
    # afcd = interpolate((alpha*pi/180.0,), cd, Gridded(Linear()))
    # af = AirfoilData(afcl, afcd)

    afcl = FLOWMath.Akima(alpha*pi/180, cl)
    afcd = FLOWMath.Akima(alpha*pi/180, cd)

    function af(alpha,re,mach)

        cl = afcl(alpha)
        cd = afcd(alpha)

        return cl, cd
    end

    return af#, alpha*pi/180, cl, cd
end

"""
    readaerodyn_BV(filename)

create airfoil lookup function with boeing vertol dynamic stall model for a file with only one reynolds number

# Inputs
* `filename::string`: file path/name to airfoil file formatted like in the test folder

# Outputs:
* `af::function`: cl, cd = af_BV(alpha,Re,M,env,V_twist,c,dt,U;solvestep=false) with alpha in rad, OWENSAero.Env, V_twist in rad/s, c chord in m, dt in sec, U Vloc in m/s, solvestep true during solve loop

"""
function readaerodyn_BV(filename) #TODO: use multiple dispatch to simplify this
    """currently only reads one Reynolds number if multiple exist"""

    alpha = Float64[]
    cl = Float64[]
    cd = Float64[]

    open(filename) do f

        # skip header
        for i = 1:13
            readline(f)
        end

        # read until EOT
        while true
            line = readline(f)
            if occursin(line, "EOT")
                break
            end
            parts = split(line)
            push!(alpha, parse(Float64,parts[1]))
            push!(cl, parse(Float64,parts[2]))
            push!(cd, parse(Float64,parts[3]))
        end
    end

    # af = AirfoilData(alpha*pi/180.0, cl, cd)

    # 1D interpolations for now.  ignoring Re dependence (which is very minor)
    # afcl = interpolate((alpha*pi/180.0,), cl, Gridded(Linear()))
    # afcd = interpolate((alpha*pi/180.0,), cd, Gridded(Linear()))
    # af = AirfoilData(afcl, afcd)

    afcl = FLOWMath.Akima(alpha*pi/180, cl)
    afcd = FLOWMath.Akima(alpha*pi/180, cd)
    cl=0.0
    function af(alpha,Re,umach,family_factor)

        cl = afcl(alpha)
        cd = afcd(alpha)

        return cl, cd, 0.0
    end

    #TODO: Get these from the airfoil data
    aoaStallPos = 10*pi/180
    aoaStallNeg = -10*pi/180
    AOA0 = 0.0
    tc = 0.12

    function af_BV(alpha,Re,M,env,V_twist,c,dt,U;solvestep=false,idx=1)
        if idx==1
            idxmin1 = length(env.alpha_last)
        else
            idxmin1 = idx-1
        end
        dalpha=alpha-env.alpha_last[idxmin1] + V_twist*dt #TODO: verify sign of motion induced velocity
        adotnorm=dalpha/dt*c/(2.0*FLOWMath.ksmax([U,0.001]))
        CL_BV, CD_BV, CM_BV, env.BV_DynamicFlagL[idxmin1], env.BV_DynamicFlagD[idxmin1] = OWENSAero.Boeing_Vertol(af,alpha,adotnorm,M,Re,aoaStallPos,aoaStallNeg,AOA0,tc,env.BV_DynamicFlagL[idxmin1],env.BV_DynamicFlagD[idxmin1], family_factor = 0.0)
        if !solvestep #don't update this while it is in the DMS independent solve loop
            env.alpha_last[idx] = alpha
            env.BV_DynamicFlagL[idx] = env.BV_DynamicFlagL[idxmin1]
            env.BV_DynamicFlagD[idx] = env.BV_DynamicFlagD[idxmin1]
        end
        return CL_BV, CD_BV #TODO: add CM
    end

    return af_BV #, alpha*pi/180, cl, cd
end

"""
    readaerodyn_BV_NEW(filename;DSModel="BV")

for a file with multiple reynolds numbers create airfoil lookup function with boeing vertol dynamic stall model and wrap interpolation

# Inputs
* `filename::string`: file path/name to airfoil file formatted like in the test folder
* `DSModel::string`: "BV" or "none"

# Outputs:
* `af::function`: cl, cd = af_BV(alpha,Re,M,env,V_twist,c,dt,U;solvestep=false) with alpha in rad, OWENSAero.Env, V_twist in rad/s, c chord in m, dt in sec, U Vloc in m/s, solvestep true during solve loop
* `af::function`: cl, cd = af(alpha,re,mach) with alpha in rad
"""
function readaerodyn_BV_NEW(filename;DSModel="BV") #TODO: use multiple dispatch to simplify this
    # Loop through the file and determine how many Reynolds numbers there are
    nRe = 0
    open(filename) do f
        for line in eachline(f)
            if occursin(line, "Reynolds")
                nRe += 1
            end
        end
    end

    # Set up arrays
    alphas = [collect(-180:1:-21);collect(-20:0.5:20);collect(21:1:180)].*pi/180
    cls = zeros(length(alphas),nRe)
    cds = zeros(length(alphas),nRe)
    #TODO: cms
    REs = zeros(nRe)
    aoaStallPosVec = zeros(nRe)
    aoaStallNegVec = zeros(nRe)
    clalphaVec = zeros(nRe)
    clcritposVec = zeros(nRe)
    clcritnegVec = zeros(nRe)

    tc = 0.0
    AOA0 = 0.0

    iRe = 0

    # alphaOrig = Float64[]
    # clOrig = Float64[]
    # cdOrig = Float64[]

    open(filename) do f
        # Extract the top level header information
        line = readline(f) #Title
        tc = parse(Float64, split(readline(f))[end]) #tc ratio
        AOA0 = parse(Float64, split(readline(f))[end])*pi/180 #zero lift aoa
        for line in eachline(f)
            if !isempty(line) && split(line)[1]=="Reynolds"
                alphaOrig = Float64[]
                clOrig = Float64[]
                cdOrig = Float64[]
                iRe += 1
                REs[iRe] = parse(Float64, split(line)[end])
                aoaStallPosVec[iRe] = parse(Float64, split(readline(f))[end])
                aoaStallNegVec[iRe] = parse(Float64, split(readline(f))[end])
                clalphaVec[iRe] = parse(Float64, split(readline(f))[end])
                clcritposVec[iRe] = parse(Float64, split(readline(f))[end])
                clcritnegVec[iRe] = parse(Float64, split(readline(f))[end])
                line = readline(f) #skip cl cd header

                line = readline(f) #this should have the data
                while !isempty(line)
                    parts = split(line)
                    push!(alphaOrig, parse(Float64,parts[1]))
                    push!(clOrig, parse(Float64,parts[2]))
                    push!(cdOrig, parse(Float64,parts[3]))
                    line = readline(f) # reload line until it hits the next blank spot, then break
                end

                # Interpolate to a uniform aoa array and store
                cls[:,iRe] = FLOWMath.akima(alphaOrig*pi/180,clOrig,alphas)
                cds[:,iRe] = FLOWMath.akima(alphaOrig*pi/180,cdOrig,alphas)
            end
        end
    end

    clspl = Dierckx.Spline2D(alphas, REs, cls)
    cdspl = Dierckx.Spline2D(alphas, REs, cds)
    aoaStallPosspl = FLOWMath.Akima(REs,aoaStallPosVec)
    aoaStallNegspl = FLOWMath.Akima(REs,aoaStallNegVec)

    # Create the base airfoil wrapper function
    function af2(alpha,Re,umach=0.0,family_factor=1.0)

        # if length(REs)>1
            cl = Dierckx.evaluate(clspl,alpha, Re)
            cd = Dierckx.evaluate(cdspl,alpha, Re)
            # cl = FLOWMath.interp2d(FLOWMath.akima, alphas, REs, cls, [alpha], [Re])[1]
            # cd = FLOWMath.interp2d(FLOWMath.akima, alphas, REs, cds, [alpha], [Re])[1]
        # else
        #     cl = FLOWMath.akima(alphaOrig*pi/180, clOrig,alpha)
        #     cd = FLOWMath.akima(alphaOrig*pi/180, cdOrig,alpha)
        # end

        return cl, cd, 0.0
    end

    # end

    function af_BV2(alpha,Re,M,env,V_twist,c,dt,U;solvestep=false,idx=1)
        if idx==1
            idxmin1 = length(env.alpha_last)
        else
            idxmin1 = idx-1
        end
        dalpha=alpha-env.alpha_last[idxmin1] + V_twist*dt #TODO: verify sign of motion induced velocity
        adotnorm=dalpha/dt*c/(2.0*FLOWMath.ksmax([U,0.001]))
        # if length(REs)>1
            aoaStallPos = aoaStallPosspl(Re)*pi/180#FLOWMath.akima(REs,aoaStallPosVec,Re)*pi/180
            aoaStallNeg = aoaStallNegspl(Re)*pi/180#FLOWMath.akima(REs,aoaStallNegVec,Re)*pi/180
        # else
        #     aoaStallPos = aoaStallPosVec[1]
        #     aoaStallNeg = aoaStallNegVec[1]
        # end
        CL_BV, CD_BV, CM_BV, env.BV_DynamicFlagL[idxmin1], env.BV_DynamicFlagD[idxmin1] = OWENSAero.Boeing_Vertol(af2,alpha,adotnorm,M,Re,aoaStallPos,aoaStallNeg,AOA0,tc,env.BV_DynamicFlagL[idxmin1],env.BV_DynamicFlagD[idxmin1], family_factor = 0.0)
        if !solvestep #don't update this while it is in the DMS independent solve loop #TODO: should we include the flags as well?
            env.alpha_last[idx] = alpha
            env.BV_DynamicFlagL[idx] = env.BV_DynamicFlagL[idxmin1]
            env.BV_DynamicFlagD[idx] = env.BV_DynamicFlagD[idxmin1]
        end
        return CL_BV, CD_BV #TODO: add CM
    end

    if DSModel == "BV"
        return af_BV2 #, alpha*pi/180, cl, cd
    else
        return af2
    end
end

# end

# filename = "airfoils/naca0015-wt.dat"
# af, a0, cl0, cd0 = readaerodyn(filename)
# alpha = linspace(-pi, pi, 2000)
# cl, cd = af(alpha)
# using PyPlot
# figure()
# plot(alpha, cl)
# plot(a0, cl0, "o")
# figure()
# plot(alpha, cd)
# plot(a0, cd0, "o")
# show()
