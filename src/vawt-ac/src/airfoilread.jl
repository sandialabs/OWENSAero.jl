# module AirfoilRead: From the repository associated with: doi:10.5194/wes-1-327-2016 Actuator cylinder theory for multiple vertical axis wind turbines

# export AirfoilData, readaerodyn

# import Interpolations: interpolate, Gridded, Linear, GriddedInterpolation
import FLOWMath
import Interpolations

_parse_airfoil_header_number(line) = parse(Float64, split(strip(line))[end])

function _read_first_reynolds_airfoil_table(filename)
    alpha = Float64[]
    cl = Float64[]
    cd = Float64[]
    cm = Float64[]

    thickness_to_chord = 0.0
    zero_lift_aoa = 0.0
    positive_stall_aoa = 10 * pi / 180
    negative_stall_aoa = -10 * pi / 180
    in_table = false

    open(filename) do f
        for line in eachline(f)
            stripped = strip(line)
            isempty(stripped) && continue

            if startswith(stripped, "Thickness to Chord Ratio")
                thickness_to_chord = _parse_airfoil_header_number(stripped)
            elseif startswith(stripped, "Zero Lift AOA")
                zero_lift_aoa = _parse_airfoil_header_number(stripped) * pi / 180
            elseif startswith(stripped, "BV Dyn. Stall Model - Positive Stall AOA")
                positive_stall_aoa = _parse_airfoil_header_number(stripped) * pi / 180
            elseif startswith(stripped, "BV Dyn. Stall Model - Negative Stall AOA")
                negative_stall_aoa = _parse_airfoil_header_number(stripped) * pi / 180
            elseif startswith(stripped, "AOA")
                in_table = true
            elseif stripped == "EOT" || startswith(stripped, "Reynolds Number")
                if in_table && !isempty(alpha)
                    break
                end
            elseif in_table
                parts = split(stripped)
                length(parts) >= 3 ||
                    throw(ArgumentError("airfoil data row in $(filename) has fewer than 3 columns"))
                push!(alpha, parse(Float64, parts[1]))
                push!(cl, parse(Float64, parts[2]))
                push!(cd, parse(Float64, parts[3]))
                push!(cm, length(parts) >= 4 ? parse(Float64, parts[4]) : 0.0)
            end
        end
    end

    isempty(alpha) && throw(ArgumentError("no airfoil data rows found in $(filename)"))
    return (
        alpha = alpha,
        cl = cl,
        cd = cd,
        cm = cm,
        thickness_to_chord = thickness_to_chord,
        zero_lift_aoa = zero_lift_aoa,
        positive_stall_aoa = positive_stall_aoa,
        negative_stall_aoa = negative_stall_aoa,
    )
end

"""
    readaerodyn(filename)

create airfoil lookup for the first Reynolds-number table in an OWENSAero
airfoil file

# Inputs
* `filename::string`: file path/name to airfoil file formatted like in the test folder

# Outputs:
* `af::function`: cl, cd = af(alpha,re,mach) with alpha in rad. Use
  `af(alpha,re,mach; return_cm=true)` to also return Cm25.

"""
function readaerodyn(filename)
    table = _read_first_reynolds_airfoil_table(filename)
    alpha = table.alpha
    cl = table.cl
    cd = table.cd
    cm = table.cm

    # af = AirfoilData(alpha*pi/180.0, cl, cd)

    # 1D interpolations for now.  ignoring Re dependence (which is very minor)
    # afcl = interpolate((alpha*pi/180.0,), cl, Gridded(Linear()))
    # afcd = interpolate((alpha*pi/180.0,), cd, Gridded(Linear()))
    # af = AirfoilData(afcl, afcd)

    afcl = FLOWMath.Akima(alpha*pi/180, cl)
    afcd = FLOWMath.Akima(alpha*pi/180, cd)
    afcm = FLOWMath.Akima(alpha*pi/180, cm)

    function af(alphain,re,mach; return_cm=false)
        # if alphain>maximum(alpha)*pi/180 || alphain<minimum(alpha)*pi/180
        #     @error "aoa is greater or less than what is defined"
        # end
            cl = afcl(alphain)
            cd = afcd(alphain)
            cm = afcm(alphain)

        return return_cm ? (cl, cd, cm) : (cl, cd)
    end

    return af#, alpha*pi/180, cl, cd
end

"""
    readaerodyn_BV(filename)

create airfoil lookup function with Boeing-Vertol dynamic stall metadata from
the first Reynolds-number table in an OWENSAero airfoil file

# Inputs
* `filename::string`: file path/name to airfoil file formatted like in the test folder

# Outputs:
* `af::function`: cl, cd = af_BV(alpha,Re,M,env,V_twist,c,dt,U;solvestep=false) with alpha in rad, OWENSAero.Env, V_twist in rad/s, c chord in m, dt in sec, U Vloc in m/s, solvestep true during solve loop. Use `return_cm=true` to also return Cm25.

"""
function readaerodyn_BV(filename) #TODO: use multiple dispatch to simplify this
    table = _read_first_reynolds_airfoil_table(filename)
    alpha = table.alpha
    cl = table.cl
    cd = table.cd
    cm = table.cm

    # af = AirfoilData(alpha*pi/180.0, cl, cd)

    # 1D interpolations for now.  ignoring Re dependence (which is very minor)
    # afcl = interpolate((alpha*pi/180.0,), cl, Gridded(Linear()))
    # afcd = interpolate((alpha*pi/180.0,), cd, Gridded(Linear()))
    # af = AirfoilData(afcl, afcd)

    afcl = FLOWMath.Akima(alpha*pi/180, cl)
    afcd = FLOWMath.Akima(alpha*pi/180, cd)
    afcm = FLOWMath.Akima(alpha*pi/180, cm)
    cl=0.0
    function af(alphain,Re,umach,family_factor)
        if alphain>maximum(alpha)*pi/180 || alphain<minimum(alpha)*pi/180
            @error "aoa is greater or less than what is defined"
        end

        cl = afcl(alphain)
        cd = afcd(alphain)
        cm = afcm(alphain)

        return cl, cd, cm
    end

    aoaStallPos = table.positive_stall_aoa
    aoaStallNeg = table.negative_stall_aoa
    AOA0 = table.zero_lift_aoa
    tc = table.thickness_to_chord

    function af_BV(alpha,Re,M,env,V_twist,c,dt,U;solvestep=false,idx=1,return_cm=false)
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
        return return_cm ? (CL_BV, CD_BV, CM_BV) : (CL_BV, CD_BV)
    end

    return af_BV #, alpha*pi/180, cl, cd
end

"""
    readaerodyn_BV_NEW(filename;DynamicStallModel="BV")

for a file with multiple reynolds numbers create airfoil lookup function with boeing vertol dynamic stall model and wrap interpolation

# Inputs
* `filename::string`: file path/name to airfoil file formatted like in the test folder
* `DynamicStallModel::string`: "BV", "LB", or "none"

# Outputs:
* `af::function`: cl, cd = af_BV(alpha,Re,M,env,V_twist,c,dt,U;solvestep=false) with alpha in rad, OWENSAero.Env, V_twist in rad/s, c chord in m, dt in sec, U Vloc in m/s, solvestep true during solve loop. Use `return_cm=true` to also return Cm25.
* `af::function`: cl, cd = af(alpha,re,mach) with alpha in rad. Use `return_cm=true` to also return Cm25.
"""
function readaerodyn_BV_NEW(filename;DynamicStallModel="BV") #TODO: use multiple dispatch to simplify this
    DynamicStallModel = _canonical_dynamic_stall_model(DynamicStallModel)

    # Loop through the file and determine how many Reynolds numbers there are
    nRe = 0
    open(filename) do f
        for line in eachline(f)
            if startswith(strip(line), "Reynolds Number")
                nRe += 1
            end
        end
    end

    # Set up arrays
    alphas = [collect(-180:1:-21);collect(-20:0.5:20);collect(21:1:180)].*pi/180
    cls = zeros(length(alphas),nRe)
    cds = zeros(length(alphas),nRe)
    cms = zeros(length(alphas),nRe)
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
                cmOrig = Float64[]
                iRe += 1
                REs[iRe] = parse(Float64, split(line)[end])
                aoaStallPosVec[iRe] = parse(Float64, split(readline(f))[end])
                aoaStallNegVec[iRe] = parse(Float64, split(readline(f))[end])
                clalphaVec[iRe] = parse(Float64, split(readline(f))[end])
                clcritposVec[iRe] = parse(Float64, split(readline(f))[end])
                clcritnegVec[iRe] = parse(Float64, split(readline(f))[end])
                line = readline(f) #skip cl cd header

                line = readline(f) #this should have the data
                while !isempty(strip(line)) && strip(line) != "EOT"
                    parts = split(line)
                    push!(alphaOrig, parse(Float64,parts[1]))
                    push!(clOrig, parse(Float64,parts[2]))
                    push!(cdOrig, parse(Float64,parts[3]))
                    push!(cmOrig, length(parts) >= 4 ? parse(Float64,parts[4]) : 0.0)
                    line = readline(f) # reload line until it hits the next blank spot, then break
                end

                # Interpolate to a uniform aoa array and store
                cls[:,iRe] = safeakima(alphaOrig*pi/180,clOrig,alphas)
                cds[:,iRe] = safeakima(alphaOrig*pi/180,cdOrig,alphas)
                cms[:,iRe] = safeakima(alphaOrig*pi/180,cmOrig,alphas)
            end
        end
    end

    if nRe == 0
        throw(ArgumentError("No Reynolds-number airfoil tables were found in $(repr(filename))."))
    elseif nRe == 1
        itpcl = Interpolations.interpolate((alphas,), cls[:,1], Interpolations.Gridded(Interpolations.Linear()))
        cl_alpha = Interpolations.extrapolate(itpcl, Interpolations.Flat())

        itpcd = Interpolations.interpolate((alphas,), cds[:,1], Interpolations.Gridded(Interpolations.Linear()))
        cd_alpha = Interpolations.extrapolate(itpcd, Interpolations.Flat())

        itpcm = Interpolations.interpolate((alphas,), cms[:,1], Interpolations.Gridded(Interpolations.Linear()))
        cm_alpha = Interpolations.extrapolate(itpcm, Interpolations.Flat())

        clspl = (alpha, Re) -> cl_alpha(alpha)
        cdspl = (alpha, Re) -> cd_alpha(alpha)
        cmspl = (alpha, Re) -> cm_alpha(alpha)
        aoaStallPosspl = Re -> aoaStallPosVec[1]
        aoaStallNegspl = Re -> aoaStallNegVec[1]
        clalphaspl = Re -> clalphaVec[1]
        clcritposspl = Re -> clcritposVec[1]
        clcritnegspl = Re -> clcritnegVec[1]
    else
        # clspl = Dierckx.Spline2D(alphas, REs, cls)
        # cdspl = Dierckx.Spline2D(alphas, REs, cds)
        # Assume cls is filled appropriately: size = (length(alphas), length(REs))
        itpcl = Interpolations.interpolate((alphas, REs), cls, Interpolations.Gridded(Interpolations.Linear()))
        clspl = Interpolations.extrapolate(itpcl, Interpolations.Flat())

        itpcd = Interpolations.interpolate((alphas, REs), cds, Interpolations.Gridded(Interpolations.Linear()))
        cdspl = Interpolations.extrapolate(itpcd, Interpolations.Flat())

        itpcm = Interpolations.interpolate((alphas, REs), cms, Interpolations.Gridded(Interpolations.Linear()))
        cmspl = Interpolations.extrapolate(itpcm, Interpolations.Flat())

        aoaStallPosspl = FLOWMath.Akima(REs,aoaStallPosVec)
        aoaStallNegspl = FLOWMath.Akima(REs,aoaStallNegVec)
        clalphaspl = FLOWMath.Akima(REs, clalphaVec)
        clcritposspl = FLOWMath.Akima(REs, clcritposVec)
        clcritnegspl = FLOWMath.Akima(REs, clcritnegVec)
    end

    # Create the base airfoil wrapper function
    function af2(alpha,Re,umach=0.0,family_factor=1.0; return_cm=false)

        # if length(REs)>1
            # cl = Dierckx.evaluate(clspl,alpha, Re)
            # cd = Dierckx.evaluate(cdspl,alpha, Re)
            cl = clspl(alpha, Re)
            cd = cdspl(alpha, Re)
            cm = cmspl(alpha, Re)
            # cl = FLOWMath.interp2d(safeakima, alphas, REs, cls, [alpha], [Re])[1]
            # cd = FLOWMath.interp2d(safeakima, alphas, REs, cds, [alpha], [Re])[1]
        # else
        #     cl = safeakima(alphaOrig*pi/180, clOrig,alpha)
        #     cd = safeakima(alphaOrig*pi/180, cdOrig,alpha)
        # end

        return return_cm ? (cl, cd, cm) : (cl, cd)
    end

    # end

    af2_with_cm(alpha, Re, umach=0.0, family_factor=1.0) =
        af2(alpha, Re, umach, family_factor; return_cm=true)

    function af_BV2(alpha,Re,M,env,V_twist,c,dt,U;solvestep=false,idx=1,return_cm=false)
        if idx==1
            idxmin1 = length(env.alpha_last)
        else
            idxmin1 = idx-1
        end
        dalpha=alpha-env.alpha_last[idxmin1] + V_twist*dt #TODO: verify sign of motion induced velocity
        adotnorm=dalpha/dt*c/(2.0*FLOWMath.ksmax([U,0.001]))
        # if length(REs)>1
            aoaStallPos = aoaStallPosspl(Re)*pi/180#safeakima(REs,aoaStallPosVec,Re)*pi/180
            aoaStallNeg = aoaStallNegspl(Re)*pi/180#safeakima(REs,aoaStallNegVec,Re)*pi/180
        # else
        #     aoaStallPos = aoaStallPosVec[1]
        #     aoaStallNeg = aoaStallNegVec[1]
        # end
        CL_BV, CD_BV, CM_BV, env.BV_DynamicFlagL[idxmin1], env.BV_DynamicFlagD[idxmin1] = OWENSAero.Boeing_Vertol(af2_with_cm,alpha,adotnorm,M,Re,aoaStallPos,aoaStallNeg,AOA0,tc,env.BV_DynamicFlagL[idxmin1],env.BV_DynamicFlagD[idxmin1], family_factor = 0.0)
        if !solvestep #don't update this while it is in the DMS independent solve loop #TODO: should we include the flags as well?
            env.alpha_last[idx] = alpha
            env.BV_DynamicFlagL[idx] = env.BV_DynamicFlagL[idxmin1]
            env.BV_DynamicFlagD[idx] = env.BV_DynamicFlagD[idxmin1]
        end
        return return_cm ? (CL_BV, CD_BV, CM_BV) : (CL_BV, CD_BV)
    end

    function af_LB2(alpha,Re,M,env,V_twist,c,dt,U;solvestep=false,idx=1,return_cm=false)
        idxmin1 = idx == 1 ? length(env.lb_state.cl_ref_last) : idx - 1
        U_safe = FLOWMath.ksmax([U, 0.001])
        ds = 2.0 * U_safe * dt / c
        alpha75 = alpha
        alpha50 = alpha
        state = env.lb_state
        CL_LB,
        CD_LB,
        CM_LB,
        cl_ref_le_last,
        cl_ref_last,
        fstat_last,
        cv_last,
        dp,
        df,
        dcnv,
        le_separation_state,
        s_lev = OWENSAero.Leishman_Beddoes(
            af2_with_cm,
            alpha75,
            alpha50,
            ds,
            M,
            Re,
            clalphaspl(Re),
            AOA0,
            clcritposspl(Re),
            clcritnegspl(Re),
            state.cl_ref_le_last[idxmin1],
            state.cl_ref_last[idxmin1],
            state.fstat_last[idxmin1],
            state.cv_last[idxmin1],
            state.dp[idxmin1],
            state.df[idxmin1],
            state.dcnv[idxmin1],
            state.le_separation_state[idxmin1],
            state.s_lev[idxmin1],
            family_factor = 0.0,
        )
        if !solvestep
            state.cl_ref_le_last[idx] = cl_ref_le_last
            state.cl_ref_last[idx] = cl_ref_last
            state.fstat_last[idx] = fstat_last
            state.cv_last[idx] = cv_last
            state.dp[idx] = dp
            state.df[idx] = df
            state.dcnv[idx] = dcnv
            state.le_separation_state[idx] = le_separation_state
            state.s_lev[idx] = s_lev
        end
        return return_cm ? (CL_LB, CD_LB, CM_LB) : (CL_LB, CD_LB)
    end

    if DynamicStallModel == "BV"
        return af_BV2 #, alpha*pi/180, cl, cd
    elseif DynamicStallModel == "LB"
        return af_LB2
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
