"""
    Boeing_Vertol(af,alpha,adotnorm,umach,Re,aoaStallPos,aoaStallNeg,AOA0,tc,BV_DynamicFlagL,BV_DynamicFlagD; family_factor = 0.0)
Boeing-Vertol Dynamic Stall Model. All angles are in rad unless explicitely stated otherwise (e.g. alpha_d)
**Arguments**
- `af::airfoil_data4D`: airfoil function callable by: CL, CD, CM = af(aoa,Re,mach,family_factor)
- `alpha::Float64`: Static Angle of Attack (at 0.75 chord)
- `adotnorm::Float64`: Normalized Change in Angle of Attack adot*c/(2*U)
- `umach::Float64`: Blade mach number
- `Re::Float64`: Blade Reynolds number
- `aoaStallPos::Float64`: Positive Stall Angle (onset)
- `aoaStallNeg::Float64`: Negative Stall Angle (onset)
- `AOA0::Float64`: Zero Lift AOA
- `tc::Float64`: Thickness to chord ratio
- `BV_DynamicFlagL::Int`: lagged dynamic stall state for lift
- `BV_DynamicFlagD::Int`: lagged dynamic stall state for drag
- `family_factor::float64`: factor indexing airfoil family, if used
"""
function Boeing_Vertol(af,alpha,adotnorm,umach,Re,aoaStallPos,aoaStallNeg,AOA0,tc,BV_DynamicFlagL,BV_DynamicFlagD; family_factor = 0.0)

    # Calculate the static stall envelope limits and other params
    k1pos = 0.5 # Default is 1.0 per CACTUS code
    k1neg = 0.5 # Default is 0.5 per CACTUS code
    diff=0.06-tc
    smachl=0.4+5.0*diff
    hmachl=0.9+2.5*diff
    gammaxl=1.4-6.0*diff
    dgammal=gammaxl/(hmachl-smachl)
    smachm=0.2
    hmachm=0.7+2.5*diff
    gammaxm=1.0-2.5*diff
    dgammam=gammaxm/(hmachm-smachm)

    # Limit reference dalpha to a maximum to keep sign of CL the same for
    # alpha and lagged alpha (considered a reasonable lag...). Note:
    # magnitude increasing and decreasing effect ratios are maintained.
    Fac=.9 # Margin to ensure that dalphaRef is never large enough to make alrefL == AOA0 (blows up linear expansion model)
    dalphaRefMax=Fac*FLOWMath.ksmin([FLOWMath.abs_smooth(aoaStallPos-AOA0,0.0001),FLOWMath.abs_smooth(aoaStallNeg-AOA0,0.0001)],300)/FLOWMath.ksmax([k1pos,k1neg],300)
    TransA=.5*dalphaRefMax # transition region for fairing lagged AOA in pure lag model

    sign_adot=sign(adotnorm)

    # Modified Boeing-Vertol approach

    # Lift
    gammal=gammaxl-(umach-smachl)*dgammal
    dalphaLRef=gammal*sqrt(FLOWMath.abs_smooth(adotnorm,0.0001))
    dalphaLRef=FLOWMath.ksmin([dalphaLRef,dalphaRefMax],300)

    if (adotnorm*(alpha-AOA0) < 0.0)
        # Magnitude of CL decreasing
        dalphaL=k1neg*dalphaLRef
        alrefL=alpha-dalphaL*sign_adot

        # Only switch DS off using lagged alpha
        if (BV_DynamicFlagL == 1 && (alrefL > aoaStallNeg && alrefL < aoaStallPos))
            BV_DynamicFlagL=0
        end

    else
        # Magnitude of CL increasing
        dalphaL=dalphaLRef*k1pos
        alrefL=alpha-dalphaL*sign_adot

        # switch DS on or off using alpha
        if (alpha <= aoaStallNeg || alpha >= aoaStallPos)
            BV_DynamicFlagL=1
        else
            BV_DynamicFlagL=0
        end
    end

    # Drag
    gammam=gammaxm-(umach-smachm)*dgammam
    if umach < smachm
        gammam=gammaxm
    end
    dalphaDRef=gammam*sqrt(FLOWMath.abs_smooth(adotnorm,0.0001))
    dalphaDRef=FLOWMath.ksmin([dalphaDRef,dalphaRefMax],300)

    if adotnorm*(alpha-AOA0) < 0.0
        # Magnitude of CL decreasing
        dalphaD=k1neg*dalphaDRef
        alLagD=alpha-dalphaD*sign_adot

        # Only switch DS off using lagged alpha
        if BV_DynamicFlagD == 1
            delN=aoaStallNeg-alLagD
            delP=alLagD-aoaStallPos
        else
            delN=0.0
            delP=0.0
        end
    else
        # Magnitude of CL increasing
        dalphaD=dalphaDRef*k1pos
        alLagD=alpha-dalphaD*sign_adot

        # switch DS on or off using alpha
        delN=aoaStallNeg-alpha
        delP=alpha-aoaStallPos
    end

    if delN > TransA || delP > TransA
        alrefD=alLagD
        BV_DynamicFlagD=1
    elseif delN > 0 && delN < TransA
        # Transition region (fairing effect...)
        alrefD=alpha+(alLagD-alpha)*delN/TransA
        BV_DynamicFlagD=1
    elseif delP > 0 && delP < TransA
        # Transition region (fairing effect...)
        alrefD=alpha+(alLagD-alpha)*delP/TransA
        BV_DynamicFlagD=1
    else
        BV_DynamicFlagD=0
    end

    # Static characteristics
    if BV_DynamicFlagL == 0 || BV_DynamicFlagD == 0
        CL, CD, CM = af(alpha,Re,umach,family_factor)
    end

    # Static or dynamic model
    if BV_DynamicFlagL == 1
        # Dynamic stall characteristics
        # Linear expansion model for linear region coeffs
        CL, _, CM = af(Main.ForwardDiff.value(alrefL),Main.ForwardDiff.value(Re),Main.ForwardDiff.value(umach),family_factor) #TODO: Verify CM calculation
        CL=CL/(alrefL-AOA0)*(alpha-AOA0)
    end

    if BV_DynamicFlagD == 1
        # Dynamic characteristics
        # Pure lag model for drag
        _, CD, _ = af(alrefD,Re,umach,family_factor)
    end

    return CL, CD, CM, BV_DynamicFlagL, BV_DynamicFlagD
end #function
