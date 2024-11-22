import OWENSAero
import FLOWMath
using Test

path,_ = splitdir(@__FILE__)
# include("$(path)/../src/OWENSAero.jl")

function runfullturb(interpolate)
    # for AeroModel in ["DMS"]#,"AC"]
    AeroModel = "DMS"
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

    if interpolate
        ntheta = 60
    else
        ntheta = 60#450
    end
    filterwindow = 3*ntheta

    ##########################################
    ######## Unsteady Method with ifw ########
    ##########################################
    ifw=false
    
    
    # d633a3bbf6f7330ff0dc1135a4c35ef8c18dacb7 good



    # e691cffc8f02b8860e5040aaf50540fd8c595eca good
    # e55cf31fc656899358d2cae9d6dee37460789858 good
    Nslices = 30
    Vinf = 13.0
    println("Running")
    xrotor = [0.0,10,10.0017885627159,11.1553560041493,12.2884687527063,13.4014193201263,14.4944950133051,15.5679780084628,16.6221454239871,17.6572693919714,18.6736171284654,19.6714510024561,20.6510286035987,21.6126028087127,22.5564218470618,23.4827293644340,24.3917644860386,25.2837618782360,26.1589518091164,27.0175602079435,27.8598087234779,28.6859147811952,29.4960916394141,30.2905484443493,31.0694902841018,31.8331182416024,32.5816294465209,33.3152171261553,34.0340706553126,34.7383756051961,35.4283137913100,36.1040633203956,36.7657986364086,37.4136905655517,38.0479063603738,38.6686097429453,39.2759609471231,39.8701167599149,40.4512305619533,41.0194523670911,41.5749288611271,42.1178034396727,42.6482162451696,43.1663042030671,43.6722010571695,44.1660374041621,44.6479407273238,45.1180354294375,45.5764428649045,46.0232813710714,46.4586662987795,46.8827100421419,47.2955220675583,47.6972089419733,48.0878743603867,48.4676191726225,48.8365414093626,49.1947363074538,49.5422963344927,49.8793112126966,50.2058679420645,50.5220508228369,50.8279414772568,51.1236188706413,51.4091593317660,51.6846365725691,51.9501217071802,52.2056832702783,52.4513872347840,52.6872970288901,52.9134735524361,53.1299751926284,53.3368578391141,53.5341748984079,53.7219773076794,53.9003135479022,54.0692296563697,54.2287692385790,54.3789734794877,54.5198811541461,54.6515286377068,54.7739499148146,54.8871765883801,54.9912378877372,55.0861606761897,55.1719694579451,55.2486863844405,55.3163312600615,55.3749215472542,55.4244723710330,55.4649965228855,55.4965044640745,55.5190043283384,55.5325019239909,55.5370007354208,55.5325019239909,55.5190043283384,55.4965044640745,55.4649965228855,55.4244723710330,55.3749215472542,55.3163312600615,55.2486863844405,55.1719694579451,55.0861606761897,54.9912378877372,54.8871765883801,54.7739499148146,54.6515286377068,54.5198811541461,54.3789734794877,54.2287692385790,54.0692296563697,53.9003135479022,53.7219773076794,53.5341748984079,53.3368578391141,53.1299751926284,52.9134735524361,52.6872970288901,52.4513872347840,52.2056832702783,51.9501217071802,51.6846365725691,51.4091593317660,51.1236188706413,50.8279414772568,50.5220508228369,50.2058679420645,49.8793112126966,49.5422963344927,49.1947363074538,48.8365414093626,48.4676191726225,48.0878743603867,47.6972089419733,47.2955220675583,46.8827100421419,46.4586662987795,46.0232813710714,45.5764428649045,45.1180354294375,44.6479407273238,44.1660374041621,43.6722010571695,43.1663042030671,42.6482162451696,42.1178034396727,41.5749288611271,41.0194523670911,40.4512305619533,39.8701167599149,39.2759609471231,38.6686097429453,38.0479063603738,37.4136905655517,36.7657986364086,36.1040633203956,35.4283137913100,34.7383756051961,34.0340706553126,33.3152171261553,32.5816294465209,31.8331182416024,31.0694902841018,30.2905484443493,29.4960916394141,28.6859147811952,27.8598087234779,27.0175602079435,26.1589518091164,25.2837618782360,24.3917644860386,23.4827293644340,22.5564218470618,21.6126028087127,20.6510286035987,19.6714510024561,18.6736171284654,17.6572693919714,16.6221454239871,15.5679780084628,14.4944950133051,13.4014193201263,12.2884687527063,11.1553560041493,10.0017885627159,8.82746863631162,7.63209307561193,6.41535329580485,5.17693519693007,5,0.0]

    zrotor = [0.0,4.47914708496519,4.48000000000000,5.04000000000000,5.60000000000000,6.16000000000000,6.72000000000000,7.28000000000000,7.84000000000000,8.40000000000000,8.96000000000000,9.52000000000000,10.0800000000000,10.6400000000000,11.2000000000000,11.7600000000000,12.3200000000000,12.8800000000000,13.4400000000000,14,14.5600000000000,15.1200000000000,15.6800000000000,16.2400000000000,16.8000000000000,17.3600000000000,17.9200000000000,18.4800000000000,19.0400000000000,19.6000000000000,20.1600000000000,20.7200000000000,21.2800000000000,21.8400000000000,22.4000000000000,22.9600000000000,23.5200000000000,24.0800000000000,24.6400000000000,25.2000000000000,25.7600000000000,26.3200000000000,26.8800000000000,27.4400000000000,28,28.5600000000000,29.1200000000000,29.6800000000000,30.2400000000000,30.8000000000000,31.3600000000000,31.9200000000000,32.4800000000000,33.0400000000000,33.6000000000000,34.1600000000000,34.7200000000000,35.2800000000000,35.8400000000000,36.4000000000000,36.9600000000000,37.5200000000000,38.0800000000000,38.6400000000000,39.2000000000000,39.7600000000000,40.3200000000000,40.8800000000000,41.4400000000000,42,42.5600000000000,43.1200000000000,43.6800000000000,44.2400000000000,44.8000000000000,45.3600000000000,45.9200000000000,46.4800000000000,47.0400000000000,47.6000000000000,48.1600000000000,48.7200000000000,49.2800000000000,49.8400000000000,50.4000000000000,50.9600000000000,51.5200000000000,52.0800000000000,52.6400000000000,53.2000000000000,53.7600000000000,54.3200000000000,54.8800000000000,55.4400000000000,56,56.5600000000000,57.1200000000000,57.6800000000000,58.2400000000000,58.8000000000000,59.3600000000000,59.9200000000000,60.4800000000000,61.0400000000000,61.6000000000000,62.1600000000000,62.7200000000000,63.2800000000000,63.8400000000000,64.4000000000000,64.9600000000000,65.5200000000000,66.0800000000000,66.6400000000000,67.2000000000000,67.7600000000000,68.3200000000000,68.8800000000000,69.4400000000000,70,70.5600000000000,71.1200000000000,71.6800000000000,72.2400000000000,72.8000000000000,73.3600000000000,73.9200000000000,74.4800000000000,75.0400000000000,75.6000000000000,76.1600000000000,76.7200000000000,77.2800000000000,77.8400000000000,78.4000000000000,78.9600000000000,79.5200000000000,80.0800000000000,80.6400000000000,81.2000000000000,81.7600000000000,82.3200000000000,82.8800000000000,83.4400000000000,84,84.5600000000000,85.1200000000000,85.6800000000000,86.2400000000000,86.8000000000000,87.3600000000000,87.9200000000000,88.4800000000000,89.0400000000000,89.6000000000000,90.1600000000000,90.7200000000000,91.2800000000000,91.8400000000000,92.4000000000000,92.9600000000000,93.5200000000000,94.0800000000000,94.6400000000000,95.2000000000000,95.7600000000000,96.3200000000000,96.8800000000000,97.4400000000000,98,98.5600000000000,99.1200000000000,99.6800000000000,100.240000000000,100.800000000000,101.360000000000,101.920000000000,102.480000000000,103.040000000000,103.600000000000,104.160000000000,104.720000000000,105.280000000000,105.840000000000,106.400000000000,106.960000000000,107.520000000000,108.080000000000,108.640000000000,109.200000000000,109.760000000000,109.838611903775,109.9]

    shapeX_raw = xrotor./150*bladelength
    shapeY_raw = zrotor./150*bladelength
    B = 3
    Radius = maximum(shapeX_raw)
    omega = Vinf/Radius*TSR
    RPM = omega/2/pi*60
    N_Rev = 0.5

    OWENSAero.setupTurb(shapeX_raw,shapeY_raw,B,chord,omega,Vinf;Nslices,ntheta,afname = "$(path)/airfoils/NACA_0015.dat",DynamicStallModel="BV",tau = [1e-5,1e-5])

    if interpolate
        mydt = 2*pi/(ntheta*omega)/3
        mytime = collect(0.0:mydt:N_Rev/RPM*60)
        n_steps = length(mytime)#ntheta*N_Rev#120
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


        for (ii,tnew) in enumerate(mytime)
            # ii = 1
            # tnew = 0.01
            for jj = 1:2
                println(jj)
                println("time $tnew of $(maximum(mytime))")
                azi = omega*tnew

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
                power2_temp = OWENSAero.AdvanceTurbineInterpolate(tnew;azi,alwaysrecalc=true)

                CP[:,ii] = CP_temp[:,end]
                Rp[:,:,ii] = Rp_temp[:,:,end]
                Tp[:,:,ii] = Tp_temp[:,:,end]
                Zp[:,:,ii] = Zp_temp[:,:,end]
                alpha[:,:,ii] = alpha_temp[:,:,end]
                cl[:,ii] = cl_temp[:,end]
                cd_af[:,ii] = cd_af_temp[:,end]
                Vloc[:,ii] = Vloc_temp[1,:,end]
                Re[:,ii] = Re_temp[:,end]
                thetavec[ii] = thetavec_temp[1,end]
                Fx_base[ii] = Fx_base_temp[end]
                Fy_base[ii] = Fy_base_temp[end]
                Fz_base[ii] = Fz_base_temp[end]
                Mx_base[ii] = Mx_base_temp[end]
                My_base[ii] = My_base_temp[end]
                Mz_base[ii] = Mz_base_temp[end]
                power[ii] = power_temp[end]
                power2[ii] = power2_temp[end]
            end
        end
    else
        CP,Rp,Tp,Zp,alpha,cl,cd_af,Vloc,Re,thetavec,nstep,Fx_base,Fy_base,Fz_base,Mx_base,My_base,Mz_base,power,power2 = OWENSAero.advanceTurb(N_Rev/RPM*60)
        thetavec = thetavec[1,:]
        Vloc = Vloc[1,:,:]
    end

    return CP,Rp,Tp,Zp,alpha,cl,cd_af,Vloc,Re,thetavec,Fx_base,Fy_base,Fz_base,Mx_base,My_base,Mz_base,power,power2
end

CP,Rp,Tp,Zp,alpha,cl,cd_af,Vloc,Re,thetavec,Fx_base,Fy_base,Fz_base,Mx_base,My_base,Mz_base,power,power2 = runfullturb(false)
CP_I,Rp_I,Tp_I,Zp_I,alpha_I,cl_I,cd_af_I,Vloc_I,Re_I,thetavec_I,Fx_base_I,Fy_base_I,Fz_base_I,Mx_base_I,My_base_I,Mz_base_I,power_I,power2_I = runfullturb(true)

# # Blade Loads
# import PyPlot
# PyPlot.pygui(true)
# PyPlot.close("all")

# PyPlot.figure()
# PyPlot.plot(thetavec_I,Rp_I[1,1,:],label="Interp")
# PyPlot.plot(thetavec,Rp[1,1,:],".-",label="orig")
# PyPlot.legend()
# PyPlot.ylabel("Rp")

# PyPlot.figure()
# PyPlot.plot(thetavec_I,Tp_I[1,1,:],label="Interp")
# PyPlot.plot(thetavec,Tp[1,1,:],".-",label="orig")
# PyPlot.legend()
# PyPlot.ylabel("Tp")

# PyPlot.figure()
# PyPlot.plot(thetavec_I,Zp_I[1,1,:],label="Interp")
# PyPlot.plot(thetavec,Zp[1,1,:],".-",label="orig")
# PyPlot.legend()
# PyPlot.ylabel("Zp")

# # Base Loads
# PyPlot.figure()
# PyPlot.plot(thetavec_I,Fx_base_I,".-",label="interp")
# PyPlot.plot(thetavec,Fx_base,label="orig")
# PyPlot.legend()
# PyPlot.ylabel("Fx_base")

# PyPlot.figure()
# PyPlot.plot(thetavec_I,Fy_base_I,".-",label="interp")
# PyPlot.plot(thetavec,Fy_base,label="orig")
# PyPlot.legend()
# PyPlot.ylabel("Fy_base")

# PyPlot.figure()
# PyPlot.plot(thetavec_I,Fz_base_I,".-",label="interp")
# PyPlot.plot(thetavec,Fz_base,label="orig")
# PyPlot.legend()
# PyPlot.ylabel("Fz_base")

# PyPlot.figure()
# PyPlot.plot(thetavec_I,Mx_base_I,".-",label="interp")
# PyPlot.plot(thetavec,Mx_base,label="orig")
# PyPlot.legend()
# PyPlot.ylabel("Mx_base")

# PyPlot.figure()
# PyPlot.plot(thetavec_I,My_base_I,".-",label="interp")
# PyPlot.plot(thetavec,My_base,label="orig")
# PyPlot.legend()
# PyPlot.ylabel("My_base")

# PyPlot.figure()
# PyPlot.plot(thetavec_I,Mz_base_I,".-",label="interp")
# PyPlot.plot(thetavec,Mz_base,label="orig")
# PyPlot.legend()
# PyPlot.ylabel("Mz_base")

# Spline the base forces and moments
mysamplethetavec = LinRange(thetavec[3],thetavec[end-1],37)
Fx_sample = FLOWMath.linear(thetavec,Fx_base,mysamplethetavec)
Fy_sample = FLOWMath.linear(thetavec,Fy_base,mysamplethetavec)
Fz_sample = FLOWMath.linear(thetavec,Fz_base,mysamplethetavec)
Mx_sample = FLOWMath.linear(thetavec,Mx_base,mysamplethetavec)
My_sample = FLOWMath.linear(thetavec,My_base,mysamplethetavec)
Mz_sample = FLOWMath.linear(thetavec,Mz_base,mysamplethetavec)
Fx_sample_I = FLOWMath.linear(thetavec_I,Fx_base_I,mysamplethetavec)
Fy_sample_I = FLOWMath.linear(thetavec_I,Fy_base_I,mysamplethetavec)
Fz_sample_I = FLOWMath.linear(thetavec_I,Fz_base_I,mysamplethetavec)
Mx_sample_I = FLOWMath.linear(thetavec_I,Mx_base_I,mysamplethetavec)
My_sample_I = FLOWMath.linear(thetavec_I,My_base_I,mysamplethetavec)
Mz_sample_I = FLOWMath.linear(thetavec_I,Mz_base_I,mysamplethetavec)

# PyPlot.figure()
# PyPlot.plot(LinRange(0,1,37),Fy_sample_I,".")
# PyPlot.plot(LinRange(0,1,37),Fy_sample)

for itest = [1:10;13:length(Fx_sample)]
    # println(itest)
    @test isapprox(Fx_sample_I[itest],Fx_sample[itest];atol=abs(Fx_sample[1]*0.01))
    @test isapprox(Fy_sample_I[itest],Fy_sample[itest];atol=abs(Fy_sample[1]*0.02))
    @test isapprox(Fz_sample_I[itest],Fz_sample[itest];atol=abs(Fz_sample[1]*0.01))
    @test isapprox(Mx_sample_I[itest],Mx_sample[itest];atol=abs(Mx_sample[1]*0.02))
    @test isapprox(My_sample_I[itest],My_sample[itest];atol=abs(My_sample[1]*0.01))
    @test isapprox(Mz_sample_I[itest],Mz_sample[itest];atol=abs(Mz_sample[1]*0.01))
end
