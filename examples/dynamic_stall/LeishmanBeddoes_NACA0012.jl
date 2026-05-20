using OWENSAero

k = 0.048
reynolds = 3.8e6
mach = 0.3
speed_of_sound = 343.0
u_inf = mach * speed_of_sound
rho = 1.225
mu = 1.7894e-5
chord = reynolds * mu / (rho * u_inf)
omega_pitch = k * 2 * u_inf / chord
dt = 0.002
tmax = 3.0

alpha_mean = 10.0 * pi / 180
alpha_amp = 10.0 * pi / 180
lift_slope = 2pi
alpha_zero_lift = 0.0
cl_critical_positive = 1.0
cl_critical_negative = -1.0
ds = 2.0 * u_inf * dt / chord

function static_airfoil(alpha, reynolds, mach, family_factor)
    cl = 2pi * alpha
    cd = 0.01 + 0.02 * alpha^2
    cm = -0.03 * alpha
    return cl, cd, cm
end

function main()
    n_samples = round(Int, tmax / dt) + 1
    alpha = zeros(n_samples)
    cl = zeros(n_samples)
    cd = zeros(n_samples)
    cm = zeros(n_samples)

    cl_ref_le_last = 0.0
    cl_ref_last = 0.0
    fstat_last = 0.0
    cv_last = 0.0
    dp = 0.0
    df = 0.0
    dcnv = 0.0
    le_separation_state = 0
    s_lev = 0.0

    for (i, t) in enumerate(0:dt:tmax)
        alpha[i] = alpha_mean + alpha_amp * sin(omega_pitch * t)
        cl[i],
        cd[i],
        cm[i],
        cl_ref_le_last,
        cl_ref_last,
        fstat_last,
        cv_last,
        dp,
        df,
        dcnv,
        le_separation_state,
        s_lev = OWENSAero.Leishman_Beddoes(
            static_airfoil,
            alpha[i],
            alpha[i],
            ds,
            mach,
            reynolds,
            lift_slope,
            alpha_zero_lift,
            cl_critical_positive,
            cl_critical_negative,
            cl_ref_le_last,
            cl_ref_last,
            fstat_last,
            cv_last,
            dp,
            df,
            dcnv,
            le_separation_state,
            s_lev,
        )
    end

    println("Leishman-Beddoes dynamic-stall example")
    println("samples: $(length(cl))")
    println("peak CL: $(maximum(cl))")
    println("peak CD: $(maximum(cd))")
    println("minimum CM: $(minimum(cm))")
    println("final leading-edge separation state: $(le_separation_state)")
end

main()
