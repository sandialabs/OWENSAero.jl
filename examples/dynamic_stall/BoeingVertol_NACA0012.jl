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
alpha_stall_pos = 10.0 * pi / 180
alpha_stall_neg = -10.0 * pi / 180
alpha_zero_lift = 0.0
thickness_to_chord = 0.12

function static_airfoil(alpha, reynolds, mach, family_factor)
    cl = 2pi * alpha
    cd = 0.01 + 0.02 * alpha^2
    cm = -0.03 * alpha
    return cl, cd, cm
end

function main()
    n_samples = round(Int, tmax / dt) + 1
    alpha = fill(alpha_mean, n_samples + 1)
    cl = zeros(n_samples)
    cd = zeros(n_samples)
    cm = zeros(n_samples)
    dynamic_flag_lift = false
    dynamic_flag_drag = false

    for (i, t) in enumerate(0:dt:tmax)
        alpha[i+1] = alpha_mean + alpha_amp * sin(omega_pitch * t)
        dalpha = alpha[i+1] - alpha[i]
        adot_norm = dalpha / dt * chord / (2.0 * max(u_inf, 0.001))
        cl[i],
        cd[i],
        cm[i],
        dynamic_flag_lift,
        dynamic_flag_drag = OWENSAero.Boeing_Vertol(
            static_airfoil,
            alpha[i+1],
            adot_norm,
            mach,
            reynolds,
            alpha_stall_pos,
            alpha_stall_neg,
            alpha_zero_lift,
            thickness_to_chord,
            dynamic_flag_lift,
            dynamic_flag_drag,
            family_factor = 0.0,
        )
    end

    println("Boeing-Vertol dynamic-stall example")
    println("samples: $(length(cl))")
    println("peak CL: $(maximum(cl))")
    println("peak CD: $(maximum(cd))")
    println("minimum CM: $(minimum(cm))")
end

main()
