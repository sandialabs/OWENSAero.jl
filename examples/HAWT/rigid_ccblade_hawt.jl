using OWENSAero

"""
    run_rigid_ccblade_hawt_example()

Run a compact rigid HAWT BEM example through the OWENSAero CCBlade adapter and
return steady BEM metrics plus one Oye dynamic-inflow state update. The geometry
is intentionally small and analytic so the example is fast enough for CI and
easy to inspect.
"""
function run_rigid_ccblade_hawt_example()
    radial_positions = [2.0, 4.0, 6.0]
    chord = [0.5, 0.45, 0.35]
    twist = deg2rad.([12.0, 8.0, 3.0])
    rotor_speed = 2.0
    inflow_speed = 8.0
    density = 1.225
    tip_radius = 7.0

    airfoil(alpha, reynolds_number, mach) =
        (2 * pi * alpha, 0.01 + 0.02 * alpha^2, -0.05)

    steady = ccbladeHAWTSolve(
        radial_positions,
        chord,
        twist,
        airfoil,
        rotor_speed,
        inflow_speed,
        density;
        num_blades = 3,
        hub_radius = 1.0,
        tip_radius,
        pitch = deg2rad(1.0),
        npts = 8,
    )

    tau1, tau2 = oyeDynamicInflowTimeConstants(
        tip_radius,
        inflow_speed,
        sum(steady.a) / length(steady.a),
        radial_positions,
    )
    reduced_induction, dynamic_induction = oyeDynamicInflowStep(
        zeros(length(radial_positions)),
        zeros(length(radial_positions)),
        steady.a,
        0.2,
        tau1,
        tau2,
    )

    return (
        radial_positions = radial_positions,
        steady = steady,
        tau1 = tau1,
        tau2 = tau2,
        reduced_induction = reduced_induction,
        dynamic_induction = dynamic_induction,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = run_rigid_ccblade_hawt_example()
    steady = result.steady
    println("CP = $(steady.CP)")
    println("CT = $(steady.CT)")
    println("CQ = $(steady.CQ)")
    println("thrust_N = $(steady.thrust)")
    println("torque_Nm = $(steady.torque)")
    println("power_W = $(steady.power)")
    println("axial_induction = $(steady.a)")
    println("oye_tau1_s = $(result.tau1)")
    println("oye_tau2_s = $(result.tau2)")
    println("oye_dynamic_induction_after_0p2s = $(result.dynamic_induction)")
end
