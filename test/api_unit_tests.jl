using Test
using OWENSAero

const API_TEST_DIR, _ = splitdir(@__FILE__)

@testset "public API bindings" begin
    exported = names(OWENSAero)
    @test :AC in exported
    @test :DMS in exported
    @test :cpValidationMetrics in exported
    @test :wholeRevolutionIndexRange in exported
    @test :wholeRevolutionMean in exported
    @test :jointDragForce in exported
    @test :towerShadowVelocity in exported
    @test :liftingStrutForce in exported
    @test :AC_steady ∉ exported
    @test isdefined(OWENSAero, :AC)
    @test !isdefined(OWENSAero, :AC_steady)
end

@testset "public constructors and defaults" begin
    ntheta = 6
    af(alpha, Re, M) = (2.0 * alpha, 0.01 + alpha^2)

    turbine = OWENSAero.Turbine(
        2.0,
        fill(2.0, ntheta),
        fill(0.2, ntheta),
        zeros(ntheta),
        zeros(ntheta),
        fill(3.0, ntheta),
        3,
        af,
        ntheta,
        false,
    )
    @test turbine isa OWENSAero.Turbine
    @test turbine.R == 2.0
    @test turbine.z == 1.0
    @test turbine.thick == 0.18
    @test turbine.B == 3
    @test turbine.ntheta == ntheta
    @test turbine.r_delta_influence == false
    @test turbine.helical_offset == zeros(1)
    @test turbine.rhoA[] == 0

    env = OWENSAero.Environment(
        1.225,
        1.7894e-5,
        fill(5.0, ntheta),
        "none",
        "DMS",
        zeros(2 * ntheta),
    )
    @test env isa OWENSAero.Environment
    @test env.V_y == zeros(ntheta)
    @test env.V_z == zeros(ntheta)
    @test env.V_twist == zeros(ntheta)
    @test env.windangle == 0.0
    @test env.DynamicStallModel == "none"
    @test env.AeroModel == "DMS"
    @test env.V_wake_old == fill(5.0, ntheta)
    @test env.gravity == [0.0, 0.0, -9.81]

    env_mixed_strings = OWENSAero.Environment(
        1.225,
        1.7894e-5,
        fill(5.0, ntheta),
        "none",
        strip(" DMS "),
        zeros(2 * ntheta),
    )
    @test env_mixed_strings.DynamicStallModel == "none"
    @test env_mixed_strings.AeroModel == "DMS"

    env_aliases = OWENSAero.Environment(
        1.225,
        1.7894e-5,
        fill(5.0, ntheta),
        "NONE",
        "ac",
        zeros(2 * ntheta),
    )
    @test env_aliases.DynamicStallModel == "none"
    @test env_aliases.AeroModel == "AC"

    env_bv = OWENSAero.Environment(
        1.225,
        1.7894e-5,
        fill(5.0, ntheta),
        "boeing-vertol",
        "dms",
        zeros(2 * ntheta),
    )
    @test env_bv.DynamicStallModel == "BV"
    @test env_bv.AeroModel == "DMS"

    env_analytic = OWENSAero.Environment(
        1.225,
        1.7894e-5,
        fill(5.0, ntheta),
        "analytic",
        "DMS",
        zeros(2 * ntheta),
    )
    @test env_analytic.DynamicStallModel == "analytic"

    @test_throws ArgumentError OWENSAero.Environment(
        1.225,
        1.7894e-5,
        fill(5.0, ntheta),
        "LB",
        "DMS",
        zeros(2 * ntheta),
    )
    @test_throws ArgumentError OWENSAero.Environment(
        1.225,
        1.7894e-5,
        fill(5.0, ntheta),
        "none",
        "BEM",
        zeros(2 * ntheta),
    )

    us = OWENSAero.UnsteadyParams(true, [0.3, 3.0], false)
    @test us isa OWENSAero.UnsteadyParams
    @test us.RPI == true
    @test us.tau == [0.3, 3.0]
    @test us.ifw == false
    @test us.IECgust == false
    @test us.nominalVinf == 1.0
end

@testset "setupTurb model option normalization" begin
    blade_z = [0.0, 0.5, 1.0]
    blade_x = [1.0, 1.0, 1.0]

    OWENSAero.setupTurb(
        blade_x,
        blade_z,
        2,
        [0.2],
        2.0,
        5.0;
        DynamicStallModel = "NONE",
        AeroModel = "dms",
        Nslices = 1,
        ntheta = 4,
        RPI = false,
    )
    @test OWENSAero.envslices[1].DynamicStallModel == "none"
    @test OWENSAero.envslices[1].AeroModel == "DMS"
    @test OWENSAero.turbslices[1].af(0.0, 1.0e5, 0.0; return_cm = true) isa Tuple{Float64,Float64,Float64}

    @test_throws ArgumentError OWENSAero.setupTurb(
        blade_x,
        blade_z,
        2,
        [0.2],
        2.0,
        5.0;
        DynamicStallModel = "LB",
        AeroModel = "DMS",
        Nslices = 1,
        ntheta = 4,
    )
    @test_throws ArgumentError OWENSAero.setupTurb(
        blade_x,
        blade_z,
        2,
        [0.2],
        2.0,
        5.0;
        DynamicStallModel = "none",
        AeroModel = "BEM",
        Nslices = 1,
        ntheta = 4,
    )
end

@testset "blade azimuth indexing wraps over all bins" begin
    @test [OWENSAero._blade_azimuth_index(6, 2, ibld, 6) for ibld = 1:3] == [6, 2, 4]
    @test [OWENSAero._blade_azimuth_index(1, 2, ibld, 6) for ibld = 1:3] == [1, 3, 5]
    @test [OWENSAero._blade_azimuth_index(30, 10, ibld, 30) for ibld = 1:3] == [30, 10, 20]
end

@testset "helical azimuth mapping helpers" begin
    shape_x = ones(4)
    shape_y = [0.0, 1.0, 0.0, -1.0]
    @test OWENSAero._helical_azimuth_offset_bins(shape_x, shape_y, 8) == [0, 1, 0, -1]
    @test OWENSAero._helical_azimuth_offset_bins([0.0, 1.0], [1.0, 1.0], 12) == [0, 2]

    @test [OWENSAero._helical_blade_index(1, 1, ibld, 2, 8) for ibld = 1:2] ==
          [2, 6]
    @test [OWENSAero._helical_blade_index(1, -1, ibld, 2, 8) for ibld = 1:2] ==
          [8, 4]
    @test OWENSAero._helical_blade_index(9, 1, 1, 2, 8) == 2

    @test_throws ArgumentError OWENSAero._helical_azimuth_offset_bins([1.0], [0.0, 1.0], 8)
    @test_throws ArgumentError OWENSAero._helical_azimuth_offset_bins([1.0], [0.0], 0)
    @test_throws ArgumentError OWENSAero._helical_blade_index(1, 0, 1, 3, 8)
end

@testset "whole-revolution averaging helpers" begin
    azimuth = collect(0.0:(pi/2):(7pi/2))
    values = collect(1.0:length(azimuth))

    one_rev_window = OWENSAero.wholeRevolutionIndexRange(azimuth)
    @test one_rev_window == 4:7
    @test azimuth[first(one_rev_window)] + 2pi == azimuth[last(one_rev_window)+1]
    @test OWENSAero.wholeRevolutionMean(values, azimuth) == 5.5

    two_rev_azimuth = collect(0.0:(pi/2):4pi)
    two_rev_values = collect(1.0:length(two_rev_azimuth))
    @test OWENSAero.wholeRevolutionIndexRange(two_rev_azimuth) == 1:8
    @test OWENSAero.wholeRevolutionIndexRange(two_rev_azimuth; revolutions = 1) == 5:8
    @test OWENSAero.wholeRevolutionMean(two_rev_values, two_rev_azimuth; revolutions = 1) ==
          6.5

    @test OWENSAero.wholeRevolutionIndexRange([0.0, pi/2, pi]; allow_partial = true) == 1:3
    @test_throws ArgumentError OWENSAero.wholeRevolutionIndexRange(Float64[])
    @test_throws ArgumentError OWENSAero.wholeRevolutionIndexRange([0.0, pi, pi/2])
    @test_throws ArgumentError OWENSAero.wholeRevolutionIndexRange([0.0, Inf, 2pi])
    @test_throws ArgumentError OWENSAero.wholeRevolutionIndexRange(azimuth; period = 0.0)
    @test_throws ArgumentError OWENSAero.wholeRevolutionIndexRange(azimuth; revolutions = 0)
    @test_throws ArgumentError OWENSAero.wholeRevolutionIndexRange(
        azimuth;
        revolutions = 1.5,
    )
    @test_throws ArgumentError OWENSAero.wholeRevolutionIndexRange([0.0, pi/2, pi])
    @test_throws ArgumentError OWENSAero.wholeRevolutionMean([1.0, 2.0], azimuth)
end

@testset "pInt periodic integration" begin
    ntheta = 24
    dtheta = 2 * pi / ntheta
    theta = collect((dtheta/2):dtheta:(2*pi-dtheta/2))

    @test OWENSAero.pInt(theta, ones(ntheta)) ≈ 2 * pi atol=1e-14
    @test OWENSAero.pInt(theta, fill(3.0, ntheta)) ≈ 6 * pi atol=1e-14
    @test OWENSAero.pInt(theta, sin.(theta)) ≈ 0.0 atol=1e-14
    @test OWENSAero.pInt(theta, cos.(theta)) ≈ 0.0 atol=1e-14
    @test OWENSAero.pInt(theta, sin.(theta) .^ 2) ≈ pi atol=1e-14
end

@testset "added-mass and buoyancy geometry helpers" begin
    chord = 0.2
    thickness = 0.04
    @test OWENSAero.added_mass_flap_volume_per_unit_span(chord) ≈ 0.031415926535897934 atol=1e-16
    @test OWENSAero.added_mass_edge_volume_per_unit_span(thickness) ≈ 1.2566370614359173e-5 atol=1e-20
    @test OWENSAero.buoyancy_section_area_per_unit_span(chord, thickness) == 0.004

    chords = [0.2, 0.4]
    thicknesses = [0.04, 0.08]
    @test OWENSAero.buoyancy_section_area_per_unit_span.(chords, thicknesses) ==
          [0.004, 0.016]
end

@testset "lumped joint drag helper" begin
    @test OWENSAero.jointDragForce(2.0, [0.0, 0.0, 0.0], 0.5) == [0.0, 0.0, 0.0]
    @test OWENSAero.jointDragForce(2.0, [3.0, 4.0, 0.0], 0.0) == [0.0, 0.0, 0.0]

    force = OWENSAero.jointDragForce(2.0, [3.0, 4.0, 0.0], 0.5)
    @test force == [-7.5, -10.0, -0.0]
    @test sqrt(sum(abs2, force)) == 0.5 * 2.0 * 0.5 * 5.0^2
    @test sum(force .* [3.0, 4.0, 0.0]) < 0.0

    @test OWENSAero.jointDragForce(2.0, [-2.0, 1.0], 0.5) == [sqrt(5.0), -sqrt(5.0) / 2]
    @test_throws ArgumentError OWENSAero.jointDragForce(-1.0, [1.0, 0.0, 0.0], 0.5)
    @test_throws ArgumentError OWENSAero.jointDragForce(1.0, [1.0, 0.0, 0.0], -0.5)
    @test_throws ArgumentError OWENSAero.jointDragForce(1.0, [1.0], 0.5)
    @test_throws ArgumentError OWENSAero.jointDragForce(1.0, [1.0, Inf, 0.0], 0.5)
end

@testset "tower shadow velocity helper" begin
    velocity = [10.0, 0.0]
    @test OWENSAero.towerShadowVelocity(
        velocity,
        [2.0, 0.0];
        tower_radius = 1.0,
        centerline_deficit = 0.2,
    ) == velocity
    @test OWENSAero.towerShadowVelocity(
        velocity,
        [-2.0, 0.0];
        active = true,
        tower_radius = 1.0,
        centerline_deficit = 0.2,
    ) == velocity
    @test OWENSAero.towerShadowVelocity(
        [0.0, 0.0],
        [2.0, 0.0];
        active = true,
        tower_radius = 1.0,
        centerline_deficit = 0.2,
    ) == [0.0, 0.0]

    centerline = OWENSAero.towerShadowVelocity(
        velocity,
        [2.0, 0.0];
        active = true,
        tower_radius = 1.0,
        centerline_deficit = 0.2,
    )
    @test centerline == [8.0, 0.0]

    lateral = OWENSAero.towerShadowVelocity(
        velocity,
        [2.0, 1.0];
        active = true,
        tower_radius = 1.0,
        centerline_deficit = 0.2,
    )
    @test lateral == [10.0 * (1 - 0.2 * exp(-0.5)), 0.0]

    expanded = OWENSAero.towerShadowVelocity(
        [0.0, 12.0, 0.0],
        [0.0, 4.0, 1.0];
        active = true,
        tower_radius = 1.0,
        wake_expansion = 0.5,
        centerline_deficit = 0.3,
    )
    expected_deficit = 0.3 * (1 / 3)^2 * exp(-0.5 / 9)
    @test expanded == [0.0, 12.0 * (1 - expected_deficit), 0.0]

    @test_throws ArgumentError OWENSAero.towerShadowVelocity([1.0], [1.0]; active = true)
    @test_throws ArgumentError OWENSAero.towerShadowVelocity([1.0, 0.0], [1.0, 0.0]; active = "yes")
    @test_throws ArgumentError OWENSAero.towerShadowVelocity([1.0, 0.0], [1.0, 0.0]; tower_radius = -1.0)
    @test_throws ArgumentError OWENSAero.towerShadowVelocity([1.0, 0.0], [1.0, 0.0]; wake_expansion = -0.1)
    @test_throws ArgumentError OWENSAero.towerShadowVelocity([1.0, 0.0], [1.0, 0.0]; centerline_deficit = 1.0)
    @test_throws ArgumentError OWENSAero.towerShadowVelocity([1.0, NaN], [1.0, 0.0])
end

@testset "lifting strut force helper" begin
    @test OWENSAero.liftingStrutForce(
        1.2,
        [0.0, 0.0],
        0.5,
        2.0,
        1.0,
        0.1,
        [0.0, 1.0],
    ) == [0.0, 0.0]
    @test OWENSAero.liftingStrutForce(
        1.2,
        [10.0, 0.0],
        0.0,
        2.0,
        1.0,
        0.1,
        [0.0, 1.0],
    ) == [0.0, 0.0]
    @test OWENSAero.liftingStrutForce(1.2, [10.0, 0.0], 0.5, 0.0, 1.0, 0.1, [0.0, 1.0]) ==
          [0.0, 0.0]
    @test OWENSAero.liftingStrutForce(1.2, [0.0, 0.0], 0.5, 2.0, 1.0, 0.1, [0.0, 0.0]) ==
          [0.0, 0.0]

    force = OWENSAero.liftingStrutForce(1.2, [10.0, 0.0], 0.5, 2.0, 1.0, 0.1, [0.0, 1.0])
    @test force isa Vector{Float64}
    @test force == [-6.0, 60.0]
    @test sum(force .* [10.0, 0.0]) == -60.0
    @test sum(force .* [0.0, 1.0]) == 60.0
    @test sum(force .* [10.0, 0.0]) < 0.0

    @test OWENSAero.liftingStrutForce(1.2, [10.0, 0.0], 0.5, 2.0, 0.0, 0.1, [0.0, 1.0]) ==
          [-6.0, 0.0]
    @test OWENSAero.liftingStrutForce(1.2, [10.0, 0.0], 0.5, 2.0, 1.0, 0.0, [0.0, 1.0]) ==
          [0.0, 60.0]

    negative_lift =
        OWENSAero.liftingStrutForce(1.2, [10.0, 0.0], 0.5, 2.0, -1.0, 0.1, [0.0, 1.0])
    @test negative_lift == [-6.0, -60.0]

    force_3d =
        OWENSAero.liftingStrutForce(2.0, [0.0, 0.0, 5.0], 0.25, 4.0, 0.5, 0.2, [1.0, 0.0, 0.0])
    @test force_3d == [12.5, 0.0, -5.0]

    oblique_2d_velocity = [3.0, 4.0]
    oblique_2d_lift = [-8.0, 6.0]
    oblique_2d_force = OWENSAero.liftingStrutForce(
        2.0,
        oblique_2d_velocity,
        0.5,
        4.0,
        0.7,
        0.2,
        oblique_2d_lift,
    )
    oblique_2d_velocity_hat = oblique_2d_velocity ./ 5.0
    oblique_2d_lift_hat = oblique_2d_lift ./ 10.0
    @test oblique_2d_force ≈ [-34.0, 13.0] atol=1e-14
    @test sum(oblique_2d_force .* oblique_2d_velocity_hat) ≈ -10.0 atol=1e-14
    @test sum(oblique_2d_force .* oblique_2d_lift_hat) ≈ 35.0 atol=1e-14

    oblique_3d_velocity = [1.0, 2.0, 2.0]
    oblique_3d_lift = [0.0, -1.0, 1.0]
    oblique_3d_force = OWENSAero.liftingStrutForce(
        2.0,
        oblique_3d_velocity,
        2.0,
        0.5,
        0.5,
        1 / 3,
        oblique_3d_lift,
    )
    oblique_3d_velocity_hat = oblique_3d_velocity ./ 3.0
    oblique_3d_lift_hat = oblique_3d_lift ./ sqrt(2.0)
    @test sum(oblique_3d_force .* oblique_3d_velocity_hat) ≈ -3.0 atol=1e-14
    @test sum(oblique_3d_force .* oblique_3d_lift_hat) ≈ 4.5 atol=1e-14
    @test oblique_3d_force ≈ 9.0 .* (0.5 .* oblique_3d_lift_hat .- (1 / 3) .* oblique_3d_velocity_hat)

    @test_throws ArgumentError OWENSAero.liftingStrutForce(-1.0, [1.0, 0.0], 1.0, 1.0, 0.0, 0.0, [0.0, 1.0])
    @test_throws ArgumentError OWENSAero.liftingStrutForce(1.0, [1.0, 0.0], -1.0, 1.0, 0.0, 0.0, [0.0, 1.0])
    @test_throws ArgumentError OWENSAero.liftingStrutForce(1.0, [1.0, 0.0], 1.0, -1.0, 0.0, 0.0, [0.0, 1.0])
    @test_throws ArgumentError OWENSAero.liftingStrutForce(1.0, [1.0, 0.0], 1.0, 1.0, NaN, 0.0, [0.0, 1.0])
    @test_throws ArgumentError OWENSAero.liftingStrutForce(1.0, [1.0, 0.0], 1.0, 1.0, 0.0, -0.1, [0.0, 1.0])
    @test_throws ArgumentError OWENSAero.liftingStrutForce(1.0, [1.0, 0.0], 1.0, 1.0, 0.0, 0.0, [0.0, 0.0])
    @test_throws ArgumentError OWENSAero.liftingStrutForce(1.0, [1.0, 0.0], 1.0, 1.0, 0.0, 0.0, [0.0, Inf])
    @test_throws ArgumentError OWENSAero.liftingStrutForce(1.0, [1.0, 0.0], 1.0, 1.0, 0.0, 0.0, [0.0, 1.0, 0.0])
    @test_throws ArgumentError OWENSAero.liftingStrutForce(1.0, [1.0, 0.0], 1.0, 1.0, 0.0, 0.0, [1.0, 0.0])
    @test_throws ArgumentError OWENSAero.liftingStrutForce(1.0, [1.0, Inf], 1.0, 1.0, 0.0, 0.0, [0.0, 1.0])
end

@testset "DMS added-mass force sign convention" begin
    ntheta = 4
    chord = 0.2
    thickness = 0.04
    rho = 1000.0
    added_mass_coeff = 1.0
    accel_flap = 0.3
    accel_edge = 0.4
    af(alpha, Re, M) = (0.0, 0.0, 0.0)

    turbine = OWENSAero.Turbine(
        1.0,
        fill(1.0, ntheta),
        0.5,
        fill(chord, ntheta),
        fill(thickness, ntheta),
        zeros(ntheta),
        zeros(ntheta),
        fill(2.0, ntheta),
        1,
        af,
        ntheta,
        false,
        zeros(ntheta),
        zeros(ntheta),
        zeros(1),
        0.0,
    )
    env = OWENSAero.Environment(
        rho,
        1.0e-3,
        fill(1.0, ntheta),
        zeros(ntheta),
        zeros(ntheta),
        zeros(ntheta),
        0.0,
        "none",
        "DMS",
        true,
        false,
        false,
        added_mass_coeff,
        false,
        zeros(2 * ntheta),
        zeros(Int, 1),
        collect(1:ntheta),
        fill(1.0, ntheta),
        zeros(Int, ntheta),
        zeros(Int, ntheta),
        zeros(ntheta),
        false,
        fill(accel_flap, ntheta),
        fill(accel_edge, ntheta),
        [0.0, 0.0, -9.81],
    )

    result = OWENSAero.streamtube(0.0, pi / 2, turbine, env; output_all = true)
    rp = result[3]
    tp = result[4]
    zp = result[5]
    m_addedmass_np = result[13]
    m_addedmass_tp = result[14]
    f_addedmass_np = result[15]
    f_addedmass_tp = result[16]

    expected_mass_np =
        rho * added_mass_coeff * OWENSAero.added_mass_flap_volume_per_unit_span(chord)
    expected_mass_tp =
        rho * added_mass_coeff * OWENSAero.added_mass_edge_volume_per_unit_span(thickness)

    @test m_addedmass_np ≈ 31.415926535897935 atol=1e-14
    @test m_addedmass_tp ≈ 0.012566370614359173 atol=1e-17
    @test m_addedmass_np ≈ expected_mass_np atol=1e-14
    @test m_addedmass_tp ≈ expected_mass_tp atol=1e-17
    @test f_addedmass_np ≈ 9.42477796076938 atol=1e-14
    @test f_addedmass_tp ≈ 0.005026548245743669 atol=1e-17
    @test rp ≈ -f_addedmass_np atol=1e-14
    @test tp ≈ f_addedmass_tp atol=1e-17
    @test zp == 0.0
end

@testset "AC added-mass force sign convention" begin
    ntheta = 4
    theta = collect(((2*pi/ntheta)/2):(2*pi/ntheta):(2*pi))
    chord = 0.2
    thickness = 0.04
    rho = 1000.0
    added_mass_coeff = 1.0
    accel_flap = 0.3
    accel_edge = 0.4
    af(alpha, Re, M) = (zero(alpha), zero(alpha), zero(alpha))

    turbine = OWENSAero.Turbine(
        1.0,
        fill(1.0, ntheta),
        0.5,
        chord,
        thickness,
        zeros(ntheta),
        zeros(ntheta),
        fill(2.0, ntheta),
        1,
        af,
        ntheta,
        false,
        zeros(ntheta),
        zeros(ntheta),
        zeros(1),
        0.0,
    )
    env = OWENSAero.Environment(
        rho,
        1.0e-3,
        fill(1.0, ntheta),
        zeros(ntheta),
        zeros(ntheta),
        zeros(ntheta),
        0.0,
        "none",
        "AC",
        true,
        false,
        false,
        added_mass_coeff,
        false,
        zeros(2 * ntheta),
        zeros(Int, 1),
        collect(1:ntheta),
        fill(1.0, ntheta),
        zeros(Int, ntheta),
        zeros(Int, ntheta),
        zeros(ntheta),
        false,
        fill(accel_flap, ntheta),
        fill(accel_edge, ntheta),
        [0.0, 0.0, -9.81],
    )

    result = OWENSAero.radialforce(zeros(ntheta), zeros(ntheta), theta, turbine, env)
    rp = result[5]
    tp = result[6]
    zp = result[7]
    m_addedmass_np = result[16]
    m_addedmass_tp = result[17]
    f_addedmass_np = result[18]
    f_addedmass_tp = result[19]

    expected_mass_np =
        rho * added_mass_coeff * OWENSAero.added_mass_flap_volume_per_unit_span(chord)
    expected_mass_tp =
        rho * added_mass_coeff * OWENSAero.added_mass_edge_volume_per_unit_span(thickness)

    @test m_addedmass_np ≈ fill(31.415926535897935, ntheta) atol=1e-14
    @test m_addedmass_tp ≈ fill(0.012566370614359173, ntheta) atol=1e-17
    @test m_addedmass_np ≈ fill(expected_mass_np, ntheta) atol=1e-14
    @test m_addedmass_tp ≈ fill(expected_mass_tp, ntheta) atol=1e-17
    @test f_addedmass_np ≈ fill(9.42477796076938, ntheta) atol=1e-14
    @test f_addedmass_tp ≈ fill(0.005026548245743669, ntheta) atol=1e-17
    @test rp ≈ -f_addedmass_np atol=1e-14
    @test tp ≈ f_addedmass_tp atol=1e-17
    @test zp == zeros(ntheta)
end

@testset "AC negative RPM power-frame convention" begin
    ntheta = 8
    theta = collect(((2*pi/ntheta)/2):(2*pi/ntheta):(2*pi))
    af(alpha, Re, M) = (2.0 .* alpha, 0.01 .+ alpha .^ 2)

    function run_ac(omega)
        turbine = OWENSAero.Turbine(
            1.5,
            fill(1.5, ntheta),
            1.0,
            0.2,
            0.1,
            zeros(ntheta),
            zeros(ntheta),
            fill(omega, ntheta),
            2,
            af,
            ntheta,
            false,
            zeros(ntheta),
            zeros(ntheta),
            zeros(1),
            0.0,
        )
        env = OWENSAero.Environment(
            1.225,
            1.7894e-5,
            fill(5.0, ntheta),
            fill(0.5, ntheta),
            zeros(ntheta),
            zeros(ntheta),
            0.0,
            "none",
            "AC",
            false,
            false,
            false,
            0.0,
            false,
            zeros(2 * ntheta),
            zeros(Int, 1),
            collect(1:ntheta),
            fill(1.0, ntheta),
            zeros(Int, ntheta),
            zeros(Int, ntheta),
            zeros(ntheta),
            false,
            zeros(ntheta),
            zeros(ntheta),
            [0.0, 0.0, -9.81],
        )

        return OWENSAero.radialforce(zeros(ntheta), zeros(ntheta), theta, turbine, env)
    end

    positive = run_ac(2.0)
    negative = run_ac(-2.0)
    expected_torque = [
        0.27967876123505137,
        4.474352791343868,
        8.957294191955247,
        8.361694418259795,
        7.755298442377075,
        9.355425618531791,
        5.854983184917381,
        0.9809290791619916,
    ]
    expected_tp = [
        0.18645250749003425,
        2.982901860895912,
        5.971529461303498,
        5.5744629455065295,
        5.170198961584717,
        6.236950412354528,
        3.903322123278254,
        0.6539527194413277,
    ]

    @test positive[4] ≈ 0.09869472822476016 atol=1e-14
    @test negative[4] ≈ 0.09869472822476016 atol=1e-14
    @test positive[15] ≈ expected_torque atol=1e-14
    @test negative[15] ≈ expected_torque atol=1e-14
    @test positive[6] ≈ expected_tp atol=1e-14
    @test negative[6] ≈ -expected_tp atol=1e-14
end

@testset "DMS and AC no tip-loss finite-span baseline" begin
    ntheta = 8
    theta = collect(((2*pi/ntheta)/2):(2*pi/ntheta):(2*pi))
    af(alpha, Re, M) = (2.0 .* alpha, 0.01 .+ alpha .^ 2)

    function dms_at_slice_z(z)
        turbine = OWENSAero.Turbine(
            1.5,
            fill(1.5, ntheta),
            z,
            fill(0.2, ntheta),
            fill(0.1, ntheta),
            fill(0.05, ntheta),
            fill(2.0, ntheta),
            2,
            af,
            ntheta,
            false,
        )
        env = OWENSAero.Environment(
            1.225,
            1.7894e-5,
            fill(5.0, ntheta),
            fill(0.5, ntheta),
            zeros(ntheta),
            zeros(ntheta),
            0.0,
            "none",
            "DMS",
            zeros(2 * ntheta),
        )
        return OWENSAero.DMS(turbine, env; w = zeros(2 * ntheta), solve = false)
    end

    function ac_at_slice_z(z)
        turbine = OWENSAero.Turbine(
            1.5,
            fill(1.5, ntheta),
            z,
            0.2,
            0.1,
            zeros(ntheta),
            zeros(ntheta),
            fill(2.0, ntheta),
            2,
            af,
            ntheta,
            false,
            zeros(ntheta),
            zeros(ntheta),
            zeros(1),
            0.0,
        )
        env = OWENSAero.Environment(
            1.225,
            1.7894e-5,
            fill(5.0, ntheta),
            fill(0.5, ntheta),
            zeros(ntheta),
            zeros(ntheta),
            0.0,
            "none",
            "AC",
            false,
            false,
            false,
            0.0,
            false,
            zeros(2 * ntheta),
            zeros(Int, 1),
            collect(1:ntheta),
            fill(1.0, ntheta),
            zeros(Int, ntheta),
            zeros(Int, ntheta),
            zeros(ntheta),
            false,
            zeros(ntheta),
            zeros(ntheta),
            [0.0, 0.0, -9.81],
        )
        return OWENSAero.radialforce(zeros(ntheta), zeros(ntheta), theta, turbine, env)
    end

    dms_midspan = dms_at_slice_z(0.0)
    dms_near_tip = dms_at_slice_z(1.0)
    ac_midspan = ac_at_slice_z(0.0)
    ac_near_tip = ac_at_slice_z(1.0)

    expected_dms_torque = [
        0.13216300752723198,
        4.155328939728862,
        8.747381434947037,
        10.518927273475686,
        11.899127897963252,
        10.482424731039607,
        5.94768269035212,
        0.8606922816284484,
    ]
    expected_dms_rp = [
        1.2736810051148109,
        6.8359693596987245,
        6.214005470597589,
        1.0483645079053878,
        0.586123014722692,
        -6.72391144149577,
        -10.101933047419232,
        -6.636463965253183,
    ]
    expected_dms_tp = [
        0.08821892236601551,
        2.7736856779962813,
        5.838884708767798,
        7.021393095350895,
        7.9426782115064585,
        6.997027616556331,
        3.970083363983335,
        0.5745128458760576,
    ]
    expected_ac_rp = [
        2.8550942523884117,
        8.151470249030385,
        7.1681706174268305,
        1.5576140367875473,
        0.025532492405510754,
        -5.968263633754219,
        -8.55357364320575,
        -4.632716960246749,
    ]
    expected_ac_tp = [
        0.18645250749003425,
        2.982901860895912,
        5.971529461303498,
        5.5744629455065295,
        5.170198961584717,
        6.236950412354528,
        3.903322123278254,
        0.6539527194413277,
    ]
    expected_ac_torque = [
        0.27967876123505137,
        4.474352791343868,
        8.957294191955247,
        8.361694418259795,
        7.755298442377075,
        9.355425618531791,
        5.854983184917381,
        0.9809290791619916,
    ]

    @test dms_midspan[1] ≈ 0.1148162791981763 atol=1e-14
    @test dms_midspan[3] ≈ expected_dms_torque atol=1e-14
    @test dms_midspan[4] ≈ expected_dms_rp atol=1e-14
    @test dms_midspan[5] ≈ expected_dms_tp atol=1e-14
    @test dms_near_tip[1] ≈ dms_midspan[1] atol=0.0
    @test dms_near_tip[3] ≈ dms_midspan[3] atol=0.0
    @test dms_near_tip[4] ≈ dms_midspan[4] atol=0.0
    @test dms_near_tip[5] ≈ dms_midspan[5] atol=0.0

    @test ac_midspan[4] ≈ 0.09869472822476016 atol=1e-14
    @test ac_midspan[5] ≈ expected_ac_rp atol=1e-14
    @test ac_midspan[6] ≈ expected_ac_tp atol=1e-14
    @test ac_midspan[15] ≈ expected_ac_torque atol=1e-14
    @test ac_near_tip[4] ≈ ac_midspan[4] atol=0.0
    @test ac_near_tip[5] ≈ ac_midspan[5] atol=0.0
    @test ac_near_tip[6] ≈ ac_midspan[6] atol=0.0
    @test ac_near_tip[15] ≈ ac_midspan[15] atol=0.0
end

@testset "CP validation metrics" begin
    metrics = OWENSAero.cpValidationMetrics(
        [3.0, 1.0, 2.0],
        [9.0, 1.0, 4.0],
        [1.0, 1.5, 2.5],
        [1.25, 2.75, 6.0],
    )

    @test metrics.n == 3
    @test metrics.model_cp_on_reference == [1.0, 2.5, 6.5]
    @test metrics.reference_tsr == [1.0, 1.5, 2.5]
    @test metrics.reference_cp == [1.25, 2.75, 6.0]
    @test metrics.rmse ≈ 0.3535533905932738 atol=1e-16
    @test metrics.mean_bias == 0.0
    @test metrics.mean_abs_error ≈ 1 / 3 atol=1e-16
    @test metrics.max_abs_error == 0.5
    @test metrics.model_peak_tsr == 3.0
    @test metrics.model_peak_cp == 9.0
    @test metrics.reference_peak_tsr == 2.5
    @test metrics.reference_peak_cp == 6.0
    @test metrics.peak_cp_error == 3.0

    @test_throws ArgumentError OWENSAero.cpValidationMetrics([1.0], [1.0], [1.0], [1.0])
    @test_throws ArgumentError OWENSAero.cpValidationMetrics(
        [1.0, 1.0],
        [1.0, 2.0],
        [1.0],
        [1.0],
    )
    @test_throws ArgumentError OWENSAero.cpValidationMetrics(
        [1.0, 2.0],
        [1.0, 2.0],
        [3.0],
        [1.0],
    )
    @test_throws ArgumentError OWENSAero.cpValidationMetrics(
        [1.0, 2.0],
        [1.0, NaN],
        [1.0],
        [1.0],
    )
end

@testset "AeroDyn airfoil readers" begin
    filename = joinpath(API_TEST_DIR, "airfoils", "NACA_0015_RE3E5.dat")

    af = OWENSAero.readaerodyn(filename)
    cl0, cd0 = af(0.0, 3.0e5, 0.0)
    clp5, cdp5 = af(5 * pi / 180, 3.0e5, 0.0)
    clm5, cdm5 = af(-5 * pi / 180, 3.0e5, 0.0)

    @test cl0 ≈ 0.0 atol=1e-14
    @test cd0 ≈ 0.0091 atol=1e-14
    @test clp5 ≈ 0.55 atol=1e-14
    @test cdp5 ≈ 0.0114 atol=1e-14
    @test clm5 ≈ -0.55 atol=1e-14
    @test cdm5 ≈ 0.0114 atol=1e-14

    af_new_static = OWENSAero.readaerodyn_BV_NEW(filename; DynamicStallModel = "NONE")
    cl_new, cd_new, cm_new = af_new_static(5 * pi / 180, 3.0e5, 0.0; return_cm = true)
    @test cl_new ≈ 0.55 atol=1e-14
    @test cd_new ≈ 0.0114 atol=1e-14
    @test cm_new == 0.0
    @test_throws ArgumentError OWENSAero.readaerodyn_BV_NEW(filename; DynamicStallModel = "LB")

    af_bv = OWENSAero.readaerodyn_BV(filename)
    env = OWENSAero.Environment(
        1.225,
        1.7894e-5,
        fill(5.0, 4),
        zeros(4),
        zeros(4),
        zeros(4),
        0.0,
        "BV",
        "DMS",
        zeros(8),
    )

    cl_bv, cd_bv =
        af_bv(5 * pi / 180, 3.0e5, 0.0, env, 0.0, 0.2, 0.1, 5.0; solvestep = true, idx = 2)
    @test cl_bv ≈ 0.55 atol=1e-14
    @test cd_bv ≈ 0.0114 atol=1e-14
    @test env.alpha_last[2] == 0.0

    cl_bv_update, cd_bv_update =
        af_bv(5 * pi / 180, 3.0e5, 0.0, env, 0.0, 0.2, 0.1, 5.0; idx = 2)
    @test cl_bv_update ≈ 0.55 atol=1e-14
    @test cd_bv_update ≈ 0.0114 atol=1e-14
    @test env.alpha_last[2] ≈ 5 * pi / 180 atol=1e-14
    @test env.BV_DynamicFlagL[2] == 0
    @test env.BV_DynamicFlagD[2] == 0

    stall_env = OWENSAero.Environment(
        1.225,
        1.7894e-5,
        fill(5.0, 4),
        zeros(4),
        zeros(4),
        zeros(4),
        0.0,
        "BV",
        "DMS",
        zeros(8),
    )
    cl_stall, cd_stall, cm_stall = af_bv(
        7 * pi / 180,
        3.0e5,
        0.0,
        stall_env,
        0.0,
        0.2,
        0.1,
        5.0;
        idx = 2,
        return_cm = true,
    )
    @test cl_stall ≈ 0.77 atol=1e-14
    @test cd_stall ≈ 0.01260229047028014 atol=1e-16
    @test cm_stall == 0.0
    @test stall_env.BV_DynamicFlagL == [1, 1, 0, 0]
    @test stall_env.BV_DynamicFlagD == [1, 1, 0, 0]
end

@testset "airfoil Cm25 propagation" begin
    filename = joinpath(API_TEST_DIR, "airfoils", "cm25_unit.dat")

    af = OWENSAero.readaerodyn(filename)
    cl, cd = af(5 * pi / 180, 3.0e5, 0.0)
    cl_cm, cd_cm, cm = af(5 * pi / 180, 3.0e5, 0.0; return_cm = true)
    @test (cl, cd) == (0.5, 0.02)
    @test (cl_cm, cd_cm, cm) == (0.5, 0.02, 0.04)

    cl_fallback, cd_fallback, cm_fallback = OWENSAero._airfoil_coefficients(
        (alpha, Re, mach) -> (2.0 * alpha, 0.01),
        0.25,
        1.0e5,
        0.0,
    )
    @test (cl_fallback, cd_fallback, cm_fallback) == (0.5, 0.01, 0.0)

    af_bv = OWENSAero.readaerodyn_BV(filename)
    env_bv = OWENSAero.Environment(
        1.225,
        1.7894e-5,
        fill(5.0, 4),
        zeros(4),
        zeros(4),
        zeros(4),
        0.0,
        "BV",
        "DMS",
        zeros(8),
    )
    cl_bv, cd_bv, cm_bv = af_bv(
        5 * pi / 180,
        3.0e5,
        0.0,
        env_bv,
        0.0,
        0.2,
        0.1,
        5.0;
        solvestep = true,
        idx = 2,
        return_cm = true,
    )
    @test cl_bv == 0.5
    @test cd_bv == 0.02
    @test cm_bv == 0.04

    ntheta = 6
    rho = 1.225
    chord = 0.2
    cm25 = 0.125
    af_with_cm(alpha, Re, mach; return_cm = false) = begin
        cl_local = zero(alpha)
        cd_local = 0.01 .+ zero(alpha)
        cm_local = cm25 .+ zero(alpha)
        return_cm ? (cl_local, cd_local, cm_local) : (cl_local, cd_local)
    end
    turbine = OWENSAero.Turbine(
        2.0,
        fill(2.0, ntheta),
        fill(chord, ntheta),
        zeros(ntheta),
        zeros(ntheta),
        fill(3.0, ntheta),
        3,
        af_with_cm,
        ntheta,
        false,
    )
    env = OWENSAero.Environment(
        rho,
        1.7894e-5,
        fill(5.0, ntheta),
        "none",
        "DMS",
        zeros(2 * ntheta),
    )

    streamtube_result =
        OWENSAero.streamtube(0.0, pi / ntheta, turbine, env; output_all = true)
    Vloc = streamtube_result[6]
    @test streamtube_result[end-1] == cm25
    @test streamtube_result[end] ≈ cm25 * 0.5 * rho * chord * Vloc^2 * chord atol=1e-14

    turbine_ac = OWENSAero.Turbine(
        2.0,
        fill(2.0, ntheta),
        chord,
        zeros(ntheta),
        zeros(ntheta),
        fill(3.0, ntheta),
        3,
        af_with_cm,
        ntheta,
        false,
    )
    theta = collect(((2*pi/ntheta)/2):(2*pi/ntheta):(2*pi))
    radialforce_result =
        OWENSAero.radialforce(zeros(ntheta), zeros(ntheta), theta, turbine_ac, env)
    Vn = radialforce_result[12]
    Vt = radialforce_result[13]
    W = sqrt.(Vn .^ 2 .+ Vt .^ 2)
    @test radialforce_result[end-1] == fill(cm25, ntheta)
    @test radialforce_result[end] ≈ cm25 .* 0.5 .* rho .* W .^ 2 .* chord .^ 2 atol=1e-14
end
