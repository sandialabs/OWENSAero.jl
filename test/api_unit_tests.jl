using Test
using OWENSAero

const API_TEST_DIR, _ = splitdir(@__FILE__)

@testset "public API bindings" begin
    exported = names(OWENSAero)
    @test :AC in exported
    @test :DMS in exported
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

    env = OWENSAero.Environment(1.225, 1.7894e-5, fill(5.0, ntheta), "none", "DMS", zeros(2 * ntheta))
    @test env isa OWENSAero.Environment
    @test env.V_y == zeros(ntheta)
    @test env.V_z == zeros(ntheta)
    @test env.V_twist == zeros(ntheta)
    @test env.windangle == 0.0
    @test env.DynamicStallModel == "none"
    @test env.AeroModel == "DMS"
    @test env.V_wake_old == fill(5.0, ntheta)
    @test env.gravity == [0.0, 0.0, -9.81]

    us = OWENSAero.UnsteadyParams(true, [0.3, 3.0], false)
    @test us isa OWENSAero.UnsteadyParams
    @test us.RPI == true
    @test us.tau == [0.3, 3.0]
    @test us.ifw == false
    @test us.IECgust == false
    @test us.nominalVinf == 1.0
end

@testset "blade azimuth indexing wraps over all bins" begin
    @test [OWENSAero._blade_azimuth_index(6, 2, ibld, 6) for ibld in 1:3] == [6, 2, 4]
    @test [OWENSAero._blade_azimuth_index(1, 2, ibld, 6) for ibld in 1:3] == [1, 3, 5]
    @test [OWENSAero._blade_azimuth_index(30, 10, ibld, 30) for ibld in 1:3] == [30, 10, 20]
end

@testset "pInt periodic integration" begin
    ntheta = 24
    dtheta = 2 * pi / ntheta
    theta = collect(dtheta / 2:dtheta:2 * pi - dtheta / 2)

    @test OWENSAero.pInt(theta, ones(ntheta)) ≈ 2 * pi atol=1e-14
    @test OWENSAero.pInt(theta, fill(3.0, ntheta)) ≈ 6 * pi atol=1e-14
    @test OWENSAero.pInt(theta, sin.(theta)) ≈ 0.0 atol=1e-14
    @test OWENSAero.pInt(theta, cos.(theta)) ≈ 0.0 atol=1e-14
    @test OWENSAero.pInt(theta, sin.(theta) .^ 2) ≈ pi atol=1e-14
end

@testset "added-mass and buoyancy geometry helpers" begin
    chord = 0.2
    thickness = 0.04
    @test OWENSAero.added_mass_flap_volume_per_unit_span(chord) ≈
          0.031415926535897934 atol=1e-16
    @test OWENSAero.added_mass_edge_volume_per_unit_span(thickness) ≈
          1.2566370614359173e-5 atol=1e-20
    @test OWENSAero.buoyancy_section_area_per_unit_span(chord, thickness) == 0.004

    chords = [0.2, 0.4]
    thicknesses = [0.04, 0.08]
    @test OWENSAero.buoyancy_section_area_per_unit_span.(chords, thicknesses) ==
          [0.004, 0.016]
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

    cl_bv, cd_bv = af_bv(5 * pi / 180, 3.0e5, 0.0, env, 0.0, 0.2, 0.1, 5.0; solvestep=true, idx=2)
    @test cl_bv ≈ 0.55 atol=1e-14
    @test cd_bv ≈ 0.0114 atol=1e-14
    @test env.alpha_last[2] == 0.0

    cl_bv_update, cd_bv_update = af_bv(5 * pi / 180, 3.0e5, 0.0, env, 0.0, 0.2, 0.1, 5.0; idx=2)
    @test cl_bv_update ≈ 0.55 atol=1e-14
    @test cd_bv_update ≈ 0.0114 atol=1e-14
    @test env.alpha_last[2] ≈ 5 * pi / 180 atol=1e-14
    @test env.BV_DynamicFlagL[2] == 0
    @test env.BV_DynamicFlagD[2] == 0
end

@testset "airfoil Cm25 propagation" begin
    filename = joinpath(API_TEST_DIR, "airfoils", "cm25_unit.dat")

    af = OWENSAero.readaerodyn(filename)
    cl, cd = af(5 * pi / 180, 3.0e5, 0.0)
    cl_cm, cd_cm, cm = af(5 * pi / 180, 3.0e5, 0.0; return_cm=true)
    @test (cl, cd) == (0.5, 0.02)
    @test (cl_cm, cd_cm, cm) == (0.5, 0.02, 0.04)

    cl_fallback, cd_fallback, cm_fallback =
        OWENSAero._airfoil_coefficients((alpha, Re, mach) -> (2.0 * alpha, 0.01), 0.25, 1.0e5, 0.0)
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
    cl_bv, cd_bv, cm_bv =
        af_bv(5 * pi / 180, 3.0e5, 0.0, env_bv, 0.0, 0.2, 0.1, 5.0; solvestep=true, idx=2, return_cm=true)
    @test cl_bv == 0.5
    @test cd_bv == 0.02
    @test cm_bv == 0.04

    ntheta = 6
    rho = 1.225
    chord = 0.2
    cm25 = 0.125
    af_with_cm(alpha, Re, mach; return_cm=false) = begin
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

    streamtube_result = OWENSAero.streamtube(0.0, pi / ntheta, turbine, env; output_all=true)
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
    theta = collect((2 * pi / ntheta) / 2:(2 * pi / ntheta):2 * pi)
    radialforce_result = OWENSAero.radialforce(zeros(ntheta), zeros(ntheta), theta, turbine_ac, env)
    Vn = radialforce_result[12]
    Vt = radialforce_result[13]
    W = sqrt.(Vn .^ 2 .+ Vt .^ 2)
    @test radialforce_result[end-1] == fill(cm25, ntheta)
    @test radialforce_result[end] ≈ cm25 .* 0.5 .* rho .* W .^ 2 .* chord .^ 2 atol=1e-14
end
