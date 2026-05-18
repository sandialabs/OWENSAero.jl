using Test
using OWENSAero

@testset "DMS scalar streamtube" begin
    af(alpha, Re, M) = (2.0 * alpha, 0.01 + alpha^2)
    turbine = OWENSAero.Turbine(1.5, 1.5, 0.2, 0.1, 0.05, 2.0, 2, af, 4, false)
    env = OWENSAero.Environment(
        1.225,
        1.7894e-5,
        [5.0],
        [0.5],
        [0.25],
        [0.0],
        pi / 7,
        "none",
        "DMS",
        zeros(8),
    )

    out = OWENSAero.streamtube(0.2, pi / 3, turbine, env; output_all=true)

    @test out isa NTuple{19, Any}
    @test out[1] ≈ 3.120186136285371 atol=1e-12
    @test out[2] ≈ 2.100708299227373 atol=1e-12
    @test out[6] ≈ 6.53305459604201 atol=1e-12
    @test out[8] ≈ 0.04943576554896863 atol=1e-12
    @test out[9] ≈ 0.431135250046925 atol=1e-12
    @test out[10] ≈ 0.86227050009385 atol=1e-12
    @test out[12] ≈ 89448.88655584514 atol=1e-8
    @test out[17] == [0.0, 0.0, 0.0]
    @test out[18] == 0.0
    @test out[19] == 0.0
end

@testset "DMS preserves caller wind frame" begin
    ntheta = 4
    af(alpha, Re, M) = (2.0 * alpha, 0.01 + alpha^2)
    turbine = OWENSAero.Turbine(
        1.5,
        fill(1.5, ntheta),
        1.0,
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
        pi / 7,
        "none",
        "DMS",
        zeros(2 * ntheta),
    )

    vx_global = copy(env.V_x)
    vy_global = copy(env.V_y)
    result_1 = OWENSAero.DMS(turbine, env; w=zeros(2 * ntheta), solve=false)
    result_2 = OWENSAero.DMS(turbine, env; w=zeros(2 * ntheta), solve=false)

    @test env.V_x == vx_global
    @test env.V_y == vy_global
    @test result_1[1] ≈ 0.08429279232669457 atol=1e-14
    @test result_2[1] ≈ result_1[1] atol=1e-14
    @test result_1[10] == 0.0
    @test result_1[12] ≈ [0.7610205044672681, 1.629220639954444, -1.2396968194374036, -0.42285018527264395] atol=1e-14
    @test result_1[16] ≈ [82105.99225005694, 29404.070092675687, 68549.89675104182, 91517.5549025191] atol=1e-8
end
