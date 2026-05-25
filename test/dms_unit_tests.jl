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

    out = OWENSAero.streamtube(0.2, pi / 3, turbine, env; output_all = true)

    @test out isa NTuple{19,Any}
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
    result_1 = OWENSAero.DMS(turbine, env; w = zeros(2 * ntheta), solve = false)
    result_2 = OWENSAero.DMS(turbine, env; w = zeros(2 * ntheta), solve = false)

    @test env.V_x == vx_global
    @test env.V_y == vy_global
    @test result_1[1] ≈ 0.08429279232669457 atol=1e-14
    @test result_2[1] ≈ result_1[1] atol=1e-14
    @test result_1[10] == 0.0
    @test result_1[12] ≈
          [0.7610205044672681, 1.629220639954444, -1.2396968194374036, -0.42285018527264395] atol=1e-14
    @test result_1[16] ≈
          [82105.99225005694, 29404.070092675687, 68549.89675104182, 91517.5549025191] atol=1e-8
end

@testset "DMS CP uses signed torque" begin
    ntheta = 8
    af(alpha, Re, M) = (0.0, 0.5)
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
        0.0,
        "none",
        "DMS",
        zeros(2 * ntheta),
    )

    result = OWENSAero.DMS(turbine, env; w = zeros(2 * ntheta), solve = false)
    torque = result[3]

    @test torque ≈ [
        -6.383478423454256,
        -3.636152142701073,
        -0.6225483502083994,
        0.5266952578259413,
        0.5628607215812882,
        -0.17690166794402218,
        -2.850013476873858,
        -5.972358432573441,
    ] atol=1e-14
    @test result[1] ≈ -0.04038508084755988 atol=1e-14
    @test abs(result[1] - 0.04512872592797231) > 0.08
end

@testset "DMS negative RPM power-frame convention" begin
    ntheta = 8
    af(alpha, Re, M) = (2.0 * alpha, 0.01 + alpha^2)

    function run_dms(omega)
        turbine = OWENSAero.Turbine(
            1.5,
            fill(1.5, ntheta),
            1.0,
            fill(0.2, ntheta),
            fill(0.1, ntheta),
            fill(0.05, ntheta),
            fill(omega, ntheta),
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

    positive = run_dms(2.0)
    negative = run_dms(-2.0)
    expected_torque = [
        0.13216300752723198,
        4.155328939728862,
        8.747381434947037,
        10.518927273475686,
        11.899127897963252,
        10.482424731039607,
        5.94768269035212,
        0.8606922816284484,
    ]
    expected_tp = [
        0.08821892236601551,
        2.7736856779962813,
        5.838884708767798,
        7.021393095350895,
        7.9426782115064585,
        6.997027616556331,
        3.970083363983335,
        0.5745128458760576,
    ]

    @test positive[1] ≈ 0.1148162791981763 atol=1e-14
    @test negative[1] ≈ 0.1148162791981763 atol=1e-14
    @test positive[3] ≈ expected_torque atol=1e-14
    @test negative[3] ≈ expected_torque atol=1e-14
    @test positive[5] ≈ expected_tp atol=1e-14
    @test negative[5] ≈ -expected_tp atol=1e-14
end

@testset "DMS wind-angle rotation invariance" begin
    ntheta = 8
    af(alpha, Re, M) = (2.0 * alpha, 0.01 + alpha^2)

    function run_dms(windangle)
        freestream = 5.0
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
            fill(freestream * cos(windangle), ntheta),
            fill(freestream * sin(windangle), ntheta),
            zeros(ntheta),
            zeros(ntheta),
            windangle,
            "none",
            "DMS",
            zeros(2 * ntheta),
        )

        return OWENSAero.DMS(turbine, env; w = zeros(2 * ntheta), solve = false)
    end

    baseline = run_dms(0.0)
    rotated = run_dms(pi / 6)

    @test baseline[1] ≈ 0.09762438364772608 atol=1e-14
    @test rotated[1] ≈ 0.09762438364772608 atol=1e-14
    @test baseline[15][1] ≈ 0.39269908169872414 atol=1e-14
    @test rotated[15][1] ≈ -0.13089969389957468 atol=1e-14
    @test rotated[15][end] ≈ 5.3668874498825625 atol=1e-14
    @test rotated[15] ≈ baseline[15] .- pi / 6 atol=1e-14
    @test rotated[4] ≈ baseline[4] atol=1e-14
    @test rotated[5] ≈ baseline[5] atol=1e-14
    @test rotated[7] ≈ baseline[7] atol=1e-14
    @test rotated[12] ≈ baseline[12] atol=1e-14
end
