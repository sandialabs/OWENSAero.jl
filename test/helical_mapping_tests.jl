using Test
using OWENSAero

@testset "helical setup and unsteady scalar diagnostics" begin
    ntheta = 8
    nslices = 3
    blade_count = 2
    rho = 1.225
    mu = 1.7894e-5
    chord = fill(0.15, nslices)
    omega = 2.0
    vinf = 4.0

    blade_z = collect(0.0:1.0:3.0)
    blade_x = ones(length(blade_z))
    blade_y = [0.0, 1.0, 0.0, -1.0]

    OWENSAero.setupTurb(
        blade_x,
        blade_z,
        blade_count,
        chord,
        omega,
        vinf;
        bld_y = blade_y,
        rho,
        mu,
        RPI = false,
        DynamicStallModel = "analytic",
        AeroModel = "DMS",
        Nslices = nslices,
        ntheta,
    )

    @test [OWENSAero.turbslices[i].helical_offset for i = 1:nslices] == [0, 1, 0]

    result = OWENSAero.advanceTurb(0.0; azi = 2pi / ntheta)
    alpha = result[5]
    cl = result[6]
    cd_af = result[7]
    vloc = result[8]
    reynolds = result[9]

    @test size(alpha) == (blade_count, nslices, 1)
    @test size(cl) == (nslices, 1)
    @test alpha isa Array{Float64,3}
    @test cl isa Matrix{Float64}

    # Slice 2 has a +1 helical-bin offset, so its scalar diagnostics should
    # describe blade 1 at the shifted blade index, matching the returned load
    # and angle arrays rather than the unshifted revolution step.
    @test cl[2, 1] ≈ 6.2 * alpha[1, 2, 1] atol = 1e-12
    @test cd_af[2, 1] ≈ 0.008 - 0.003 * cl[2, 1] + 0.01 * cl[2, 1]^2 atol = 1e-12
    @test reynolds[2, 1] ≈ rho * vloc[1, 2, 1] * chord[2] / mu atol = 1e-8

    unshifted_cl = 6.2 * alpha[2, 2, 1]
    @test abs(cl[2, 1] - unshifted_cl) > 1e-4
end
