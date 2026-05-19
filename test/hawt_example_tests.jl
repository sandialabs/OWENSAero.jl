using Test

include(joinpath(@__DIR__, "..", "examples", "HAWT", "rigid_ccblade_hawt.jl"))

@testset "rigid CCBlade HAWT example" begin
    result = run_rigid_ccblade_hawt_example()
    steady = result.steady

    @test result.radial_positions == [2.0, 4.0, 6.0]
    @test steady.CP ≈ 0.1774837007705498 atol=1e-14
    @test steady.CT ≈ 0.22007046759243654 atol=1e-14
    @test steady.CQ ≈ 0.10141925758317133 atol=1e-14
    @test steady.thrust ≈ 1327.9868847915095 atol=1e-10
    @test steady.torque ≈ 4284.01010759976 atol=1e-10
    @test steady.power ≈ 8568.02021519952 atol=1e-10
    @test steady.a ≈ [0.12148242946405534, 0.0866902801884321, 0.1057684801932912] atol=1e-14
    @test result.tau1 ≈ 1.1140577308129298 atol=1e-15
    @test result.tau2 ≈ [0.41083720807733964, 0.3399012872582306, 0.22167475255971567] atol=1e-15
    @test result.reduced_induction ≈
          [0.007985381710187242, 0.005698395899075594, 0.006952459635322927] atol=1e-16
    @test result.dynamic_induction ≈
          [0.029801800796282796, 0.024563777147845293, 0.04015039372941627] atol=1e-16
end
