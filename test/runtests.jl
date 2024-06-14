
using Test

println("Start Tests")

include("simple_example.jl")

# include("dyn_stall_tests.jl")

@testset "full turbine" begin
    include("full_turb.jl")
end

@testset "undersampling" begin
    include("full_turb_undersampling.jl")
end
