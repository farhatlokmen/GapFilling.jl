using GapFilling
using Test

@testset "GapFilling.jl" begin
    # Write your tests here.

    x = 2
    y = 3
    @test GapFilling.addValues(x,y) == 5
end
