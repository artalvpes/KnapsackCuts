using Test

include("toy/run.jl")
include("GeneralizedAssignment/src/run.jl")

@testset "Solving a GAP instance by BC" begin
    sol, cost, total_cuts = run_gap([])
    @test cost == 1931
    @test sol == [
        5, 3, 2, 4, 2, 1, 4, 4, 3, 1, 4, 2, 2, 3, 5, 5, 5, 3, 1, 5, 2, 1, 5, 3, 1, 2,
        5, 1, 2, 3, 1, 3, 3, 5, 3, 1, 1, 3, 4, 1, 5, 5, 1, 4, 3, 2, 2, 5, 5, 4, 3, 1,
        5, 5, 4, 2, 4, 1, 2, 2, 3, 1, 4, 5, 4, 4, 4, 1, 3, 2, 4, 2, 5, 5, 2, 3, 1, 3,
        1, 1, 2, 4, 3, 4, 1, 3, 4, 2, 2, 5, 2, 2, 4, 4, 4, 3, 5, 1, 3, 5
    ]
    @test total_cuts == 111
    run_gap(["-i","GeneralizedAssignment/data/gapD-10-100.txt", "-u", "6348"])
end