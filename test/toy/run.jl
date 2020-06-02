#
# Test routines for the module KnapsackCuts
#
# by Artur Pessoa (2020)
#

# include("../../src/KnapsackCuts.jl")

using KnapsackCuts

const SepParams = Dict(
    "minimum_cut_violation" => 1e-4,
    "separation_point_precision" => 1e-7,
    "coefficient_tolerance" => 1e-7
)

# Note: Instead of using a MIP to convert the Knapsack cut coefficients in integers as in Avella
# et al (2010), this code runs an extension of the Euclides' algorithm for fractional numbers to
# obtain the greatest common fractional divisor. In this context, "coefficient_tolerance" is the
# maximum remainder that is considered as zero and the inverse of "maximum_right-hand_side" is
# the minimum remainder that allows the algorithm to proceed. Any remainder between the two values
# is interpreted as a numerical error, and the algorithm is aborted.

function run_toy(inst::Int)
    # read all the input file to the vector "numbers"
    fname = "data/toy$(inst).txt"
    f = open(fname, "r")
    numbers = map(y -> parse(Float64, y), split(read(f, String)))
    close(f)

    # get the data from "numbers"
    n = Int(numbers[1])
    b = Int(numbers[2])
    a = Int.(numbers[3:(n + 2)])
    _x_ = numbers[(n + 3):(2 * n + 2)]

    # show the data
    @show n
    @show b
    @show a
    @show _x_

    # call the separation algorithm
    return separate_knapsack_cut(n, b, a, _x_, SepParams)
end

println("=====================================================")
cut = run_toy(1)
@show cut
println("-----------------------------------------------------")
cut = run_toy(2)
@show cut
println("-----------------------------------------------------")
cut = run_toy(3)
@show cut
println("=====================================================")
