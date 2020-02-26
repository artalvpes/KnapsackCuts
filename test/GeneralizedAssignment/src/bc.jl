# include("../../../src/KnapsackCuts.jl")

using KnapsackCuts

using JuMP
using CPLEX
import MathProgBase

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

# Solve the Generalized Assignment Problem instance in the file named
# filename and return the optimal solution, as a job-indexed array of machines, and its cost
function solve_gap(filename::String, upcutoff::Float64)::Tuple{Array{Int}, Int}
    # read the instance data
    f = open(filename, "r")
    numbers = map(y -> parse(Int, y), split(read(f, String)))
    close(f)
    m = numbers[1]
    n = numbers[2]
    c = transpose(reshape(numbers[3:(n * m + 2)], n, m))
    a = transpose(reshape(numbers[(n * m + 3):(2 * n * m + 2)], n, m))
    b = numbers[(2 * n * m + 3):(2 * n * m + m + 2)]

    # set the separation parameters
    minimum_cut_violation = SepParams["minimum_cut_violation"]
    minimum_gap_to_separate = (0.05 / n)

    # create the model object
    M = Model(solver=CplexSolver(CPX_PARAM_MIPDISPLAY=4,
            CPX_PARAM_SCRIND=0, CPX_PARAM_THREADS=1,
            CPX_PARAM_EPGAP=2e-6, CPX_PARAM_CUTUP=upcutoff))

    # set the model formulation
    @variable(M, x[i = 0:m, j = 1:n], Bin)
    @objective(M, Min, sum(c[i, j] * x[i, j] for i = 1:m, j = 1:n))
    @constraint(M, assign[j = 1:n], sum(x[i, j] for i = 1:m) == 1)
    @constraint(M, cap[i = 1:m],  sum(a[i, j] * x[i, j] for j = 1:n) <= b[i])

    # define the cut separation callback function
    function separate(cb)
        # get the optimal solution of the current relaxation
        _x_ = getvalue(x)
        # @show _x_

        # check if the current gap is sufficiently large to justify the separation effort
        # @show cbgetnodeobjval(cb), MathProgBase.cbgetobj(cb)
        gap = 1.0 - cbgetnodeobjval(cb) / min(upcutoff, MathProgBase.cbgetobj(cb))
        if gap < minimum_gap_to_separate
            return
        end

        # separate one exact knapsack cut for each machine
        num_cuts = 0
        sum_viol = 0.0
        for i in 1:m
            cut = separate_knapsack_cut(n, b[i], [a[i, j] for j in 1:n],
                [_x_[i, j] for j in 1:n], SepParams
            )
            violation = 0.0
            if !isempty(cut.coeffs)
                violation = sum(cut.coeffs[k] * _x_[i, cut.indices[k]]
                    for k in 1:length(cut.coeffs)) - cut.rhs
            end

            if violation >= minimum_cut_violation
                @usercut(cb, sum(cut.coeffs[k] * x[i, cut.indices[k]]
                    for k in 1:length(cut.coeffs)) <= cut.rhs
                )
                sum_viol += violation / cut.rhs
                num_cuts += 1
            end
        end

        if num_cuts > 0
            # println("Found $num_cuts cuts with avg violation $(sum_viol / num_cuts)")

            # tell CPLEX to keep calling this routine in this node
            unsafe_store!(cb.userinteraction_p, convert(Cint,2), 1)
        end
    end

    # set the callback cut separation function and the output files, and solve
    addcutcallback(M, separate)
    writeLP(M, "gap.lp", genericnames=false)
    JuMP.build(M)
    cpxM = getrawsolver(M)
    CPLEX.set_logfile(cpxM.env, "gap.log")
    status = solve(M)

    # get the solution and return its objective value
    sol = zeros(Int, n)
    cost = 0
    if status == :Optimal
        _x_ = getvalue(x)
        cost = round(Int, getobjectivevalue(M))
        for i = 1:m, j in 1:n
            if _x_[i, j] > 0.5
                sol[j] = i
            end
        end
    end
    return (sol, cost)
end
