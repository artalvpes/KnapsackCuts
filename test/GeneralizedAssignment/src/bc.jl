# include("../../../src/KnapsackCuts.jl")

using KnapsackCuts

using JuMP
using CPLEX
import MathProgBase

const SepParams = Dict(
    "minimum_cut_violation" => 1e-4,
    "separation_point_precision" => 1e-7,
    "coefficient_tolerance" => 1e-7,
)

# Note: Instead of using a MIP to convert the Knapsack cut coefficients in
# integers as in Avella et al (2010), this code runs an extension of the
# Euclides' algorithm for fractional numbers to obtain the greatest common
# fractional divisor. In this context, "coefficient_tolerance" is the
# maximum remainder that is considered as zero and the inverse of
# "maximum_right-hand_side" is the minimum remainder that allows the algorithm
# to proceed. Any remainder between the two values is interpreted as a
# numerical error, and the algorithm is aborted.

# Solve the Generalized Assignment Problem instance in the file named
# filename and return the optimal solution, as a job-indexed array of machines,
# and its cost
function solve_gap(filename::String, upcutoff::Float64)::Tuple{Array{Int},Int}
    # read the instance data
    f = open(filename, "r")
    numbers = map(y -> parse(Int, y), split(read(f, String)))
    close(f)
    m = numbers[1]
    n = numbers[2]
    c = transpose(reshape(numbers[3:(n*m+2)], n, m))
    a = transpose(reshape(numbers[(n*m+3):(2*n*m+2)], n, m))
    b = numbers[(2*n*m+3):(2*n*m+m+2)]

    # set the separation parameters
    minimum_cut_violation = SepParams["minimum_cut_violation"]
    minimum_gap_to_separate = (0.05 / n)

    # create the model object
    M = Model(CPLEX.Optimizer)
    set_optimizer_attribute(M, "CPXPARAM_MIP_Display", 4)
    set_optimizer_attribute(M, "CPXPARAM_ScreenOutput", 1)
    set_optimizer_attribute(M, "CPXPARAM_Threads", 1)
    set_optimizer_attribute(M, "CPXPARAM_MIP_Tolerances_MIPGap", 2e-6)
    set_optimizer_attribute(M, "CPXPARAM_MIP_Tolerances_UpperCutoff", upcutoff)

    # set the model formulation
    @variable(M, x[i = 1:m, j = 1:n], Bin)
    # indices = [(i, j) for i in 1:m, j in 1:n]
    n_vars = num_variables(M)
    @objective(M, Min, sum(c[i, j] * x[i, j] for i in 1:m, j in 1:n))
    @constraint(M, assign[j = 1:n], sum(x[i, j] for i in 1:m) == 1)
    @constraint(M, cap[i = 1:m], sum(a[i, j] * x[i, j] for j in 1:n) <= b[i])
    callback_count = 0

    # define the cut separation callback function
    function separate(cb_data, context_id)
        callback_count += 1
        if context_id != CPX_CALLBACKCONTEXT_RELAXATION
            return
        end
        CPLEX.load_callback_variable_primal(cb_data, context_id)

        # get the optimal solution of the current relaxation
        _x_cplex = Vector{Cdouble}(undef, n_vars)
        obj_node = Ref{Cdouble}()
        obj_p = Ref{Cdouble}()
        # Get relaxed obj
        ret = CPXcallbackgetrelaxationpoint(cb_data, _x_cplex, 0, n_vars - 1, obj_node)
        # Get primal
        ret = CPXcallbackgetincumbent(cb_data, _x_cplex, 0, n_vars - 1, obj_p)

        # check if the current gap is sufficiently large to justify the separation effort
        gap = 1.0 - obj_node[] / min(upcutoff, obj_p[])
        if gap < minimum_gap_to_separate
            return
        end

        sol_MOI = Dict()
        for i in 1:m, j in 1:n
            sol_MOI[i, j] = callback_value(cb_data, x[i, j])
        end

        # separate one exact knapsack cut for each machine
        num_cuts = 0
        sum_viol = 0.0
        for i in 1:m
            cut = separate_knapsack_cut(
                n,
                b[i],
                [a[i, j] for j in 1:n],
                [sol_MOI[i, j] for j in 1:n],
                SepParams,
            )
            violation = 0.0
            if !isempty(cut.coeffs)
                violation =
                    sum(
                        cut.coeffs[k] * sol_MOI[i, cut.indices[k]]
                        for k in 1:length(cut.coeffs)
                    ) - cut.rhs
            end

            if violation >= minimum_cut_violation
                con = @build_constraint(
                    sum(
                        cut.coeffs[k] * x[i, cut.indices[k]] for k in 1:length(cut.coeffs)
                    ) <= cut.rhs
                )
                MOI.submit(M, MOI.UserCut(cb_data), con)
                sum_viol += violation / cut.rhs
                num_cuts += 1
            end
        end

        # TODO: Find a found a way of doing this with the current API
        # if num_cuts > 0
        #     println("Found $num_cuts cuts with avg violation $(sum_viol / num_cuts)")

        #     # tell CPLEX to keep calling this routine in this node
        #     unsafe_store!(cb_data.userinteraction_p, convert(Cint, 2), 1)
        # end
    end

    # set the callback cut separation function and the output files, and solve
    MOI.set(M, CPLEX.CallbackFunction(), separate)
    # writeLP(M, "gap.lp", genericnames = false)
    # View solver output
    unset_silent(M)
    optimize!(M)
    status = termination_status(M)

    # get the solution and return its objective value
    sol = zeros(Int, n)
    cost = 0
    if status == MOI.OPTIMAL
        _x_ = value.(x)
        cost = round(Int, objective_value(M))
        for i in 1:m, j in 1:n
            if _x_[i, j] > 0.5
                sol[j] = i
            end
        end
    end
    return (sol, cost)
end
