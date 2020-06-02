module KnapsackCuts

using JuMP, CPLEX, PisingerKnapsack

#
# Julia implementation of the Knapsack cut separation proposed in (Avella et al, 2010).
#
# Code by Artur Pessoa (2020)
#
# AVELLA ET AL, 2010: A computational study of exact knapsack separation for the generalized
# assignment problem. Computational Optimization and Applications, 45. pp 543-555.

export separate_knapsack_cut, CutData

struct CutData
    coeffs::Array{Float64}
    indices::Array{Int}
    rhs::Float64
end

function build_initial_solutions(b::Int, a::Array{Int}, _x_::Array{Float64}, nz::Array{Int}
    )::Array{Array{Int}}
    # sort the nonzero items by nonincreasing order of _x_ / a
    ratios = [(_x_[j] / a[j], j) for j in nz]
    sort!(ratios, rev = true)
    ord = map(x -> x[2], ratios)

    # build the first solution greedly
    lastCovered = 0
    load = 0
    s = Int[]
    for k in 1:length(nz)
        load += a[ord[k]]
        if load > b
            break
        end
        push!(s, ord[k])
        lastCovered = k
    end
    sort!(s)
    S = [s]

    # build more solutions until all items are covered
    for l in (length(s) + 1):length(ord)
        if l > lastCovered
            load = a[ord[l]]
            s = [ord[l]]
            lastCovered = l
            for k in 1:length(nz)
                if k != l
                    load += a[ord[k]]
                    if load > b
                        break
                    end
                    push!(s, ord[k])
                    if k > l
                        lastCovered = k
                    end
                end
            end
            sort!(s)
            push!(S, s)
        end
    end

    return S
end

# This function uses an extension of the Euclides' algorithm to find a commom integer
# denominator for an array of real coefficients so the all numerators are also integer,
# respecting the error tolerance for the conversion.
function rationalize(coeffs::Array{Float64}, tolerance::Float64)::Tuple{Array{Int},Int}

    # make a sorted copy of the coefficients with the corresponding indices
    work = sort(vcat(1.0, coeffs), rev = true)

    # keep dividing each coefficient by the immediately smaller one while the second is non-zero
    while work[2] > tolerance
        # divide each coefficient but the smallest one
        for j in 1:(length(work) - 1)
            # check the stopping condition
            if work[j + 1] <= tolerance
                break
            end

            # replace the largest of work[j] and work[j + 1] by the remainder
            work[j] -= floor(work[j] / work[j + 1]) * work[j + 1]

            # treat negative errors
            if work[j] > (work[j + 1] / 2)
                work[j] = work[j + 1] - work[j]
            end
        end

        # sort again the coefficients
        sort!(work, rev = true)
    end

    # Here, the greatest common divisor is work[1]: build the result
    return ([round(Cint, coeff / work[1]) for coeff in coeffs], round(Cint, 1.0 / work[1]))
end

function separate_knapsack_cut(
    n::Int, b::Int, a::Array{Int}, _x_::Array{Float64}, params::Dict{String, Float64}
)::CutData

    # get used parameters
    ϵ = params["separation_point_precision"]
    tolerance = params["coefficient_tolerance"]
    min_violation = params["minimum_cut_violation"]

    # find the nonzero components in the sep point
    nz = [j for j in 1:n if _x_[j] > ϵ]
    a_nz = Cint[a[j] for j in nz]

    # if there are js in nz with a[j] > b, return a cut setting them to zero
    big = [j for j in nz if a[j] > b]
    if !isempty(big)
        big = [j for j in 1:n if a[j] > b]   # lift this cut for null _x_[j]
        return CutData([1.0 for j in big], big, 0.0)
    end

    # build the initial knapsack solutions
    S = build_initial_solutions(b, a, _x_, nz)

    # create the model object
    sep = Model(solver=CplexSolver(CPX_PARAM_SCRIND=0, CPX_PARAM_THREADS=1))

    # set the initial model formulation
    @variable(sep, α[j in nz] >= 0.0)
    @objective(sep, Max, sum(_x_[j] * α[j] for j in nz))
    @constraint(sep, sol[i in 1:length(S)], sum(α[j] for j in S[i]) <= 1.0)

    # keep solving the LP and adding new violated constraints while any
    coeffs = Int[]
    rhs = 0
    converged = false
    prev_S = zeros(Cint, length(nz))
    while !converged
        # solve the separation LP, check the violation and get the solution
        status = solve(sep)
        (status != :Optimal) && error("KnapsackCuts: could not solve the separation LP")
        if getobjectivevalue(sep) < (1.0 + min_violation)
            return CutData(Float64[], Int[], 0.0)
        end
        _α_ = [getvalue(α[j]) for j in nz]
        #println(M)
        #@show _α_
        #@show nz

        # find rational coefficients for the cut
        coeffs, rhs = rationalize(_α_, tolerance)

        # run the separation algorithm for the current (rationalized) solution
        z, S′ = minknap(coeffs, a_nz, Cint(b))
        if z > rhs
            if S′ == prev_S
                @warn "KnapsackCuts: separation failed due to numerical issues!"
                converged = true
            else
                @constraint(sep, sum(α[nz[k]] for k in 1:length(nz) if S′[k] == 1) <= 1.0)
                prev_S .= S′
            end
        else
            converged = true
        end
    end

    # apply the sequential lifting to the cut (except for infeasible items)
    all_coeffs = fill(0, n)
    for k in 1:length(nz)
        all_coeffs[nz[k]] = coeffs[k]
    end
    for j in 1:n
        if (_x_[j] <= ϵ) && (a[j] <= b)
            z, S′ = minknap(all_coeffs, a, b - a[j])
            all_coeffs[j] = rhs - z
        end
    end

    # return the lifted cut
    return CutData([Float64(all_coeffs[j]) for j in 1:n if all_coeffs[j] > 0],
        [j for j in 1:n if all_coeffs[j] > 0],
        Float64(rhs)
    )
end

end # module
