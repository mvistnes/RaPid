using PowerSystems
using JuMP
using Printf
include("N-1_SCOPF.jl")

function make_first_stage_problem(system::System, optimizer; time_limit_sec = 60, voll = nothing, contingencies = nothing, prob = nothing, unit_commit::Bool=false)
    M = -1000
    voll = isnothing(voll) ? make_voll(system) : voll
    opfm = opfmodel(system, optimizer, time_limit_sec, voll)
    opfm = init_var_dc_SCOPF!(opfm) |> set_ref_angle! |> add_c_bus! |> add_c_branch! |> add_obj!
    if unit_commit
        opfm = add_unit_commit!(opfm)
    end
    @variable(opfm.mod, θ >= M)
    # @objective(opfm.mod, Min, θ)
    set_objective_function(opfm.mod, objective_function(opfm.mod) + θ)
    MOI.set(opfm.mod, MOI.LazyConstraintCallback(), benders_cb)
    solve_model!(opfm.mod)
    return opfm
end

function solve_subproblem(x, contingency)
    opfm = opfmodel(system, optimizer, time_limit_sec, voll)
        # add pgu and pgd here and make pg0 a parameter!!!
    opfm = init_var_dc_SCOPF!(opfm, max_shed) |> add_c_bus! |> add_c_branch! |> add_obj!
    if unit_commit
        opfm = add_unit_commit!(opfm)
    end
    @constraint(opfm, opfm[:pf0][contingency] == 0)
    solve_model!(opfm.mod)
    return (
        obj = objective_value(opfm.mod), 
        # pgu = value.(opfm[:pgu]), 
        # λ = dual.(opfm[:pfcc_lim]), 
        lsc = value.(opfm[:lsc]), 
        π = dual.(opfm[:pfc_lim])
    )
end

MAXIMUM_ITERATIONS = 100
ABSOLUTE_OPTIMALITY_GAP = 1e-6
k=0
""" Callback for Benders decomposition. """
function benders_cb(cb_data, k=0)
    global k += 1
    x_k = callback_value.(cb_data, x)
    θ_k = callback_value(cb_data, θ)
    lower_bound = c_1' * x_k + θ_k
    ret = [solve_subproblem(x_k,c) for c in contingencies]
    upper_bound = c_1' * x_k + c_2' * ret.y
    gap = (upper_bound - lower_bound) / upper_bound
    print_iteration(_k, lower_bound, upper_bound, gap)
    if gap < ABSOLUTE_OPTIMALITY_GAP
        println("Terminating with the optimal solution")
        return
    end
    cut = @build_constraint(θ >= ret.obj + -ret.π' * A_1 * (x .- x_k))
    MOI.submit(model, MOI.LazyConstraint(cb_data), cut)
    return
end


function print_iteration(k, args...)
    f(x) = Printf.@sprintf("%12.4e", x)
    println(lpad(k, 9), " ", join(f.(args), " "))
    return
end