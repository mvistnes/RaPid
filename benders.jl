using PowerSystems
using JuMP
using Printf
import MathOptInterface
const MOI = MathOptInterface
include("N-1_SCOPF.jl")

function benders_scopf(system::System, optimizer; 
        time_limit_sec::Int64 = 600,
        voll = nothing, 
        contingencies = nothing, 
        unit_commit::Bool = false,
        max_shed::Float64 = 0.1,
        ratio::Float64= 0.5, 
        circuit_breakers::Bool=false,
        short_term_limit_multi::Float64 = 1.5,
        ramp_minutes::Int64 = 10,
        prob = nothing)
    M = -1000
    opfm = p_scopf(system, optimizer, time_limit_sec=time_limit_sec, voll=voll, contingencies=contingencies, 
        unit_commit=unit_commit, max_shed=max_shed, ratio=ratio, circuit_breakers=circuit_breakers,
        short_term_limit_multi=short_term_limit_multi)
    @variable(opfm.mod, θ >= M)
    # @objective(opfm.mod, Min, θ)
    set_objective_function(opfm.mod, objective_function(opfm.mod) + θ)
    
    ABSOLUTE_OPTIMALITY_GAP = 1e-6
    k = 0
    """ Callback for Benders decomposition. """
    function benders_cb(cb_data)
        k += 1
        pg0_k = callback_value.(cb_data, pg0)
        ls0_k = callback_value.(cb_data, ls0)
        θ_k = callback_value(cb_data, θ)
        lower_bound = voll' * pg0_k + θ_k
        ret = [solve_subproblem(pg0_k, ls0_k, c) for c in contingencies]
        upper_bound = voll' * pg0_k + prob * voll' * ret.y
        gap = (upper_bound - lower_bound) / upper_bound
        print_iteration(k, lower_bound, upper_bound, gap)
        if gap < ABSOLUTE_OPTIMALITY_GAP
            println("Terminating with the optimal solution")
            return
        end
        cut = @build_constraint(θ >= ret.obj + -ret.π' * A_1 * (pg0 .- pg0_k))
        MOI.submit(model, MOI.LazyConstraint(cb_data), cut)
        return
    end
    MOI.set(opfm.mod, MOI.LazyConstraintCallback(), benders_cb)
    solve_model!(opfm.mod)
    return opfm
end

function solve_subproblem(pg0, ls0, contingency, u = nothing)
    opfm = c_scopf(system, optimizer, pg0, ls0, u, time_limit_sec=time_limit_sec, voll=voll, 
        contingencies=contingency, unit_commit=unit_commit, max_shed=max_shed, 
        short_term_limit_multi=short_term_limit_multi, ramp_minutes=ramp_minutes, prob=prob)
    solve_model!(opfm.mod)
    return (
        obj = objective_value(opfm.mod), 
        pgu = value.(opfm[:pgu]), 
        λ = dual.(opfm[:pfcc_lim]), 
        lsc = value.(opfm[:lsc]), 
        π = dual.(opfm[:pfcc_lim])
    )
end

function print_iteration(k, args...)
    f(x) = Printf.@sprintf("%12.4e", x)
    println(lpad(k, 9), " ", join(f.(args), " "))
    return
end