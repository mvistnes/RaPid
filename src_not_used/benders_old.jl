# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

using PowerSystems
using JuMP
using Printf
using Gurobi
import MathOptInterface
const MOI = MathOptInterface
# include("N-1_SCOPF.jl")
# include("short_long_SCOPF.jl")

const ABSOLUTE_OPTIMALITY_GAP = 1e-6
k::Int64 = 0
global _opfm::OPFmodel
global _unit_commit::Bool
global _max_shed::Float64
global _max_curtail::Float64
global _ratio::Float64
global _circuit_breakers::Bool
global _short_term_limit_multi::Float64
global _ramp_minutes::Int64
global _repair_time

function benders_scopf(system::System, optimizer; 
        voll = nothing, 
        contingencies = nothing,
        prob = nothing, 
        time_limit_sec::Int64 = 600,
        unit_commit::Bool = false,
        max_shed::Float64 = 0.1,
        max_curtail::Float64 = 1.0,
        ratio::Float64= 0.5, 
        circuit_breakers::Bool=false,
        short_term_limit_multi::Float64 = 1.5,
        ramp_minutes::Int64 = 10,
        repair_time::Float64 = 1.0)
    system = set_renewable_prod!(system, ratio)
    voll = isnothing(voll) ? make_voll(system) : voll
    contingencies = isnothing(contingencies) ? get_name.(branches(system)) : contingencies
    prob = isnothing(prob) ? make_prob(contingencies) : prob
    
    # Set global variables
    M = -1000
    global _opfm = opfmodel(system, optimizer, time_limit_sec, voll, contingencies, prob)
    global _unit_commit=unit_commit
    global _max_shed=max_shed
    global _max_curtail=max_curtail
    global _ratio=ratio 
    global _circuit_breakers=circuit_breakers
    global _short_term_limit_multi=short_term_limit_multi
    global _ramp_minutes=ramp_minutes
    global _repair_time=repair_time
    
    _opfm = p_scopf(_opfm, unit_commit=_unit_commit, max_shed=_max_shed, max_curtail=_max_curtail, ratio=_ratio, 
        circuit_breakers=_circuit_breakers, short_term_limit_multi=_short_term_limit_multi, ramp_minutes=_ramp_minutes)
    @variable(_opfm.mod, θ >= M)
    @objective(_opfm.mod, Min, θ)
    # set_objective_function(_opfm.mod, objective_function(_opfm.mod) + θ)
    MOI.set(_opfm.mod, MOI.LazyConstraintCallback(), benders_cb)
    _opfm.mod = solve_model!(_opfm.mod)
    return _opfm
end

""" Callback for Benders decomposition. """
function benders_cb(cb_data)
    global k += 1
    pg0_k = callback_value.(cb_data, _opfm.mod[:pg0])
    ls0_k = callback_value.(cb_data, _opfm.mod[:ls0])
    lsc_k = callback_value.(cb_data, _opfm.mod[:lsc])
    θ_k = callback_value(cb_data, _opfm.mod[:θ])

    # ret = [solve_subproblem(pg0_k, ls0_k, lsc_k, c) for c in _opfm.contingencies]
    ret = solve_subproblem(pg0_k, ls0_k, lsc_k)
    upper_bound = sum(_opfm.cost.data' * pg0_k) + sum(_opfm.voll.data' * ls0_k) + sum(_opfm.prob.data * _opfm.voll.data' * lsc_k) +
        sum(_opfm.prob.data * _opfm.cost.data' * ret.pgu) + sum(_opfm.prob.data * _opfm.voll.data' * ret.lscc)
    gap = 1 - (θ_k / upper_bound)

    print_iteration(k, θ_k, upper_bound, gap)
    if gap < ABSOLUTE_OPTIMALITY_GAP
        println("Terminating with the optimal solution")
        return
    end

    p_lim = [get_active_power_limits(x).max for x in get_ctrl_generation(_opfm.sys)]
    cut = @build_constraint(_opfm.mod[:θ] >= ret.obj + ret.λ' * p_lim * (_opfm.mod[:pg0] .- pg0_k) + ret.π' * _opfm.voll.data * (_opfm.mod[:ls0] .- ls0_k))
    MOI.submit(model, MOI.LazyConstraint(cb_data), cut)
    @info "Adding the cut $(cut)"
    return
end
    
function solve_subproblem(pg0, ls0, lsc, u = nothing)
    opfm = c_scopf(_opfm.sys, Gurobi.Optimizer, _opfm.voll, _opfm.contingencies, _opfm.prob, pg0, ls0, lsc, u, 
        unit_commit=_unit_commit, max_shed=_max_shed, short_term_limit_multi=_short_term_limit_multi, 
        ramp_minutes=_ramp_minutes, repair_time=_repair_time)
    solve_model!(opfm.mod)
    return (
        obj = objective_value(opfm.mod), 
        pgu = value.(opfm[:pgu]), 
        λ = dual.(opfm[:pfcc_lim]),
        lscc = value.(opfm[:lscc]), 
        π = dual.(opfm[:pfcc_lim])
    )
end

function print_iteration(k, args...)
    f(x) = Printf.@sprintf("%12.4e", x)
    println(lpad(k, 9), " ", join(f.(args), " "))
    return
end
