using PowerSystems
using JuMP
using Ipopt
using GLPK

""" OPF model type """
mutable struct OPFmodel
    sys::PowerSystems.System
    mod::JuMP.Model
    # beta::Union{Nothing, JuMP.Containers.DenseAxisArray} # Should convert to sparse some time in the future
    voll::Union{Nothing, JuMP.Containers.DenseAxisArray}
    contingencies::Union{Nothing, Vector{String}}
    prob::Union{Nothing, JuMP.Containers.DenseAxisArray}
end

""" Constructor for OPFmodel """
function opfmodel(sys::System, optimizer, time_limit_sec, voll=nothing, contingencies=nothing, prob=nothing)
    mod = Model(optimizer)
    if GLPK.Optimizer == optimizer 
        set_optimizer_attribute(mod, "msg_lev", GLPK.GLP_MSG_ON)
    end
    set_time_limit_sec(mod, time_limit_sec)
    
    # β = JuMP.Containers.DenseAxisArray{Int64}(undef, get_name.(nodes(sys)), get_name.(branches(sys)))
    # for l in branches(sys)
    #     β[get_name(l.arc.from), get_name(l)] = 1
    #     β[get_name(l.arc.to), get_name(l)] = -1
    # end

    # @assert length(voll) == length(get_nonctrl_generation(sys)) # should be implemented
    # @assert isempty(prob) || length(contingencies) == length(prob)

    return OPFmodel(sys, mod, voll, contingencies, prob)
end

"""Find value of β, β = 1 if from-bus is the bus, β = -1 if to-bus is the bus, 0 else"""
beta(bus::Bus, branch::Branch) = bus == branch.arc.from ? 1 : bus == branch.arc.to ? -1 : 0
beta(bus::String, branch::Branch) = bus == branch.arc.from.name ? 1 : bus == branch.arc.to.name ? -1 : 0
function beta(sys::System, branch::String) 
    l = get_component(ACBranch, sys, branch)
    return [beta(b, l) for b in nodes(sys)]
end

a(line::String,contingency::String) = line != contingency ? 1 : 0

get_generator_cost(gen) = (get_cost ∘ get_variable ∘ get_operation_cost)(gen) |> get_value
get_value(x::Vector{Tuple{Float64, Float64}}) = x[1][2]
get_value(x::Tuple{Float64, Float64}) = x[2]

""" An iterator to a type of power system component """
gens_t(sys::System) = get_components(ThermalGen, sys)
gens_h(sys::System) = get_components(HydroGen, sys)
branches(sys::System) = get_components(ACBranch, sys) # Includes both Line and Phase Shifting Transformer
nodes(sys::System) = get_components(Bus, sys)
demands(sys::System) = get_components(StaticLoad, sys)
renewables(sys::System) = get_components(RenewableGen, sys) # Renewable modelled as negative demand
get_ctrl_generation(sys::System) = Iterators.flatten((gens_t(sys), gens_h(sys))) # An iterator of all controllable generators
get_nonctrl_generation(sys::System) = Iterators.flatten((demands(sys), renewables(sys))) # An iterator of all non-controllable load and generation

# it_name(::Type{T}, mos::Model) where {T <: Component} = get_name.(get_components(T,mod))

make_named_array(value_func, list) = JuMP.Containers.DenseAxisArray(
    [value_func(x) for x in list], get_name.(list) 
)

""" An array of the Value Of Lost Load for the demand and renewables """
make_voll(sys::System) = JuMP.Containers.DenseAxisArray(
        [rand(1000:3000, length(demands(sys))); rand(1:30, length(renewables(sys)))], 
        [get_name.(demands(sys)); get_name.(renewables(sys))]
    )

""" An array of the outage probability of the contingencies """
make_prob(contingencies::Vector{String}) = JuMP.Containers.DenseAxisArray(
        (rand(length(contingencies)).*(0.5-0.02).+0.02)./8760, 
        contingencies
    )

function set_renewable_prod!(system::System, ratio=0.5)
    for g in get_components(RenewableGen, system)
        set_active_power!(g, get_max_active_power(g)*ratio)
    end
end

""" Initialize variables for dc SCOPF """
function init_var_dc_SCOPF!(opfm::OPFmodel, max_shed=1.0)
    p_lim = make_named_array(get_active_power_limits, get_ctrl_generation(opfm.sys))
    demand = make_named_array(get_active_power, get_nonctrl_generation(opfm.sys))
    @variables(opfm.mod, begin
        0 <= pg0[g in get_name.(get_ctrl_generation(opfm.sys))] <= p_lim[g].max     # active power variables for the generators
            # restricted to positive numbers (can be changed to minimum active power) and less than maximum active power
        pf0[l in get_name.(branches(opfm.sys))]      # power flow on branches in base case
        va0[b in get_name.(nodes(opfm.sys))]         # voltage angle at a node in base case
        0 <= ls0[d in get_name.(get_nonctrl_generation(opfm.sys))] <= demand[d] * max_shed  # demand curtailment variables
            # restrict load shedding to max load per node
    end)
    @info "Variables added: pg0, pf0, va0, ls0"
    return opfm
end
""" Initialize variables for dc P-SCOPF """
function init_var_dc_P_SCOPF!(opfm::OPFmodel, max_shed=1.0)
    opfm = init_var_dc_SCOPF!(opfm, max_shed)    
    # power flow on branches in contingencies
    @variable(opfm.mod, pfc[l in get_name.(branches(opfm.sys)), c in opfm.contingencies])
    # voltage angle at a node in contingencies
    @variable(opfm.mod, vac[b in get_name.(nodes(opfm.sys)), c in opfm.contingencies])
    @info "Variables added: pfc, vac"
    return opfm
end
""" Initialize variables for dc PC-SCOPF """
function init_var_dc_PC_SCOPF!(opfm::OPFmodel, max_shed=1.0)
    opfm = init_var_dc_P_SCOPF!(opfm, max_shed)
    @variables(opfm.mod, begin
        # pgc[g in get_name.(get_ctrl_generation(opfm.sys)), c in opfm.contingencies] >= 0    # active power variables for the generators in contingencies 
        pgu[g in get_name.(get_ctrl_generation(opfm.sys)), c in opfm.contingencies] >= 0    # active power variables for the generators in contingencies ramp up 
        pgd[g in get_name.(get_ctrl_generation(opfm.sys)), c in opfm.contingencies] >= 0       # and ramp down
        pfcc[l in get_name.(branches(opfm.sys)), c in opfm.contingencies]         # power flow on branches in in contingencies after corrective actions
        vacc[b in get_name.(nodes(opfm.sys)), c in opfm.contingencies]            # voltage angle at a node in in contingencies after corrective actions
        lsc[d in get_name.(get_nonctrl_generation(opfm.sys)), c in opfm.contingencies] >= 0 # load curtailment variables in in contingencies
        # cbc[l in get_name.(branches(opfm.sys)), c in opfm.contingencies], Bin      # circuit breakers on branches in in contingencies before 
        # cbcc[l in get_name.(branches(opfm.sys)), c in opfm.contingencies], Bin         # and after corrective actions
    end)
    @info "Variables added: pgu, pgd, pfcc, vacc, lsc"
    return opfm
end

""" Objective with base case generation and load shedding """
function add_obj!(opfm::OPFmodel)
    @objective(opfm.mod, Min, 
        sum(opfm.mod[:pg0][get_name(g)] * get_generator_cost(g) for g in gens_t(opfm.sys)) + 
        sum(opfm.mod[:pg0][get_name(g)] * 5 for g in gens_h(opfm.sys)) + 
        sum(opfm.voll[d] * opfm.mod[:ls0][d] for d in get_name.(get_nonctrl_generation(opfm.sys)))
    )
    @info "Objective added: base case generation and load shedding"
    return opfm
end
""" Objective with base case and contingency generation and load shedding """
function add_obj_corrective!(opfm::OPFmodel)
    @objective(opfm.mod, Min, 
        sum(get_generator_cost(g) * (opfm.mod[:pg0][get_name(g)]  + 
            # sum(opfm.prob[c] * (opfm.mod[:pgc][get_name(g),c] - opfm.mod[:pg0][get_name(g)]) for c in contingencies))
            sum(opfm.prob[c] * (opfm.mod[:pgu][get_name(g),c] #=+ opfm.mod[:pgd][get_name(g),c]*0.1=#) for c in opfm.contingencies))
            for g in gens_t(opfm.sys)
        ) + 
        # sum(5 * (opfm.mod[:pg0][get_name(g)] + sum(opfm.prob[c] * (opfm.mod[:pgc][get_name(g),c] - opfm.mod[:pg0][get_name(g)])
        sum(5 * (opfm.mod[:pg0][get_name(g)] + sum(opfm.prob[c] * (opfm.mod[:pgu][get_name(g),c] #=+ opfm.mod[:pgd][get_name(g),c]*0.1=#)
            for c in opfm.contingencies) ) for g in gens_h(opfm.sys)
        ) +
        sum(opfm.voll[d] * (opfm.mod[:ls0][d] + sum(opfm.prob[c] * opfm.mod[:lsc][d,c] for c in opfm.contingencies)) 
            for d in get_name.(get_nonctrl_generation(opfm.sys))
        )
    )
    @info "Objective added: base case and contingency generation and load shedding"
    return opfm
end

""" Set voltage angle at reference bus """
function set_ref_angle!(opfm::OPFmodel)
    slack = find_slack(opfm.sys)
    @constraint(opfm.mod, opfm.mod[:va0][get_name(slack)] == 0)
    @info "Voltage angle at the reference bus set to 0"
    return opfm
end

""" Incerted power at each bus for the base case """
function add_c_bus!(opfm::OPFmodel, slack = nothing)
    @expression(opfm.mod, ctrl[b = get_name.(nodes(opfm.sys))], 
        sum(g.bus.name == b ? opfm.mod[:pg0][get_name(g)] : 0 for g in get_ctrl_generation(opfm.sys)))
    @expression(opfm.mod, nonctrl[b = get_name.(nodes(opfm.sys))], 
        sum((d.bus.name == b ? get_active_power(d) - opfm.mod[:ls0][get_name(d)] : 0 for d in demands(opfm.sys)), init = 0.0) + 
        sum((d.bus.name == b ? -get_active_power(d) + opfm.mod[:ls0][get_name(d)] : 0 for d in renewables(opfm.sys)), init = 0.0))
    @constraint(opfm.mod, inj_p[n = get_name.(nodes(opfm.sys))], ctrl[n] .- sum(beta(n, l) * opfm.mod[:pf0][get_name(l)] for l in branches(opfm.sys)) .== nonctrl[n])
    # @constraint(opfm.mod, va_lim, -π/2 .<= opfm.mod[:va0] .<= π/2) # Not really needed, could be implemented with spesific line angle limits
    slack = isnothing(slack) ? find_slack(opfm.sys) : slack
    @constraint(opfm.mod, opfm.mod[:va0][get_name(slack)] == 0) # Set voltage angle at reference bus
    @info "Constraints for power balance at each bus: Base case"
    return opfm
end
""" Incerted power at each bus for the base case and short-term contingencies """
function add_c_bus_cont!(opfm::OPFmodel, slack = nothing)
    slack = isnothing(slack) ? find_slack(opfm.sys) : slack
    opfm = add_c_bus!(opfm, slack)
    @expression(opfm.mod, 
        n_pfc[n = get_name.(nodes(opfm.sys)), c = opfm.contingencies], 
        sum(beta(n,l) * opfm.mod[:pfc][get_name(l),c] for l in branches(opfm.sys))
    )
    @constraint(opfm.mod, inj_pc[n = get_name.(nodes(opfm.sys)), c = opfm.contingencies], opfm.mod[:ctrl][n] .- n_pfc[n,c] .== opfm.mod[:nonctrl][n])
    # @constraint(opfm.mod, vac_lim, -π/2 .<= opfm.mod[:vac] .<= π/2) # Not really needed
    @constraint(opfm.mod, ref_vac[c = opfm.contingencies], opfm.mod[:vac][get_name(slack),c] == 0) # Set voltage angle at reference bus
    @info "- After contingency, before corrective actions"
    return opfm
end
""" Incerted power at each bus for the base case and short-term and long-term contingencies """
function add_c_bus_ccont!(opfm::OPFmodel, slack = nothing)
    slack = isnothing(slack) ? find_slack(opfm.sys) : slack
    opfm = add_c_bus_cont!(opfm, slack)
    @constraint(opfm.mod, inj_pcc[n = get_name.(nodes(opfm.sys)), c = opfm.contingencies],
        # sum(get_name(get_bus(g)) == n ? pgc[get_name(g),c] : 0 for g in get_ctrl_generation(opfm)) -
        sum(get_name(get_bus(g)) == n ? opfm.mod[:pg0][get_name(g)] + opfm.mod[:pgu][get_name(g),c] - opfm.mod[:pgd][get_name(g),c] : 0 for g in get_ctrl_generation(opfm.sys)) -
        sum(beta(n,l) * opfm.mod[:pfcc][get_name(l),c] for l in branches(opfm.sys)) == 
        sum(get_name(get_bus(d)) == n ? get_active_power(d) - opfm.mod[:lsc][get_name(d),c] : 0 for d in demands(opfm.sys)) +
        sum(get_name(get_bus(d)) == n ? -get_active_power(d) + opfm.mod[:lsc][get_name(d),c] : 0 for d in renewables(opfm.sys))
    )
    # @constraint(opfm.mod, vacc_lim, -π/2 .<= opfm.mod[:vacc] .<= π/2) # Not really needed
    @constraint(opfm.mod, ref_vacc[c = opfm.contingencies], opfm.mod[:vacc][get_name(slack),c] == 0) # Set voltage angle at reference bus
    @info "- After contingency and corrective actions"
    return opfm
end

""" Power flow on branch and branch limits for the base case """
function add_c_branch!(opfm::OPFmodel, branch_rating = nothing)
    branch_rating = isnothing(branch_rating) ? make_named_array(get_rate, branches(opfm.sys)) : branch_rating
    @constraint(opfm.mod, pf0_lim[l = get_name.(branches(opfm.sys))], 
        -branch_rating[l] <= opfm.mod[:pf0][l] <= branch_rating[l]
    )
    x = make_named_array(get_x, branches(opfm.sys))
    @constraint(opfm.mod, pb0[l = get_name.(branches(opfm.sys))],
        opfm.mod[:pf0][l] - sum(beta(opfm.sys,l) .* opfm.mod[:va0]) / x[l] == 0
    )
    @info "Power flow on lines: Base case"
    return opfm
end
""" Power flow on branch and branch limits for the base case and short-term contingency """
function add_c_branch_cont!(opfm::OPFmodel, short_term_limit_multi::Float64 = 1.0, branch_rating = nothing)
    branch_rating = isnothing(branch_rating) ? make_named_array(get_rate, branches(opfm.sys)) : branch_rating
    @constraint(opfm.mod, pfc_lim[l = get_name.(branches(opfm.sys)), c = opfm.contingencies], 
        -branch_rating[l] .* a(l,c) .* short_term_limit_multi .<= opfm.mod[:pfc][l,c] .<= branch_rating[l] .* a(l,c) .* short_term_limit_multi
    )
    opfm = add_c_branch!(opfm, branch_rating)
    x = make_named_array(get_x, branches(opfm.sys))
    @constraint(opfm.mod, pbc[l = get_name.(branches(opfm.sys)), c = opfm.contingencies],
        opfm.mod[:pfc][l,c] .- a(l,c) .* sum(beta(opfm.sys,l) .* opfm.mod[:vac][:,c]) ./ x[l] .== 0
    )
    @info "- After contingency, before corrective actions"
    return opfm
end
""" Power flow on branch and branch limits for the base case and short-term and long-term contingency """
function add_c_branch_ccont!(opfm::OPFmodel, short_term_limit_multi::Float64 = 1.0, branch_rating = nothing)
    branch_rating = isnothing(branch_rating) ? make_named_array(get_rate, branches(opfm.sys)) : branch_rating
    @constraint(opfm.mod, pfcc_lim[l = get_name.(branches(opfm.sys)), c = opfm.contingencies], 
        -branch_rating[l] .* a(l,c) .<= opfm.mod[:pfcc][l,c] .<= branch_rating[l] .* a(l,c)
    )
    opfm = add_c_branch_cont!(opfm, short_term_limit_multi)
    x = make_named_array(get_x, branches(opfm.sys))
    @constraint(opfm.mod, pbcc[l = get_name.(branches(opfm.sys)), c = opfm.contingencies],
        opfm.mod[:pfcc][l,c] .- a(l,c) .* sum(beta(opfm.sys,l) .* opfm.mod[:vacc][:,c]) ./ x[l] .== 0
    )
    @info "- After contingency and corrective actions"
    return opfm
end

""" Restrict active power generation to a minimum level """
function add_min_P_gen!(opfm::OPFmodel)
    for g in get_ctrl_generation(opfm.sys)
        set_lower_bound(opfm.mod.pg0[get_name(g)], get_active_power_limits(g).min) 
    end
    @info "Minimum power generation limits added"
    return opfm
end

""" Restrict active power generation ramp to min and max values """
function add_lim_ramp_P_gen!(opfm::OPFmodel, ramp_minutes)
    p_lim = make_named_array(get_active_power_limits, get_ctrl_generation(opfm.sys))
    for g in get_ctrl_generation(opfm.sys)
        for c in opfm.contingencies
            # @constraint(opfm.mod, 0 <= opfm.mod[:pg0][get_name(g)] <= get_active_power_limits(g).max)
            # @constraint(opfm.mod, 0 <= opfm.mod[:pgc][get_name(g),c] <= get_active_power_limits(g).max)
            # @constraint(opfm.mod, -get_ramp_limits(g).down * ramp_minutes <= 
            #     opfm.mod[:pgc][get_name(g),c] - opfm.mod[:pg0][get_name(g)] <= get_ramp_limits(g).up * ramp_minutes
            # )
            set_upper_bound(opfm.mod[:pgu][get_name(g),c], get_ramp_limits(g).up * ramp_minutes)
            set_upper_bound(opfm.mod[:pgd][get_name(g),c], get_ramp_limits(g).down * ramp_minutes)
        end
    end
    @constraint(opfm.mod, pg_lim[g = get_name.(get_ctrl_generation(opfm.sys)), c = opfm.contingencies], 
        0 <= opfm.mod[:pg0][g] + (opfm.mod[:pgu][g,c] - opfm.mod[:pgd][g,c]) <= p_lim[g].max
    )
    @info "Ramping restrictions on generators added"
    return opfm
end

""" Restrict load shedding to max shed amount of load per node """
function add_lim_load_shed!(opfm::OPFmodel, max_shed::Float64 = 1.0)
    p = make_named_array(get_active_power, demands(opfm.sys))
    @constraint(opfm.mod, load_shed[d in get_name.(demands(opfm.sys)), c in opfm.contingencies], 
        opfm.mod[:ls0][d] + opfm.mod[:lsc][d,c] <= p[d] * max_shed
    )
    @info "Restriction of load shedding to $(max_shed) of load per node"
    return opfm
end
""" Restrict load shedding to renewable production per node """
function add_lim_renewable_shed!(opfm::OPFmodel)
    p = make_named_array(get_active_power, renewables(opfm.sys))
    @constraint(opfm.mod, renew_shed[d = get_name.(renewables(opfm.sys)), c = opfm.contingencies], 
        opfm.mod[:ls0][d] + opfm.mod[:lsc][d,c] <= p[d]
    )
    @info "Restriction of renewable shedding"
    return opfm
end
""" Restrict load shedding to load (or renewable production) per node """
add_lim_nonctrl_shed!(opfm::OPFmodel, max_shed::Float64 = 1.0) = add_lim_load_shed!(opfm, max_shed) |> add_lim_renewable_shed!

""" Add unit commitment """
function add_unit_commit!(opfm::OPFmodel)
    @variable(opfm.mod, u[g in get_name.(get_ctrl_generation(opfm.sys))], Bin)
    delete_lower_bound.(opfm.mod[:pg0])
    delete_upper_bound.(opfm.mod[:pg0])
    p_lim = make_named_array(get_active_power_limits, get_ctrl_generation(opfm.sys))
    @constraint(opfm.mod, ucu, opfm.mod[:pg0] .<= getindex.(p_lim,2) .* u) 
    @constraint(opfm.mod, ucd, opfm.mod[:pg0] .>= getindex.(p_lim,1) .* u) 
    set_objective_function(opfm.mod, objective_function(opfm.mod) + 
        sum(u[get_name(g)] * (get_cost ∘ get_variable ∘ get_operation_cost)(g)[end] for g in get_ctrl_generation(opfm.sys)))
    @info "Unit commitment added to the base case"
    return opfm
end

""" Return the (first) slack bus in the system. """
function find_slack(sys::System)
    for x in nodes(sys)
        x.bustype == BusTypes.REF && return x
    end
    @warn "No slack bus found!"
end

""" Run optimizer to solve the model and check for optimality """
function solve_model!(model::Model)
    optimize!(model)
    if termination_status(model) != MOI.OPTIMAL 
        @warn "Model not optimally solved with status $(termination_status(model))!"
    else
        @info "Model solved in $(solve_time(model)) seconds with an objective value of $(objective_value(model))"
    end
    return model
end

""" Run a SCOPF of a power system """
function scopf(system::System, optimizer; 
        time_limit_sec = 60, 
        voll = nothing, 
        unit_commit::Bool=false,
        max_shed = 0.1,
        ratio = 0.5)
    voll = isnothing(voll) ? make_voll(system) : voll
    set_renewable_prod!(system, ratio)

    opfm = opfmodel(system, optimizer, time_limit_sec, voll)
    opfm = init_var_dc_SCOPF!(opfm, max_shed) |> add_c_bus! |> add_c_branch! |> add_obj!
    if unit_commit
        opfm = add_unit_commit!(opfm)
    end
    return opfm
end

""" Run a P-SCOPF of a power system """
function p_scopf(system::System, optimizer; 
        time_limit_sec = 60, 
        voll = nothing, 
        contingencies = nothing, 
        unit_commit::Bool=false,
        max_shed = 0.1,
        ratio = 0.5,
        short_term_limit_multi = 1.5)
    voll = isnothing(voll) ? make_voll(system) : voll
    contingencies = isnothing(contingencies) ? get_name.(branches(system)) : contingencies
    set_renewable_prod!(system, ratio)

    opfm = opfmodel(system, optimizer, time_limit_sec, voll, contingencies)
    opfm = init_var_dc_P_SCOPF!(opfm, max_shed) |> add_c_bus_cont! |> add_obj!
    opfm = add_c_branch_cont!(opfm, short_term_limit_multi)
    if unit_commit
        opfm = add_unit_commit!(opfm)
    end
    return opfm
end

""" Run a PC-SCOPF of a power system """
function pc_scopf(system::System, optimizer; 
        time_limit_sec = 60,
        voll = nothing, 
        contingencies = nothing, 
        unit_commit::Bool=false,
        max_shed = 0.1,
        ratio = 0.5, 
        short_term_limit_multi = 1.5,
        ramp_minutes = 10,
        prob = nothing)
    voll = isnothing(voll) ? make_voll(system) : voll
    contingencies = isnothing(contingencies) ? get_name.(branches(system)) : contingencies
    prob = isnothing(prob) ? make_prob(contingencies) : prob
    set_renewable_prod!(system, ratio)
    
    opfm = opfmodel(system, optimizer, time_limit_sec, voll, contingencies, prob)
    opfm = init_var_dc_PC_SCOPF!(opfm, max_shed) |> add_c_bus_ccont! |> add_obj_corrective!
    # opfm = init_var_dc_PC_SCOPF!(opfm, max_shed) |> add_c_bus_ccont! |> add_obj!
    opfm = add_lim_nonctrl_shed!(opfm, max_shed)
    opfm = add_c_branch_ccont!(opfm, short_term_limit_multi) 
    opfm = add_lim_ramp_P_gen!(opfm, ramp_minutes)
    if unit_commit
        opfm = add_unit_commit!(opfm)
    end

    return opfm
end

