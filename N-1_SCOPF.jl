# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
module SCOPF

using PowerSystems
using JuMP
using Ipopt
using GLPK

@enum OPF SC=0 PSC=1 PCSC=2 # SCOPF, P-SCOPF, SC-SCOPF

""" OPF model type """
mutable struct OPFmodel
    sys::PowerSystems.System
    mod::JuMP.Model
    # beta::Union{Nothing, JuMP.Containers.DenseAxisArray} # Should convert to sparse some time in the future
    voll::JuMP.Containers.DenseAxisArray
    cost::JuMP.Containers.DenseAxisArray
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
    
    cost = JuMP.Containers.DenseAxisArray(
        [[get_generator_cost(g)[2] for g in gens_t(sys)]; [5 for _ in gens_h(sys)]],
        get_name.(get_ctrl_generation(sys))
    )

    # β = JuMP.Containers.DenseAxisArray{Int64}(undef, get_name.(nodes(sys)), get_name.(branches(sys)))
    # for l in branches(sys)
    #     β[get_name(l.arc.from), get_name(l)] = 1
    #     β[get_name(l.arc.to), get_name(l)] = -1
    # end

    # @assert length(voll) == length(get_nonctrl_generation(sys)) # should be implemented
    # @assert isempty(prob) || length(contingencies) == length(prob)

    return OPFmodel(sys, mod, voll, cost, contingencies, prob)
end

"""Find value of β, β = 1 if from-bus is the bus, β = -1 if to-bus is the bus, 0 else"""
beta(bus::Bus, branch::Branch) = bus == branch.arc.from ? 1 : bus == branch.arc.to ? -1 : 0
beta(bus::String, branch::Branch) = bus == branch.arc.from.name ? 1 : bus == branch.arc.to.name ? -1 : 0
function beta(sys::System, branch::String) 
    l = get_component(ACBranch, sys, branch)
    return [beta(b, l) for b in nodes(sys)]
end

a(line::String,contingency::String) = line != contingency ? 1 : 0

get_generator_cost(gen) = get_operation_cost(gen) |> get_variable |> get_cost |> get_value
get_value(x::Vector{Tuple{Float64, Float64}}) = x[1]
get_value(x::Tuple{Float64, Float64}) = x

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

""" Make a DenseAxisArray using the list and function for the value of each element """
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

""" Set the renewable production to a ratio of maximum active power """
function set_renewable_prod!(system::System, ratio::Float64=0.5)
    for g in get_components(RenewableGen, system)
        set_active_power!(g, get_max_active_power(g)*ratio)
    end
    return system
end

""" Initialize variables for dc SCOPF """
function init_var_dc_SCOPF!(opfm::OPFmodel)
    p_lim = make_named_array(get_active_power_limits, get_ctrl_generation(opfm.sys))
    @variables(opfm.mod, begin
        0 <= pg0[g in get_name.(get_ctrl_generation(opfm.sys))] <= p_lim[g].max     # active power variables for the generators
            # restricted to positive numbers (can be changed to minimum active power) and less than maximum active power
        pf0[l in get_name.(branches(opfm.sys))]      # power flow on branches in base case
        va0[b in get_name.(nodes(opfm.sys))]         # voltage angle at a node in base case
        0 <= ls0[d in get_name.(get_nonctrl_generation(opfm.sys))]  # demand curtailment variables
    end)
    @info "Variables added: pg0, pf0, va0, ls0"
    return opfm
end
""" Initialize variables for dc P-SCOPF """
function init_var_dc_P_SCOPF!(opfm::OPFmodel)
    opfm = init_var_dc_SCOPF!(opfm)    
    @variables(opfm.mod, begin
        pfc[l in get_name.(branches(opfm.sys)), c in opfm.contingencies] # power flow on branches in contingencies
        vac[b in get_name.(nodes(opfm.sys)), c in opfm.contingencies] # voltage angle at a node in contingencies
        0 <= lsc[d in get_name.(get_nonctrl_generation(opfm.sys)), c in opfm.contingencies] # demand curtailment variables
    end)
    @info "Variables added: pfc, vac, lsc"
    return opfm
end
""" Initialize variables for dc PC-SCOPF """
function init_var_dc_PC_SCOPF!(opfm::OPFmodel)
    opfm = init_var_dc_P_SCOPF!(opfm)
    @variables(opfm.mod, begin
        0 <= pgu[g in get_name.(get_ctrl_generation(opfm.sys)), c in opfm.contingencies]    # active power variables for the generators in contingencies ramp up 
        0 <= pgd[g in get_name.(get_ctrl_generation(opfm.sys)), c in opfm.contingencies]       # and ramp down
        pfcc[l in get_name.(branches(opfm.sys)), c in opfm.contingencies]         # power flow on branches in in contingencies after corrective actions
        vacc[b in get_name.(nodes(opfm.sys)), c in opfm.contingencies]            # voltage angle at a node in in contingencies after corrective actions
        0 <= lscc[d in get_name.(get_nonctrl_generation(opfm.sys)), c in opfm.contingencies] # load curtailment variables in in contingencies
    end)
    @info "Variables added: pgu, pgd, pfcc, vacc, lscc"
    return opfm
end

""" Objective with base case generation and load shedding """
function add_obj!(opfm::OPFmodel)
    @objective(opfm.mod, Min, opfm.cost.data' * opfm.mod[:pg0] + opfm.voll.data' * opfm.mod[:ls0])
    @info "Objective added: base case generation and load shedding"
    return opfm
end
""" Objective with base case and contingency generation and load shedding """
function add_obj_cont!(opfm::OPFmodel, ramp_minutes)
    opfm = add_obj!(opfm)
    set_objective_function(opfm.mod, objective_function(opfm.mod) +
        sum(opfm.voll[d] * sum(opfm.prob[c] * opfm.mod[:lsc][d,c] * ramp_minutes / 60 for c in opfm.contingencies)
            for d in get_name.(get_nonctrl_generation(opfm.sys))
        )
    )
    @info "- contingency load shedding"
    return opfm
end
""" Objective with base case and contingency generation and load shedding """
function add_obj_ccont!(opfm::OPFmodel, ramp_minutes, repair_time = 1.0)
    opfm = add_obj_cont!(opfm, ramp_minutes)
    set_objective_function(opfm.mod, objective_function(opfm.mod) +
        sum(opfm.cost[g] * sum(opfm.prob[c] * repair_time * (opfm.mod[:pgu][g,c] #=+ opfm.mod[:pgd][g,c]*0.1=#) 
            for c in opfm.contingencies) for g in get_name.(get_ctrl_generation(opfm.sys))
        ) +
        sum(opfm.voll[d] * (sum(opfm.prob[c] * opfm.mod[:lscc][d,c] for c in opfm.contingencies)) 
            for d in get_name.(get_nonctrl_generation(opfm.sys))
        )
    )
    @info "- corrective generation and load shedding"
    return opfm
end

""" Add unit commitment to thermal generation (not hydro) """
function add_unit_commit!(opfm::OPFmodel)
    @variable(opfm.mod, u[g in get_name.(gens_t(opfm.sys))], Bin)
    delete_lower_bound.(opfm.mod[:pg0])
    delete_upper_bound.(opfm.mod[:pg0])
    p_lim = make_named_array(get_active_power_limits, gens_t(opfm.sys))
    @constraint(opfm.mod, ucu[g in get_name.(gens_t(opfm.sys))], opfm.mod[:pg0][g] .<= getindex.(p_lim,2)[g] .* u[g]) 
    @constraint(opfm.mod, ucd[g in get_name.(gens_t(opfm.sys))], opfm.mod[:pg0][g] .>= getindex.(p_lim,1)[g] .* u[g]) 
    set_objective_function(opfm.mod, objective_function(opfm.mod) + 
        sum(u[get_name(g)] * get_generator_cost(g)[end] for g in gens_t(opfm.sys)))
    @info "Unit commitment added to the base case"
    return opfm
end
""" Add unit commitment """
function add_unit_commit_ccont!(opfm::OPFmodel)
    opfm = add_unit_commit!(opfm)
    delete.(opfm.mod, opfm.mod[:pg_lim])
    unregister(opfm.mod, :pg_lim)
    p_lim = make_named_array(get_active_power_limits, get_ctrl_generation(opfm.sys))
    @constraint(opfm.mod, pg_lim_d[g = get_name.(gens_t(opfm.sys)), c = opfm.contingencies], 
        0 <= opfm.mod[:pg0][g] + (opfm.mod[:pgu][g,c] - opfm.mod[:pgd][g,c])
    )
    @constraint(opfm.mod, pg_lim_u[g = get_name.(gens_t(opfm.sys)), c = opfm.contingencies], 
        opfm.mod[:pg0][g] + (opfm.mod[:pgu][g,c] - opfm.mod[:pgd][g,c]) <= p_lim[g].max * opfm.mod[:u][g]
    )
    @constraint(opfm.mod, pg_lim[g = get_name.(gens_h(opfm.sys)), c = opfm.contingencies], 
        0 <= opfm.mod[:pg0][g] + (opfm.mod[:pgu][g,c] - opfm.mod[:pgd][g,c]) <= p_lim[g].max
    )
    @info "- Constricted to contingencies"
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
    @constraint(opfm.mod, inj_p[n = get_name.(nodes(opfm.sys))], ctrl[n] .- 
        sum(beta(n, l) * opfm.mod[:pf0][get_name(l)] for l in branches(opfm.sys)) .== 
        sum((d.bus.name == n ? get_active_power(d) - opfm.mod[:ls0][get_name(d)] : 0 for d in demands(opfm.sys)), init = 0.0) + 
        sum((d.bus.name == n ? -get_active_power(d) + opfm.mod[:ls0][get_name(d)] : 0 for d in renewables(opfm.sys)), init = 0.0)
    )
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
    @constraint(opfm.mod, inj_pc[n = get_name.(nodes(opfm.sys)), c = opfm.contingencies], opfm.mod[:ctrl][n] .- 
        sum(beta(n,l) * opfm.mod[:pfc][get_name(l),c] for l in branches(opfm.sys)) .== 
        sum(d.bus.name == n ? get_active_power(d) - opfm.mod[:lsc][get_name(d),c] : 0 for d in demands(opfm.sys)) +
        sum(d.bus.name == n ? -get_active_power(d) + opfm.mod[:lsc][get_name(d),c] : 0 for d in renewables(opfm.sys))
    )
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
        sum(get_name(get_bus(g)) == n ? opfm.mod[:pg0][get_name(g)] + opfm.mod[:pgu][get_name(g),c] - 
            opfm.mod[:pgd][get_name(g),c] : 0 for g in get_ctrl_generation(opfm.sys)) -
        sum(beta(n,l) * opfm.mod[:pfcc][get_name(l),c] for l in branches(opfm.sys)) == 
        sum(d.bus.name == n ? get_active_power(d) - opfm.mod[:lscc][get_name(d),c] : 0 for d in demands(opfm.sys)) +
        sum(d.bus.name == n ? -get_active_power(d) + opfm.mod[:lscc][get_name(d),c] : 0 for d in renewables(opfm.sys))
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
    opfm = add_c_branch_cont!(opfm, short_term_limit_multi, branch_rating)
    x = make_named_array(get_x, branches(opfm.sys))
    @constraint(opfm.mod, pbcc[l = get_name.(branches(opfm.sys)), c = opfm.contingencies],
        opfm.mod[:pfcc][l,c] .- a(l,c) .* sum(beta(opfm.sys,l) .* opfm.mod[:vacc][:,c]) ./ x[l] .== 0
    )
    @info "- After contingency and corrective actions"
    return opfm
end

""" Add circuit breakers, then power flow on branch and branch limits for the base case and short-term contingency"""
function add_circuit_breakers_cont!(opfm::OPFmodel, short_term_limit_multi::Float64 = 1.0, branch_rating = nothing)
    @variable(opfm.mod, cbc[l in get_name.(branches(opfm.sys)), c in opfm.contingencies], Bin) # circuit breakers on branches in in contingencies before 
    branch_rating = isnothing(branch_rating) ? make_named_array(get_rate, branches(opfm.sys)) : branch_rating
    @constraint(opfm.mod, pfc_lim_l[l = get_name.(branches(opfm.sys)), c = opfm.contingencies], 
        -branch_rating[l] .* a(l,c) .* cbc[l,c] .* short_term_limit_multi .<= opfm.mod[:pfc][l,c]
    )
    @constraint(opfm.mod, pfc_lim_u[l = get_name.(branches(opfm.sys)), c = opfm.contingencies], 
        opfm.mod[:pfc][l,c] .<= branch_rating[l] .* a(l,c) .* cbc[l,c] .* short_term_limit_multi
    )
    opfm = add_c_branch!(opfm, branch_rating)
    x = make_named_array(get_x, branches(opfm.sys))
    @constraint(opfm.mod, pbc[l = get_name.(branches(opfm.sys)), c = opfm.contingencies],
        opfm.mod[:pfc][l,c] .- a(l,c) .* cbc[l,c] .* sum(beta(opfm.sys,l) .* opfm.mod[:vac][:,c]) ./ x[l] .== 0
    )
    @info "- After contingency, before corrective actions, with circuit breakers"
    return opfm
end
""" Add circuit breakers, then power flow on branch and branch limits for the base case and short-term and long-term contingency """
function add_circuit_breakers_ccont!(opfm::OPFmodel, short_term_limit_multi::Float64 = 1.0, branch_rating = nothing)
    @variable(opfm.mod, cbcc[l in get_name.(branches(opfm.sys)), c in opfm.contingencies], Bin) # and after corrective actions
    branch_rating = isnothing(branch_rating) ? make_named_array(get_rate, branches(opfm.sys)) : branch_rating
    @constraint(opfm.mod, pfcc_lim_l[l = get_name.(branches(opfm.sys)), c = opfm.contingencies], 
        -branch_rating[l] .* a(l,c) .* cbcc[l,c] .<= opfm.mod[:pfcc][l,c]
    )
    @constraint(opfm.mod, pfcc_lim_u[l = get_name.(branches(opfm.sys)), c = opfm.contingencies], 
        opfm.mod[:pfcc][l,c] .<= branch_rating[l] .* a(l,c) .* cbcc[l,c]
    )
    opfm = add_circuit_breakers_cont!(opfm, short_term_limit_multi, branch_rating)
    x = make_named_array(get_x, branches(opfm.sys))
    @constraint(opfm.mod, pbcc[l = get_name.(branches(opfm.sys)), c = opfm.contingencies],
        opfm.mod[:pfcc][l,c] .- a(l,c) .* cbcc[l,c] .* sum(beta(opfm.sys,l) .* opfm.mod[:vacc][:,c]) ./ x[l] .== 0
    )
    @info "- After contingency and corrective actions, with circuit breakers"
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
    @constraint(opfm.mod, load_shed[d in get_name.(demands(opfm.sys))], 
        opfm.mod[:ls0][d] <= p[d] * max_shed
    )
    @info "Restriction of load shedding to $(max_shed) of load per node: Base case"
    return opfm
end
function add_lim_load_shed_cont!(opfm::OPFmodel, max_shed::Float64 = 1.0)
    p = make_named_array(get_active_power, demands(opfm.sys))
    @constraint(opfm.mod, load_shed_cont[d in get_name.(demands(opfm.sys)), c in opfm.contingencies], 
        opfm.mod[:lsc][d,c] <= p[d] * max_shed
    )
    @info "- After contingency, before corrective actions"
    return opfm
end
function add_lim_load_shed_ccont!(opfm::OPFmodel, max_shed::Float64 = 1.0)
    p = make_named_array(get_active_power, demands(opfm.sys))
    @constraint(opfm.mod, load_shed_ccont[d in get_name.(demands(opfm.sys)), c in opfm.contingencies], 
        opfm.mod[:lscc][d,c] <= p[d] * max_shed
    )
    @info "- After contingency and corrective actions"
    return opfm
end
""" Restrict load shedding to renewable production per node """
function add_lim_renewable_shed!(opfm::OPFmodel, max_curtail::Float64 = 1.0)
    p = make_named_array(get_active_power, renewables(opfm.sys))
    @constraint(opfm.mod, renew_shed[d = get_name.(renewables(opfm.sys))], 
        opfm.mod[:ls0][d] <= p[d] * max_curtail
    )
    @info "Restriction of renewable shedding: Base case"
    return opfm
end
function add_lim_renewable_shed_cont!(opfm::OPFmodel, max_curtail::Float64 = 1.0)
    p = make_named_array(get_active_power, renewables(opfm.sys))
    @constraint(opfm.mod, renew_shed_cont[d = get_name.(renewables(opfm.sys)), c = opfm.contingencies], 
        opfm.mod[:lsc][d,c] <= p[d] * max_curtail
    )
    @info "- After contingency, before corrective actions"
    return opfm
end
function add_lim_renewable_shed_ccont!(opfm::OPFmodel, max_curtail::Float64 = 1.0)
    p = make_named_array(get_active_power, renewables(opfm.sys))
    @constraint(opfm.mod, renew_shed_ccont[d = get_name.(renewables(opfm.sys)), c = opfm.contingencies], 
        opfm.mod[:lscc][d,c] <= p[d] * max_curtail
    )
    @info "- After contingency and corrective actions"
    return opfm
end
""" Restrict load shedding to load (or renewable production) per node """
add_lim_nonctrl_shed!(opfm::OPFmodel, max_shed::Float64 = 1.0, max_curtail::Float64 = 1.0) = 
    add_lim_renewable_shed!(add_lim_load_shed!(opfm, max_shed), max_curtail)
add_lim_nonctrl_shed_cont!(opfm::OPFmodel, max_shed::Float64 = 1.0, max_curtail::Float64 = 1.0) = 
    add_lim_load_shed_cont!(add_lim_renewable_shed_cont!(add_lim_nonctrl_shed!(opfm, max_shed, max_curtail), max_curtail), max_shed)
add_lim_nonctrl_shed_ccont!(opfm::OPFmodel, max_shed::Float64 = 1.0, max_curtail::Float64 = 1.0) = 
    add_lim_load_shed_ccont!(add_lim_renewable_shed_ccont!(add_lim_nonctrl_shed_cont!(opfm, max_shed, max_curtail), max_curtail), max_shed)

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
function scopf(type::OPF, system::System, optimizer; 
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
    voll = isnothing(voll) ? make_voll(system) : voll
    contingencies = (type != SC::OPF && isnothing(contingencies)) ? get_name.(branches(system)) : contingencies
    prob = (type == PCSC::OPF && isnothing(prob)) ? make_prob(contingencies) : prob
    set_renewable_prod!(system, ratio)

    opfm = opfmodel(system, optimizer, time_limit_sec, voll, contingencies, prob)
    if type == SC::OPF
        opfm = scopf(opfm, unit_commit=unit_commit, max_shed=max_shed, max_curtail=max_curtail, ratio=ratio)
    elseif type == PSC::OPF
        opfm = p_scopf(opfm, unit_commit=unit_commit, max_shed=max_shed, max_curtail=max_curtail, ratio=ratio, 
        circuit_breakers=circuit_breakers, short_term_limit_multi=short_term_limit_multi, ramp_minutes=ramp_minutes)
    elseif type == PCSC::OPF
        opfm = pc_scopf(opfm, unit_commit=unit_commit, max_shed=max_shed, max_curtail=max_curtail, ratio=ratio, 
        circuit_breakers=circuit_breakers, short_term_limit_multi=short_term_limit_multi, ramp_minutes=ramp_minutes,
        repair_time=repair_time)
    end
    return opfm
end

function scopf(opfm::OPFmodel; 
        unit_commit::Bool = false,
        max_shed::Float64 = 0.1,
        max_curtail::Float64 = 1.0,
        ratio::Float64 = 0.5)
    opfm = init_var_dc_SCOPF!(opfm) |> add_c_bus! |> add_c_branch! |> add_obj!
    opfm = add_lim_nonctrl_shed!(opfm, max_shed, max_curtail)
    if unit_commit
        opfm = add_unit_commit!(opfm)
    end
    return opfm
end

""" Run a P-SCOPF of a power system """
function p_scopf(opfm::OPFmodel; 
        unit_commit::Bool=false,
        max_shed::Float64 = 0.1,
        max_curtail::Float64 = 1.0,
        ratio::Float64 = 0.5,
        circuit_breakers::Bool = false,
        short_term_limit_multi::Float64 = 1.5,
        ramp_minutes::Int64 = 10)
    opfm = init_var_dc_P_SCOPF!(opfm) |> add_c_bus_cont! 
    opfm = add_obj_cont!(opfm, ramp_minutes)
    opfm = add_lim_nonctrl_shed_cont!(opfm, max_shed, max_curtail)
    if circuit_breakers
        opfm = add_circuit_breakers_cont!(opfm, short_term_limit_multi)
    else
        opfm = add_c_branch_cont!(opfm, short_term_limit_multi)
    end
    if unit_commit
        opfm = add_unit_commit!(opfm)
    end
    return opfm
end

""" Run a PC-SCOPF of a power system """
function pc_scopf(opfm::OPFmodel; 
        unit_commit::Bool = false,
        max_shed::Float64 = 0.1,
        max_curtail::Float64 = 1.0,
        ratio::Float64= 0.5, 
        circuit_breakers::Bool=false,
        short_term_limit_multi::Float64 = 1.5,
        ramp_minutes::Int64 = 10,
        repair_time::Float64 = 1.0)
    opfm = init_var_dc_PC_SCOPF!(opfm) |> add_c_bus_ccont! 
    opfm = add_obj_ccont!(opfm, ramp_minutes, repair_time)
    # opfm = init_var_dc_PC_SCOPF!(opfm, max_shed) |> add_c_bus_ccont! |> add_obj!
    opfm = add_lim_nonctrl_shed_ccont!(opfm, max_shed, max_curtail)
    if circuit_breakers
        opfm = add_circuit_breakers_ccont!(opfm, short_term_limit_multi)
    else
        opfm = add_c_branch_ccont!(opfm, short_term_limit_multi)
    end
    opfm = add_lim_ramp_P_gen!(opfm, ramp_minutes)
    if unit_commit
        opfm = add_unit_commit_ccont!(opfm)
    end

    return opfm
end

end # module