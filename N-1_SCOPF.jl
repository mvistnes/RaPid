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

""" Initialize variables for dc SCOPF """
function init_var_dc_SCOPF!(OPFm::OPFmodel)
    p_lim = JuMP.Containers.DenseAxisArray(
        [get_active_power_limits(g).max for g in get_ctrl_generation(OPFm.sys)],
        get_name.(get_ctrl_generation(OPFm.sys))
    )
    demand = JuMP.Containers.DenseAxisArray(
        [get_active_power(d) for d in demands(OPFm.sys)],
        get_name.(get_nonctrl_generation(OPFm.sys))
    )
    @variables(OPFm.mod, begin
        0 <= pg0[g in get_name.(get_ctrl_generation(OPFm.sys))] <= p_lim[g]     # active power variables for the generators
            # restricted to positive numbers (can be changed to minimum active power) and less than maximum active power
        pf0[l in get_name.(branches(OPFm.sys))]      # power flow on branches in base case
        va0[b in get_name.(nodes(OPFm.sys))]         # voltage angle at a node in base case
        0 <= ls0[d in get_name.(get_nonctrl_generation(OPFm.sys))] <= demand[d]    # demand curtailment variables
            # restrict load shedding to max load per node
    end)
    @info "Variables added: pg0, pf0, va0, ls0"
    return OPFm
end
""" Initialize variables for dc P-SCOPF """
function init_var_dc_P_SCOPF!(OPFm::OPFmodel)
    OPFm = init_var_dc_SCOPF!(OPFm)
    # power flow on branches in contingencies
    @variable(OPFm.mod, pfc[l in get_name.(branches(OPFm.sys)), c in OPFm.contingencies])
    # voltage angle at a node in contingencies
    @variable(OPFm.mod, vac[b in get_name.(nodes(OPFm.sys)), c in OPFm.contingencies])
    @info "Variables added: pfc, vac"
    return OPFm
end
""" Initialize variables for dc PC-SCOPF """
function init_var_dc_PC_SCOPF!(OPFm::OPFmodel)
    OPFm = init_var_dc_P_SCOPF!(OPFm)
    @variables(OPFm.mod, begin
        # pgc[g in get_name.(get_ctrl_generation(OPFm.sys)), c in OPFm.contingencies] >= 0    # active power variables for the generators in contingencies 
        pgu[g in get_name.(get_ctrl_generation(OPFm.sys)), c in OPFm.contingencies] >= 0    # active power variables for the generators in contingencies ramp up 
        pgd[g in get_name.(get_ctrl_generation(OPFm.sys)), c in OPFm.contingencies] >= 0       # and ramp down
        pfcc[l in get_name.(branches(OPFm.sys)), c in OPFm.contingencies]         # power flow on branches in in contingencies after corrective actions
        vacc[b in get_name.(nodes(OPFm.sys)), c in OPFm.contingencies]            # voltage angle at a node in in contingencies after corrective actions
        lsc[d in get_name.(get_nonctrl_generation(OPFm.sys)), c in OPFm.contingencies] >= 0 # load curtailment variables in in contingencies
        # cbc[l in get_name.(branches(OPFm.sys)), c in OPFm.contingencies], Bin      # circuit breakers on branches in in contingencies before 
        # cbcc[l in get_name.(branches(OPFm.sys)), c in OPFm.contingencies], Bin         # and after corrective actions
    end)
    @info "Variables added: pgu, pgd, pfcc, vacc, lsc"
    return OPFm
end

""" Objective with base case generation and load shedding """
function add_obj!(OPFm::OPFmodel)
    @objective(OPFm.mod, Min, 
        sum(OPFm.mod[:pg0][get_name(g)] * (get_cost ∘ get_variable ∘ get_operation_cost)(g)[2] for g in gens_t(OPFm.sys)) + 
        sum(OPFm.mod[:pg0][get_name(g)] * 5 for g in gens_h(OPFm.sys)) + 
        sum(OPFm.voll[d] * OPFm.mod[:ls0][d] for d in get_name.(get_nonctrl_generation(OPFm.sys)))
    )
    @info "Objective added: base case generation and load shedding"
    return OPFm
end
""" Objective with base case and contingency generation and load shedding """
function add_obj_corrective!(OPFm::OPFmodel)
    @objective(OPFm.mod, Min, 
        sum((OPFm.mod[:pg0][get_name(g)] * (get_cost ∘ get_variable ∘ get_operation_cost)(g)[2] + 
            # sum(OPFm.prob[c] * (OPFm.mod[:pgc][get_name(g),c] - OPFm.mod[:pg0][get_name(g)]) for c in contingencies))
            sum(OPFm.prob[c] * (OPFm.mod[:pgu][get_name(g),c] #=+ OPFm.mod[:pgd][get_name(g),c]*0.1=#) for c in OPFm.contingencies))
            for g in gens_t(OPFm.sys)
        ) + 
        # sum(5 * (OPFm.mod[:pg0][get_name(g)] + sum(OPFm.prob[c] * (OPFm.mod[:pgc][get_name(g),c] - OPFm.mod[:pg0][get_name(g)])
        sum(5 * (OPFm.mod[:pg0][get_name(g)] + sum(OPFm.prob[c] * (OPFm.mod[:pgu][get_name(g),c] #=+ OPFm.mod[:pgd][get_name(g),c]*0.1=#)
            for c in OPFm.contingencies) ) for g in gens_h(OPFm.sys)
        ) +
        sum(OPFm.voll[d] * (OPFm.mod[:ls0][d] + sum(OPFm.prob[c] * OPFm.mod[:lsc][d,c] for c in OPFm.contingencies)) 
            for d in get_name.(get_nonctrl_generation(OPFm.sys)))
    )
    @info "Objective added: base case and contingency generation and load shedding"
    return OPFm
end

""" Set voltage angle at reference bus """
function set_ref_angle!(OPFm::OPFmodel)
    slack = find_slack(OPFm.sys)
    @constraint(OPFm.mod, OPFm.mod[:va0][get_name(slack)] == 0)
    @info "Voltage angle at the reference bus set to 0"
    return OPFm
end

""" Incerted power at each bus for the base case """
function add_c_bus!(OPFm::OPFmodel)
    @expression(OPFm.mod, ctrl[b = get_name.(nodes(OPFm.sys))], 
        sum(g.bus.name == b ? OPFm.mod[:pg0][get_name(g)] : 0 for g in get_ctrl_generation(OPFm.sys)))
    @expression(OPFm.mod, nonctrl[b = get_name.(nodes(OPFm.sys))], 
        sum((d.bus.name == b ? get_active_power(d) - OPFm.mod[:ls0][get_name(d)] : 0 for d in demands(OPFm.sys)), init = 0.0) + 
        sum((d.bus.name == b ? -get_active_power(d) + OPFm.mod[:ls0][get_name(d)] : 0 for d in renewables(OPFm.sys)), init = 0.0))
    @constraint(OPFm.mod, inj_p[n = get_name.(nodes(OPFm.sys))], ctrl[n] .- sum(beta(n, l) * OPFm.mod[:pf0][get_name(l)] for l in branches(OPFm.sys)) .== nonctrl[n])
    # @constraint(OPFm.mod, va_lim, -π/2 .<= OPFm.mod[:va0] .<= π/2) # Not really needed, could be implemented with spesific line angle limits
    @info "Contraints for power balance at each bus: Base case"
    return OPFm
end
""" Incerted power at each bus for the base case and short-term contingencies """
function add_c_bus_cont!(OPFm::OPFmodel)
    OPFm = add_c_bus!(OPFm)
    @expression(OPFm.mod, 
        n_pfc[n = get_name.(nodes(OPFm.sys)), c = OPFm.contingencies], 
        sum(beta(n,l) * OPFm.mod[:pfc][get_name(l),c] for l in branches(OPFm.sys))
    )
    @constraint(OPFm.mod, inj_pc[n = get_name.(nodes(OPFm.sys)), c = OPFm.contingencies], OPFm.mod[:ctrl][n] .- n_pfc[n,c] .== OPFm.mod[:nonctrl][n])
    # @constraint(OPFm.mod, vac_lim, -π/2 .<= OPFm.mod[:vac] .<= π/2) # Not really needed
    @info "- After contingency, before corrective actions"
    return OPFm
end
""" Incerted power at each bus for the base case and short-term and long-term contingencies """
function add_c_bus_ccont!(OPFm::OPFmodel)
    OPFm = add_c_bus_cont!(OPFm)
    @constraint(OPFm.mod, inj_pcc[n = get_name.(nodes(OPFm.sys)), c = OPFm.contingencies],
        # sum(get_name(get_bus(g)) == n ? pgc[get_name(g),c] : 0 for g in get_ctrl_generation(OPFm)) -
        sum(get_name(get_bus(g)) == n ? OPFm.mod[:pg0][get_name(g)] + OPFm.mod[:pgu][get_name(g),c] - OPFm.mod[:pgd][get_name(g),c] : 0 for g in get_ctrl_generation(OPFm.sys)) -
        sum(beta(n,l) * OPFm.mod[:pfcc][get_name(l),c] for l in branches(OPFm.sys)) == 
        sum(get_name(get_bus(d)) == n ? get_active_power(d) - OPFm.mod[:lsc][get_name(d),c] : 0 for d in demands(OPFm.sys)) +
        sum(get_name(get_bus(d)) == n ? -get_active_power(d) + OPFm.mod[:lsc][get_name(d),c] : 0 for d in renewables(OPFm.sys))
    )
    # @constraint(OPFm.mod, vacc_lim, -π/2 .<= OPFm.mod[:vacc] .<= π/2) # Not really needed
    @info "- After contingency and corrective actions"
    return OPFm
end

""" Power flow on branch and branch limits for the base case """
function add_c_branch!(OPFm::OPFmodel, branch_rating = nothing)
    branch_rating = isnothing(branch_rating) ? get_ratings(OPFm.sys) : branch_rating
    @constraint(OPFm.mod, pf0_lim[l = get_name.(branches(OPFm.sys))], 
        -branch_rating[l] <= OPFm.mod[:pf0][l] <= branch_rating[l]
    )
    x = JuMP.Containers.DenseAxisArray([get_x(l) for l in branches(OPFm.sys)], get_name.(branches(OPFm.sys)))
    @constraint(OPFm.mod, pb0[l = get_name.(branches(OPFm.sys))],
        OPFm.mod[:pf0][l] - sum(beta(OPFm.sys,l) .* OPFm.mod[:va0]) / x[l] == 0
    )
    @info "Power flow on lines: Base case"
    return OPFm
end
""" Power flow on branch and branch limits for the base case and short-term contingency """
function add_c_branch_cont!(OPFm::OPFmodel, short_term_limit_multi::Float64 = 1.0, branch_rating = nothing)
    branch_rating = isnothing(branch_rating) ? get_ratings(OPFm.sys) : branch_rating
    @constraint(OPFm.mod, pfc_lim[l = get_name.(branches(OPFm.sys)), c = OPFm.contingencies], 
        -branch_rating[l] .* short_term_limit_multi .<= a(l,c) .* OPFm.mod[:pfc][l,c] .<= branch_rating[l] .* short_term_limit_multi
    )
    OPFm = add_c_branch!(OPFm, branch_rating)
    x = JuMP.Containers.DenseAxisArray([get_x(l) for l in branches(OPFm.sys)], get_name.(branches(OPFm.sys)))
    @constraint(OPFm.mod, pbc[l = get_name.(branches(OPFm.sys)), c = OPFm.contingencies],
        OPFm.mod[:pfc][l,c] .- a(l,c) .* sum(beta(OPFm.sys,l) .* OPFm.mod[:vac][:,c]) ./ x[l] .== 0
    )
    @info "- After contingency, before corrective actions"
    return OPFm
end
""" Power flow on branch and branch limits for the base case and short-term and long-term contingency """
function add_c_branch_ccont!(OPFm::OPFmodel, short_term_limit_multi::Float64 = 1.0, branch_rating = nothing)
    branch_rating = isnothing(branch_rating) ? get_ratings(OPFm.sys) : branch_rating
    @constraint(OPFm.mod, pfcc_lim[l = get_name.(branches(OPFm.sys)), c = OPFm.contingencies], 
        -branch_rating[l] .<= a(l,c) .* OPFm.mod[:pfcc][l,c] .<= branch_rating[l]
    )
    OPFm = add_c_branch_cont!(OPFm, short_term_limit_multi)
    x = JuMP.Containers.DenseAxisArray([get_x(l) for l in branches(OPFm.sys)], get_name.(branches(OPFm.sys)))
    @constraint(OPFm.mod, pbcc[l = get_name.(branches(OPFm.sys)), c = OPFm.contingencies],
        OPFm.mod[:pfcc][l,c] .- a(l,c) .* sum(beta(OPFm.sys,l) .* OPFm.mod[:vacc][:,c]) ./ x[l] .== 0
    )
    @info "- After contingency and corrective actions"
    return OPFm
end
get_ratings(sys::System) = JuMP.Containers.DenseAxisArray(
            [get_rate(l) for l in branches(sys)], get_name.(branches(sys)))

""" Restrict active power generation to a minimum level """
function add_min_P_gen!(OPFm::OPFmodel)
    for g in get_ctrl_generation(OPFm.sys)
        set_lower_bound(OPFm.mod.pg0[get_name(g)], get_active_power_limits(g).min) 
    end
    @info "Minimum power generation limits added"
    return OPFm
end

""" Restrict active power generation ramp to min and max values """
function add_lim_ramp_P_gen!(OPFm::OPFmodel, ramp_minutes)
    p_lim = JuMP.Containers.DenseAxisArray(
        [get_active_power_limits(g).max for g in get_ctrl_generation(OPFm.sys)],
        get_name.(get_ctrl_generation(OPFm.sys))
    )
    for g in get_ctrl_generation(OPFm.sys)
        for c in OPFm.contingencies
            # @constraint(OPFm.mod, 0 <= OPFm.mod[:pg0][get_name(g)] <= get_active_power_limits(g).max)
            # @constraint(OPFm.mod, 0 <= OPFm.mod[:pgc][get_name(g),c] <= get_active_power_limits(g).max)
            # @constraint(OPFm.mod, -get_ramp_limits(g).down * ramp_minutes <= 
            #     OPFm.mod[:pgc][get_name(g),c] - OPFm.mod[:pg0][get_name(g)] <= get_ramp_limits(g).up * ramp_minutes
            # )
            set_upper_bound(OPFm.mod[:pgu][get_name(g),c], get_ramp_limits(g).up * ramp_minutes)
            set_upper_bound(OPFm.mod[:pgd][get_name(g),c], get_ramp_limits(g).down * ramp_minutes)
        end
    end
    @constraint(OPFm.mod, pg_lim[g = get_name.(get_ctrl_generation(OPFm.sys)), c = OPFm.contingencies], 
        0 <= OPFm.mod[:pg0][g] + (OPFm.mod[:pgu][g,c] - OPFm.mod[:pgd][g,c]) <= p_lim[g]
    )
    @info "Ramping restrictions on generators added"
    return OPFm
end

""" Restrict load shedding to max shed amount of load per node """
function add_lim_load_shed!(OPFm::OPFmodel, max_shed::Float64 = 1.0)
    p = JuMP.Containers.DenseAxisArray(
        [get_active_power(d) for d in demands(OPFm.sys)],
        get_name.(demands(OPFm.sys))
    )
    @constraint(OPFm.mod, load_shed[d in get_name.(demands(OPFm.sys)), c in OPFm.contingencies], 
        OPFm.mod[:ls0][d] + OPFm.mod[:lsc][d,c] <= p[d] * max_shed
    )
    @info "Restriction of load shedding to $(max_shed) of load per node"
    return OPFm
end
""" Restrict load shedding to renewable production per node """
function add_lim_renewable_shed!(OPFm::OPFmodel)
    p = JuMP.Containers.DenseAxisArray(
        [get_active_power(d) for d in renewables(OPFm.sys)],
        get_name.(renewables(OPFm.sys))
    )
    @constraint(OPFm.mod, renew_shed[d = get_name.(renewables(OPFm.sys)), c = OPFm.contingencies], 
        OPFm.mod[:ls0][d] + OPFm.mod[:lsc][d,c] <= p[d]
    )
    @info "Restriction of renewable shedding"
    return OPFm
end
""" Restrict load shedding to load (or renewable production) per node """
add_lim_nonctrl_shed!(OPFm::OPFmodel, max_shed::Float64 = 1.0) = add_lim_load_shed!(OPFm, max_shed) |> add_lim_renewable_shed!

""" Add unit commitment """
function add_unit_commit!(OPFm::OPFmodel)
    @variable(OPFm.mod, u[g in get_name.(get_ctrl_generation(OPFm.sys))], Bin)
    delete_lower_bound.(OPFm.mod[:pg0])
    delete_upper_bound.(OPFm.mod[:pg0])
    p_lim = JuMP.Containers.DenseAxisArray(
        [get_active_power_limits(g) for g in get_ctrl_generation(OPFm.sys)],
        get_name.(get_ctrl_generation(OPFm.sys))
    )
    @constraint(OPFm.mod, ucu, OPFm.mod[:pg0] .<= getindex.(p_lim,2) .* u) 
    @constraint(OPFm.mod, ucd, OPFm.mod[:pg0] .>= getindex.(p_lim,1) .* u) 
    set_objective_function(OPFm.mod, objective_function(OPFm.mod) + 
        sum(u[get_name(g)] * (get_cost ∘ get_variable ∘ get_operation_cost)(g)[end] for g in get_ctrl_generation(OPFm.sys)))
    @info "Unit commitment added to the base case"
    return OPFm
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
    termination_status(model) != MOI.OPTIMAL && @warn "Model not optimally solved with status $(termination_status(model))!"
    @info "Model solved in $(solve_time(model)) seconds"
    return model
end

""" Run a SCOPF of a power system """
function scopf(system::System, optimizer; time_limit_sec = 60, voll = nothing, unit_commit::Bool=false)
    voll = isnothing(voll) ? make_voll(system) : voll
    opfm = opfmodel(system, optimizer, time_limit_sec, voll)
    opfm = init_var_dc_SCOPF!(opfm) |> set_ref_angle! |> add_c_bus! |> add_c_branch! |> add_obj!
    if unit_commit
        opfm = add_unit_commit!(opfm)
    end
    return opfm
end

""" Run a P-SCOPF of a power system """
function p_scopf(system::System, optimizer; time_limit_sec = 60, voll = nothing, contingencies = nothing, unit_commit::Bool=false)
    voll = isnothing(voll) ? make_voll(system) : voll
    contingencies = isnothing(contingencies) ? get_name.(branches(system)) : contingencies
    opfm = opfmodel(system, optimizer, time_limit_sec, voll, contingencies)
    opfm = init_var_dc_P_SCOPF!(opfm) |> set_ref_angle! |> add_c_bus_cont! |> add_c_branch_cont! |> add_obj!
    if unit_commit
        opfm = add_unit_commit!(opfm)
    end
    return opfm
end

""" Run a PC-SCOPF of a power system """
function pc_scopf(system::System, optimizer; time_limit_sec = 60, voll = nothing, contingencies = nothing, prob = nothing, unit_commit::Bool=false)
    ramp_minutes = 10
    short_term_limit_multi = 1.5
    max_shed = 0.1

    voll = isnothing(voll) ? make_voll(system) : voll
    contingencies = isnothing(contingencies) ? get_name.(branches(system)) : contingencies
    prob = isnothing(prob) ? make_prob(contingencies) : prob
    opfm = opfmodel(system, optimizer, time_limit_sec, voll, contingencies, prob)
    opfm = init_var_dc_PC_SCOPF!(opfm) |> set_ref_angle! |> add_c_bus_ccont! |> add_obj_corrective!
    opfm = add_lim_nonctrl_shed!(opfm, max_shed)
    opfm = add_c_branch_ccont!(opfm, short_term_limit_multi) 
    opfm = add_lim_ramp_P_gen!(opfm, ramp_minutes)
    if unit_commit
        opfm = add_unit_commit!(opfm)
    end

    return opfm
end

