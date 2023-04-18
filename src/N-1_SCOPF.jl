# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

using PowerSystems
using JuMP
# nclude("utils.jl")
# include("SCOPF_ext.jl")

@enum OPF SC=0 PSC=1 PCSC=2 # SCOPF, P-SCOPF, PC-SCOPF

""" Run a SCOPF of a power system """
function scopf(type::OPF, system::System, optimizer; 
            voll = nothing, 
            cost_renewables = nothing, 
            contingencies = nothing, 
            prob = nothing,
            time_limit_sec::Int64 = 600,
            unit_commit::Bool = false,
            max_shed::Float64 = 0.1,
            max_curtail::Float64 = 1.0,
            renewable_prod::Float64= 0.5, 
            circuit_breakers::Bool=false,
            short_term_limit_multi::Float64 = 1.5,
            ramp_minutes::Int64 = 10,
            ramp_mult::Real = 10
        )
    voll = isnothing(voll) ? make_voll(system) : voll
    cost_renewables = isnothing(cost_renewables) ? make_cost_renewables(system) : cost_renewables
    contingencies = (type != SC::OPF && isnothing(contingencies)) ? get_name.(get_branches(system)) : contingencies
    prob = (type != SC::OPF && isnothing(prob)) ? make_prob(contingencies) : prob
    set_renewable_prod!(system, renewable_prod)

    opfm = opfmodel(system, optimizer, time_limit_sec, voll, cost_renewables, contingencies, prob)
    if type == SC::OPF
        return scopf(opfm, unit_commit=unit_commit, max_shed=max_shed, max_curtail=max_curtail)
    elseif type == PSC::OPF
        return p_scopf(opfm, unit_commit=unit_commit, max_shed=max_shed, max_curtail=max_curtail, 
        circuit_breakers=circuit_breakers, short_term_limit_multi=short_term_limit_multi)
    elseif type == PCSC::OPF
        return pc_scopf(opfm, unit_commit=unit_commit, max_shed=max_shed, max_curtail=max_curtail, 
        circuit_breakers=circuit_breakers, short_term_limit_multi=short_term_limit_multi,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes)
    end
end

""" Run a SCOPF of a power system, base case only """
function scopf(opfm::OPFmodel; 
            unit_commit::Bool = false,
            max_shed::Float64 = 0.1,
            max_curtail::Float64 = 1.0
        )
    init_var_dc_SCOPF!(opfm)
    add_c_bus!(opfm)
    add_c_branch!(opfm)
    add_c_dc_branch!(opfm)
    add_obj!(opfm)
    add_lim_load_shed!(opfm, max_shed=max_shed)
    listPr = make_named_array(get_active_power, get_renewables(opfm.sys)) 
    !isempty(listPr) && add_lim_renewable_shed!(opfm, listPr, max_curtail=max_curtail)
    if unit_commit
        add_unit_commit!(opfm)
    end
    return opfm
end

""" Run a P-SCOPF of a power system, preventive actions allowed for securing contingencies """
function p_scopf(opfm::OPFmodel; 
            unit_commit::Bool=false,
            max_shed::Float64 = 0.1,
            max_curtail::Float64 = 1.0,
            circuit_breakers::Bool = false,
            short_term_limit_multi::Float64 = 1.5
        )
    init_var_dc_P_SCOPF!(opfm) 
    islands = get_all_islands(opfm, find_slack(opfm.sys)[1]) 
    add_c_bus_cont!(opfm, islands)
    add_obj_cont!(opfm)
    add_lim_load_shed_cont!(opfm, islands, max_shed=max_shed)
    listPr = make_named_array(get_active_power, get_renewables(opfm.sys)) 
    !isempty(listPr) && add_lim_renewable_shed_cont!(opfm, islands, listPr, max_curtail=max_curtail)
    if circuit_breakers
        add_circuit_breakers_cont!(opfm, short_term_limit_multi)
    else
        add_c_branch_cont!(opfm, short_term_limit_multi=short_term_limit_multi)
        add_c_dc_branch_cont!(opfm)
    end
    if unit_commit
        add_unit_commit!(opfm)
    end
    return opfm
end

""" Run a PC-SCOPF of a power system, preventive and corrective actions allowed for securing contingencies """
function pc_scopf(opfm::OPFmodel; 
            unit_commit::Bool = false,
            max_shed::Float64 = 0.1,
            max_curtail::Float64 = 1.0,
            circuit_breakers::Bool=false,
            short_term_limit_multi::Float64 = 1.5,
            ramp_mult::Real = 10,
            ramp_minutes::Int64 = 10
        )
    init_var_dc_PC_SCOPF!(opfm)
    islands = get_all_islands(opfm, find_slack(opfm.sys)[1]) 
    add_c_bus_ccont!(opfm, islands)
    add_obj_ccont!(opfm, ramp_mult)
    # init_var_dc_PC_SCOPF!(opfm, max_shed) |> add_c_bus_ccont! |> add_obj!
    add_lim_load_shed_ccont!(opfm, islands, max_shed=max_shed)
    listPr = make_named_array(get_active_power, get_renewables(opfm.sys)) 
    !isempty(listPr) && add_lim_renewable_shed_ccont!(opfm, islands, listPr, max_curtail=max_curtail)
    if circuit_breakers
        add_circuit_breakers_ccont!(opfm, short_term_limit_multi)
    else
        add_c_branch_ccont!(opfm, short_term_limit_multi=short_term_limit_multi)
        add_c_dc_branch_ccont!(opfm)
    end
    add_lim_ramp_P_gen!(opfm, ramp_minutes)
    if unit_commit
        add_unit_commit_ccont!(opfm)
    end
    return opfm
end

""" Initialize variables for dc SCOPF """
function init_var_dc_SCOPF!(opfm::OPFmodel)
    p_lim = make_named_array(get_active_power_limits, get_ctrl_generation(opfm.sys))
    @variables(opfm.mod, begin
            0 <= pg0[g in get_name.(get_ctrl_generation(opfm.sys))] <= p_lim[g].max     
                # active power variables for the generators
                # restricted to positive numbers (can be changed to minimum active power) and less than maximum active power
            pf0[l in get_name.(get_branches(opfm.sys))]
                # power flow on AC branches in base case
            pfdc0[l in get_name.(get_dc_branches(opfm.sys))]
                # power flow on DC branches
            va0[b in get_name.(get_nodes(opfm.sys))]
                # voltage angle at a node in base case
            0 <= ls0[d in get_name.(get_demands(opfm.sys))]
                # demand curtailment variables
            0 <= pr0[d in get_name.(get_renewables(opfm.sys))]
                # renewable curtailment variables
        end)
    @info "Variables added: pg0, pf0, va0, ls0"
end
""" Initialize variables for dc P-SCOPF """
function init_var_dc_P_SCOPF!(opfm::OPFmodel)
    init_var_dc_SCOPF!(opfm)    
    @variables(opfm.mod, begin
            pfc[l in get_name.(get_branches(opfm.sys)), c in opfm.contingencies]
                # power flow on AC branches in contingencies
            pfdcc[l in get_name.(get_dc_branches(opfm.sys)), c in opfm.contingencies]
            # power flow on DC branches in contingencies
            vac[b in get_name.(get_nodes(opfm.sys)), c in opfm.contingencies]
                # voltage angle at a node in contingencies
            0 <= lsc[d in get_name.(get_demands(opfm.sys)), c in opfm.contingencies]
                # demand curtailment variables
            0 <= prc[d in get_name.(get_renewables(opfm.sys)), c in opfm.contingencies]
                # renewable curtailment variables
        end)
    @info "Variables added: pfc, vac, lsc"
end
""" Initialize variables for dc PC-SCOPF """
function init_var_dc_PC_SCOPF!(opfm::OPFmodel)
    init_var_dc_P_SCOPF!(opfm)
    @variables(opfm.mod, begin
            0 <= pgu[g in get_name.(get_ctrl_generation(opfm.sys)), c in opfm.contingencies]
                # active power variables for the generators in contingencies ramp up 
            0 <= pgd[g in get_name.(get_ctrl_generation(opfm.sys)), c in opfm.contingencies]
                # and ramp down
            pfcc[l in get_name.(get_branches(opfm.sys)), c in opfm.contingencies]
                # power flow on branches in in contingencies after corrective actions
            pfdccc[l in get_name.(get_dc_branches(opfm.sys)), c in opfm.contingencies]
                # power flow on DC branches in contingencies
            vacc[b in get_name.(get_nodes(opfm.sys)), c in opfm.contingencies]
                # voltage angle at a node in in contingencies after corrective actions
            0 <= lscc[d in get_name.(get_demands(opfm.sys)), c in opfm.contingencies]
                # load curtailment variables in in contingencies
            0 <= prcc[d in get_name.(get_renewables(opfm.sys)), c in opfm.contingencies]
                # renewable curtailment variables in in contingencies
        end)
    @info "Variables added: pgu, pgd, pfcc, vacc, lscc"
end

""" Objective with base case generation and load shedding """
function add_obj!(opfm::OPFmodel)
    @objective(opfm.mod, Min, opfm.cost_gen.data' * opfm.mod[:pg0] + opfm.cost_renewables.data' * opfm.mod[:pr0] + opfm.voll.data' * opfm.mod[:ls0])
    @info "Objective added: base case generation and load shedding"
end
""" Objective with base case and contingency generation and load shedding """
function add_obj_cont!(opfm::OPFmodel)
    add_obj!(opfm)
    @objective(opfm.mod, Min, objective_function(opfm.mod) + 
            sum(opfm.prob[c] * (
                sum(opfm.voll.data' * opfm.mod[:lsc][:,c]) +
                sum(opfm.cost_renewables.data' * opfm.mod[:prc][:,c]))
                    for c in opfm.contingencies
                )
        )
    @info "- contingency load shedding"
end
""" Objective with base case and contingency generation and load shedding """
function add_obj_ccont!(opfm::OPFmodel, ramp_mult = 10)
    add_obj_cont!(opfm)
    @objective(opfm.mod, Min, objective_function(opfm.mod) + 
            sum(opfm.prob[c] * (
                sum(opfm.cost_gen.data' * ramp_mult * opfm.mod[:pgu][:,c]) + #+ opfm.mod[:pgd][g,c]*0.1) 
                sum(opfm.voll.data' * opfm.mod[:lscc][:,c]) +
                sum(opfm.cost_renewables.data' * opfm.mod[:prcc][:,c]))
                    for c in opfm.contingencies
            )
        )
    @info "- corrective generation and load shedding"
end

""" Set voltage angle at reference bus """
function set_ref_angle!(opfm::OPFmodel)
    i, slack = find_slack(opfm.sys)
    @constraint(opfm.mod, opfm.mod[:va0][get_name(slack)] == 0)
    @info "Voltage angle at the reference bus set to 0"
end

""" Incerted power at each bus for the base case """
function add_c_bus!(opfm::OPFmodel, slack::Tuple{Integer, Bus} = find_slack(opfm.sys), 
        listGen = make_list(opfm, get_ctrl_generation), 
        listD = make_list(opfm, get_demands), listRe = make_list(opfm, get_renewables))
    @expression(opfm.mod, ctrl[n = get_name.(get_nodes(opfm.sys))], 
        sum(opfm.mod[:pg0][g.name] for g in listGen[n]))
    @constraint(opfm.mod, inj_p[n = get_name.(get_nodes(opfm.sys))], ctrl[n] .- 
            sum(beta(n, l) * opfm.mod[:pf0][get_name(l)] for l in get_branches(opfm.sys)) .-
            sum(beta(n, l) * opfm.mod[:pfdc0][get_name(l)] for l in get_dc_branches(opfm.sys)) .== 
            sum((get_active_power(d) - opfm.mod[:ls0][get_name(d)] for d in listD[n]), init = 0.0) + 
            sum((-get_active_power(d) + opfm.mod[:pr0][get_name(d)] for d in listRe[n]), init = 0.0)
        )
    # @constraint(opfm.mod, va_lim, -π/2 .<= opfm.mod[:va0] .<= π/2) # Not really needed, could be implemented with spesific line angle limits
    @constraint(opfm.mod, opfm.mod[:va0][get_name(slack[2])] == 0) # Set voltage angle at reference bus
    @info "Constraints for power balance at each bus: Base case"
end
""" Incerted power at each bus for the base case and short-term contingencies """
function add_c_bus_cont!(opfm::OPFmodel, islands::Vector, slack::Tuple{Integer, Bus} = find_slack(opfm.sys), 
        listGen = make_list(opfm, get_ctrl_generation), 
        listD = make_list(opfm, get_demands), listRe = make_list(opfm, get_renewables))
    add_c_bus!(opfm, slack, listGen, listD, listRe)
    for (i_c,c) in zip(islands, opfm.contingencies)
        isempty(i_c) && continue
        for n in get_name.(get_sorted_nodes(opfm.sys))[i_c]
            @constraint(opfm.mod, opfm.mod[:ctrl][n] .- 
                sum(beta(n,l) * opfm.mod[:pfc][get_name(l),c] for l in get_branches(opfm.sys)) .- 
                sum(beta(n,l) * opfm.mod[:pfdcc][get_name(l),c] for l in get_dc_branches(opfm.sys)) .== 
                sum((get_active_power(d) - opfm.mod[:lsc][get_name(d),c] for d in listD[n]), init = 0.0) +
                sum((-get_active_power(d) + opfm.mod[:prc][get_name(d),c] for d in listRe[n]), init = 0.0)
            )
        end
    end
    # @constraint(opfm.mod, vac_lim, -π/2 .<= opfm.mod[:vac] .<= π/2) # Not really needed
    @constraint(opfm.mod, ref_vac[c = opfm.contingencies], opfm.mod[:vac][get_name(slack[2]),c] == 0) # Set voltage angle at reference bus
    @info "- After contingency, before corrective actions"
end
""" Incerted power at each bus for the base case and short-term and long-term contingencies """
function add_c_bus_ccont!(opfm::OPFmodel, islands::Vector, slack::Tuple{Integer, Bus} = find_slack(opfm.sys), 
        listGen = make_list(opfm, get_ctrl_generation), 
        listD = make_list(opfm, get_demands), listRe = make_list(opfm, get_renewables))
    add_c_bus_cont!(opfm, islands, slack, listGen, listD, listRe)
    for (i_c,c) in zip(islands, opfm.contingencies)
        isempty(i_c) && continue
        for n in get_name.(get_sorted_nodes(opfm.sys))[i_c]
            @constraint(opfm.mod,
                    sum(opfm.mod[:pg0][get_name(g)] + opfm.mod[:pgu][get_name(g),c] - 
                        opfm.mod[:pgd][get_name(g),c] for g in listGen[n]) -
                    sum(beta(n,l) * opfm.mod[:pfcc][get_name(l),c] for l in get_branches(opfm.sys)) -
                    sum(beta(n,l) * opfm.mod[:pfdccc][get_name(l),c] for l in get_dc_branches(opfm.sys)) == 
                    sum(get_active_power(d) - opfm.mod[:lscc][get_name(d),c] for d in listD[n]) +
                    sum(-get_active_power(d) + opfm.mod[:prcc][get_name(d),c] for d in listRe[n])
                )
        end
    end
    # @constraint(opfm.mod, vacc_lim, -π/2 .<= opfm.mod[:vacc] .<= π/2) # Not really needed
    @constraint(opfm.mod, ref_vacc[c = opfm.contingencies], opfm.mod[:vacc][get_name(slack[2]),c] == 0) # Set voltage angle at reference bus
    @info "- After contingency and corrective actions"
end

""" Power flow on branch and branch limits for the base case """
function add_c_branch!(
            opfm::OPFmodel, 
            x::JuMP.Containers.DenseAxisArray = make_named_array(get_x, get_branches(opfm.sys)), 
            branch_rating::JuMP.Containers.DenseAxisArray = make_named_array(get_rate, get_branches(opfm.sys))
        )
    @constraint(opfm.mod, pf0_lim[l = get_name.(get_branches(opfm.sys))], 
            -branch_rating[l] <= opfm.mod[:pf0][l] <= branch_rating[l]
        )
    @constraint(opfm.mod, pb0[l = get_name.(get_branches(opfm.sys))],
            opfm.mod[:pf0][l] - sum(beta(opfm.sys,l) .* opfm.mod[:va0]) / x[l] == 0
        )
    @info "Power flow on lines: Base case"
end
""" Power flow on branch and branch limits for the base case and short-term contingency """
function add_c_branch_cont!(
            opfm::OPFmodel, 
            x::JuMP.Containers.DenseAxisArray = make_named_array(get_x, get_branches(opfm.sys)), 
            branch_rating::JuMP.Containers.DenseAxisArray = make_named_array(get_rate, get_branches(opfm.sys)); 
            short_term_limit_multi::Float64 = 1.0
        )
    @constraint(opfm.mod, pfc_lim[l = get_name.(get_branches(opfm.sys)), c = opfm.contingencies], 
            -branch_rating[l] .* a(l,c) .* short_term_limit_multi .<= opfm.mod[:pfc][l,c] .<= branch_rating[l] .* a(l,c) .* short_term_limit_multi
        )
    add_c_branch!(opfm, x, branch_rating)
    @constraint(opfm.mod, pbc[l = get_name.(get_branches(opfm.sys)), c = opfm.contingencies],
            opfm.mod[:pfc][l,c] .- a(l,c) .* sum(beta(opfm.sys,l) .* opfm.mod[:vac][:,c]) ./ x[l] .== 0
        )
    @info "- After contingency, before corrective actions"
end
""" Power flow on branch and branch limits for the base case and short-term and long-term contingency """
function add_c_branch_ccont!(
            opfm::OPFmodel, 
            x::JuMP.Containers.DenseAxisArray = make_named_array(get_x, get_branches(opfm.sys)), 
            branch_rating::JuMP.Containers.DenseAxisArray = make_named_array(get_rate, get_branches(opfm.sys)); 
            short_term_limit_multi::Float64 = 1.0
        )
    @constraint(opfm.mod, pfcc_lim[l = get_name.(get_branches(opfm.sys)), c = opfm.contingencies], 
            -branch_rating[l] .* a(l,c) .<= opfm.mod[:pfcc][l,c] .<= branch_rating[l] .* a(l,c)
        )
    add_c_branch_cont!(opfm, x, branch_rating, short_term_limit_multi=short_term_limit_multi)
    @constraint(opfm.mod, pbcc[l = get_name.(get_branches(opfm.sys)), c = opfm.contingencies],
            opfm.mod[:pfcc][l,c] .- a(l,c) .* sum(beta(opfm.sys,l) .* opfm.mod[:vacc][:,c]) ./ x[l] .== 0
        )
    @info "- After contingency and corrective actions"
end

""" Power flow on DC branch and branch limits for the base case """
function add_c_dc_branch!(
            opfm::OPFmodel, 
            branch_rating::JuMP.Containers.DenseAxisArray = make_named_array(get_active_power_limits_from, get_dc_branches(opfm.sys))
        )
    @constraint(opfm.mod, pfdc0_lim[l = get_name.(get_dc_branches(opfm.sys))], 
            branch_rating[l].min <= opfm.mod[:pfdc0][l] <= branch_rating[l].max
        )
    @info "Power flow on dc lines: Base case"
end   
""" Power flow on DC branch and branch limits for the base case and short-term contingency """
function add_c_dc_branch_cont!(
            opfm::OPFmodel, 
            branch_rating::JuMP.Containers.DenseAxisArray = make_named_array(get_active_power_limits_from, get_dc_branches(opfm.sys))
        )
    @constraint(opfm.mod, pfdcc_lim[l = get_name.(get_dc_branches(opfm.sys)), c = opfm.contingencies], 
            branch_rating[l].min <= opfm.mod[:pfdcc][l,c] <= branch_rating[l].max
        )
    add_c_dc_branch!(opfm, branch_rating)
    @info "- After contingency, before corrective actions"
end  
""" Power flow on DC branch and branch limits for the base case and short-term and long-term contingency """
function add_c_dc_branch_ccont!(
            opfm::OPFmodel, 
            branch_rating::JuMP.Containers.DenseAxisArray = make_named_array(get_active_power_limits_from, get_dc_branches(opfm.sys))
        )
    @constraint(opfm.mod, pfdccc_lim[l = get_name.(get_dc_branches(opfm.sys)), c = opfm.contingencies], 
            branch_rating[l].min <= opfm.mod[:pfdccc][l,c] <= branch_rating[l].max
        )
    add_c_dc_branch_cont!(opfm, branch_rating)
    @info "- After contingency and corrective actions"
end   
    
""" Restrict active power generation to a minimum level """
function add_min_P_gen!(opfm::OPFmodel)
    for g in get_ctrl_generation(opfm.sys)
        set_lower_bound(opfm.mod.pg0[get_name(g)], get_active_power_limits(g).min) 
    end
    @info "Minimum power generation limits added"
end

""" Restrict active power generation ramp to min and max values """
function add_lim_ramp_P_gen!(opfm::OPFmodel, ramp_minutes)
    p_lim = make_named_array(get_active_power_limits, get_ctrl_generation(opfm.sys))
    for g in get_ctrl_generation(opfm.sys), c in opfm.contingencies
        set_upper_bound(opfm.mod[:pgu][get_name(g),c], get_ramp_limits(g).up * ramp_minutes)
        set_upper_bound(opfm.mod[:pgd][get_name(g),c], get_ramp_limits(g).down * ramp_minutes)
    end
    @constraint(opfm.mod, pg_lim[g = get_name.(get_ctrl_generation(opfm.sys)), c = opfm.contingencies], 
            0 <= opfm.mod[:pg0][g] + (opfm.mod[:pgu][g,c] - opfm.mod[:pgd][g,c]) <= p_lim[g].max
        )
    @info "Ramping restrictions on generators added"
end

""" Restrict load shedding to max shed amount of load per node """
function add_lim_load_shed!(opfm::OPFmodel, listPd = make_named_array(get_active_power, get_demands(opfm.sys)); 
        max_shed::Float64 = 1.0)
    @constraint(opfm.mod, load_shed[d in get_name.(get_demands(opfm.sys))], 
            opfm.mod[:ls0][d] <= listPd[d] * max_shed
        )
    @info "Restriction of load shedding to $(max_shed) of load per node: Base case"
end
function add_lim_load_shed_cont!(opfm::OPFmodel, islands::Vector, 
        listPd = make_named_array(get_active_power, get_demands(opfm.sys)), 
        listD = make_list(opfm, get_demands); max_shed::Float64 = 1.0)
    add_lim_load_shed!(opfm, listPd, max_shed=max_shed)
    for (i_c,c) in zip(islands, opfm.contingencies)
        isempty(i_c) && continue
        for (i,n) in enumerate(get_name.(get_sorted_nodes(opfm.sys)))
            if i ∈ i_c
                @constraint(opfm.mod, [d in get_name.(listD[n])], opfm.mod[:lsc][d,c] <= listPd[d] * max_shed)
            else
                @constraint(opfm.mod, [d in get_name.(listD[n])], opfm.mod[:lsc][d,c] == listPd[d])
            end
        end
    end
    @info "- After contingency, before corrective actions"
end
function add_lim_load_shed_ccont!(opfm::OPFmodel, islands::Vector, 
        listPd = make_named_array(get_active_power, get_demands(opfm.sys)), 
        listD = make_list(opfm, get_demands); max_shed::Float64 = 1.0)
    add_lim_load_shed_cont!(opfm, islands, listPd, listD, max_shed=max_shed)
    for (i_c,c) in zip(islands, opfm.contingencies)
        isempty(i_c) && continue
        for (i,n) in enumerate(get_name.(get_sorted_nodes(opfm.sys)))
            if i ∈ i_c
                @constraint(opfm.mod, [d in get_name.(listD[n])], opfm.mod[:lscc][d,c] <= listPd[d] * max_shed)
            else
                @constraint(opfm.mod, [d in get_name.(listD[n])], opfm.mod[:lscc][d,c] == listPd[d])
            end
        end
    end
    @info "- After contingency and corrective actions"
end

""" Restrict load shedding to renewable production per node """
function add_lim_renewable_shed!(opfm::OPFmodel, listPr = make_named_array(get_active_power, get_renewables(opfm.sys)); 
        max_curtail::Float64 = 1.0)
    @constraint(opfm.mod, renew_shed[d = get_name.(get_renewables(opfm.sys))], 
            opfm.mod[:pr0][d] <= listPr[d] * max_curtail
        )
    @info "Restriction of renewable shedding: Base case"
end
function add_lim_renewable_shed_cont!(opfm::OPFmodel, islands::Vector, 
        listPr = make_named_array(get_active_power, get_renewables(opfm.sys)), 
        listD = make_list(opfm, get_renewables); max_curtail::Float64 = 1.0)
    add_lim_renewable_shed!(opfm, listPr, max_curtail=max_curtail)
    for (i_c,c) in zip(islands, opfm.contingencies)
        isempty(i_c) && continue
        for (i,n) in enumerate(get_name.(get_sorted_nodes(opfm.sys)))
            if i ∈ i_c
                @constraint(opfm.mod, [d = get_name.(listD[n])], opfm.mod[:prc][d,c] <= listPr[d] * max_curtail)
            else
                @constraint(opfm.mod, [d = get_name.(listD[n])], opfm.mod[:prc][d,c] == listPr[d])
            end
        end
    end
    @info "- After contingency, before corrective actions"
end
function add_lim_renewable_shed_ccont!(opfm::OPFmodel, islands::Vector, 
        listPr = make_named_array(get_active_power, get_renewables(opfm.sys)), 
        listD = make_list(opfm, get_renewables); max_curtail::Float64 = 1.0)
    add_lim_renewable_shed_cont!(opfm, islands, listPr, listD, max_curtail=max_curtail)
    for (i_c,c) in zip(islands, opfm.contingencies)
        isempty(i_c) && continue
        for (i,n) in enumerate(get_name.(get_sorted_nodes(opfm.sys)))
            if i ∈ i_c
                @constraint(opfm.mod, [d = get_name.(listD[n])], opfm.mod[:prcc][d,c] <= listPr[d] * max_curtail)
            else
                @constraint(opfm.mod, [d = get_name.(listD[n])], opfm.mod[:prcc][d,c] == listPr[d])
            end
        end
    end
    @info "- After contingency and corrective actions"
end

