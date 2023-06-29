# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

using PowerSystems
using JuMP
# include("utils.jl")
# include("SCOPF_ext.jl")

@enum OPF SC=0 PSC=1 PCSC=2 # SCOPF, P-SCOPF, PC-SCOPF

""" Run a SCOPF of a power system """
function scopf(type::OPF, system::System, optimizer; 
            voll = nothing, 
            contingencies = nothing, 
            prob = nothing,
            time_limit_sec::Integer = 600,
            unit_commit::Bool = false,
            max_shed::Real = 1.0,
            max_curtail::Real = 1.0,
            renewable_prod::Real= 0.5, 
            circuit_breakers::Bool=false,
            short_term_limit_multi::Real = 1.5,
            ramp_minutes::Real = 10,
            ramp_mult::Real = 10,
            debug=false
        )
    contingencies = isnothing(contingencies) ? sort_components!(get_branches(system)) : contingencies
    prob = isnothing(prob) ? make_prob(contingencies) : prob
    # set_renewable_prod!(system, renewable_prod)
    # set_active_power_demand!(system)

    opfm = isnothing(voll) ? opfmodel(system, optimizer, time_limit_sec, debug=debug) : opfmodel(system, optimizer, time_limit_sec, voll, contingencies, prob, debug=debug)
    idx = get_nodes_idx(opfm.nodes)
    list = make_list(opfm, idx, opfm.nodes)
    if type == SC::OPF
        return scopf(opfm, list, idx, unit_commit=unit_commit, max_shed=max_shed, max_curtail=max_curtail)
    elseif type == PSC::OPF
        return p_scopf(opfm, list, idx, unit_commit=unit_commit, max_shed=max_shed, max_curtail=max_curtail, 
        circuit_breakers=circuit_breakers, short_term_limit_multi=short_term_limit_multi)
    elseif type == PCSC::OPF
        return pc_scopf(opfm, list, idx, unit_commit=unit_commit, max_shed=max_shed, max_curtail=max_curtail, 
        circuit_breakers=circuit_breakers, short_term_limit_multi=short_term_limit_multi,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes)
    end
end

""" Run a SCOPF of a power system, base case only """
function scopf(opfm::OPFmodel, list::Vector{CTypes}, idx::Dict{<:Any, <:Int}; 
            unit_commit::Bool = false,
            max_shed::Real = 0.1,
            max_curtail::Real = 1.0
        )
    init_var_dc_SCOPF!(opfm)
    slack = find_slack(opfm.nodes)
    add_c_bus!(opfm, slack[1], list)
    add_c_branch!(opfm, idx)
    add_c_dc_branch!(opfm)
    add_obj!(opfm)
    add_lim_load_shed!(opfm, max_shed=max_shed)
    listPr = get_active_power.(opfm.renewables)
    !isempty(listPr) && add_lim_renewable_shed!(opfm, listPr, max_curtail=max_curtail)
    if unit_commit
        add_unit_commit!(opfm)
    end
    return opfm
end

""" Run a P-SCOPF of a power system, preventive actions allowed for securing contingencies """
function p_scopf(opfm::OPFmodel, list::Vector{CTypes}, idx::Dict{<:Any, <:Int}; 
            unit_commit::Bool=false,
            max_shed::Real = 0.1,
            max_curtail::Real = 1.0,
            circuit_breakers::Bool = false,
            short_term_limit_multi::Real = 1.5
        )
    init_var_dc_P_SCOPF!(opfm) 
    slack = find_slack(opfm.nodes)
    islands = get_all_islands(opfm, slack[1]) 
    add_c_bus_cont!(opfm, islands, slack[1], list)
    add_obj_cont!(opfm)
    add_lim_P_gen!(opfm, islands, list)
    add_lim_load_shed_cont!(opfm, islands, list, max_shed=max_shed)
    listPr = get_active_power.(opfm.renewables)
    !isempty(listPr) && add_lim_renewable_shed_cont!(opfm, islands, list, listPr, max_curtail=max_curtail)
    if circuit_breakers
        add_circuit_breakers_cont!(opfm, short_term_limit_multi)
    else
        add_c_branch_cont!(opfm, idx, short_term_limit_multi=short_term_limit_multi)
    end
    if unit_commit
        add_unit_commit!(opfm)
    end
    return opfm
end

""" Run a PC-SCOPF of a power system, preventive and corrective actions allowed for securing contingencies """
function pc_scopf(opfm::OPFmodel, list::Vector{CTypes}, idx::Dict{<:Any, <:Int}; 
            unit_commit::Bool = false,
            max_shed::Real = 0.1,
            max_curtail::Real = 1.0,
            circuit_breakers::Bool=false,
            short_term_limit_multi::Real = 1.5,
            ramp_mult::Real = 10,
            ramp_minutes::Real = 10
        )
    init_var_dc_PC_SCOPF!(opfm)
    slack = find_slack(opfm.nodes)
    islands = get_all_islands(opfm, slack[1]) 
    add_c_bus_ccont!(opfm, islands, slack[1], list)
    add_obj_ccont!(opfm, ramp_mult)
    add_lim_P_gen!(opfm, islands, list)
    # init_var_dc_PC_SCOPF!(opfm, max_shed) |> add_c_bus_ccont! |> add_obj!
    add_lim_load_shed_ccont!(opfm, islands, list, max_shed=max_shed)
    listPr = get_active_power.(opfm.renewables)
    !isempty(listPr) && add_lim_renewable_shed_ccont!(opfm, islands, list, listPr, max_curtail=max_curtail)
    if circuit_breakers
        add_circuit_breakers_ccont!(opfm, short_term_limit_multi)
    else
        add_c_branch_ccont!(opfm, idx, short_term_limit_multi=short_term_limit_multi)
        add_c_dc_branch_ccont!(opfm)
    end
    add_lim_ramp_P_gen!(opfm, islands, list, ramp_minutes)
    if unit_commit
        add_unit_commit_ccont!(opfm)
    end
    return opfm
end

""" Initialize variables for dc SCOPF """
function init_var_dc_SCOPF!(opfm::OPFmodel)
    p_lim = get_active_power_limits.(opfm.ctrl_generation)
    @variables(opfm.mod, begin
            0 <= pg0[g in 1:length(opfm.ctrl_generation)] <= p_lim[g].max     
                # active power variables for the generators
                # restricted to positive numbers (can be changed to minimum active power) and less than maximum active power
            pf0[l in 1:length(opfm.branches)]
                # power flow on AC branches in base case
            pfdc0[l in 1:length(opfm.dc_branches)]
                # power flow on DC branches
            va0[b in 1:length(opfm.nodes)]
                # voltage angle at a node in base case
            0 <= ls0[d in 1:length(opfm.demands)]
                # demand curtailment variables
            0 <= pr0[d in 1:length(opfm.renewables)]
                # renewable curtailment variables
        end)
    @info "Variables added: pg0, pf0, va0, ls0"
end
""" Initialize variables for dc P-SCOPF """
function init_var_dc_P_SCOPF!(opfm::OPFmodel)
    init_var_dc_SCOPF!(opfm)    
    @variables(opfm.mod, begin
            0 <= pgc[g in 1:length(opfm.ctrl_generation), c in 1:length(opfm.contingencies)]
            pfc[l in 1:length(opfm.branches), c in 1:length(opfm.contingencies)]
                # power flow on AC branches in contingencies
            vac[b in 1:length(opfm.nodes), c in 1:length(opfm.contingencies)]
                # voltage angle at a node in contingencies
            0 <= lsc[d in 1:length(opfm.demands), c in 1:length(opfm.contingencies)]
                # demand curtailment variables
            0 <= prc[d in 1:length(opfm.renewables), c in 1:length(opfm.contingencies)]
                # renewable curtailment variables
        end)
    @info "Variables added: pfc, vac, lsc"
end
""" Initialize variables for dc PC-SCOPF """
function init_var_dc_PC_SCOPF!(opfm::OPFmodel)
    init_var_dc_P_SCOPF!(opfm)
    @variables(opfm.mod, begin
            0 <= pgu[g in 1:length(opfm.ctrl_generation), c in 1:length(opfm.contingencies)]
                # active power variables for the generators in contingencies ramp up 
            0 <= pgd[g in 1:length(opfm.ctrl_generation), c in 1:length(opfm.contingencies)]
                # and ramp down
            pfcc[l in 1:length(opfm.branches), c in 1:length(opfm.contingencies)]
                # power flow on branches in in contingencies after corrective actions
            pfdccc[l in 1:length(opfm.dc_branches), c in 1:length(opfm.contingencies)]
                # power flow on DC branches in contingencies
            vacc[b in 1:length(opfm.nodes), c in 1:length(opfm.contingencies)]
                # voltage angle at a node in in contingencies after corrective actions
            0 <= lscc[d in 1:length(opfm.demands), c in 1:length(opfm.contingencies)]
                # load curtailment variables in in contingencies
            0 <= prcc[d in 1:length(opfm.renewables), c in 1:length(opfm.contingencies)]
                # renewable curtailment variables in in contingencies
        end)
    @info "Variables added: pgu, pgd, pfcc, vacc, lscc"
end

""" Objective with base case generation and load shedding """
function add_obj!(opfm::OPFmodel)
    @objective(opfm.mod, Min, sum(x * opfm.mod[:pg0][i] for (i,x) in enumerate(opfm.cost_ctrl_gen))+ 
        # opfm.cost_renewables' * opfm.mod[:pr0] + 
        sum(x * opfm.mod[:ls0][i] for (i,x) in enumerate(opfm.voll)))
    @info "Objective added: base case generation and load shedding"
end
""" Objective with base case and contingency generation and load shedding """
function add_obj_cont!(opfm::OPFmodel)
    add_obj!(opfm)
    @objective(opfm.mod, Min, objective_function(opfm.mod) + 
            sum(opfm.prob[c] * (
                sum(x * opfm.mod[:lsc][i,c] for (i,x) in enumerate(opfm.voll)) # +
                # sum(opfm.voll' * opfm.mod[:lsc][:,c]) # +
                # sum(opfm.cost_renewables' * opfm.mod[:prc][:,c])
                ) for c in 1:length(opfm.contingencies)
            )
        )
    @info "- contingency load shedding"
end
""" Objective with base case and contingency generation and load shedding """
function add_obj_ccont!(opfm::OPFmodel, ramp_mult = 10)
    add_obj_cont!(opfm)
    @objective(opfm.mod, Min, objective_function(opfm.mod) + 
            sum(opfm.prob[c] * (
                sum(x * ramp_mult * (opfm.mod[:pgu][i,c] .+ opfm.mod[:pgd][i,c]) for (i,x) in enumerate(opfm.cost_ctrl_gen)) + 
                sum(x * opfm.mod[:lscc][i,c] for (i,x) in enumerate(opfm.voll)) # +
                # sum(opfm.cost_ctrl_gen' * ramp_mult * (opfm.mod[:pgu][:,c] .+ opfm.mod[:pgd][:,c])) + 
                # sum(opfm.voll' * opfm.mod[:lscc][:,c]) # +
                # sum(opfm.cost_renewables' * opfm.mod[:prcc][:,c])
                ) for c in 1:length(opfm.contingencies)
            )
        )
    @info "- corrective generation and load shedding"
end

""" Set voltage angle at reference bus """
function set_ref_angle!(opfm::OPFmodel)
    i, slack = find_slack(opfm.sys)
    @constraint(opfm.mod, opfm.mod[:va0][i] == 0)
    @info "Voltage angle at the reference bus set to 0"
end

""" Incerted power at each bus for the base case """
function add_c_bus!(opfm::OPFmodel, slack::Integer, list::Vector{CTypes})
    @constraint(opfm.mod, inj_p[n = 1:length(opfm.nodes)], 
            sum(beta(opfm.nodes[n], opfm.branches[l]) * opfm.mod[:pf0][l] for l in list[n].branches) == 
            sum(beta(opfm.nodes[n], opfm.dc_branches[l]) * opfm.mod[:pfdc0][l] for l in list[n].dc_branches) +
            sum(opfm.mod[:pg0][g] for g in list[n].ctrl_generation) + 
            sum((get_active_power(opfm.renewables[d]) - opfm.mod[:pr0][d] for d in list[n].renewables), init = 0.0) - 
            sum((get_active_power(opfm.demands[d]) - opfm.mod[:ls0][d] for d in list[n].demands), init = 0.0)
        )
    # @constraint(opfm.mod, va_lim, -π/2 .<= opfm.mod[:va0] .<= π/2) # Not really needed, could be implemented with spesific line angle limits
    @constraint(opfm.mod, opfm.mod[:va0][slack[1]] == 0) # Set voltage angle at reference bus
    @info "Constraints for power balance at each bus: Base case"
end
""" Incerted power at each bus for the base case and short-term contingencies """
function add_c_bus_cont!(opfm::OPFmodel, islands::Vector, slack::Integer, list::Vector{CTypes})
    add_c_bus!(opfm, slack, list)
    for (c,island) in enumerate(islands)
        isempty(island) && continue
        for n in island
            @constraint(opfm.mod, 
                    sum(beta(opfm.nodes[n], opfm.branches[l]) * opfm.mod[:pfc][l,c] for l in list[n].branches) == 
                    sum(beta(opfm.nodes[n], opfm.dc_branches[l]) * opfm.mod[:pfdc][l,c] for l in list[n].dc_branches) + 
                    sum(opfm.mod[:pg0][g] - opfm.mod[:pgc][g] for g in list[n].ctrl_generation) + 
                    sum((get_active_power(opfm.renewables[d]) - opfm.mod[:prc][d,c] for d in list[n].renewables), init = 0.0) -
                    sum((get_active_power(opfm.demands[d]) - opfm.mod[:lsc][d,c] for d in list[n].demands), init = 0.0)
                )
        end
    end
    # @constraint(opfm.mod, vac_lim, -π/2 .<= opfm.mod[:vac] .<= π/2) # Not really needed
    @constraint(opfm.mod, ref_vac[l = 1:length(opfm.contingencies)], opfm.mod[:vac][slack[1],l] .== 0) # Set voltage angle at reference bus
    # @constraint(opfm.mod, sum(opfm.mod[:lsc]) == sum(opfm.mod[:prc]) + sum(opfm.mod[:pgc]))
    @info "- After contingency, before corrective actions"
end
""" Incerted power at each bus for the base case and short-term and long-term contingencies """
function add_c_bus_ccont!(opfm::OPFmodel, islands::Vector, slack::Integer, list::Vector{CTypes})
    add_c_bus_cont!(opfm, islands, slack, list)
    for (c,island) in enumerate(islands)
        isempty(island) && continue
        for n in island
            @constraint(opfm.mod,
                    sum(beta(opfm.nodes[n], opfm.branches[l]) * opfm.mod[:pfcc][l,c] for l in list[n].branches) == 
                    sum(beta(opfm.nodes[n], opfm.dc_branches[l]) * opfm.mod[:pfdccc][l,c] for l in list[n].dc_branches) +
                    sum(opfm.mod[:pg0][g] + opfm.mod[:pgu][g,c] - opfm.mod[:pgd][g,c] for g in list[n].ctrl_generation) +
                    sum(get_active_power(opfm.renewables[d]) - opfm.mod[:prcc][d,c] for d in list[n].renewables) -
                    sum(get_active_power(opfm.demands[d]) - opfm.mod[:lscc][d,c] for d in list[n].demands)
                )
        end
    end
    # @constraint(opfm.mod, vacc_lim, -π/2 .<= opfm.mod[:vacc] .<= π/2) # Not really needed
    @constraint(opfm.mod, ref_vacc[l = 1:length(opfm.contingencies)], opfm.mod[:vacc][slack[1],l] .== 0) # Set voltage angle at reference bus
    # @constraint(opfm.mod, sum(opfm.mod[:pgu]) + sum(opfm.mod[:lscc]) == sum(opfm.mod[:pgd]) + sum(opfm.mod[:prcc]))
    @info "- After contingency and corrective actions"
end

""" Power flow on branch and branch limits for the base case """
function add_c_branch!(
            opfm::OPFmodel, idx,
            x::AbstractVector{<:Real} = get_x.(opfm.branches), 
            branch_rating::AbstractVector{<:Real} = get_rate.(opfm.branches)
        )
    @constraint(opfm.mod, pf0_lim_n, -branch_rating .<= opfm.mod[:pf0])
    @constraint(opfm.mod, pf0_lim_p, opfm.mod[:pf0] .<= branch_rating)
    @constraint(opfm.mod, pb0[l = 1:length(opfm.branches)],
            opfm.mod[:pf0][l] - (opfm.mod[:va0][idx[opfm.branches[l].arc.from.number]] - 
                opfm.mod[:va0][idx[opfm.branches[l].arc.to.number]]) / x[l] == 0
        )
    @info "Power flow on lines: Base case"
end
""" Power flow on branch and branch limits for the base case and short-term contingency """
function add_c_branch_cont!(
            opfm::OPFmodel, idx,
            x::AbstractVector{<:Real} = get_x.(opfm.branches), 
            branch_rating::AbstractVector{<:Real} = get_rate.(opfm.branches); 
            short_term_limit_multi::Real = 1.0
        )
    @constraint(opfm.mod, pfc_lim_n[l = 1:length(opfm.branches), c = 1:length(opfm.contingencies)], 
            -branch_rating[l] * !isequal(opfm.branches[l],opfm.contingencies[c]) * short_term_limit_multi <= opfm.mod[:pfc][l,c]
        )
    @constraint(opfm.mod, pfc_lim_p[l = 1:length(opfm.branches), c = 1:length(opfm.contingencies)], 
            opfm.mod[:pfc][l,c] <= !isequal(opfm.branches[l],opfm.contingencies[c]) * branch_rating[l] * short_term_limit_multi
        )
    add_c_branch!(opfm, idx, x, branch_rating)
    @constraint(opfm.mod, pbc[l = 1:length(opfm.branches), c = 1:length(opfm.contingencies)],
            opfm.mod[:pfc][l,c] - !isequal(opfm.branches[l],opfm.contingencies[c]) * (opfm.mod[:vac][idx[opfm.branches[l].arc.from.number],c] - 
            opfm.mod[:vac][idx[opfm.branches[l].arc.to.number],c]) / x[l] == 0
        )
    @info "- After contingency, before corrective actions"
end
""" Power flow on branch and branch limits for the base case and short-term and long-term contingency """
function add_c_branch_ccont!(
            opfm::OPFmodel, idx, 
            x::AbstractVector{<:Real} = get_x.(opfm.branches), 
            branch_rating::AbstractVector = get_rate.(opfm.branches); 
            short_term_limit_multi::Real = 1.0
        )
    @constraint(opfm.mod, pfcc_lim_n[l = 1:length(opfm.branches), c = 1:length(opfm.contingencies)], 
            -branch_rating[l] * !isequal(opfm.branches[l],opfm.contingencies[c]) <= opfm.mod[:pfcc][l,c]
        )
    @constraint(opfm.mod, pfcc_lim_p[l = 1:length(opfm.branches), c = 1:length(opfm.contingencies)], 
            opfm.mod[:pfcc][l,c] <= !isequal(opfm.branches[l],opfm.contingencies[c]) * branch_rating[l]
        )
    add_c_branch_cont!(opfm, idx, x, branch_rating, short_term_limit_multi=short_term_limit_multi)
    @constraint(opfm.mod, pbcc[l = 1:length(opfm.branches), c = 1:length(opfm.contingencies)],
            opfm.mod[:pfcc][l,c] - !isequal(opfm.branches[l],opfm.contingencies[c]) * (opfm.mod[:vacc][idx[opfm.branches[l].arc.from.number],c] - 
            opfm.mod[:vacc][idx[opfm.branches[l].arc.to.number],c]) / x[l] == 0
        )
    @info "- After contingency and corrective actions"
end

""" Power flow on DC branch and branch limits for the base case """
function add_c_dc_branch!(
            opfm::OPFmodel, 
            branch_rating::AbstractVector = get_active_power_limits_from.(opfm.dc_branches)
        )
    @constraint(opfm.mod, pfdc0_lim_n[l = 1:length(opfm.dc_branches)], branch_rating[l].min <= opfm.mod[:pfdc0][l])
    @constraint(opfm.mod, pfdc0_lim_p[l = 1:length(opfm.dc_branches)], opfm.mod[:pfdc0][l] <= branch_rating[l].max)
    @info "Power flow on dc lines: Base case"
end 
""" Power flow on DC branch and branch limits for the base case and short-term and long-term contingency """
function add_c_dc_branch_ccont!(
            opfm::OPFmodel, 
            branch_rating::AbstractVector = get_active_power_limits_from.(opfm.dc_branches)
        )
    @constraint(opfm.mod, pfdccc_lim_n[l = 1:length(opfm.dc_branches), c = 1:length(opfm.contingencies)], 
            branch_rating[l].min <= opfm.mod[:pfdccc][l,c]
        )
    @constraint(opfm.mod, pfdccc_lim_p[l = 1:length(opfm.dc_branches), c = 1:length(opfm.contingencies)], 
            opfm.mod[:pfdccc][l,c] <= branch_rating[l].max
        )
    add_c_dc_branch!(opfm, branch_rating)
    @info "- After contingency and corrective actions"
end   
    
""" Restrict active power generation to a minimum level """
function add_min_P_gen!(opfm::OPFmodel)
    for (i,g) in enumerate(opfm.ctrl_generation)
        set_lower_bound(opfm.mod.pg0[i], get_active_power_limits(g).min) 
    end
    @info "Minimum power generation limits added"
end

""" Restrict active power generation in contingency"""
function add_lim_P_gen!(opfm::OPFmodel, islands::Vector, list::Vector{CTypes})
    for (c,island) in enumerate(islands)
        if length(island) < length(opfm.nodes)
            for n in island
                @constraint(opfm.mod, [g in list[n].ctrl_generation], 0 <= opfm.mod[:pg0][g] - opfm.mod[:pgc][g,c])
            end
            for n in setdiff(1:length(opfm.nodes), island)
                @constraint(opfm.mod, [g in list[n].ctrl_generation], opfm.mod[:pg0][g] == opfm.mod[:pgc][g,c])
            end
        else
            @constraint(opfm.mod, [g in 1:length(opfm.ctrl_generation)], 0 .<= opfm.mod[:pg0] .- opfm.mod[:pgc][g,c])
        end
    end
    @info "Restriction power generation limits added"
end

""" Restrict active power generation ramp to min and max values """
function add_lim_ramp_P_gen!(opfm::OPFmodel, islands::Vector, list::Vector{CTypes}, ramp_minutes)
    p_lim = getindex.(get_active_power_limits.(opfm.ctrl_generation), :max)
    (rampup, rampdown) = split_pair(get_ramp_limits.(opfm.ctrl_generation))
    for (c,island) in enumerate(islands)
        if length(island) < length(opfm.nodes)
            for n in island
                @constraint(opfm.mod, [g in list[n].ctrl_generation], opfm.mod[:pgu][g,c] .<= rampup * ramp_minutes)
                @constraint(opfm.mod, [g in list[n].ctrl_generation], opfm.mod[:pgd][g,c] .<= rampdown * ramp_minutes)
                @constraint(opfm.mod, [g in list[n].ctrl_generation], 
                    0 <= opfm.mod[:pg0][g] + opfm.mod[:pgu][g,c] - opfm.mod[:pgd][g,c])
                @constraint(opfm.mod, [g in list[n].ctrl_generation], 
                    opfm.mod[:pg0][g] + opfm.mod[:pgu][g,c] - opfm.mod[:pgd][g,c] <= p_lim[g])
            end
            for n in setdiff(1:length(opfm.nodes), island)
                @constraint(opfm.mod, [g in list[n].ctrl_generation], 
                    opfm.mod[:pg0][g] == opfm.mod[:pgd][g,c] - opfm.mod[:pgu][g,c])
            end
        else
            @constraint(opfm.mod, [g in 1:length(opfm.ctrl_generation)], opfm.mod[:pgu][g,c] <= rampup[g] * ramp_minutes)
            @constraint(opfm.mod, [g in 1:length(opfm.ctrl_generation)], opfm.mod[:pgd][g,c] <= rampdown[g] * ramp_minutes)
            @constraint(opfm.mod, [g in 1:length(opfm.ctrl_generation)], 
                0 <= opfm.mod[:pg0][g] + (opfm.mod[:pgu][g,c] - opfm.mod[:pgd][g,c]))
            @constraint(opfm.mod, [g in 1:length(opfm.ctrl_generation)], 
                opfm.mod[:pg0][g] + (opfm.mod[:pgu][g,c] - opfm.mod[:pgd][g,c]) <= p_lim[g])
        end
    end
    @info "Ramping restrictions on generators added"
end

""" Restrict load shedding to max shed amount of load per node """
function add_lim_load_shed!(opfm::OPFmodel, listPd = get_active_power.(opfm.demands); 
        max_shed::Real = 1.0)
    @constraint(opfm.mod, load_shed, 
            opfm.mod[:ls0] .<= listPd .* max_shed
        )
    @info "Restriction of load shedding to $(max_shed) of load per node: Base case"
end
function add_lim_load_shed_cont!(opfm::OPFmodel, islands::Vector, 
        list::Vector{CTypes}, 
        listPd = get_active_power.(opfm.demands); max_shed::Real = 1.0)
    add_lim_load_shed!(opfm, listPd, max_shed=max_shed)
    for (c,island) in enumerate(islands)
        if length(island) < length(opfm.nodes)
            for n in island
                @constraint(opfm.mod, [d in list[n].demands], opfm.mod[:lsc][d,c] <= listPd[d] * max_shed)
            end
            for n in setdiff(1:length(opfm.nodes), island)
                @constraint(opfm.mod, [d in list[n].demands], opfm.mod[:lsc][d,c] == listPd[d])
            end
        else
            @constraint(opfm.mod, [d in 1:length(opfm.demands)], opfm.mod[:lsc][d,c] <= listPd[d] * max_shed)
        end
    end
    @info "- After contingency, before corrective actions"
end
function add_lim_load_shed_ccont!(opfm::OPFmodel, islands::Vector, 
        list::Vector{CTypes}, 
        listPd = get_active_power.(opfm.demands); max_shed::Real = 1.0)
    add_lim_load_shed_cont!(opfm, islands, list, listPd, max_shed=max_shed)
    for (c,island) in enumerate(islands)
        if length(island) < length(opfm.nodes)
            for n in island
                @constraint(opfm.mod, [d in list[n].demands], opfm.mod[:lscc][d,c] <= listPd[d] * max_shed)
            end
            for n in setdiff(1:length(opfm.nodes), island)
                @constraint(opfm.mod, [d in list[n].demands], opfm.mod[:lscc][d,c] == listPd[d])
            end
        else
            @constraint(opfm.mod, [d in 1:length(opfm.demands)], opfm.mod[:lscc][d,c] <= listPd[d] * max_shed)
        end
    end
    @info "- After contingency and corrective actions"
end

""" Restrict load shedding to renewable production per node """
function add_lim_renewable_shed!(opfm::OPFmodel, listPr = get_active_power.(opfm.renewables); 
        max_curtail::Real = 1.0)
    @constraint(opfm.mod, renew_shed, 
            opfm.mod[:pr0] .<= listPr .* max_curtail
        )
    @info "Restriction of renewable shedding: Base case"
end
function add_lim_renewable_shed_cont!(opfm::OPFmodel, islands::Vector, 
        list::Vector{CTypes}, listPr = get_active_power.(opfm.renewables); max_curtail::Real = 1.0)
    add_lim_renewable_shed!(opfm, listPr, max_curtail=max_curtail)
    for (c,island) in enumerate(islands)
        if length(island) < length(opfm.nodes)
            for n in island
                @constraint(opfm.mod, [d in list[n].renewables], opfm.mod[:prc][d,c] <= listPr[d] * max_curtail)
            end
            for n in setdiff(1:length(opfm.nodes), island)
                @constraint(opfm.mod, [d in list[n].renewables], opfm.mod[:prc][d,c] == listPr[d])
            end
        else
            @constraint(opfm.mod, [d in 1:length(opfm.renewables)], opfm.mod[:prc][d,c] <= listPr[d] * max_curtail)
        end
    end
    @info "- After contingency, before corrective actions"
end
function add_lim_renewable_shed_ccont!(opfm::OPFmodel, islands::Vector, 
        list::Vector{CTypes}, listPr = get_active_power.(opfm.renewables); max_curtail::Real = 1.0)
    add_lim_renewable_shed_cont!(opfm, islands, list, listPr, max_curtail=max_curtail)
    for (c,island) in enumerate(islands)
        if length(island) < length(opfm.nodes)
            for n in island
                @constraint(opfm.mod, [d in list[n].renewables], opfm.mod[:prcc][d,c] <= listPr[d] * max_curtail)
            end
            for n in setdiff(1:length(opfm.nodes), island)
                @constraint(opfm.mod, [d in list[n].renewables], opfm.mod[:prcc][d,c] == listPr[d])
            end
        else
            @constraint(opfm.mod, [d in 1:length(opfm.demands)], opfm.mod[:prcc][d,c] <= listPr[d] * max_curtail)
        end
    end
    @info "- After contingency and corrective actions"
end

