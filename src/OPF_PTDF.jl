# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

using PowerSystems
using JuMP

""" Run an OPF of a power system """
function opf(type::OPF, system::System, optimizer; 
            voll = nothing, 
            contingencies = nothing, 
            prob = nothing,
            time_limit_sec::Int64 = 600,
            unit_commit::Bool = false,
            ramp_minutes = 10, 
            ramp_mult = 10, 
            max_shed::Float64 = 1.0,
            max_curtail::Float64 = 1.0,
            short_term_limit_multi::Real = 1.5, long_term_limit_multi::Real = 1.2, 
            debug=false
        )
    contingencies = isnothing(contingencies) ? get_name.(sort_components!(get_branches(system))) : contingencies
    prob = isnothing(prob) ? make_prob(contingencies) : prob
    # set_renewable_prod!(system, renewable_prod)
    # set_active_power_demand!(system)

    opfm = isnothing(voll) ? opfmodel(system, optimizer, time_limit_sec, debug=debug) : opfmodel(system, optimizer, time_limit_sec, voll, contingencies, prob, debug=debug)
    idx = get_nodes_idx(opfm.nodes)
    list = make_list(opfm, idx, opfm.nodes)
    pf = DCPowerFlow(opfm.nodes, opfm.branches, idx)
    Pc = Dict{Int, NTuple{3, Any}}() # Holds the short term variables for contingencies
    Pcc = Dict{Int, NTuple{5, Any}}() # Holds the long term variables for contingencies

    (pg_lim_min, pg_lim_max) = split_pair(get_active_power_limits.(opfm.ctrl_generation))
    pr_lim = get_active_power.(opfm.renewables)
    (dc_lim_min, dc_lim_max) = split_pair(get_active_power_limits_from.(opfm.dc_branches))
    pd_lim = get_active_power.(opfm.demands)
    (rampup, rampdown) = split_pair(get_ramp_limits.(opfm.ctrl_generation))
    branch_rating = get_rate.(opfm.branches)

    @variables(opfm.mod, begin
            p0[n in 1:length(opfm.nodes)]
            0 <= pg0[g in 1:length(opfm.ctrl_generation)] <= pg_lim_max[g]
                # active power variables for the generators
            pfdc0[l in 1:length(opfm.dc_branches)]
                # power flow on DC branches
            0 <= ls0[d in 1:length(opfm.demands)]
                # demand curtailment variables
            0 <= pr0[d in 1:length(opfm.renewables)]
                # renewable curtailment variables
        end)
    
    @objective(opfm.mod, Min, opfm.cost_ctrl_gen' * pg0 + opfm.voll' * ls0)

    # @expression(opfm.mod, inj_p[n = 1:length(opfm.nodes)], 
    #     sum(beta(opfm.nodes[n], opfm.dc_branches[l]) * pfdc0[l] for l in list[n].dc_branches) +
    #     sum(pg0[g] for g in list[n].ctrl_generation) + 
    #     sum((get_active_power(opfm.renewables[d]) - pr0[d] for d in list[n].renewables), init = 0.0) - 
    #     sum((get_active_power(opfm.demands[d]) - ls0[d] for d in list[n].demands), init = 0.0)
    # )
    @expression(opfm.mod, inj_p, -p0)
    for n = 1:length(opfm.nodes)
        add_to_expression!.(inj_p[n], sum((pg0[g] for g in list[n].ctrl_generation), init = 0.0))
        add_to_expression!.(inj_p[n], sum((get_active_power(opfm.renewables[d]) for d in list[n].renewables), init = 0.0))
        add_to_expression!.(inj_p[n], -1, sum((get_active_power(opfm.demands[d]) for d in list[n].demands), init = 0.0))
    end
    @expression(opfm.mod, inj_p0, inj_p)
    for n = 1:length(opfm.nodes)
        add_to_expression!.(inj_p0[n], sum((beta(opfm.nodes[n], opfm.dc_branches[l]) * pfdc0[l] for l in list[n].dc_branches), init = 0.0))
        add_to_expression!.(inj_p0[n], -1, sum((pr0[d] for d in list[n].renewables), init = 0.0))
        add_to_expression!.(inj_p0[n], sum((ls0[d] for d in list[n].demands), init = 0.0))
    end
    @constraint(opfm.mod, inj_p0 .== 0)

    @constraint(opfm.mod, branch_lim_n, -branch_rating .<= pf.ϕ * p0)
    @constraint(opfm.mod, branch_lim_p, pf.ϕ * p0 .<= branch_rating)
    @constraint(opfm.mod, power_balance, sum(pg0, init=0.0) + sum(ls0, init=0.0) + sum(get_active_power.(opfm.renewables), init=0.0) == 
        sum(get_active_power.(opfm.demands), init=0.0) + sum(pr0, init=0.0))

    branch_rating_dc = get_active_power_limits_from.(opfm.dc_branches)
    @constraint(opfm.mod, pfdc0_lim_n[l = 1:length(opfm.dc_branches)], branch_rating_dc[l].min <= pfdc0[l])
    @constraint(opfm.mod, pfdc0_lim_p[l = 1:length(opfm.dc_branches)], pfdc0[l] <= branch_rating_dc[l].max)


    listPd = get_active_power.(opfm.demands)
    @constraint(opfm.mod, load_shed, ls0 .<= listPd .* max_shed)

    listPr = get_active_power.(opfm.renewables)
    !isempty(listPr) && @constraint(opfm.mod, renew_shed, pr0 .<= listPr .* max_curtail)
    if unit_commit
        add_unit_commit!(opfm)
    end

    if type != SC::OPF
        ptdf = copy(pf.ϕ)
        for (c,cont) in enumerate(get_bus_idx.(opfm.contingencies, [idx]))
            @info "Contingency $(cont[1])-$(cont[2]) is added"
            islands, island, island_b, ptdf = find_system_state(pf, cont, findfirst(x -> x == opfm.contingencies[c], opfm.branches), branch_rating, short_term_limit_multi, long_term_limit_multi, ptdf)
            add_short_term_contingencies(opfm, Pc, islands, island, ptdf, list, pr_lim, pd_lim, max_shed, branch_rating * short_term_limit_multi, cont, c)
            add_long_term_contingencies(opfm, Pcc, islands, island, ptdf, list, pr_lim, pd_lim, max_shed, branch_rating * long_term_limit_multi, ramp_mult, ramp_minutes, rampup, rampdown, dc_lim_min, dc_lim_max, pg_lim_max, cont, c)
        end
    end

    return opfm, Pc, Pcc
end

function add_short_term_contingencies(opfm, Pc, islands, island, ptdf, list, pr_lim, pd_lim, max_shed, branch_rating, cont, c)
    pc = JuMP.@variable(opfm.mod, [n in 1:length(opfm.nodes)], base_name = @sprintf("pc%s", c))
    pgc = JuMP.@variable(opfm.mod, [g in 1:length(opfm.ctrl_generation)], base_name = @sprintf("pgc%s", c), lower_bound = 0)
    prc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.renewables)], base_name = @sprintf( "prc%s", c), lower_bound = 0)
    lsc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.demands)], base_name = @sprintf( "lsc%s", c), lower_bound = 0)
    Pc[c] = (pgc, prc, lsc)
    
    obj = objective_function(opfm.mod)
    add_to_expression!.(obj, opfm.prob[c], sum(opfm.voll' * lsc))
    set_objective_function(opfm.mod, obj)

    # Add new constraints that limit the corrective variables within operating limits
    JuMP.@constraint(opfm.mod, sum(lsc) == sum(prc) + sum(pgc))
    if isempty(islands)
        JuMP.@constraint(opfm.mod, 0 .<= opfm.mod[:pg0] .- pgc)
        JuMP.@constraint(opfm.mod, prc .<= pr_lim .* max_shed)
        JuMP.@constraint(opfm.mod, lsc .<= pd_lim .* max_shed)
    else
        for n in islands[island]
            JuMP.@constraint(opfm.mod, [g = list[n].ctrl_generation], 
                0 <= opfm.mod[:pg0][g] - pgc[g])
            JuMP.@constraint(opfm.mod, [g = list[n].renewables], 
                prc[g] <= pr_lim[g] * max_shed)
            JuMP.@constraint(opfm.mod, [d = list[n].demands], 
                lsc[d] <= pd_lim[d] * max_shed)
        end
        for in_vec in islands[1:end .!= island]
            for n in in_vec
                JuMP.@constraint(opfm.mod, [g = list[n].ctrl_generation], 
                    opfm.mod[:pg0][g] == pgc[g])
                JuMP.@constraint(opfm.mod, [g = list[n].renewables], 
                    prc[g] == pr_lim[g])
                JuMP.@constraint(opfm.mod, [d = list[n].demands], 
                    lsc[d] == pd_lim[d])
            end
        end
    end
    inj_pc = @expression(opfm.mod, opfm.mod[:inj_p] .+ opfm.mod[:p0] .- pc)
    for n = 1:length(opfm.nodes)
        add_to_expression!.(inj_pc[n], -1, sum((pgc[g] for g in list[n].ctrl_generation), init = 0.0))
        add_to_expression!.(inj_pc[n], sum((beta(opfm.nodes[n], opfm.dc_branches[l]) * opfm.mod[:pfdc0][l] for l in list[n].dc_branches), init = 0.0))
        add_to_expression!.(inj_pc[n], sum((prc[d] for d in list[n].renewables), init = 0.0)) 
        add_to_expression!.(inj_pc[n], -1, sum((lsc[d] for d in list[n].demands), init = 0.0))
    end
    @constraint(opfm.mod, inj_pc .== 0)
    @constraint(opfm.mod, -branch_rating .<= ptdf * pc)
    @constraint(opfm.mod, ptdf * pc .<= branch_rating)
    # unregister(opfm.mod, :inj_p)
end

function add_long_term_contingencies(opfm, Pcc, islands, island, ptdf, list, pr_lim, pd_lim, max_shed, branch_rating, ramp_mult, ramp_minutes, rampup, rampdown, dc_lim_min, dc_lim_max, pg_lim_max, cont, c)
    pcc = JuMP.@variable(opfm.mod, [n in 1:length(opfm.nodes)], base_name = @sprintf("pcc%s", c))
    pgu = JuMP.@variable(opfm.mod, [g in 1:length(opfm.ctrl_generation)], base_name = @sprintf( "pgu%s", c), lower_bound = 0)
        # active power variables for the generators in contingencies ramp up 
    pgd = JuMP.@variable(opfm.mod, [g in 1:length(opfm.ctrl_generation)], base_name = @sprintf( "pgd%s", c), lower_bound = 0)
            # and ramp down
    pfdccc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.dc_branches)], base_name = @sprintf( "pfdccc%s", c))
    prcc = JuMP.@variable(opfm.mod, [r in 1:length(opfm.renewables)], base_name = @sprintf( "prcc%s", c), lower_bound = 0)
    lscc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.demands)], base_name = @sprintf( "lscc%s", c), lower_bound = 0)
            # load curtailment variables in in contingencies
    Pcc[c] = (pgu, pgd, pfdccc, prcc, lscc)

    # Extend the objective with the corrective variables
    obj = objective_function(opfm.mod)
    add_to_expression!.(obj, opfm.prob[c],
        # (1.0 - p_failure) * (sum(opfm.voll' * lscc) + sum(60 * (pgu + pgd)) # + # TODO: remove 60 and uncomment next lines for non-4-area analysis!!!!
        (sum(opfm.voll' * lscc) # + sum(opfm.cost_ctrl_gen' * ramp_mult * (pgu + pgd))
        # (sum(opfm.voll' * lscc) +
        # sum(opfm.cost_ctrl_gen' * ramp_mult * pgu) # +
        # sum(opfm.cost_ctrl_gen' * pgd)
        ) )
    add_to_expression!.(obj, opfm.prob[c], (sum(opfm.cost_ctrl_gen' * ramp_mult * pgu)) )
    add_to_expression!.(obj, opfm.prob[c], (sum(opfm.cost_ctrl_gen' * ramp_mult * pgd)) )
    set_objective_function(opfm.mod, obj)

    # Add new constraints that limit the corrective variables within operating limits
    JuMP.@constraint(opfm.mod, sum(pgu) + sum(lscc) == sum(pgd) + sum(prcc))
    if isempty(islands)
        JuMP.@constraint(opfm.mod, pgu .<= rampup * ramp_minutes)
        JuMP.@constraint(opfm.mod, pgd .<= rampdown * ramp_minutes)
        JuMP.@constraint(opfm.mod, 0 .<= opfm.mod[:pg0] .+ pgu .- pgd)
        JuMP.@constraint(opfm.mod, opfm.mod[:pg0] .+ pgu .- pgd .<= pg_lim_max)
        JuMP.@constraint(opfm.mod, dc_lim_min .<= pfdccc)
        JuMP.@constraint(opfm.mod, pfdccc .<= dc_lim_max)
        JuMP.@constraint(opfm.mod, prcc .<= pr_lim .* max_shed)
        JuMP.@constraint(opfm.mod, lscc .<= pd_lim .* max_shed)
    else
        for n in islands[island]
            JuMP.@constraint(opfm.mod, [g = list[n].ctrl_generation], 
                pgu[g] <= rampup[g] * ramp_minutes)
            JuMP.@constraint(opfm.mod, [g = list[n].ctrl_generation], 
                pgd[g] <= rampdown[g] * ramp_minutes)
            JuMP.@constraint(opfm.mod, [g = list[n].ctrl_generation], 
                0 <= opfm.mod[:pg0][g] + pgu[g] - pgd[g])
            JuMP.@constraint(opfm.mod, [g = list[n].ctrl_generation], 
                opfm.mod[:pg0][g] + pgu[g] - pgd[g] <= pg_lim_max[g])
            JuMP.@constraint(opfm.mod, [g = list[n].dc_branches], 
                dc_lim_min[g] <= pfdccc[g])
            JuMP.@constraint(opfm.mod, [g = list[n].dc_branches], 
                pfdccc[g] <= dc_lim_max[g])
            JuMP.@constraint(opfm.mod, [g = list[n].renewables], 
                prcc[g] <= pr_lim[g] * max_shed)
            JuMP.@constraint(opfm.mod, [d = list[n].demands], 
                lscc[d] <= pd_lim[d] * max_shed)
        end
        for in_vec in islands[1:end .!= island]
            for n in in_vec
                JuMP.@constraint(opfm.mod, [g = list[n].ctrl_generation], 
                    opfm.mod[:pg0][g] == pgd[g] - pgu[g])
                JuMP.@constraint(opfm.mod, [g = list[n].dc_branches], 
                    pfdccc[g] == 0)
                JuMP.@constraint(opfm.mod, [g = list[n].renewables], 
                    prcc[g] == pr_lim[g])
                JuMP.@constraint(opfm.mod, [d = list[n].demands], 
                    lscc[d] == pd_lim[d])
            end
        end
    end
    inj_pcc = @expression(opfm.mod, opfm.mod[:inj_p] .+ opfm.mod[:p0] .- pcc)
    for n = 1:length(opfm.nodes)
        add_to_expression!.(inj_pcc[n], sum((pgu[g] - pgd[g] for g in list[n].ctrl_generation), init = 0.0))
        add_to_expression!.(inj_pcc[n], sum((beta(opfm.nodes[n], opfm.dc_branches[l]) * opfm.mod[:pfdccc][l] for l in list[n].dc_branches), init = 0.0))
        add_to_expression!.(inj_pcc[n], sum((prcc[d] for d in list[n].renewables), init = 0.0)) 
        add_to_expression!.(inj_pcc[n], -1, sum((lscc[d] for d in list[n].demands), init = 0.0))
    end
    @constraint(opfm.mod, inj_pcc .== 0)
    @constraint(opfm.mod, -branch_rating .<= ptdf * pcc)
    @constraint(opfm.mod, ptdf * pcc .<= branch_rating)
end