# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

using PowerSystems
using JuMP

""" Run an OPF of a power system """
function opf(system::System, optimizer; 
            voll = nothing, 
            contingencies = nothing, 
            prob = nothing,
            time_limit_sec::Int64 = 600,
            unit_commit::Bool = false,
            max_shed::Float64 = 0.1,
            max_curtail::Float64 = 1.0,
            renewable_prod::Float64= 0.5
        )
    contingencies = isnothing(contingencies) ? get_name.(get_sorted_branches(system)) : contingencies
    prob = isnothing(prob) ? make_prob(contingencies) : prob
    set_renewable_prod!(system, renewable_prod)
    set_active_power_demand!(system)

    opfm = isnothing(voll) ? opfmodel(system, optimizer, time_limit_sec) : opfmodel(system, optimizer, time_limit_sec, voll, contingencies, prob)
    idx = get_nodes_idx(opfm.nodes)
    list = make_list(opfm, idx, opfm.nodes)
    isf = get_isf(opfm.branches, opfm.nodes)

    p_lim = get_active_power_limits.(opfm.ctrl_generation)
    @variables(opfm.mod, begin
            0 <= pg0[g in 1:length(opfm.ctrl_generation)] <= p_lim[g].max     
                # active power variables for the generators
            pfdc0[l in 1:length(opfm.dc_branches)]
                # power flow on DC branches
            0 <= ls0[d in 1:length(opfm.demands)]
                # demand curtailment variables
            0 <= pr0[d in 1:length(opfm.renewables)]
                # renewable curtailment variables
        end)
    
    branch_rating = get_rate.(opfm.branches)
    @expression(opfm.mod, inj_p[n = 1:length(opfm.nodes)], 
        sum(beta(opfm.nodes[n], opfm.dc_branches[l]) * pfdc0[l] for l in list[n].dc_branches) +
        sum(pg0[g] for g in list[n].ctrl_generation) + 
        sum((get_active_power(opfm.renewables[d]) - pr0[d] for d in list[n].renewables), init = 0.0) - 
        sum((get_active_power(opfm.demands[d]) - ls0[d] for d in list[n].demands), init = 0.0)
    )
    
    # @expression(opfm.mod, inj_p[n = 1:length(opfm.nodes)], 0) 
    # for n = 1:length(opfm.nodes)
    #     for l in list[n].dc_branches
    #         add_to_expression!(inj_p[n], beta(opfm.nodes[n], opfm.dc_branches[l]) * pfdc0[l])
    #     end
    #     for g in list[n].ctrl_generation
    #         add_to_expression!(inj_p[n], pg0[g])
    #     end
    #     for d in list[n].renewables
    #         add_to_expression!(inj_p[n], get_active_power(opfm.renewables[d]) - pr0[d])
    #     end
    #     for d in list[n].demands
    #         add_to_expression!(inj_p[n], -get_active_power(opfm.demands[d]) + ls0[d])
    #     end
    # end

    @constraint(opfm.mod, branch_lim, -branch_rating .<= isf * inj_p .<= branch_rating)
    @constraint(opfm.mod, sum(pg0, init=0.0) + sum(ls0, init=0.0) + sum(get_active_power.(opfm.renewables), init=0.0) == 
        sum(get_active_power.(opfm.demands), init=0.0) + sum(pr0, init=0.0))

    branch_rating = get_active_power_limits_from.(opfm.dc_branches)
    @constraint(opfm.mod, pfdc0_lim[l = 1:length(opfm.dc_branches)], 
            branch_rating[l].min <= pfdc0[l] <= branch_rating[l].max
        )

    @objective(opfm.mod, Min, opfm.cost_ctrl_gen' * pg0 + opfm.voll' * ls0)

    listPd = get_active_power.(opfm.demands)
    @constraint(opfm.mod, load_shed, ls0 .<= listPd .* max_shed)

    listPr = get_active_power.(opfm.renewables)
    !isempty(listPr) && @constraint(opfm.mod, renew_shed, pr0 .<= listPr .* max_curtail)
    if unit_commit
        add_unit_commit!(opfm)
    end
    return opfm
end
