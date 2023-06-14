# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022 

using PowerSystems
import JuMP
import Printf
import Gurobi
import MathOptInterface
const MOI = MathOptInterface
import LinearAlgebra
import SparseArrays

""" 
Solve the optimization model using Benders decomposition.

Function creates extra production constraints on generators in the system
based on the power transfer distribution factors of the generators in the
system. 
"""
function run_benders(type::OPF, system::System, voll=nothing, prob=nothing, contingencies=nothing; 
        ramp_minutes = 10, ramp_mult = 10, max_shed = 1.0, lim = 1e-6, 
        branch_short_term_limit_multi::Real = 1.5, branch_long_term_limit_multi::Real = 1.2, 
        gen_short_term_limit_multi::Real = 2.0, gen_long_term_limit_multi::Real = 1.0, p_failure = 0.0, branch_c = contingencies[1], rate_c = get_x(contingencies[1]))
    total_time = time()
    opfm = scopf(SC, system, Gurobi.Optimizer, voll=voll, prob=prob, contingencies=contingencies, ramp_minutes=ramp_minutes, max_shed=max_shed, 
        renewable_prod = 1.0)
    JuMP.@constraint(opfm.mod, opfm.mod[:pf0][findfirst(x -> x == branch_c, opfm.branches)] == rate_c)
    MOI.set(opfm.mod, MOI.Silent(), true) # supress output from the solver
    solve_model!(opfm.mod)
    termination_status(opfm.mod) != MOI.OPTIMAL && return opfm, [], [], [], []
    total_solve_time = solve_time(opfm.mod)
    @debug "lower_bound = $(objective_value(opfm.mod))"

    # Set variables
    idx = get_nodes_idx(opfm.nodes)
    list = make_list(opfm, idx, opfm.nodes)
    # cgen = connectivitymatrix(system, length(nodes), idx)
    linerating = get_rate.(opfm.branches)
    (pg_lim_min, pg_lim_max) = split_pair(get_active_power_limits.(opfm.ctrl_generation))
    pr_lim = get_active_power.(opfm.renewables)
    (dc_lim_min, dc_lim_max) = split_pair(get_active_power_limits_from.(opfm.dc_branches))
    pd_lim = get_active_power.(opfm.demands)
    (rampup, rampdown) = split_pair(get_ramp_limits.(opfm.ctrl_generation))
    Pc = Dict{Int, NTuple{3, Any}}() # Holds the short term variables for contingencies
    Pcc = Dict{Int, NTuple{5, Any}}() # Holds the long term variables for contingencies
    Pccx = Dict{Int, NTuple{2, Any}}() # Holds the long term variables for contingencies, no ramp up allowed
    # overloads = Dict{Int, Vector{Tuple{Int64, Float64}}}()
    island = Vector{Vector{Int64}}() 
    island_b = Vector{Vector{Int64}}()

    # Get initial state
    pf = SCOPF.DCPowerFlow(opfm.mod, opfm.nodes, opfm.branches, idx)
    ptdf = copy(pf.ϕ)
    Pg = get_controllable(opfm, idx)
    Pd = get_Pd(opfm, idx) # Fixed injection
    Pᵢ = calc_Pᵢ(pf)
    ΔPc = zeros(length(Pg))
    ΔPcc = zeros(length(Pg))
    ΔPccx = zeros(length(Pg))
    # print_contingency_overflow(opfm, pf, Pcc, 1.0)

    it = get_branch_bus_idx(opfm.branches, opfm.contingencies, idx)
    next = iterate(it)
    cut_added = 0
    iterations = 0
    while next !== nothing || cut_added > 0 # loops until no new cuts are added for the contingencies
        if next === nothing
            next = iterate(it)
            cut_added = 0
            iterations += 1
            @info "Iteration $iterations"
            # print_contingency_overflow(opfm, pf, Pcc, 1.0)
            if iterations >= length(opfm.branches)^2
                @error "Reached $(iterations) iterations without a stable solution."
                return opfm, pf, Pc, Pcc, Pccx
            end
        end

        # Calculate the power flow with the new outage and find if there are any overloads
        (c, cont), state = next
        islands, island, island_b, ptdf = find_overload(pf, cont, c, linerating, branch_short_term_limit_multi, branch_long_term_limit_multi, ptdf)
        overloads_c = filter_overload(ptdf * (Pᵢ .+ get_ΔPc(length(opfm.nodes), list, Pc, c)), linerating * branch_short_term_limit_multi)
        if type == PCSC::OPF
            overloads_cc = filter_overload(ptdf * (Pᵢ .+ get_ΔPcc(opfm, length(opfm.nodes), list, Pcc, c)), linerating * branch_long_term_limit_multi)
            if p_failure > 0.0
                overloads_ccx = filter_overload(ptdf * (Pᵢ .+ get_ΔPccx(opfm, length(opfm.nodes), list, Pccx, c)), linerating * branch_long_term_limit_multi)
            end
        end

        # Cannot change the model before all data is exctracted
        if !isempty(overloads_c)
            # overloads[c] = overload
            # add_contingency(opfm, opfm.branches[c].name, ramp_minutes=ramp_minutes, ramp_mult=ramp_mult, 
            #     max_shed=max_shed, lim=lim, branch_short_term_limit_multi=branch_short_term_limit_multi)
            cut_added = set_Pc(opfm, Pc, ΔPc, list, ptdf, Pg, overloads_c, islands, island, pr_lim, pd_lim, max_shed, cont, c, cut_added, lim)
        end
        if type == PCSC::OPF 
            if !isempty(overloads_cc) 
                cut_added = set_Pcc(opfm, Pcc, ΔPcc, list, ptdf, Pg, overloads_cc, islands, island, pr_lim, pd_lim, max_shed, ramp_mult, ramp_minutes, rampup, rampdown, dc_lim_min, dc_lim_max, pg_lim_max, p_failure, cont, c, cut_added, lim)
            end
            if p_failure > 0.0 && !isempty(overloads_ccx)
                cut_added = set_Pccx(opfm, Pccx, ΔPccx, list, ptdf, Pg, overloads_ccx, islands, island, pr_lim, pd_lim, max_shed, ramp_mult, ramp_minutes, rampup, rampdown, dc_lim_min, dc_lim_max, pg_lim_max, p_failure, cont, c, cut_added, lim)
            end
        end
        if !isempty(islands) && get(Pc, c, 0) == 0
            _ = get_Pc(opfm, Pc, islands, island, list, pr_lim, pd_lim, max_shed, cont, c)
            _ = get_Pcc(opfm, Pcc, islands, island, list, pr_lim, pd_lim, max_shed, ramp_mult, ramp_minutes, rampup, rampdown, dc_lim_min, dc_lim_max, pg_lim_max, p_failure, cont, c)
            _ = get_Pccx(opfm, Pccx, islands, island, list, pr_lim, pd_lim, max_shed, ramp_mult, ramp_minutes, rampup, rampdown, dc_lim_min, dc_lim_max, pg_lim_max, p_failure, cont, c)
            @info "Island: Contingency line $(cont[1])-$(cont[2])-i_$(c)"
        end
        if cut_added > 1
            Pᵢ, Pg = update_model!(opfm, pf, Pd, total_solve_time)
            cut_added = 1
        end
        termination_status(opfm.mod) != MOI.OPTIMAL && return opfm, pf, Pc, Pcc, Pccx
        next = iterate(it, state)
    end
    @printf "END: Total solve time %.4f. Total time %.4f.\n" total_solve_time (time() - total_time)
    return opfm, pf, Pc, Pcc, Pccx
end         

""" Solve model and update the power flow object """
function update_model!(opfm, pf, Pd, total_solve_time)
    solve_model!(opfm.mod)
    total_solve_time += solve_time(opfm.mod)

    Pᵢ = calc_Pᵢ(pf)
    Pg = Pᵢ - Pd
    set_θ!(pf, opfm.mod)
    calc_Pline!(pf)
    return Pᵢ, Pg
end

function find_overload(pf, cont, c, linerating, branch_short_term_limit_multi, branch_long_term_limit_multi, ptdf)
    islands = Vector{Vector{Int64}}() 
    island = 0
    island_b = Vector{Int64}()
    # lmt = min(branch_short_term_limit_multi, branch_long_term_limit_multi)
    try
        # flow = calculate_line_flows(pf, cont, c)
        # overloads = filter_overload(flow, linerating * lmt)
        # isempty(overloads) && return false, islands, island, island_b, ptdf
        ptdf = get_isf(pf, cont, c) # ptdf calculation is more computational expensive than line flow
    catch DivideError
        islands = island_detection(pf.B, cont[1], cont[2])
        island = find_ref_island(islands, pf.slack)
        island_b = find_island_branches(islands[island], pf.DA, c)
        fill!(ptdf, zero(eltype(ptdf)))
        ptdf[island_b, islands[island]] = get_isf(pf.DA, pf.B, cont, c, pf.slack, islands[island], island_b)
    end
    
    return islands, island, island_b, ptdf
end

""" Return the short term power injection change at each node. """
function get_ΔPc(numnodes::Integer, list, Pc, c)
    ΔP = zeros(numnodes)
    x = get(Pc, c, 0)
    if x != 0
        for (i, n) in enumerate(list)     
            for g in n.ctrl_generation
                ΔP[i] -= value(x[1][g]) # pgc
            end
            for g in n.renewables
                ΔP[i] -= value(x[2][g]) # prc
            end
            for g in n.demands
                ΔP[i] += value(x[3][g]) # lsc
            end
        end
    end
    return ΔP
end

function get_Pc(opfm::OPFmodel, Pc, islands, island, list, pr_lim, pd_lim, max_shed, cont, c)
    x = get(Pc, c, 0)
    if x == 0
        pgc = JuMP.@variable(opfm.mod, [g in 1:length(opfm.ctrl_generation)], base_name = "pgc", lower_bound = 0)
        prc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.renewables)], base_name = "prc", lower_bound = 0)
        lsc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.demands)], base_name = "lsc", lower_bound = 0)
        Pc[c] = (pgc, prc, lsc)

        @objective(opfm.mod, Min, objective_function(opfm.mod) + opfm.prob[c] * sum(opfm.voll' * lsc))

        # Add new constraints that limit the corrective variables within operating limits
        JuMP.@constraint(opfm.mod, sum(lsc) == sum(prc) + sum(pgc))
        if isempty(islands)
            JuMP.@constraint(opfm.mod, 0 .<= opfm.mod[:pg0] .- pgc)
            JuMP.@constraint(opfm.mod, 0 .<= prc .<= pr_lim .* max_shed)
            JuMP.@constraint(opfm.mod, 0 .<= lsc .<= pd_lim .* max_shed)
        else
            for in in islands[island]
                JuMP.@constraint(opfm.mod, [g = list[in].ctrl_generation], 
                    0 <= opfm.mod[:pg0][g] - pgc[g])
                JuMP.@constraint(opfm.mod, [g = list[in].renewables], 
                    0 <= prc[g] <= pr_lim[g] * max_shed)
                JuMP.@constraint(opfm.mod, [d = list[in].demands], 
                    0 <= lsc[d] <= pd_lim[d] * max_shed)
            end
            for in_vec in islands[1:end .!= island]
                for in in in_vec
                    JuMP.@constraint(opfm.mod, [g = list[in].ctrl_generation], 
                        opfm.mod[:pg0][g] == pgc[g])
                    JuMP.@constraint(opfm.mod, [g = list[in].renewables], 
                        prc[g] == pr_lim[g])
                    JuMP.@constraint(opfm.mod, [d = list[in].demands], 
                        lsc[d] == pd_lim[d])
                end
            end
        end
    else # If the contingency is run before, a set of variables already exist
        pgc = x[1]
        prc = x[2]
        lsc = x[3]
    end
    return pgc, prc, lsc
end

function set_Pc(opfm::OPFmodel, Pc, ΔPc, list, ptdf, Pg, overloads, islands, island, pr_lim, pd_lim, max_shed, cont, c, cut_added, lim)
    if !isempty(overloads)
        pgc, prc, lsc = get_Pc(opfm, Pc, islands, island, list, pr_lim, pd_lim, max_shed, cont, c)
        for (i,ol) in overloads
            expr = JuMP.@expression(opfm.mod, sum(
                    (ptdf[i, in] * (Pg[in] + ΔPc[in] - 
                    sum((opfm.mod[:pg0][ctrl] - pgc[ctrl] for ctrl in sublist.ctrl_generation), init=0.0) -
                    sum((beta(n, opfm.dc_branches[d]) * opfm.mod[:pfdc0][d] for d in sublist.dc_branches), init=0.0) +
                    sum((prc[r] for r in sublist.renewables), init=0.0) -
                    sum((lsc[d] for d in sublist.demands), init=0.0)
                    ) for (in, (n, sublist)) in enumerate(zip(opfm.nodes, list))), init=0.0
                ) )
            
            @info @sprintf "Pre: Contingency line %d-%d-i_%d; overload on %s of %.2f" cont[1] cont[2] c opfm.branches[i].name ol
            @debug "Cut added: $(sprint_expr(expr,lim))\n"
            # set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
            if ol < 0
                pre_cut = JuMP.@constraint(opfm.mod, expr <= ol)
            else
                pre_cut = JuMP.@constraint(opfm.mod, expr >= ol)
            end
            cut_added = 2
        end
    end
    return cut_added
end

""" Return the long term power injection change at each node. """
function get_ΔPcc(opfm::OPFmodel, numnodes::Integer, list, Pcc, c)
    ΔP = zeros(numnodes)
    x = get(Pcc, c, 0)
    if x != 0
        for (i, n) in enumerate(list)
            for g in n.ctrl_generation
                ΔP[i] += (value(x[1][g]) - value(x[2][g])) # pgu - pgd
            end
            for g in n.dc_branches
                ΔP[i] += beta(opfm.nodes[i], opfm.dc_branches[g]) * value(x[3][g]) # pfdccc
            end
            for g in n.renewables
                ΔP[i] -= value(x[4][g]) # prcc
            end
            for g in n.demands
                ΔP[i] += value(x[5][g]) # lscc
            end
        end
    end
    return ΔP
end

function get_Pcc(opfm::OPFmodel, Pcc, islands, island, list, pr_lim, pd_lim, max_shed, ramp_mult, ramp_minutes, rampup, rampdown, dc_lim_min, dc_lim_max, pg_lim_max, p_failure, cont, c)
    x = get(Pcc, c, 0)
    if x == 0

        # Add corrective variables
        pgu = JuMP.@variable(opfm.mod, [g in 1:length(opfm.ctrl_generation)], base_name = "pgu", lower_bound = 0)
            # active power variables for the generators in contingencies ramp up 
        pgd = JuMP.@variable(opfm.mod, [g in 1:length(opfm.ctrl_generation)], base_name = "pgd", lower_bound = 0)
                # and ramp down
        pfdccc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.dc_branches)], base_name = "pfdccc")
        prcc = JuMP.@variable(opfm.mod, [r in 1:length(opfm.renewables)], base_name = "prcc", lower_bound = 0)
        lscc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.demands)], base_name = "lscc", lower_bound = 0)
                # load curtailment variables in in contingencies
        Pcc[c] = (pgu, pgd, pfdccc, prcc, lscc)

        # Extend the objective with the corrective variables
        @objective(opfm.mod, Min, objective_function(opfm.mod) + opfm.prob[c] * 
            (1.0 - p_failure) * (sum(opfm.voll' * lscc) + sum(60 * (pgu + pgd)) # + # TODO: remove 60 and uncomment next lines for non-4-area analysis!!!!
            # (sum(opfm.voll' * lscc) +
            # sum(opfm.cost_ctrl_gen' * ramp_mult * pgu) # +
            # sum(opfm.cost_ctrl_gen' * pgd)
            ) )

        # Add new constraints that limit the corrective variables within operating limits
        JuMP.@constraint(opfm.mod, sum(pgu) + sum(lscc) == sum(pgd) + sum(prcc))
        if isempty(islands)
            JuMP.@constraint(opfm.mod, pgu .<= rampup * ramp_minutes)
            JuMP.@constraint(opfm.mod, pgd .<= rampdown * ramp_minutes)
            JuMP.@constraint(opfm.mod, 0 .<= opfm.mod[:pg0] .+ pgu .- pgd .<= pg_lim_max)
            JuMP.@constraint(opfm.mod, dc_lim_min .<= pfdccc .<= dc_lim_max)
            JuMP.@constraint(opfm.mod, prcc .<= pr_lim .* max_shed)
            JuMP.@constraint(opfm.mod, lscc .<= pd_lim .* max_shed)
        else
            for in in islands[island]
                JuMP.@constraint(opfm.mod,  [g = list[in].ctrl_generation], 
                    pgu[g] <= rampup[g] * ramp_minutes)
                JuMP.@constraint(opfm.mod,  [g = list[in].ctrl_generation], 
                    pgd[g] <= rampdown[g] * ramp_minutes)
                JuMP.@constraint(opfm.mod, [g = list[in].ctrl_generation], 
                    0 <= opfm.mod[:pg0][g] + pgu[g] - pgd[g] <= pg_lim_max[g])
                JuMP.@constraint(opfm.mod, [g = list[in].dc_branches], 
                    dc_lim_min[g] <= pfdccc[g] <= dc_lim_max[g])
                JuMP.@constraint(opfm.mod, [g = list[in].renewables], 
                    prcc[g] <= pr_lim[g] * max_shed)
                JuMP.@constraint(opfm.mod, [d = list[in].demands], 
                    lscc[d] <= pd_lim[d] * max_shed)
            end
            for in_vec in islands[1:end .!= island]
                for in in in_vec
                    JuMP.@constraint(opfm.mod, [g = list[in].ctrl_generation], 
                        opfm.mod[:pg0][g] == pgd[g] - pgu[g])
                    JuMP.@constraint(opfm.mod, [g = list[in].dc_branches], 
                        pfdccc[g] == 0)
                    JuMP.@constraint(opfm.mod, [g = list[in].renewables], 
                        prcc[g] == pr_lim[g])
                    JuMP.@constraint(opfm.mod, [d = list[in].demands], 
                        lscc[d] == pd_lim[d])
                end
            end
        end            

    else # If the contingency is run before, a set of corrective variables already exist
        pgu = x[1]
        pgd = x[2]
        pfdccc = x[3]
        prcc = x[4]
        lscc = x[5]
    end
    return pgu, pgd, pfdccc, prcc, lscc
end

function set_Pcc(opfm::OPFmodel, Pcc, ΔPcc, list, ptdf, Pg, overloads, islands, island, pr_lim, pd_lim, max_shed, ramp_mult, ramp_minutes, rampup, rampdown, dc_lim_min, dc_lim_max, pg_lim_max, p_failure, cont, c, cut_added, lim)
    if !isempty(overloads)
        pgu, pgd, pfdccc, prcc, lscc = get_Pcc(opfm, Pcc, islands, island, list, pr_lim, pd_lim, max_shed, ramp_mult, ramp_minutes, rampup, rampdown, dc_lim_min, dc_lim_max, pg_lim_max, p_failure, cont, c)

        # sort!(overloads, rev = true, by = x -> abs(x[2]))
        # (i,ol) = first(overloads)
        for (i,ol) in overloads
            # Finding and adding the Benders cut
            expr = JuMP.@expression(opfm.mod, sum((ptdf[i, in] * (
                    Pg[in] + ΔPcc[in] -
                    sum((opfm.mod[:pg0][ctrl] + pgu[ctrl] - pgd[ctrl] 
                        for ctrl in sublist.ctrl_generation), init=0.0) -
                    sum((beta(n, opfm.dc_branches[d]) * pfdccc[d] for d in sublist.dc_branches), init=0.0) +
                    sum((prcc[r] for r in sublist.renewables), init=0.0) -
                    sum((lscc[d] for d in sublist.demands), init=0.0)
                    ) for (in, (n, sublist)) in enumerate(zip(opfm.nodes, list))), init=0.0
                ) )
            
            @info @sprintf "Corr: Contingency line %d-%d-i_%d; overload on %s of %.4f" cont[1] cont[2] c opfm.branches[i].name ol
            @debug "Cut added: $(sprint_expr(expr,lim))\n"
            # set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
            if ol < 0
                corr_cut = JuMP.@constraint(opfm.mod, expr <= ol)
            else
                corr_cut = JuMP.@constraint(opfm.mod, expr >= ol)
            end
            cut_added = 3
        end
    end
    return cut_added
end

""" Return the long term power injection change at each node. """
function get_ΔPccx(opfm::OPFmodel, numnodes::Integer, list, Pccx, c)
    ΔP = zeros(numnodes)
    x = get(Pccx, c, 0)
    if x != 0
        for (i, n) in enumerate(list)
            for g in n.ctrl_generation
                ΔP[i] -= value(x[1][g]) # - pgd
            end
            for g in n.demands
                ΔP[i] += value(x[2][g]) # lscc
            end
        end
    end
    return ΔP
end

function get_Pccx(opfm::OPFmodel, Pccx, islands, island, list, pr_lim, pd_lim, max_shed, ramp_mult, ramp_minutes, rampup, rampdown, dc_lim_min, dc_lim_max, pg_lim_max, p_failure, cont, c)
    x = get(Pccx, c, 0)
    if x == 0

        # Add corrective variables
        pgd = JuMP.@variable(opfm.mod, [g in 1:length(opfm.ctrl_generation)], base_name = "pgd", lower_bound = 0)
                # and ramp down
        lscc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.demands)], base_name = "lscc", lower_bound = 0)
                # load curtailment variables in in contingencies
        Pccx[c] = (pgd, lscc)

        # Extend the objective with the corrective variables
        @objective(opfm.mod, Min, objective_function(opfm.mod) + opfm.prob[c] * 
            (sum(opfm.voll' * lscc) + sum(60 * pgd)) * p_failure
            )

        # Add new constraints that limit the corrective variables within operating limits
        JuMP.@constraint(opfm.mod, sum(lscc) == sum(pgd))
        if isempty(islands)
            JuMP.@constraint(opfm.mod, pgd .<= rampdown * ramp_minutes)
            JuMP.@constraint(opfm.mod, 0 .<= opfm.mod[:pg0] .- pgd)
            JuMP.@constraint(opfm.mod, lscc .<= pd_lim .* max_shed)
        else
            for in in islands[island]
                JuMP.@constraint(opfm.mod,  [g = list[in].ctrl_generation], 
                    pgd[g] <= rampdown[g] * ramp_minutes)
                JuMP.@constraint(opfm.mod, [g = list[in].ctrl_generation], 
                    0 <= opfm.mod[:pg0][g] - pgd[g])
                JuMP.@constraint(opfm.mod, [d = list[in].demands], 
                    lscc[d] <= pd_lim[d] * max_shed)
            end
            for in_vec in islands[1:end .!= island]
                for in in in_vec
                    JuMP.@constraint(opfm.mod, [g = list[in].ctrl_generation], 
                        opfm.mod[:pg0][g] == pgd[g])
                    JuMP.@constraint(opfm.mod, [d = list[in].demands], 
                        lscc[d] == pd_lim[d])
                end
            end
        end            

    else # If the contingency is run before, a set of corrective variables already exist
        pgd = x[1]
        lscc = x[2]
    end
    return pgd, lscc
end

function set_Pccx(opfm::OPFmodel, Pccx, ΔPccx, list, ptdf, Pg, overloads, islands, island, pr_lim, pd_lim, max_shed, ramp_mult, ramp_minutes, rampup, rampdown, dc_lim_min, dc_lim_max, pg_lim_max, p_failure, cont, c, cut_added, lim)
    if !isempty(overloads)
        pgd, lscc = get_Pccx(opfm, Pccx, islands, island, list, pr_lim, pd_lim, max_shed, ramp_mult, ramp_minutes, rampup, rampdown, dc_lim_min, dc_lim_max, pg_lim_max, p_failure, cont, c)
        # sort!(overloads, rev = true, by = x -> abs(x[2]))
        # (i,ol) = first(overloads)
        for (i,ol) in overloads
            # Finding and adding the Benders cut
            expr = JuMP.@expression(opfm.mod, sum((ptdf[i, in] * (
                    Pg[in] + ΔPccx[in] -
                    sum((opfm.mod[:pg0][ctrl] - pgd[ctrl] 
                        for ctrl in sublist.ctrl_generation), init=0.0)  -
                    sum((lscc[d] for d in sublist.demands), init=0.0)
                    ) for (in, (n, sublist)) in enumerate(zip(opfm.nodes, list))), init=0.0
                ) )
            
            @info @sprintf "Corr.x: Contingency line %d-%d-i_%d; overload on %s of %.4f" cont[1] cont[2] c opfm.branches[i].name ol
            @debug "Cut added: $(sprint_expr(expr,lim))\n"
            # set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
            if ol < 0
                corr_cut = JuMP.@constraint(opfm.mod, expr <= ol)
            else
                corr_cut = JuMP.@constraint(opfm.mod, expr >= ol)
            end
            cut_added = 3
        end
    end
    return cut_added
end

" An AffExpr nicely formatted to a string "
sprint_expr(expr::AffExpr, lim = 1e-6) = 
    join(Printf.@sprintf("%s%5.2f %s ", (x[2] > 0 ? "+" : "-"), abs(x[2]), x[1]) 
            for x in expr.terms if abs(x[2]) > lim) * 
        Printf.@sprintf("<= %s%.2f", (expr.constant > 0 ? "-" : " "), abs(expr.constant)
    )

function print_benders_results(opfm::OPFmodel, Pc::Dict{Int, NTuple{3, Any}} = Dict{Int, NTuple{3, Any}}(), 
        Pcc::Dict{Int, NTuple{5, Any}} = Dict{Int, NTuple{5, Any}}(), 
        Pccx::Dict{Int, NTuple{2, Any}} = Dict{Int, NTuple{2, Any}}(), lim::Real = 1e-6)
    function print_c(itr, val, symb, i_g, lim)
        for (i,c) in itr
            if JuMP.value(c[val][i_g]) > lim
                @printf("          c %12s: %s: %.3f\n", opfm.branches[i].name, symb, JuMP.value(c[val][i_g]))
            end
        end
    end
    for (i_g,g) in enumerate(opfm.ctrl_generation)
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(opfm.mod[:pg0][i_g]), get_active_power_limits(g).max)
        print_c(Pc, 1, "pgc", i_g, lim)
        print_c(Pcc, 1, "pgu", i_g, lim)
        print_c(Pcc, 2, "pgd", i_g, lim)
        print_c(Pccx, 1, "pgdx", i_g, lim)
    end
    for (i_g,g) in enumerate(opfm.dc_branches)
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(opfm.mod[:pfdc0][i_g]), get_active_power_limits(g).max)
        print_c(Pcc, 3, "pfdccc", i_g, lim)
    end
    for (i_g,g) in enumerate(opfm.renewables)
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(opfm.mod[:pr0][i_g]), get_active_power_limits(g).max)
        print_c(Pc, 2, "prc", i_g, lim)
        print_c(Pcc, 4, "prcc", i_g, lim)
    end
    for (i_g,g) in enumerate(opfm.demands)
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(opfm.mod[:ls0][i_g]), get_active_power(g))
        print_c(Pc, 3, "lsc", i_g, lim)
        print_c(Pcc, 5, "lscc", i_g, lim)
        print_c(Pccx, 2, "lsccx", i_g, lim)
    end
end

function add_contingency(opfm, c_name; ramp_minutes = 10, ramp_mult = 10, max_shed = 0.1, lim = 1e-6, branch_short_term_limit_multi::Real = 1.5)

    pgu = JuMP.@variable(opfm.mod, [g in get_name.(get_ctrl_generation(opfm.sys))], base_name = "pgu", lower_bound = 0)
    # active power variables for the generators in contingencies ramp up 
    pgd = JuMP.@variable(opfm.mod, [g in get_name.(get_ctrl_generation(opfm.sys))], base_name = "pgd", lower_bound = 0)
            # and ramp down
    pfdccc = JuMP.@variable(opfm.mod, [d in get_name.(get_dc_branches(opfm.sys))], base_name = "pfdccc")
    prcc = JuMP.@variable(opfm.mod, [r in get_name.(get_renewables(opfm.sys))], base_name = "prcc", lower_bound = 0)
    lscc = JuMP.@variable(opfm.mod, [d in get_name.(get_demands(opfm.sys))], base_name = "lscc", lower_bound = 0)
            # load curtailment variables in in contingencies
    pfcc = JuMP.@variable(opfm.mod, [l in get_name.(get_branches(opfm.sys))])         # and after corrective actions
    vacc = JuMP.@variable(opfm.mod, [b in get_name.(get_nodes(opfm.sys))])        # and after corrective actions
    for g in get_ctrl_generation(opfm.sys)
        set_upper_bound(pgu[get_name(g)], get_ramp_limits(g).up * ramp_minutes)
        set_upper_bound(pgd[get_name(g)], get_ramp_limits(g).down * ramp_minutes)
    end

    @objective(opfm.mod, Min, objective_function(opfm.mod) + opfm.prob[c_name] * 
        (sum(opfm.voll[d] * (# ramp_minutes / 60 * opfm.mod[:lsc][d,c] +
            lscc[d]) for d in get_name.(get_demands(opfm.sys))
            ) +
        sum(opfm.voll[r] * prcc[r] for r in get_name.(get_renewables(opfm.sys))) +
        sum(opfm.cost[g] * ramp_mult * pgu[g] for g in get_name.(get_ctrl_generation(opfm.sys)))
        ) )

    @constraint(opfm.mod, [n = get_name.(get_nodes(opfm.sys))],
        sum(g.bus.name == n ? opfm.mod[:pg0][get_name(g)] + pgu[get_name(g)] - pgd[get_name(g)] : 0 for g in get_ctrl_generation(opfm.sys)) -
        sum(beta(n,l) * pfcc[get_name(l)] for l in get_branches(opfm.sys)) -
        sum(beta(n,l) * pfdccc[get_name(l)] for l in get_dc_branches(opfm.sys)) == 
        sum(d.bus.name == n ? get_active_power(d) - lscc[get_name(d)] : 0 for d in get_demands(opfm.sys)) +
        sum(d.bus.name == n ? -get_active_power(d) + prcc[get_name(d)] : 0 for d in get_renewables(opfm.sys))
    )
    i, slack = find_slack(opfm.sys)
    @constraint(opfm.mod, vacc[get_name(slack)] == 0)

    branch_rating = make_named_array(get_rate, get_branches(opfm.sys))
    @constraint(opfm.mod, [l = get_name.(get_branches(opfm.sys))], 
        -branch_rating[l] .<= pfcc[l] .<= branch_rating[l]
    )
    x = make_named_array(get_x, get_branches(opfm.sys))
    @constraint(opfm.mod, [l = get_name.(get_branches(opfm.sys))],
        pfcc[l] .- sum(beta(opfm.sys,l) .* vacc[:]) ./ x[l] .== 0
    )
    branch_rating = make_named_array(get_active_power_limits_from, get_dc_branches(opfm.sys))
    @constraint(opfm.mod, [l = get_name.(get_dc_branches(opfm.sys)), c = opfm.contingencies], 
        branch_rating[l].min <= opfm.mod[:pfdccc][l,c] <= branch_rating[l].max
    )

    for g in get_ctrl_generation(opfm.sys)
        g_name = get_name(g)
        @constraint(opfm.mod, 0 <= opfm.mod[:pg0][g_name] + (pgu[g_name] - pgd[g_name]) <= 
            get_active_power_limits(g).max)
    end

    for l in get_demands(opfm.sys)
        @constraint(opfm.mod, lscc[get_name(l)] <= get_active_power(l) * max_shed)
    end
    for l in get_renewables(opfm.sys)
        @constraint(opfm.mod, prcc[get_name(l)] <= get_active_power(l))
    end
    @info "Added constraints for contingency on line $c_name"
end
