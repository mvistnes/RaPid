# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022 

using PowerSystems
import JuMP
import Printf
import Gurobi
import MathOptInterface
const MOI = MathOptInterface
import LinearAlgebra
import SparseArrays

""" Benders type """
mutable struct Benders
    idx::Dict{<:Int, <:Int}
    list::Vector{CTypes}

    linerating::Vector
    pg_lim_min::Vector
    pg_lim_max::Vector
    pr_lim::Vector
    dc_lim_min::Vector
    dc_lim_max::Vector
    pd_lim::Vector
    rampup::Vector
    rampdown::Vector

    Pc::Dict{Int, NTuple{3, Any}} # Holds the short term variables for contingencies
    Pcc::Dict{Int, NTuple{5, Any}} # Holds the long term variables for contingencies
    Pccx::Dict{Int, NTuple{2, Any}} # Holds the long term variables for contingencies, no ramp up allowed

    Pᵢ::Vector{<:Real}
    Pg::Vector{<:Real}
    Pd::Vector{<:Real}
end

""" Constructor for OPFmodel """
function benders(opfm::OPFmodel)
    idx = get_nodes_idx(opfm.nodes)
    list = make_list(opfm, idx, opfm.nodes)
    # cgen = connectivitymatrix(system, length(nodes), idx)
    linerating = get_rate.(opfm.branches)
    (pg_lim_min, pg_lim_max) = split_pair(get_active_power_limits.(opfm.ctrl_generation))
    pr_lim = get_active_power.(opfm.renewables)
    (dc_lim_min, dc_lim_max) = split_pair(get_active_power_limits_from.(opfm.dc_branches))
    pd_lim = get_active_power.(opfm.demands)
    (rampup, rampdown) = split_pair(get_ramp_limits.(opfm.ctrl_generation))
    Pc = Dict{Int, NTuple{3, Any}}() 
    Pcc = Dict{Int, NTuple{5, Any}}() 
    Pccx = Dict{Int, NTuple{2, Any}}() 
    Pᵢ = get_net_Pᵢ(opfm, idx)
    Pg = get_controllable(opfm, idx)
    Pd = get_Pd(opfm, idx) # Fixed injection
    @assert isapprox(Pg, (Pᵢ - Pd); atol=1e-5) string(Pg - (Pᵢ - Pd))
    return Benders(idx, list, linerating, pg_lim_min, pg_lim_max, pr_lim, dc_lim_min, dc_lim_max, pd_lim, rampup, rampdown, Pc, Pcc, Pccx, Pᵢ, Pg, Pd)
end

""" 
Solve the optimization model using Benders decomposition.

Function creates extra production constraints on generators in the system
based on the power transfer distribution factors of the generators in the
system. 
"""
function run_benders(type::OPF, system::System, voll=nothing, prob=nothing, contingencies=nothing; 
        ramp_minutes = 10, ramp_mult = 10, max_shed = 1.0, lim = 1e-6, 
        branch_short_term_limit_multi::Real = 1.5, branch_long_term_limit_multi::Real = 1.2, 
        gen_short_term_limit_multi::Real = 2.0, gen_long_term_limit_multi::Real = 1.0, p_failure = 0.0, branch_c = nothing, rate_c = 0.0, debug = false)
    total_time = time()
    # opfm = scopf(SC, system, Gurobi.Optimizer, voll=voll, prob=prob, contingencies=contingencies, ramp_minutes=ramp_minutes, max_shed=max_shed, 
    #     renewable_prod = 1.0, debug = debug)
    opfm, Pc_ptdf, Pcc_ptdf = SCOPF.opf(SC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed)
    # if !isnothing(branch_c) 
    #     @info "Flow constraint on branch $branch_c is $rate_c."
    #     JuMP.@constraint(opfm.mod, opfm.mod[:pf0][findfirst(x -> x == branch_c, opfm.branches)] == rate_c)
    # end
    MOI.set(opfm.mod, MOI.Silent(), true) # supress output from the solver
    solve_model!(opfm.mod)
    termination_status(opfm.mod) != MOI.OPTIMAL && return opfm, [], [], [], []
    total_solve_time = solve_time(opfm.mod)
    @debug "lower_bound = $(objective_value(opfm.mod))"

    # Set variables
    bd = benders(opfm)
    # overloads = Dict{Int, Vector{Tuple{Int64, Float64}}}()
    island = Vector{Vector{Int64}}() 
    island_b = Vector{Vector{Int64}}()

    pf = SCOPF.DCPowerFlow(opfm.nodes, opfm.branches, bd.Pᵢ, bd.idx)
    ptdf = copy(pf.ϕ)
    ΔPc = zeros(length(bd.Pg))
    ΔPcc = zeros(length(bd.Pg))
    ΔPccx = zeros(length(bd.Pg))
    # print_contingency_overflow(opfm, pf, Pcc, 1.0)

    it = enumerate(get_bus_idx.(opfm.contingencies, [bd.idx]))
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
                return opfm, pf, bd.Pc, bd.Pcc, bd.Pccx
            end
        end

        (c, cont), state = next
        islands, island, island_b, ptdf = find_system_state(pf, cont, findfirst(x -> x == opfm.contingencies[c], opfm.branches), 
            bd.linerating, branch_short_term_limit_multi, branch_long_term_limit_multi, ptdf)
        if !isempty(islands) && get(bd.Pc, c, 0) == 0
            _ = get_Pc(opfm, bd, islands, island, max_shed, cont, c)
            if type == PCSC::OPF
                _ = get_Pcc(opfm, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c)
                if p_failure > 0.0
                    _ = get_Pccx(opfm, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c)
                end
            end
            @info "Island: Contingency line $(cont[1])-$(cont[2])-i_$(c)"
            bd.Pᵢ, bd.Pg, total_solve_time = update_model!(opfm, pf, bd.Pd, total_solve_time, bd.idx)
        end

        # Calculate the power flow with the new outage and find if there are any overloads
        ΔPc = get_ΔPc(opfm, bd.list, ΔPc, bd.Pc, c)
        overloads_c = filter_overload(ptdf * (bd.Pᵢ .+ ΔPc), bd.linerating * branch_short_term_limit_multi)
        if type == PCSC::OPF
            ΔPcc = get_ΔPcc(opfm, bd.list, ΔPcc, bd.Pcc, c)
            overloads_cc = filter_overload(ptdf * (bd.Pᵢ .+ ΔPcc), bd.linerating * branch_long_term_limit_multi)
            if p_failure > 0.0
                ΔPccx = get_ΔPccx(opfm, bd.list, ΔPccx, bd.Pccx, c)
                overloads_ccx = filter_overload(ptdf * (bd.Pᵢ .+ ΔPccx), bd.linerating * branch_long_term_limit_multi)
            end
        end

        # Cannot change the model before all data is exctracted!
        if !isempty(overloads_c)
            # overloads[c] = overload
            cut_added = set_Pc(opfm, bd, ΔPc, ptdf, overloads_c, islands, island, max_shed, cont, c, cut_added, lim)
        end
        if type == PCSC::OPF 
            if !isempty(overloads_cc) 
                cut_added = set_Pcc(opfm, bd, ΔPcc, ptdf, overloads_cc, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c, cut_added, lim)
            end
            if p_failure > 0.0 && !isempty(overloads_ccx)
                cut_added = set_Pccx(opfm, bd, ΔPccx, ptdf, overloads_ccx, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c, cut_added, lim)
            end
        end
        if cut_added > 1
            bd.Pᵢ, bd.Pg, total_solve_time = update_model!(opfm, pf, bd.Pd, total_solve_time, bd.idx)
            cut_added = 1
        end
        termination_status(opfm.mod) != MOI.OPTIMAL && return opfm, pf, bd.Pc, bd.Pcc, bd.Pccx
        next = iterate(it, state)
    end
    @printf "END: Total solve time %.4f. Total time %.4f.\n" total_solve_time (time() - total_time)
    return opfm, pf, bd.Pc, bd.Pcc, bd.Pccx
end         

""" Solve model and update the power flow object """
function update_model!(opfm, pf, Pd, total_solve_time, idx)
    solve_model!(opfm.mod)
    total_solve_time += solve_time(opfm.mod)

    # Pᵢ = calc_Pᵢ(pf)
    Pᵢ = get_net_Pᵢ(opfm, idx)
    Pg = Pᵢ - Pd
    calc_θ!(pf, Pᵢ)
    calc_Pline!(pf)
    return Pᵢ, Pg, total_solve_time
end

function find_system_state(pf, cont, c, linerating, branch_short_term_limit_multi, branch_long_term_limit_multi, ptdf)
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
function get_ΔPc(opfm::OPFmodel, list, ΔP, Pc, c)
    fill!(ΔP, zero(eltype(ΔP)))
    x = get(Pc, c, 0)
    if x != 0
        pgc = get_value(opfm.mod, x[1])
        prc = get_value(opfm.mod, x[2])
        lsc = get_value(opfm.mod, x[3])
        for (i, n) in enumerate(list)
            for g in n.ctrl_generation
                ΔP[i] -= pgc[g]
            end
            for g in n.renewables
                ΔP[i] -= prc[g]
            end
            for g in n.demands
                ΔP[i] += lsc[g]
            end
        end
    end
    return ΔP
end

function get_Pc(opfm::OPFmodel, bd, islands, island, max_shed, cont, c)
    x = get(bd.Pc, c, 0)
    if x == 0
        pgc = JuMP.@variable(opfm.mod, [g in 1:length(opfm.ctrl_generation)], base_name = @sprintf("pgc%s", c), lower_bound = 0)
        prc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.renewables)], base_name = @sprintf( "prc%s", c), lower_bound = 0)
        lsc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.demands)], base_name = @sprintf( "lsc%s", c), lower_bound = 0)
        bd.Pc[c] = (pgc, prc, lsc)

        set_objective_function(opfm.mod, objective_function(opfm.mod) + opfm.prob[c] * sum(opfm.voll' * lsc))

        # Add new constraints that limit the corrective variables within operating limits
        JuMP.@constraint(opfm.mod, sum(lsc) == sum(prc) + sum(pgc))
        if isempty(islands)
            JuMP.@constraint(opfm.mod, 0 .<= opfm.mod[:pg0] .- pgc)
            JuMP.@constraint(opfm.mod, prc .<= bd.pr_lim .* max_shed)
            JuMP.@constraint(opfm.mod, lsc .<= bd.pd_lim .* max_shed)
        else
            for n in islands[island]
                JuMP.@constraint(opfm.mod, [g = bd.list[n].ctrl_generation], 
                    0 <= opfm.mod[:pg0][g] - pgc[g])
                JuMP.@constraint(opfm.mod, [g = bd.list[n].renewables], 
                    prc[g] <= bd.pr_lim[g] * max_shed)
                JuMP.@constraint(opfm.mod, [d = bd.list[n].demands], 
                    lsc[d] <= bd.pd_lim[d] * max_shed)
            end
            for in_vec in islands[1:end .!= island]
                for n in in_vec
                    JuMP.@constraint(opfm.mod, [g = bd.list[n].ctrl_generation], 
                        opfm.mod[:pg0][g] == pgc[g])
                    JuMP.@constraint(opfm.mod, [g = bd.list[n].renewables], 
                        prc[g] == bd.pr_lim[g])
                    JuMP.@constraint(opfm.mod, [d = bd.list[n].demands], 
                        lsc[d] == bd.pd_lim[d])
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

function set_Pc(opfm::OPFmodel, bd, ΔPc, ptdf, overloads, islands, island, max_shed, cont, c, cut_added, lim)
    if !isempty(overloads)
        pgc, prc, lsc = get_Pc(opfm, bd, islands, island, max_shed, cont, c)
        for (i,ol) in overloads
            expr = JuMP.@expression(opfm.mod, sum(
                    (ptdf[i, inode] * (bd.Pg[inode] + ΔPc[inode] - 
                    sum((opfm.mod[:pg0][ctrl] - pgc[ctrl] for ctrl in sublist.ctrl_generation), init=0.0) -
                    sum((beta(sublist.node, opfm.dc_branches[d]) * opfm.mod[:pfdc0][d] for d in sublist.dc_branches), init=0.0) +
                    sum((prc[r] for r in sublist.renewables), init=0.0) -
                    sum((lsc[d] for d in sublist.demands), init=0.0)
                    ) for (inode, sublist) in enumerate(bd.list)), init=0.0
                ) )
            
            @info @sprintf "Pre: Contingency line %d-%d-i_%d; overload on %s of %.4f" cont[1] cont[2] c opfm.branches[i].name ol
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
function get_ΔPcc(opfm::OPFmodel, list, ΔP, Pcc, c)
    fill!(ΔP, zero(eltype(ΔP)))
    x = get(Pcc, c, 0)
    if x != 0
        pgu = get_value(opfm.mod, x[1])
        pgd = get_value(opfm.mod, x[2])
        pfdccc = get_value(opfm.mod, x[3])
        prcc = get_value(opfm.mod, x[4])
        lscc = get_value(opfm.mod, x[5])
        for (i, n) in enumerate(list)
            for g in n.ctrl_generation
                ΔP[i] += (pgu[g] - pgd[g])
            end
            for g in n.dc_branches
                ΔP[i] += beta(opfm.nodes[i], opfm.dc_branches[g]) * pfdccc[g]
            end
            for g in n.renewables
                ΔP[i] -= prcc[g]
            end
            for g in n.demands
                ΔP[i] += lscc[g]
            end
        end
    end
    return ΔP
end

function get_Pcc(opfm::OPFmodel, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c)
    x = get(bd.Pcc, c, 0)
    if x == 0

        # Add corrective variables
        pgu = JuMP.@variable(opfm.mod, [g in 1:length(opfm.ctrl_generation)], base_name = @sprintf( "pgu%s", c), lower_bound = 0)
            # active power variables for the generators in contingencies ramp up 
        pgd = JuMP.@variable(opfm.mod, [g in 1:length(opfm.ctrl_generation)], base_name = @sprintf( "pgd%s", c), lower_bound = 0)
                # and ramp down
        pfdccc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.dc_branches)], base_name = @sprintf( "pfdccc%s", c))
        prcc = JuMP.@variable(opfm.mod, [r in 1:length(opfm.renewables)], base_name = @sprintf( "prcc%s", c), lower_bound = 0)
        lscc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.demands)], base_name = @sprintf( "lscc%s", c), lower_bound = 0)
                # load curtailment variables in in contingencies
        bd.Pcc[c] = (pgu, pgd, pfdccc, prcc, lscc)

        # Extend the objective with the corrective variables
        set_objective_function(opfm.mod, objective_function(opfm.mod) + opfm.prob[c] * 
            # (1.0 - p_failure) * (sum(opfm.voll' * lscc) + sum(60 * (pgu + pgd)) # + # TODO: remove 60 and uncomment next lines for non-4-area analysis!!!!
            (1.0 - p_failure) * (sum(opfm.voll' * lscc) + sum(opfm.cost_ctrl_gen' * ramp_mult * (pgu + pgd))
            # (sum(opfm.voll' * lscc) +
            # sum(opfm.cost_ctrl_gen' * ramp_mult * pgu) # +
            # sum(opfm.cost_ctrl_gen' * pgd)
            ) )

        # Add new constraints that limit the corrective variables within operating limits
        JuMP.@constraint(opfm.mod, sum(pgu) + sum(lscc) == sum(pgd) + sum(prcc))
        if isempty(islands)
            JuMP.@constraint(opfm.mod, pgu .<= bd.rampup * ramp_minutes)
            JuMP.@constraint(opfm.mod, pgd .<= bd.rampdown * ramp_minutes)
            JuMP.@constraint(opfm.mod, 0 .<= opfm.mod[:pg0] .+ pgu .- pgd)
            JuMP.@constraint(opfm.mod, opfm.mod[:pg0] .+ pgu .- pgd .<= bd.pg_lim_max)
            JuMP.@constraint(opfm.mod, bd.dc_lim_min .<= pfdccc)
            JuMP.@constraint(opfm.mod, pfdccc .<= bd.dc_lim_max)
            JuMP.@constraint(opfm.mod, prcc .<= bd.pr_lim .* max_shed)
            JuMP.@constraint(opfm.mod, lscc .<= bd.pd_lim .* max_shed)
        else
            for n in islands[island]
                JuMP.@constraint(opfm.mod,  [g = bd.list[n].ctrl_generation], 
                    pgu[g] <= bd.rampup[g] * ramp_minutes)
                JuMP.@constraint(opfm.mod,  [g = bd.list[n].ctrl_generation], 
                    pgd[g] <= bd.rampdown[g] * ramp_minutes)
                JuMP.@constraint(opfm.mod, [g = bd.list[n].ctrl_generation], 
                    0 <= opfm.mod[:pg0][g] + pgu[g] - pgd[g])
                JuMP.@constraint(opfm.mod, [g = bd.list[n].ctrl_generation], 
                    opfm.mod[:pg0][g] + pgu[g] - pgd[g] <= bd.pg_lim_max[g])
                JuMP.@constraint(opfm.mod, [g = bd.list[n].dc_branches], 
                    bd.dc_lim_min[g] <= pfdccc[g])
                JuMP.@constraint(opfm.mod, [g = bd.list[n].dc_branches], 
                    pfdccc[g] <= bd.dc_lim_max[g])
                JuMP.@constraint(opfm.mod, [g = bd.list[n].renewables], 
                    prcc[g] <= bd.pr_lim[g] * max_shed)
                JuMP.@constraint(opfm.mod, [d = bd.list[n].demands], 
                    lscc[d] <= bd.pd_lim[d] * max_shed)
            end
            for in_vec in islands[1:end .!= island]
                for n in in_vec
                    JuMP.@constraint(opfm.mod, [g = bd.list[n].ctrl_generation], 
                        opfm.mod[:pg0][g] == pgd[g] - pgu[g])
                    JuMP.@constraint(opfm.mod, [g = bd.list[n].dc_branches], 
                        pfdccc[g] == 0)
                    JuMP.@constraint(opfm.mod, [g = bd.list[n].renewables], 
                        prcc[g] == bd.pr_lim[g])
                    JuMP.@constraint(opfm.mod, [d = bd.list[n].demands], 
                        lscc[d] == bd.pd_lim[d])
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

function set_Pcc(opfm::OPFmodel, bd, ΔPcc, ptdf, overloads, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c, cut_added, lim)
    if !isempty(overloads)
        pgu, pgd, pfdccc, prcc, lscc = get_Pcc(opfm, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c)

        # sort!(overloads, rev = true, by = x -> abs(x[2]))
        # (i,ol) = first(overloads)
        for (i,ol) in overloads
            # Finding and adding the Benders cut
            expr = JuMP.@expression(opfm.mod, sum((ptdf[i, inode] * (
                    bd.Pg[inode] + ΔPcc[inode] -
                    sum((opfm.mod[:pg0][ctrl] + pgu[ctrl] - pgd[ctrl] 
                        for ctrl in sublist.ctrl_generation), init=0.0) -
                    sum((beta(sublist.node, opfm.dc_branches[d]) * pfdccc[d] for d in sublist.dc_branches), init=0.0) +
                    sum((prcc[r] for r in sublist.renewables), init=0.0) -
                    sum((lscc[d] for d in sublist.demands), init=0.0)
                    ) for (inode, sublist) in enumerate(bd.list)), init=0.0
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
function get_ΔPccx(opfm::OPFmodel, list, ΔP, Pccx, c)
    fill!(ΔP, zero(eltype(ΔP)))
    x = get(Pccx, c, 0)
    if x != 0
        pgd = get_value(opfm.mod, x[1])
        lscc = get_value(opfm.mod, x[2])
        for (i, n) in enumerate(list)
            for g in n.ctrl_generation
                ΔP[i] -= pgd[g]
            end
            for g in n.demands
                ΔP[i] += lscc[g]
            end
        end
    end
    return ΔP
end

function get_Pccx(opfm::OPFmodel, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c)
    x = get(bd.Pccx, c, 0)
    if x == 0

        # Add corrective variables
        pgd = JuMP.@variable(opfm.mod, [g in 1:length(opfm.ctrl_generation)], base_name = @sprintf( "pgdx%s", c), lower_bound = 0)
                # and ramp down
        lscc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.demands)], base_name = @sprintf( "lsccx%s", c), lower_bound = 0)
                # load curtailment variables in in contingencies
        bd.Pccx[c] = (pgd, lscc)

        # Extend the objective with the corrective variables
        set_objective_function(opfm.mod, objective_function(opfm.mod) + opfm.prob[c] * 
            (sum(opfm.voll' * lscc) + sum(60 * pgd)) * p_failure
            )

        # Add new constraints that limit the corrective variables within operating limits
        JuMP.@constraint(opfm.mod, sum(lscc) == sum(pgd))
        if isempty(islands)
            JuMP.@constraint(opfm.mod, pgd .<= bd.rampdown * ramp_minutes)
            JuMP.@constraint(opfm.mod, 0 .<= opfm.mod[:pg0] .- pgd)
            JuMP.@constraint(opfm.mod, lscc .<= bd.pd_lim .* max_shed)
        else
            for n in islands[island]
                JuMP.@constraint(opfm.mod,  [g = bd.list[n].ctrl_generation], 
                    pgd[g] <= bd.rampdown[g] * ramp_minutes)
                JuMP.@constraint(opfm.mod, [g = bd.list[n].ctrl_generation], 
                    0 <= opfm.mod[:pg0][g] - pgd[g])
                JuMP.@constraint(opfm.mod, [d = bd.list[n].demands], 
                    lscc[d] <= bd.pd_lim[d] * max_shed)
            end
            for in_vec in islands[1:end .!= island]
                for n in in_vec
                    JuMP.@constraint(opfm.mod, [g = bd.list[n].ctrl_generation], 
                        opfm.mod[:pg0][g] == pgd[g])
                    JuMP.@constraint(opfm.mod, [d = bd.list[n].demands], 
                        lscc[d] == bd.pd_lim[d])
                end
            end
        end            

    else # If the contingency is run before, a set of corrective variables already exist
        pgd = x[1]
        lscc = x[2]
    end
    return pgd, lscc
end

function set_Pccx(opfm::OPFmodel, bd, ΔPccx, ptdf, overloads, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c, cut_added, lim)
    if !isempty(overloads)
        pgd, lscc = get_Pccx(opfm, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c)
        # sort!(overloads, rev = true, by = x -> abs(x[2]))
        # (i,ol) = first(overloads)
        for (i,ol) in overloads
            # Finding and adding the Benders cut
            expr = JuMP.@expression(opfm.mod, sum((ptdf[i, inode] * (
                    bd.Pg[inode] + ΔPccx[inode] -
                    sum((opfm.mod[:pg0][ctrl] - pgd[ctrl] 
                        for ctrl in sublist.ctrl_generation), init=0.0)  -
                    sum((lscc[d] for d in sublist.demands), init=0.0)
                    ) for (inode, sublist) in enumerate(bd.list)), init=0.0
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
        for i in 1:length(opfm.contingencies)
            c = get(itr, i, 0)
            if c != 0 && JuMP.value(c[val][i_g]) > lim
                @printf("          c %12s: %s: %.3f\n", opfm.contingencies[i].name, symb, JuMP.value(c[val][i_g]))
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

# function add_contingency(opfm, c_name; ramp_minutes = 10, ramp_mult = 10, max_shed = 0.1, lim = 1e-6, branch_short_term_limit_multi::Real = 1.5)

#     pgu = JuMP.@variable(opfm.mod, [g in get_name.(get_ctrl_generation(opfm.sys))], base_name = "pgu", lower_bound = 0)
#     # active power variables for the generators in contingencies ramp up 
#     pgd = JuMP.@variable(opfm.mod, [g in get_name.(get_ctrl_generation(opfm.sys))], base_name = "pgd", lower_bound = 0)
#             # and ramp down
#     pfdccc = JuMP.@variable(opfm.mod, [d in get_name.(get_dc_branches(opfm.sys))], base_name = "pfdccc")
#     prcc = JuMP.@variable(opfm.mod, [r in get_name.(get_renewables(opfm.sys))], base_name = "prcc", lower_bound = 0)
#     lscc = JuMP.@variable(opfm.mod, [d in get_name.(get_demands(opfm.sys))], base_name = "lscc", lower_bound = 0)
#             # load curtailment variables in in contingencies
#     pfcc = JuMP.@variable(opfm.mod, [l in get_name.(get_branches(opfm.sys))])         # and after corrective actions
#     vacc = JuMP.@variable(opfm.mod, [b in get_name.(get_nodes(opfm.sys))])        # and after corrective actions
#     for g in get_ctrl_generation(opfm.sys)
#         set_upper_bound(pgu[get_name(g)], get_ramp_limits(g).up * ramp_minutes)
#         set_upper_bound(pgd[get_name(g)], get_ramp_limits(g).down * ramp_minutes)
#     end

#     set_objective_function(opfm.mod, objective_function(opfm.mod) + opfm.prob[c_name] * 
#         (sum(opfm.voll[d] * (# ramp_minutes / 60 * opfm.mod[:lsc][d,c] +
#             lscc[d]) for d in get_name.(get_demands(opfm.sys))
#             ) +
#         sum(opfm.voll[r] * prcc[r] for r in get_name.(get_renewables(opfm.sys))) +
#         sum(opfm.cost[g] * ramp_mult * pgu[g] for g in get_name.(get_ctrl_generation(opfm.sys)))
#         ) )

#     @constraint(opfm.mod, [n = get_name.(get_nodes(opfm.sys))],
#         sum(g.bus.name == n ? opfm.mod[:pg0][get_name(g)] + pgu[get_name(g)] - pgd[get_name(g)] : 0 for g in get_ctrl_generation(opfm.sys)) -
#         sum(beta(n,l) * pfcc[get_name(l)] for l in get_branches(opfm.sys)) -
#         sum(beta(n,l) * pfdccc[get_name(l)] for l in get_dc_branches(opfm.sys)) == 
#         sum(d.bus.name == n ? get_active_power(d) - lscc[get_name(d)] : 0 for d in get_demands(opfm.sys)) +
#         sum(d.bus.name == n ? -get_active_power(d) + prcc[get_name(d)] : 0 for d in get_renewables(opfm.sys))
#     )
#     i, slack = find_slack(opfm.sys)
#     @constraint(opfm.mod, vacc[get_name(slack)] == 0)

#     branch_rating = make_named_array(get_rate, get_branches(opfm.sys))
#     @constraint(opfm.mod, [l = get_name.(get_branches(opfm.sys))], 
#         -branch_rating[l] .<= pfcc[l] .<= branch_rating[l]
#     )
#     x = make_named_array(get_x, get_branches(opfm.sys))
#     @constraint(opfm.mod, [l = get_name.(get_branches(opfm.sys))],
#         pfcc[l] .- sum(beta(opfm.sys,l) .* vacc[:]) ./ x[l] .== 0
#     )
#     branch_rating = make_named_array(get_active_power_limits_from, get_dc_branches(opfm.sys))
#     @constraint(opfm.mod, [l = get_name.(get_dc_branches(opfm.sys)), c = opfm.contingencies], 
#         branch_rating[l].min <= opfm.mod[:pfdccc][l,c] <= branch_rating[l].max
#     )

#     for g in get_ctrl_generation(opfm.sys)
#         g_name = get_name(g)
#         @constraint(opfm.mod, 0 <= opfm.mod[:pg0][g_name] + (pgu[g_name] - pgd[g_name]) <= 
#             get_active_power_limits(g).max)
#     end

#     for l in get_demands(opfm.sys)
#         @constraint(opfm.mod, lscc[get_name(l)] <= get_active_power(l) * max_shed)
#     end
#     for l in get_renewables(opfm.sys)
#         @constraint(opfm.mod, prcc[get_name(l)] <= get_active_power(l))
#     end
#     @info "Added constraints for contingency on line $c_name"
# end
