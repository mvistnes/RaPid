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
function run_benders(type::OPF, system::System, voll, prob; ramp_minutes = 10, ramp_mult = 10, max_shed = 0.1, lim = 1e-6, short_term_limit_multi::Float64 = 1.5)
    # set_rate!.(get_branches(system), get_rate.(get_branches(system))*0.9)
    opfm = scopf(SC, system, Gurobi.Optimizer, voll=voll, ramp_minutes=ramp_minutes, max_shed=max_shed, short_term_limit_multi=short_term_limit_multi,renewable_prod = 1.0)
    MOI.set(opfm.mod, MOI.Silent(), true) # supress output from the solver
    solve_model!(opfm.mod)
    total_solve_time = solve_time(opfm.mod)
    @debug "lower_bound = $(objective_value(opfm.mod))"
    if type == SC::OPF
        short_term_limit_multi = 1.0
    end

    # Set variables
    nodes = get_sorted_nodes(opfm.sys)
    branches = get_sorted_branches(opfm.sys)
    idx = get_nodes_idx(nodes)
    list_ctrl = make_list(opfm, get_ctrl_generation)
    list_renewables = make_list(opfm, get_renewables)
    list_demands = make_list(opfm, get_demands)
    # cgen = connectivitymatrix(opfm.sys, length(nodes), idx)
    ctrl_generation = get_name.(get_ctrl_generation(opfm.sys))
    renewables = get_name.(get_renewables(opfm.sys))
    demands = get_name.(get_demands(opfm.sys))
    linerating = get_rate.(branches)
    pg_lim = make_named_array(get_active_power_limits, get_ctrl_generation(opfm.sys))
    pd_lim = make_named_array(get_active_power, get_nonctrl_generation(opfm.sys))
    ΔP = Dict{Int, NTuple{4, Any}}() # Holds the ramp up and down for contingencies that need
    overloads = []
    island = Vector{Vector{Int64}}() 
    island_b = Vector{Vector{Int64}}()

    # Get initial state
    pf = SCOPF.DCPowerFlow(nodes, branches, idx, get_net_Pᵢ(opfm, nodes, idx))
    ptdf = pf.ϕ
    P = get_Pᵢ(opfm, nodes)
    # print_contingency_overflow(opfm, pf, ΔP, 1.0)

    it = enumerate(get_bus_idx.(branches, [idx]))
    next = iterate(it)
    cut_added = 0
    iterations = 0
    while next !== nothing || cut_added > 0 # loops until no new cuts are added for the contingencies
        if next === nothing
            next = iterate(it)
            cut_added = 0
            iterations += 1
            @info "Iteration $iterations"
            # print_contingency_overflow(opfm, pf, ΔP, 1.0)
            if iterations >= length(branches)^2
                @error "Reached $(iterations) iterations without a stable solution."
                return opfm, pf, ΔP
            end
        end

        # Calculate the power flow with the new outage and find if there are any overloads
        (c, cont), state = next
        island = Vector{Vector{Int64}}() 
        island_b = Vector{Vector{Int64}}()
        # try
            try
                flow = calculate_line_flows(pf, cont, c)
                overloads = filter_overload(flow, linerating)
            catch DivideError
                island, island_b = handle_islands(pf, cont, c)
                try
                    flow = calculate_line_flows(pf, cont, c, island, island_b)
                    if isempty(flow)
                        @debug "Reference node isolated due to contingency on line $(cont[1])-$(cont[2])-i_$c."
                        overloads = [] # the reference node is not connected to any other nodes
                    else
                        @debug "Islands forming due to contingency on line $(cont[1])-$(cont[2])-i_$c."
                        overloads = filter_overload(flow, view(linerating, island_b))
                    end
                catch DivideError
                    overloads = [1]
                end
            end
            if !isempty(overloads) # ptdf calculation is more expensive than line flow
                δP = get_ΔP(opfm, length(nodes), idx, ΔP, c)
                # if overloads == [1] 
                if isempty(island)
                    ptdf = get_isf(pf.DA, pf.B, cont, c, pf.slack)
                else
                    fill!(ptdf, zero(eltype(ptdf)))
                    DA = pf.DA[island_b, island]
                    cn = searchsortedfirst(island_b, c)
                    ptdf[island_b, island] = get_isf(DA, pf.B[island, island], Tuple(DA[cn,:].nzind), cn,
                        searchsortedfirst(island, pf.slack) # searchsorted(island, pf.slack)
                    )
                end
                overloads = filter_overload(ptdf * (pf.Pᵢ .+ δP), linerating)
                # elseif isempty(island)
                #     ptdf = get_isf(pf, cont, c)
                #     overloads = filter_overload(ptdf * (pf.Pᵢ .+ δP), linerating)
                #     # print(c)
                #     # println(pf.Pᵢ)
                # else
                #     fill!(ptdf, zero(eltype(ptdf)))
                #     ptdf[island_b, island] = get_isf(pf.X, pf.B, pf.DA, cont, c, island)
                #     flow = ptdf * (pf.Pᵢ .+ δP)
                #     overloads = filter_overload(flow, linerating)
                # end
            end
        # catch e
        #     overloads = []
        #     @warn "Contingency on line $(cont[1])-$(cont[2])-i_$c resulted in $e."
        # end
        if !isempty(overloads)
            # add_contingency(opfm, prob, branches[c].name, ramp_minutes=ramp_minutes, ramp_mult=ramp_mult, 
            #     max_shed=max_shed, lim=lim, short_term_limit_multi=short_term_limit_multi)
            overloads_st = filter_overload(overloads, linerating * (short_term_limit_multi - 1.0))
            
            # Add preventive actions
            for (i,ol) in overloads_st
                expr = JuMP.@expression(opfm.mod, sum(
                        (ptdf[i, in] * (P[in] - 
                        sum((opfm.mod[:pg0][ctrl.name] for ctrl in list_ctrl[n.name]), init=0.0) +
                        sum((opfm.mod[:ls0][r.name] for r in list_renewables[n.name]), init=0.0) +
                        sum((opfm.mod[:ls0][d.name] for d in list_demands[n.name]), init=0.0))
                        for (in, n) in enumerate(nodes) if abs(ptdf[i, in]) > lim), init=0.0
                    ) )
                
                @info @sprintf "Pre: Contingency line %d-%d-i_%d; overload on %s of %.2f" cont[1] cont[2] c branches[i].name ol
                @debug "Cut added: $(sprint_expr(expr,lim))\n"
                # set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
                if ol < 0
                    pre_cut = JuMP.@constraint(opfm.mod, expr <= ol)
                else
                    pre_cut = JuMP.@constraint(opfm.mod, expr >= ol)
                end
                cut_added = 2
            end
            if type == PCSC::OPF

                # Add corrective variables
                x = get(ΔP, c, 0)
                if x == 0
                    pgu = JuMP.@variable(opfm.mod, [g in ctrl_generation], base_name = "pgu", lower_bound = 0)
                        # active power variables for the generators in contingencies ramp up 
                    pgd = JuMP.@variable(opfm.mod, [g in ctrl_generation], base_name = "pgd", lower_bound = 0)
                            # and ramp down
                    prcc = JuMP.@variable(opfm.mod, [r in renewables], base_name = "prcc", lower_bound = 0)
                    lscc = JuMP.@variable(opfm.mod, [d in demands], base_name = "lscc", lower_bound = 0)
                            # load curtailment variables in in contingencies
                    for g in get_ctrl_generation(system)
                        set_upper_bound(pgu[get_name(g)], get_ramp_limits(g).up * ramp_minutes)
                        set_upper_bound(pgd[get_name(g)], get_ramp_limits(g).down * ramp_minutes)
                    end
                    ΔP[c] = (pgu, pgd, prcc, lscc)
                else # If the contingency is run before, a set of corrective variables already exist
                    pgu = x[1]
                    pgd = x[2]
                    prcc = x[3]
                    lscc = x[4]
                end

                # sort!(overloads, rev = true, by = x -> abs(x[2]))
                # (i,ol) = first(overloads)
                for (i,ol) in overloads
                    # Finding and adding the Benders cut
                    expr = JuMP.@expression(opfm.mod, sum((ptdf[i, in] * (P[in] - 
                            sum((opfm.mod[:pg0][ctrl.name] + pgu[ctrl.name] - pgd[ctrl.name] 
                                for ctrl in list_ctrl[n.name]), init=0.0) -
                            sum((opfm.mod[:ls0][r.name] + prcc[r.name] 
                                for r in list_renewables[n.name]), init=0.0) +
                            sum((opfm.mod[:ls0][d.name] + lscc[d.name] 
                                for d in list_demands[n.name]), init=0.0))
                            for (in, n) in enumerate(nodes) if abs(ptdf[i, in]) > lim), init=0.0
                        ) )
                    
                    @info @sprintf "Corr: Contingency line %d-%d-i_%d; overload on %s of %.4f" cont[1] cont[2] c branches[i].name ol
                    @debug "Cut added: $(sprint_expr(expr,lim))\n"
                    # set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
                    if ol < 0
                        corr_cut = JuMP.@constraint(opfm.mod, expr <= ol)
                    else
                        corr_cut = JuMP.@constraint(opfm.mod, expr >= ol)
                    end
                    cut_added = 3
                end

                # Extend the objective with the corrective variables
                c_name = branches[c].name
                @objective(opfm.mod, Min, objective_function(opfm.mod) + prob[c_name] * 
                    (sum(opfm.voll[d] * (# ramp_minutes / 60 * opfm.mod[:lsc][d,c] +
                        lscc[d]) for d in demands
                        ) +
                    sum(opfm.voll[r] * prcc[r] for r in renewables) +
                    sum(opfm.cost[g] * ramp_mult * pgu[g] for g in ctrl_generation)
                    ) )

                # Add new constraints that limit the corrective variables within operating limits
                JuMP.@constraint(opfm.mod, sum(pgu) + sum(lscc) == sum(pgd) + sum(prcc))
                JuMP.@constraint(opfm.mod, [g = ctrl_generation], 
                    0 <= opfm.mod[:pg0][g] + pgu[g] - pgd[g] <= pg_lim[g].max)
                JuMP.@constraint(opfm.mod, [r = renewables], 
                    0 <= opfm.mod[:ls0][r] + prcc[r] <= pd_lim[r] * max_shed)
                JuMP.@constraint(opfm.mod, [d = demands], 
                    0 <= opfm.mod[:ls0][d] + lscc[d] <= pd_lim[d] * max_shed)
            end
            if cut_added > 1
                solve_model!(opfm.mod)
                total_solve_time += solve_time(opfm.mod)
                if termination_status(opfm.mod) != MOI.OPTIMAL 
                    @info "Removing latest cut."
                    (@isdefined pre_cut) && delete(opfm.mod, pre_cut)
                    (@isdefined corr_cut) && delete(opfm.mod, corr_cut)
                    solve_model!(opfm.mod)
                    if termination_status(opfm.mod) != MOI.OPTIMAL
                        throw(ErrorException("Infeasible model!"))
                    end
                end
                cut_added = 1
                P = update_model!(opfm, pf, nodes, idx)
            end
        end
        next = iterate(it, state)
    end
    @info "Total solve time $total_solve_time"
    return opfm, pf, ΔP
end         

""" Solve model and update the power flow object """
function update_model!(opfm, pf, nodes, idx)
    P = get_Pᵢ(opfm, nodes, idx)
    pf.Pᵢ = get_net_Pᵢ(opfm, nodes, idx, P)
    calculate_line_flows!(pf)
    run_pf!(pf)
    return P
end

""" Return the power injection change at each node. """
function get_ΔP(opfm::OPFmodel, numnodes::Integer, idx::Dict{<:Any, <:Int}, ΔP, c)
    δP = zeros(numnodes)
    x = get(ΔP, c, 0)
    if x != 0
        pgu = JuMP.value.(x[1])
        pgd = JuMP.value.(x[2])
        prcc = JuMP.value.(x[3])
        lscc = JuMP.value.(x[4])
        
        for g in get_ctrl_generation(opfm.sys)
            δP[idx[g.bus.number]] += pgu[get_name(g)] - pgd[get_name(g)]
        end
        for g in get_renewables(opfm.sys)
            δP[idx[g.bus.number]] += prcc[get_name(g)]
        end
        for g in get_demands(opfm.sys)
            δP[idx[g.bus.number]] += lscc[get_name(g)]
        end
    end
    return δP
end

" An AffExpr nicely formatted to a string "
sprint_expr(expr::AffExpr, lim = 1e-6) = 
    join(Printf.@sprintf("%s%5.2f %s ", (x[2] > 0 ? "+" : "-"), abs(x[2]), x[1]) 
            for x in expr.terms if abs(x[2]) > lim) * 
        Printf.@sprintf("<= %s%.2f", (expr.constant > 0 ? "-" : "+"), abs(expr.constant)
    )

function print_benders_results(opfm::OPFmodel, ΔP)
    for g in get_name.(get_ctrl_generation(opfm.sys))
        @printf("%s: %.3f\n", g, JuMP.value(opfm.mod[:pg0][g]))
        for (i,c) in ΔP
            if JuMP.value(c[1][g]) > 0.0001 || JuMP.value(c[2][g]) > 0.0001
                @printf("   c%s: pgu: %.3f  pgd: %.3f\n", i, JuMP.value(c[1][g]), JuMP.value(c[2][g]))
            end
        end
    end
    for g in get_name.(get_renewables(opfm.sys))
        @printf("%s: %.3f\n", g, JuMP.value(opfm.mod[:ls0][g]))
        for (i,c) in ΔP
            if JuMP.value(c[3][g]) > 0.0001
                @printf("   c%s: prcc: %.3f\n", i, JuMP.value(c[3][g]))
            end
        end
    end
    for g in get_name.(get_demands(opfm.sys))
        @printf("%s: %.3f\n", g, JuMP.value(opfm.mod[:ls0][g]))
        for (i,c) in ΔP
            if JuMP.value(c[4][g]) > 0.0001
                @printf("   c%s: lscc: %.3f\n", i, JuMP.value(c[4][g]))
            end
        end
    end
end

function findall_overloads(pf, idx, branches, linerating)
    overloads = []
    ol = []
    for (c, cont) in enumerate(get_bus_idx.(branches, [idx]))
        try
            flow = calculate_line_flows(pf, cont, c)
            ol = filter_overload(flow, linerating)
        catch DivideError
            @info "Islands forming due to contingency on line $c $(cont[1])-$(cont[2])."
            island, island_b = handle_islands(pf, cont, c)
            flow = calculate_line_flows(pf, c, cont, island, island_b)
            ol = filter_overload(flow, view(linerating, island_b))
        end
        if !isempty(ol)
            push!(overloads, sort!(ol, rev = true, by = x -> abs(x[2])))
        end
    end
    return sort!(overloads, rev = true, by = x -> abs(x[1][2]))
end


function add_contingency(opfm, prob, c_name; ramp_minutes = 10, ramp_mult = 10, max_shed = 0.1, lim = 1e-6, short_term_limit_multi::Float64 = 1.5)

    pgu = JuMP.@variable(opfm.mod, [g in get_name.(get_ctrl_generation(opfm.sys))], base_name = "pgu", lower_bound = 0)
    # active power variables for the generators in contingencies ramp up 
    pgd = JuMP.@variable(opfm.mod, [g in get_name.(get_ctrl_generation(opfm.sys))], base_name = "pgd", lower_bound = 0)
            # and ramp down
    prcc = JuMP.@variable(opfm.mod, [r in get_name.(get_renewables(opfm.sys))], base_name = "prcc", lower_bound = 0)
    lscc = JuMP.@variable(opfm.mod, [d in get_name.(get_demands(opfm.sys))], base_name = "lscc", lower_bound = 0)
            # load curtailment variables in in contingencies
    pfcc = JuMP.@variable(opfm.mod, [l in get_name.(get_branches(opfm.sys))])         # and after corrective actions
    vacc = JuMP.@variable(opfm.mod, [b in get_name.(get_nodes(opfm.sys))])        # and after corrective actions
    for g in get_ctrl_generation(opfm.sys)
        set_upper_bound(pgu[get_name(g)], get_ramp_limits(g).up * ramp_minutes)
        set_upper_bound(pgd[get_name(g)], get_ramp_limits(g).down * ramp_minutes)
    end

    @objective(opfm.mod, Min, objective_function(opfm.mod) + prob[c_name] * 
        (sum(opfm.voll[d] * (# ramp_minutes / 60 * opfm.mod[:lsc][d,c] +
            lscc[d]) for d in get_name.(get_demands(opfm.sys))
            ) +
        sum(opfm.voll[r] * prcc[r] for r in get_name.(get_renewables(opfm.sys))) +
        sum(opfm.cost[g] * ramp_mult * pgu[g] for g in get_name.(get_ctrl_generation(opfm.sys)))
        ) )

    @constraint(opfm.mod, [n = get_name.(get_nodes(opfm.sys))],
        sum(g.bus.name == n ? opfm.mod[:pg0][get_name(g)] + pgu[get_name(g)] - pgd[get_name(g)] : 0 for g in get_ctrl_generation(opfm.sys)) -
        sum(beta(n,l) * pfcc[get_name(l)] for l in get_branches(opfm.sys)) == 
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
