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
    opfm = scopf(SC, system, Gurobi.Optimizer, voll=voll, ramp_minutes=ramp_minutes, max_shed=max_shed, short_term_limit_multi=short_term_limit_multi)
    solve_model!(opfm.mod)
    @debug "lower_bound = $(objective_value(opfm.mod))"

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
    linerating = get_rate.(branches) * short_term_limit_multi
    pg_lim = make_named_array(get_active_power_limits, get_ctrl_generation(opfm.sys))
    pd_lim = make_named_array(get_active_power, get_nonctrl_generation(opfm.sys))
    ΔP = Dict{Int, NTuple{4, Any}}() # Holds the ramp up and down for contingencies that need
    overloads = []

    # Get initial state
    pf = SCOPF.DCPowerFlow(nodes, branches, idx, get_net_Pᵢ(opfm, nodes, idx))
    ptdf = pf.ϕ

    it = enumerate(get_bus_idx.(branches, [idx]))
    next = iterate(it)
    cut_added = false
    iterations = 0
    while next !== nothing || cut_added # loops until no new cuts are added for the contingencies
        if next === nothing
            next = iterate(it)
            cut_added = false
            iterations += 1
            if iterations >= 10 
                @error "Reached $(iterations) iterations without a stable solution."
                return opfm, pf, ΔP
            end
        end
        (c, cont), state = next
        δP = get_ΔP(opfm, length(nodes), idx, ΔP, c)
        try
            ptdf = get_isf(pf, cont[1], cont[2], c)
            overloads = filter_overload(ptdf * (pf.Pᵢ .+ δP), linerating)
        catch DivideError
            @info "Islands forming due to contingency on line $(cont[1])-$(cont[2])."
            flow = handle_islands(system, pf, branches, δP, (c, cont))
            if isempty(flow)
                overloads = [] # the refernece node is not connected to any other nodes
            else
                overloads = filter_overload(flow, linerating[1:end .!= c])
            end
        end
        if length(overloads) > 0
            P = get_Pᵢ(opfm, nodes)
            if type == PCSC::OPF
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
                else
                    pgu = x[1]
                    pgd = x[2]
                    prcc = x[3]
                    lscc = x[4]
                end

                sort!(overloads, rev = true, by = x -> abs(x[2]))
                (i,ol) = first(overloads)
                # for (i,ol) in overloads
                expr = JuMP.@expression(opfm.mod, sum((ptdf[i, in] * (P[in] - 
                        sum((opfm.mod[:pg0][ctrl.name] + pgu[ctrl.name] - pgd[ctrl.name] 
                            for ctrl in list_ctrl[n.name]), init=0.0) -
                        sum((opfm.mod[:ls0][r.name] + prcc[r.name] 
                            for r in list_renewables[n.name]), init=0.0) +
                        sum((opfm.mod[:ls0][d.name] + lscc[d.name] 
                            for d in list_demands[n.name]), init=0.0))
                        for (in, n) in enumerate(nodes) if abs(ptdf[i, in]) > lim), init=0.0
                    ) )
                
                if expr != 0.0
                    @info "Contingency on line $(cont[1])-$(cont[2]) resulted in overload on $(branches[i].name) of $(ol)"
                    @debug "Cut added: $(sprint_expr(expr,lim))\n"
                    # set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
                    if ol < 0
                        JuMP.@constraint(opfm.mod, expr <= ol)
                    else
                        JuMP.@constraint(opfm.mod, expr >= ol)
                    end
                end
                # end
                c_name = branches[c].name
                @objective(opfm.mod, Min, objective_function(opfm.mod) + prob[c_name] * 
                    (sum(opfm.voll[d] * (# ramp_minutes / 60 * opfm.mod[:lsc][d,c] + # extract these from the for loop in N-1 too 
                        lscc[d]) for d in demands
                        ) +
                    sum(opfm.voll[r] * prcc[r] for r in renewables) +
                    sum(opfm.cost[g] * ramp_mult * pgu[g] for g in ctrl_generation)
                    ) )
                JuMP.@constraint(opfm.mod, sum(pgu) + sum(lscc) == sum(pgd) + sum(prcc))
                JuMP.@constraint(opfm.mod, [g = ctrl_generation], 
                    0 <= opfm.mod[:pg0][g] + pgu[g] - pgd[g] <= pg_lim[g].max)
                JuMP.@constraint(opfm.mod, [r = renewables], 
                    0 <= opfm.mod[:ls0][r] + prcc[r] <= pd_lim[r] * max_shed)
                JuMP.@constraint(opfm.mod, [d = demands], 
                    0 <= opfm.mod[:ls0][d] + lscc[d] <= pd_lim[d] * max_shed)
            else
                for (i,ol) in overloads
                    expr = JuMP.@expression(opfm.mod, sum(
                            (ptdf[i, in] * (P[in] - 
                            sum((opfm.mod[:pg0][ctrl.name] for ctrl in list_ctrl[n.name]), init=0.0) +
                            sum((opfm.mod[:ls0][r.name] for r in list_renewables[n.name]), init=0.0) +
                            sum((opfm.mod[:ls0][d.name] for d in list_demands[n.name]), init=0.0))
                            for (in, n) in enumerate(nodes) if abs(ptdf[i, in]) > lim), init=0.0
                        ) )
                    
                    if expr != 0.0
                        @info "Contingency on line $(cont[1])-$(cont[2]) resulted in overload on $(branches[i].name) of $(ol)"
                        @debug "Cut added: $(sprint_expr(expr,lim))\n"
                        # set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
                        if ol < 0
                            JuMP.@constraint(opfm.mod, expr <= ol)
                        else
                            JuMP.@constraint(opfm.mod, expr >= ol)
                        end
                    end
                end
            end

            MOI.set(opfm.mod, MOI.Silent(), true) # supress output from the solver
            opfm.mod = solve_model!(opfm.mod)
            termination_status(opfm.mod) != MOI.OPTIMAL && return opfm, pf, ΔP
            cut_added = true
            pf.Pᵢ = get_net_Pᵢ(opfm, nodes, idx)
            calculate_line_flows!(pf)
        end
        next = iterate(it, state)
    end
    run_pf!(pf) # the only update of the voltage angles
    return opfm, pf, ΔP
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
        # δP[getindex.([idx], get_number.(get_bus.(get_ctrl_generation(opfm.sys))))] += pgu.data - pgd.data
        # δP[getindex.([idx], get_number.(get_bus.(get_renewables(opfm.sys))))] += prcc.data
        # δP[getindex.([idx], get_number.(get_bus.(get_demands(opfm.sys))))] += lscc.data
        
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

function handle_islands(system::System, pf::DCPowerFlow, branches::AbstractVector{<:ACBranch},
        δP::AbstractVector{<:Real}, contingency::Tuple{Integer, Tuple{Integer, Integer}})
    islands = get_islands(system, branches[1:end .!= contingency[1]])
    for island in islands
        if pf.slack ∈ get_number.(island[1]) # Nodes in island 1
            island_idx = get_nodes_idx(island[1])
            drop_id = sort!(collect(values(island_idx)))
            drop_cont = [j for (j,b) in enumerate(branches) if b.name != branches[contingency[1]].name && b.name ∈ get_name.(island[2])]
            if length(drop_id) > 1 && length(drop_cont) > 1
                ptdf = get_isf(view(pf.X, drop_id, drop_id), view(pf.B, drop_id, drop_id), view(pf.DA, drop_cont, drop_id), 
                    contingency[2][1], contingency[2][1], contingency[1])
                return ptdf * getindex((pf.Pᵢ .+ δP), drop_id)
            end
            return []
        end
    end
    return []
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