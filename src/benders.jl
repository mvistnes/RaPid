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
function run_benders(system::System, voll, contingencies, prob, lim = 1e-6)

    opfm = scopf(SC, system, Gurobi.Optimizer, voll=voll)
    solve_model!(opfm.mod)
    lower_bound = objective_value(opfm.mod)

    # Set global variables
    nodes = get_sorted_nodes(opfm.sys)
    branches = get_sorted_branches(opfm.sys)
    idx = get_nodes_idx(nodes)
    list_ctrl = make_list(opfm, get_ctrl_generation)
    contanal = iterative_cont_anal(opfm.sys, nodes, branches, contingencies, idx)

    it = enumerate(contingencies)
    next = iterate(it)
    cut_added = false
    while next !== nothing || cut_added # loops until no new cuts are added for the contingencies
        if next === nothing
            next = iterate(it)
            cut_added = false
        end
        (c, cont), state = next
        overloads = get_overload(contanal, c, get_net_Pᵢ(opfm, nodes))
        overloads = [(i,ol) for (i,ol) in enumerate(overloads) if abs(ol) > lim]
        if length(overloads) > 0
            # sort!(overloads, rev = true, by = x -> abs(x[2]))
            (i,ol) = first(overloads)
            Pᵢ = get_Pᵢ(opfm, nodes)
            ptdf = get_cont_ptdf(contanal, c)
            
            expr = sum(
                    ptdf[i, in] * 
                    (Pᵢ[in] - sum((opfm.mod[:pg0][g.name] for g in list_ctrl[n.name]), init=0.0)) 
                    for (in, n) in enumerate(nodes) if abs(ptdf[i, in]) > lim
                )
            @info "Contingency on $(cont) resulted in overload on $(branches[i].name) of $(ol) \nCut added: $(sprint_expr(expr,lim))\n"
            set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
            if ol < 0
                JuMP.@constraint(opfm.mod, expr <= ol)
            else
                JuMP.@constraint(opfm.mod, expr >= ol)
            end

            MOI.set(opfm.mod, MOI.Silent(), true) # supress output from the solver
            opfm.mod = solve_model!(opfm.mod)
            termination_status(opfm.mod) != MOI.OPTIMAL && return
            cut_added = true
        else
            next = iterate(it, state)
        end
    end
    return opfm, contanal
end

function run_benders2(system::System, voll, prob, ramp_minutes = 10, repair_time = 4, lim = 1e-6)

    opfm = scopf(SC, system, Gurobi.Optimizer, voll=voll)
    solve_model!(opfm.mod)
    # lower_bound = objective_value(opfm.mod)

    # Set variables
    nodes = get_sorted_nodes(opfm.sys)
    branches = get_sorted_branches(opfm.sys)
    idx = get_nodes_idx(nodes)
    list_ctrl = make_list(opfm, get_ctrl_generation)
    list_nonctrl = make_list(opfm, get_nonctrl_generation)
    # cgen = connectivitymatrix(opfm.sys, length(nodes), idx)
    ctrl_generation = get_name.(get_ctrl_generation(opfm.sys))
    non_ctrl_generation = get_name.(get_nonctrl_generation(opfm.sys))
    linerating = get_rate.(branches)
    ΔP = Dict() # Holds the ramp up and down for contingencies that need
    overloads = []

    # Get initial state
    pf = SCOPF.DCPowerFlow(nodes, branches, idx, get_net_Pᵢ(opfm, nodes, idx))

    it = enumerate(get_bus_idx.(branches, [idx]))
    next = iterate(it)
    cut_added = false
    while next !== nothing || cut_added # loops until no new cuts are added for the contingencies
        if next === nothing
            next = iterate(it)
            cut_added = false
        end
        (c, cont), state = next
        try
            overloads = get_overload(pf, c, cont, (pf.Pᵢ .+ get_ΔP(opfm, length(nodes), idx, ΔP, c)), linerating)
            overloads = [(i,ol) for (i,ol) in enumerate(overloads) if abs(ol) > lim]
        catch DivideError
            @warn "Contingency on line $(cont[1])-$(cont[2]) resulted in islands forming"
            overloads = []
        end
        if length(overloads) > 0
            # sort!(overloads, rev = true, by = x -> abs(x[2]))
            (i,ol) = first(overloads)
            P = get_Pᵢ(opfm, nodes)
            ptdf = get_isf(pf, cont[1], cont[2], c)
            
            pgu = JuMP.@variable(opfm.mod, [g in ctrl_generation], lower_bound = 0)
                # active power variables for the generators in contingencies ramp up 
            pgd = JuMP.@variable(opfm.mod, [g in ctrl_generation], lower_bound = 0)
                    # and ramp down
            ΔP[c] = (pgu, pgd)
            for g in get_ctrl_generation(system)
                set_upper_bound(pgu[get_name(g)], get_ramp_limits(g).up * ramp_minutes)
                set_upper_bound(pgd[get_name(g)], get_ramp_limits(g).down * ramp_minutes)
            end
                # pfcc[l in get_name.(get_branches(opfm.sys))]
                #     # power flow on get_branches in in contingencies after corrective actions
                # vacc[b in get_name.(get_nodes(opfm.sys))]
                #     # voltage angle at a node in in contingencies after corrective actions
            #lscc = JuMP.@variable(opfm.mod, [d in non_ctrl_generation], lower_bound = 0)
                    # load curtailment variables in in contingencies
            expr = JuMP.@expression(opfm.mod, sum(
                    (ptdf[i, in] * (P[in] - 
                    sum((opfm.mod[:pg0][ctrl.name] + pgu[ctrl.name] - pgd[ctrl.name] 
                        for ctrl in list_ctrl[n.name]), init=0.0) #+
                    #sum((opfm.mod[:ls0][nonctrl.name] + lscc[nonctrl.name] 
                    #=    for nonctrl in list_nonctrl[n.name]), init=0.0)=#)
                    for (in, n) in enumerate(nodes) if abs(ptdf[i, in]) > lim), init=0.0
                ) )
            if expr != 0.0
                @info "Contingency on line $(cont[1])-$(cont[2]) resulted in overload on $(branches[i].name) of $(ol) \nCut added: $(sprint_expr(expr,lim))\n"
                # set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
                if ol < 0
                    JuMP.@constraint(opfm.mod, expr <= ol)
                else
                    JuMP.@constraint(opfm.mod, expr >= ol)
                end
                c_name = branches[c].name
                @objective(opfm.mod, Min, objective_function(opfm.mod) + prob[c_name] * 
                    #(sum(opfm.voll[d] * (# ramp_minutes / 60 * opfm.mod[:lsc][d,c] + # extract these from the for loop in N-1 too 
                    #    repair_time * lscc[d]) for d in non_ctrl_generation
                    #    ) +
                    sum(opfm.cost[g] * repair_time * pgu[g]
                            for g in ctrl_generation
                        )#)
                    )
                JuMP.@constraint(opfm.mod, sum(pgu) == sum(pgd))
            end

            MOI.set(opfm.mod, MOI.Silent(), true) # supress output from the solver
            opfm.mod = solve_model!(opfm.mod)
            termination_status(opfm.mod) != MOI.OPTIMAL && return opfm, pf
            cut_added = true
            pf.Pᵢ = get_net_Pᵢ(opfm, nodes, idx)
            run_pf!(pf)
            calc_Pline!(pf)
        else
            next = iterate(it, state)
        end
    end
    return opfm, pf, ΔP
end         

""" Return the power injected at each node. """
function get_ΔP(opfm::OPFmodel, numnodes::Integer, idx::Dict{<:Any, <:Int}, ΔP, c)
    Pᵢ = zeros(numnodes)
    x = get(ΔP, c, 0)
    if x != 0
        pgu = JuMP.value.(x[1])
        pgd = JuMP.value.(x[2])
        for g in get_ctrl_generation(opfm.sys)
            Pᵢ[idx[g.bus.number]] += pgu[get_name(g)] - pgd[get_name(g)]
        end
    end
    return Pᵢ
end


" An AffExpr nicely formatted to a string "
sprint_expr(expr::AffExpr, lim = 1e-6) = 
    join(Printf.@sprintf("%s%5.2f %s ", (x[2] > 0 ? "+" : "-"), abs(x[2]), x[1]) 
            for x in expr.terms if abs(x[2]) > lim) * 
        Printf.@sprintf("<= %s%.2f", (expr.constant > 0 ? "-" : "+"), abs(expr.constant)
    )

function calc_all_line_flows(
            DA::AbstractMatrix{<:Real}, 
            X::AbstractMatrix{<:Real}, 
            B::AbstractMatrix{<:Real}, 
            δ::AbstractVector{<:Real}, 
            branches::AbstractVector{<:Branch}, 
            idx::Dict{<:Any, <:Integer},
            slack::Integer
        )
    branches_idx = get_bus_idx.(branches, [idx])
    duplex = ones(length(branches))
    b1 = 0
    for (i,b2) in enumerate(branches_idx)
        if b1 == b2
            if duplex[i] < 1.0
                x = _distribute!(branches_idx, i)
                duplex[i-x:i] .= 1.0 / x
            else
                duplex[i-1:i] .= 0.5
            end
        end
        b1 = b2
    end
    return [
            calc_Pline(DA, get_changed_angles(X, B, δ, x[1], x[2], slack, duplex[i]), i) 
            for (i,x) in enumerate(branches_idx)
        ]
end

function _distribute!(branches_idx::Tuple{<:Integer, <:Integer}, i::Integer)
    x = 0
    branch = branches_idx[i]
    while i > 1
        i -= 1
        if branches_idx[i] == branch
            x += 1
        else
            break
        end
    end
    return x
end