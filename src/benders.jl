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
    list_gen = make_list(opfm, get_ctrl_generation)
    contanal = iterative_cont_anal(opfm.sys, nodes, branches, contingencies)

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
                    (Pᵢ[in] - sum((opfm.mod[:pg0][g.name] for g in list_gen[n.name]), init=0.0)) 
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

function run_benders2(system::System, voll, prob, lim = 1e-6)

    opfm = scopf(SC, system, Gurobi.Optimizer, voll=voll)
    solve_model!(opfm.mod)
    # lower_bound = objective_value(opfm.mod)

    # Set global variables
    nodes = get_sorted_nodes(opfm.sys)
    branches = get_sorted_branches(opfm.sys)
    # cgen = connectivitymatrix(opfm.sys, length(nodes), idx)
    list_gen = make_list(opfm, get_ctrl_generation)
    idx = get_nodes_idx(nodes)
    imml = IMML(nodes, branches, idx, branches)
    lu = LinearAlgebra.factorize(imml.B)
    # bbus = get_bus_idx(branches, idx)
    # lodf = get_lodf(bbus[1], bbus[2], get_x.(branches), imml.A, imml.X)
    slack = find_slack(nodes)

    # Get initial state
    Pᵢ = get_net_Pᵢ(opfm, nodes, idx)
    angles = run_pf(lu, Pᵢ)
    Pl0 = calc_Pline(imml.DA, angles)

    it = enumerate(get_bus_idx.(branches, [idx]))
    next = iterate(it)
    cut_added = false
    while next !== nothing || cut_added # loops until no new cuts are added for the contingencies
        if next === nothing
            next = iterate(it)
            cut_added = false
        end
        (c, cont), state = next
        overloads = get_overload(Pl0, angles, imml, c, cont)
        overloads = [(i,ol) for (i,ol) in enumerate(overloads) if abs(ol) > lim]
        # Skjekk om det er øyer i nettet
        if length(overloads) > 0
            # sort!(overloads, rev = true, by = x -> abs(x[2]))
            (i,ol) = first(overloads)
            P = get_Pᵢ(opfm, nodes)
            ptdf = get_isf(imml.DA, imml.X, imml.B, cont[1], cont[2], i, slack[1], change)
            
            pgu = @variable(opfm.mod, [g in get_name.(get_ctrl_generation(opfm.sys))], lower_bound = 0)
                # active power variables for the generators in contingencies ramp up 
            pgd = @variable(opfm.mod, [g in get_name.(get_ctrl_generation(opfm.sys))], lower_bound = 0)
                    # and ramp down
                # pfcc[l in get_name.(get_branches(opfm.sys))]
                #     # power flow on get_branches in in contingencies after corrective actions
                # vacc[b in get_name.(get_nodes(opfm.sys))]
                #     # voltage angle at a node in in contingencies after corrective actions
            lscc = @variable(opfm.mod, [d in get_name.(get_nonctrl_generation(opfm.sys))], lower_bound = 0)
                    # load curtailment variables in in contingencies
            expr = sum(
                    (ptdf[i, in] * (P[in] - 
                    sum((opfm.mod[:pg0][g.name] - pgu[g.name] + pgd[g.name] for g in list_gen[n.name]), init=0.0) +
                    opfm.mod[:ls0][in] + lscc[in]) for (in, n) in enumerate(nodes) if abs(ptdf[i, in]) > lim), init=0.0
                )
            if expr != 0.0
                # @info "Contingency on $(cont) resulted in overload on $(branches[i].name) of $(ol) \nCut added: $(sprint_expr(expr,lim))\n"
                set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
                if ol < 0
                    JuMP.@constraint(opfm.mod, expr <= ol)
                else
                    JuMP.@constraint(opfm.mod, expr >= ol)
                end
                add_to_expression!(objective_function(opfm.mod), prob[c] * 
                    sum(opfm.voll[d] * (ramp_minutes / 60 * opfm.mod[:lsc][d,c] + # extract these from the for loop in N-1 too 
                        repair_time * lscc[d])
                            for d in get_name.(get_nonctrl_generation(opfm.sys))
                        ) +
                    sum(opfm.cost[g] * repair_time * pgu[g]
                            for g in get_name.(get_ctrl_generation(opfm.sys))
                        )
                    )
            end

            MOI.set(opfm.mod, MOI.Silent(), true) # supress output from the solver
            opfm.mod = solve_model!(opfm.mod)
            termination_status(opfm.mod) != MOI.OPTIMAL && return
            cut_added = true
            Pᵢ = get_net_Pᵢ(opfm, nodes, idx)
            angles = run_pf(lu, Pᵢ)
            Pl0 = calc_Pline(imml.DA, angles)
        else
            next = iterate(it, state)
        end
    end
    return opfm, imml
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