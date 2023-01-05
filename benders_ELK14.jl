# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022 
# Based on code by Sigurd Hofsmo Jakobsen, SINTEF Energy Research, 2022

using PowerSystems
import JuMP
import Printf
import Gurobi
import MathOptInterface
const MOI = MathOptInterface
import LinearAlgebra
import SparseArrays
# include("N-1_SCOPF.jl")
# include("short_long_SCOPF.jl")
# include("utils.jl")

""" 
Solve the optimization model using Benders decomposition.

Function creates extra production constraints on generators in the system
based on the power transfer distribution factors of the generators in the
system. 
"""
function run_benders(system::System, voll, contingencies, prob;
            time_limit_sec::Integer = 600,
            unit_commit::Bool = false,
            max_shed::Float64 = 0.1,
            max_curtail::Float64 = 1.0,
            ratio::Float64= 0.5, 
            circuit_breakers::Bool=false,
            short_term_limit_multi::Float64 = 1.5,
            ramp_minutes::Integer = 10,
            repair_time::Float64 = 1.0
        )

    opfm = scopf(SC, system, Gurobi.Optimizer, voll=voll)
    solve_model!(opfm.mod)
    lower_bound = objective_value(opfm.mod)

    submod = d_scopf(system, Gurobi.Optimizer, value.(opfm.mod[:pg0]),
        time_limit_sec=time_limit_sec, max_shed=max_shed, max_curtail=max_curtail,
        ramp_minutes=ramp_minutes)

    # Set global variables
    nodes = get_sorted_nodes(opfm.sys)
    branches = get_sorted_branches(opfm.sys)
    idx = get_nodes_idx(nodes)
    list_gen = make_list(opfm, get_ctrl_generation)
    contanal = iterative_cont_anal(opfm.sys, nodes, branches, idx, contingencies)

    it = enumerate(contingencies)
    next = iterate(it)
    while next !== nothing
        (c, cont), state = next
        overloads = get_overload(contanal, c, get_net_Pᵢ(opfm, nodes))
        overloads = [(i,ol) for (i,ol) in enumerate(overloads) if abs(ol) > 1e-6]
        if length(overloads) > 0
            Pᵢ = get_Pᵢ(opfm, nodes)
            ptdf = get_cont_ptdf(contanal, c)
            sort!(overloads, rev = true, by = x -> abs(x[2]))
            (i,ol) = first(overloads)

            @constraint(submod, overload, submod[:pfcc][cont] == 0)
            # submod = Model(Gurobi.Optimizer)
            # @variable(submod, 0 <= pgu[g in get_name.(get_ctrl_generation(opfm.sys))])
            # @variable(submod, 0 <= pgd[g in get_name.(get_ctrl_generation(opfm.sys))])
            # @objective(submod, Min, sum(pgu) + sum(pgd))
            # @constraint(submod, p_balance, sum(pgu) - sum(pgd) == 0)
            # expr = sum(
            #             ptdf[i, idx[n.number]] * 
            #             sum((submod[:pgu][g.name] - submod[:pgd][g.name] for g in list_gen[n.name]), init=0.0) 
            #             for n in nodes if abs(ptdf[i, idx[n.number]]) > 1e-6
            #         ) 
            # if ol < 0
            #     @constraint(submod, overload, expr >= -ol)
            # else
            #     @constraint(submod, overload, expr <= -ol)
            # end
            # MOI.set(submod, MOI.Silent(), true) # supress output from the solver
            solve_model!(submod)

            duals = get_low_dual.(submod[:pgu]) .- get_low_dual.(submod[:pgd])
            expr = JuMP.AffExpr(objective_value(submod))
            JuMP.add_to_expression!(expr, sum(
                    (sum((opfm.mod[:pg0][g.name] for g in list_gen[n.name]), init=0.0) - Pᵢ[ni]) * (d > 0.0 ? 1.0 : -1.0)
                    for (ni,(n,d)) in enumerate(zip(nodes, duals))
                ))
            @info "Contingency on $(cont) resulted in overload on $(branches[i].name) of $(abs(ol)) \nCut added: $(sprint_expr(expr))\n"
            set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
            JuMP.@constraint(opfm.mod, expr <= 0)

            MOI.set(opfm.mod, MOI.Silent(), true) # supress output from the solver
            opfm.mod = solve_model!(opfm.mod)
            termination_status(opfm.mod) != MOI.OPTIMAL && return opfm, contanal
            delete(submod, overload)
            unregister(submod, :overload)
        else
            next = iterate(it, state)
        end
    end
    return opfm, contanal
end

mutable struct IterativeDCContAnal
    lodf::Array{<:Real, 3}
    contingencies::Vector{String}
    linerating::Vector{<:Real}
end

""" Note about LODF / PTDF:
In our current implementation we will store the PTDF for the base case (i.e. no outages) in LODF[:,:,1], while 
for all contingencies we will store the N first order contingencies at indices LODF[:,:,i] for i in 2:N. 
    Second/third order contingencies are not yet covered but will likely be indexed in a similar fashion. 

It is worth noting the terminology use of LODF / PTDF. The LODF for e.g.  branch 1 is the same as the PTDF when 
branch 1 is outaged. As such, it may not make sense that, in the code below, the LODF[:,:,1] calculated as the 
PTDF of the base case is not accurate use of this terminology, but rather a simplification for us to avoid 
having to create a separate PTDF variable for the base case (and thus increase the number of variables...).

    """
function iterative_cont_anal(sys::System, nodes::Vector{Bus}, branches::Vector{<:Branch}, 
        idx::Dict{<:Any, <:Integer}, contingencies::Vector{String})
    
    numnodes = length(nodes)
    lodf = zeros(length(branches), numnodes, length(contingencies)+1)
    i_slack, slack = find_slack(nodes)
    lodf[:,:,1] = get_ptdf(branches, numnodes, idx, i_slack)
    
    for (i,cont) in enumerate(contingencies)
        branches_cont = [b for b in branches if b.name != cont]
        islands = get_islands(sys, branches_cont)
        if length(islands) == 1
            drop_cont = [j for (j,b) in enumerate(branches) if b.name != cont]
            lodf[drop_cont, :, i+1] = get_ptdf(branches_cont, numnodes, idx, i_slack)
                # The B-matrix can be altered instead of rebuilt from the ground up
        else
            # only the island with the slack bus is assumed stable
            for island in islands
                if slack.number ∈ get_number.(island[1]) # Nodes in island 1
                    island_idx = get_nodes_idx(island[1])
                    drop_id = sort!(collect(values(island_idx)))
                    drop_cont = [j for (j,b) in enumerate(branches) if b.name != cont && b.name ∈ get_name.(island[2])]
                    if length(drop_id) != 0 && length(drop_cont) != 0
                        lodf[drop_cont, drop_id, i+1] = get_ptdf(island[2], length(island[1]), island_idx, i_slack)
                    end
                    break
                end
            end
        end
    end
    return IterativeDCContAnal(lodf, contingencies, get_rate.(branches))
end

" Get the overload of all lines "
get_overload(contanal::IterativeDCContAnal, cont::Integer, Pᵢ::Vector{<:Real}) = 
    find_overload.(
            calculate_line_flows(get_cont_ptdf(contanal, cont), Pᵢ), 
            contanal.linerating
        )

get_overload(
        A::AbstractMatrix{<:Real}, 
        D::AbstractMatrix{<:Real}, 
        B::AbstractMatrix{<:Real}, 
        B_inv::AbstractMatrix{<:Real}, 
        δ::AbstractVector{<:Real}, 
        cont::Tuple{Integer, Integer}, 
        slack::Integer
    ) = 
        find_overload.(calc_pl.(A, D, [get_changed_angles(
                    B_inv[cont[1]], B_inv[cont[2]], B[cont[1], cont[2]], 
                    δ, cont[1], cont[2], slack)]
                ), 
                contanal.linerating
            )

" Get the PTDF corresponding to the contingency "
function get_cont_ptdf(contanal::IterativeDCContAnal, cont::Integer)
    @assert 0 <= cont <= length(contanal.contingencies)  # "No entry in LODF matrix for the given contingency. "
    # cont+1 as first element in LODF is PTDF for base case
    return contanal.lodf[:,:, cont+1]  
end

""" Find island(s) in the system returned in a nested Vector.
Each element of the Vector consists of two lists, nodes and branches, in that island. """
function get_islands(sys::System, branches::Vector{<:Branch})
    islands = Vector()
    visited_nodes = Vector{Bus}([branches[1].arc.from]) # start node on island 1 marked as visited
    visited_branches, new_nodes = find_connected(branches, first(visited_nodes))
        # all nodes connected are set as neighouring nodes not visited,
        # via visited branches

    bus_numbers = length(sys.bus_numbers)
    bus_number = 0
    while true

        # Visit new nodes until there are no neighouring nodes connected
        while !isempty(new_nodes)
            node = pop!(new_nodes)
            push!(visited_nodes, node)
            bn = find_connected(branches, node)
            union!(visited_branches, bn[1])
            union!(new_nodes, setdiff(bn[2], visited_nodes))
        end

        push!(islands, (sort_nodes!(visited_nodes), sort_branches!(visited_branches)))
        bus_number += length(visited_nodes)

        bus_numbers == bus_number && break # all nodes are visited 
        bus_numbers < bus_number && @error "More nodes counted, $(bus_number), than nodes in the system, $(bus_numbers)!"

        visited_nodes = Vector([setdiff(get_components(Bus, sys), visited_nodes) |> first])
                    # a random node not visited yet starts a new island
        visited_branches, new_nodes = find_connected(branches, first(visited_nodes))
    end
    return islands
end

function get_islands(nodes::Vector{Integer}, branches::Vector{<:Tuple{Integer, Integer}})
    islands = Vector()
    visited_nodes = Vector{Integer}([branches[1][1]]) # start node on island 1 marked as visited
    visited_branches, new_nodes = find_connected(branches, first(visited_nodes))
        # all nodes connected are set as neighouring nodes not visited,
        # via visited branches

    bus_numbers = length(nodes)
    bus_number = 0
    while true

        # Visit new nodes until there are no neighouring nodes connected
        while !isempty(new_nodes)
            node = pop!(new_nodes)
            push!(visited_nodes, node)
            bn = find_connected(branches, node)
            union!(visited_branches, bn[1])
            union!(new_nodes, setdiff(bn[2], visited_nodes))
        end

        push!(islands, (sort_nodes!(visited_nodes), sort_branches!(visited_branches)))
        bus_number += length(visited_nodes)

        bus_numbers == bus_number && break # all nodes are visited 
        bus_numbers < bus_number && @error "More nodes counted, $(bus_number), than nodes in the system, $(bus_numbers)!"

        visited_nodes = Vector([setdiff(nodes, visited_nodes) |> first])
                    # a random node not visited yet starts a new island
        visited_branches, new_nodes = find_connected(branches, first(visited_nodes))
    end
    return islands
end

" Find nodes connected to a node "
function find_connected(branches::Vector{<:Branch}, node::Bus)
    new_branches = Vector{Branch}()
    new_nodes = Vector{Bus}()
    for branch in branches
        if branch.arc.from.number === node.number
            push!(new_branches, branch)
            push!(new_nodes, branch.arc.to)
        elseif branch.arc.to.number === node.number
            push!(new_branches, branch)
            push!(new_nodes, branch.arc.from)
        end
    end
    return new_branches, new_nodes
end

function find_connected(branches::Vector{<:Tuple{Integer, Integer}}, node::Integer)
    new_branches = Vector{Tuple{Int, Int}}()
    new_nodes = Vector{Int}()
    for branch in branches
        if branch[1] == node
            push!(new_branches, branch)
            push!(new_nodes, branch[2])
        elseif branch[2] == node
            push!(new_branches, branch)
            push!(new_nodes, branch[1])
        end
    end
    return new_branches, new_nodes
end

" An AffExpr nicely formatted to a string "
sprint_expr(expr::AffExpr) = join(Printf.@sprintf("%s%5.2f %s ", (x[2] > 0 ? "+" : "-"), abs(x[2]), x[1]) for x in expr.terms) * 
    Printf.@sprintf("<= %s%.2f", (expr.constant > 0 ? "-" : "+"), abs(expr.constant))
