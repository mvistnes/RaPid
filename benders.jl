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
function run_benders(opfm::OPFmodel, voll, contingencies, prob;
            time_limit_sec::Int64 = 600,
            unit_commit::Bool = false,
            max_shed::Float64 = 0.1,
            max_curtail::Float64 = 1.0,
            ratio::Float64= 0.5, 
            circuit_breakers::Bool=false,
            short_term_limit_multi::Float64 = 1.5,
            ramp_minutes::Int64 = 10,
            repair_time::Float64 = 1.0
        )

    # Set global variables
    nodes = get_sorted_nodes(opfm.sys)
    branches = get_sorted_branches(opfm.sys)
    drop_idx = get_nodes_idx([x for x in nodes if x.bustype != BusTypes.REF])
    drop_id = sort!(collect(keys(drop_idx)))
    list_gen = make_list(opfm, get_ctrl_generation)
    contanal = iterative_cont_anal(opfm.sys, nodes, branches, drop_idx, contingencies)

    it = enumerate(contingencies)
    next = iterate(it)
    while next !== nothing
        (c, cont), state = next
        Pᵢ = get_Pᵢ(opfm, nodes)
        ptdf = get_cont_ptdf(contanal, c)
        overloads = get_overload(contanal, c, get_net_Pᵢ(opfm, nodes)[drop_id])
        overloads = [(i,ol) for (i,ol) in enumerate(overloads) if abs(ol) > 0.001]
        if length(overloads) > 0
            sort!(overloads, rev = true, by = x -> abs(x[2]))
            (i,ol) = first(overloads)
            submod = Model(Gurobi.Optimizer)
            @variable(submod, 0 <= pgu[g in get_name.(get_ctrl_generation(opfm.sys))])
            @variable(submod, 0 <= pgd[g in get_name.(get_ctrl_generation(opfm.sys))])
            @objective(submod, Min, sum(pgu) + sum(pgd))
            @constraint(submod, sum(pgu) - sum(pgd) == 0)
            expr = sum(
                        ptdf[i, drop_idx[n.number]] * 
                        sum((submod[:pgu][g.name] - submod[:pgd][g.name] for g in list_gen[n.name]), init=0.0) 
                        for n in nodes if n.bustype != BusTypes.REF
                    ) 
            if ol > 0
                @constraint(submod, expr <= ol + contanal.linerating[i])
            else
                @constraint(submod, expr >= abs(ol))
            end
            solve_model!(submod)

            expr = JuMP.AffExpr(objective_value(submod))
            JuMP.add_to_expression!(expr, sum(
                    (sum((opfm.mod[:pg0][g.name] for g in list_gen[n.name]), init=0.0) - Pᵢ[n.number]) *
                    (n.bustype != BusTypes.REF && abs(ptdf[i, drop_idx[n.number]]) > 0.001 ? -1 : 1) for n in nodes 
                ))
            @info "Contingency on $(cont) resulted in overload on $(branches[i].name) of $(ol) \nCut added: $(sprint_expr(expr))\n"
            set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
            JuMP.@constraint(opfm.mod, expr <= 0)
            MOI.set(opfm.mod, MOI.Silent(), true) # supress output from the solver
            opfm.mod = solve_model!(opfm.mod)
            termination_status(opfm.mod) != MOI.OPTIMAL && return
        else
            next = iterate(it, state)
        end
    end
    return opfm
end

mutable struct IterativeDCContAnal
    lodf::Array{Float64, 3}
    contingencies::Vector{String}
    linerating::Vector{Float64}
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
function iterative_cont_anal(sys::System, nodes::Vector{Bus}, branches::Vector{<: Branch}, drop_idx::Dict{<: Any, <: Int}, contingencies::Vector{String})
    lodf = zeros(
            length(branches), 
            length(nodes)-1, 
            length(contingencies)+1
        )
    i_slack, slack = find_slack(nodes)
    lodf[:,:,1] = get_ptdf(nodes, branches, drop_idx, slack.number)
    
    for (i,cont) in enumerate(contingencies)
        branches_cont = [b for b in branches if b.name != cont]
        islands = get_islands(sys, branches_cont)
        if length(islands) == 1
            drop_cont = [j for (j,b) in enumerate(branches) if b.name != cont]
            lodf[drop_cont, :, i+1] = get_ptdf(nodes, branches_cont, drop_idx, slack.number)
                # The B-matrix can be altered instead of rebuilt from the ground up
        elseif length(islands) == 2
            if slack.number ∈ get_number.(islands[1][1]) # Nodes in island 1
                island_idx = get_nodes_idx([x for x in islands[1][1] if x.bustype != BusTypes.REF])
                drop_id = sort!(collect(values(island_idx)))
                drop_cont = [j for (j,b) in enumerate(branches) if b.name != cont && b.name ∈ get_name.(islands[1][2])]
                lodf[drop_cont, drop_id, i+1] = get_ptdf(islands[1][1], islands[1][2], island_idx, slack.number)
            else
                island_idx = get_nodes_idx([x for x in islands[2][1] if x.bustype != BusTypes.REF])
                drop_id = sort!(collect(values(island_idx)))
                drop_cont = [j for (j,b) in enumerate(branches) if b.name != cont && b.name ∈ get_name.(islands[2][2])]
                lodf[drop_cont, drop_id, i+1] = get_ptdf(islands[2][1], islands[2][2], island_idx, slack.number)
            end
        else
            @warn "More than two islands is not implemented!"
        end
    end
    return IterativeDCContAnal(lodf, contingencies, get_rate.(branches))
end

" Get the overload of all lines "
get_overload(contanal::IterativeDCContAnal, cont::Int, Pᵢ::Vector{Float64}) = 
    find_overload.(
            calculate_line_flows(get_cont_ptdf(contanal, cont), Pᵢ), 
            contanal.linerating
        )

" Get the PTDF corresponding to the contingency "
function get_cont_ptdf(contanal::IterativeDCContAnal, cont::Int)
    @assert 0 <= cont <= length(contanal.contingencies)  # "No entry in LODF matrix for the given contingency. "
    # cont+1 as first element in LODF is PTDF for base case
    return contanal.lodf[:,:, cont+1]  
end

" Make the PTDF matrix for using the input nodes and branches "
function get_ptdf(nodes::Vector{Bus}, branches::Vector{<: Branch}, idx::Dict{<: Any, <: Int}, slack::Int)
    B = buildB(nodes, branches, idx, slack)
    isempty(B) && return B
    return get_ptdf(LinearAlgebra.lu(B), nodes, branches, idx, slack)
end
function get_ptdf(B, nodes::Vector{Bus}, branches::Vector{<: Branch}, idx::Dict{<: Any, <: Int}, slack::Int)
    A = zeros(Float64, size(branches,1), size(nodes,1)-1) # Container for the distribution factors
    for (i, branch) in enumerate(branches)
        branch.x == 0 && continue
        ΔP = zeros(Float64, size(nodes,1)-1)
        (f, t) = get_bus_id(branch) # (f)rom and (t)o bus at this branch
        if f != slack
            ΔP[idx[f]] = 1 / branch.x
        end
        if t != slack # ∈ keys(idx)
            ΔP[idx[t]] = -1 / branch.x
        end
        A[i, :] = B \ ΔP # append factors to matrix
    end
    return A
end

"""
Builds an admittance matrix with the line series reactance of the lines.
idx must not include the slack bus.
ToDo: Only need the upper half triagonal matrix as it's symetric.
"""
function buildB(nodes::Vector{Bus}, branches::Vector{<: Branch}, idx::Dict{<: Any, <: Int}, slack::Int)
    B = SparseArrays.spzeros(size(nodes,1)-1, size(nodes,1)-1)
    for branch in branches
        (f, t) = get_bus_id(branch)
        if f != slack # not BusTypes.REF
            B[idx[f],idx[f]] += 1 / branch.x
            if t != slack
                B[idx[t],idx[t]] += 1 / branch.x
                B[idx[f],idx[t]] = -1 / branch.x
                B[idx[t],idx[f]] = -1 / branch.x
            end
        elseif t != slack
            B[idx[t],idx[t]] += 1 / branch.x
        end
    end
    return B
end

""" Find island(s) in the system returned in a nested Vector.
Each element of the Vector consists of two lists, nodes and branches, in that island. """
function get_islands(sys::System, branches::Vector{<: Branch})
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

        push!(islands, (visited_nodes, visited_branches))
        bus_number += length(visited_nodes)
        bus_numbers == bus_number && break # all nodes are visited 
        bus_numbers < bus_number && @error "More nodes counted, $(bus_number), than nodes in the system, $(bus_numbers)!"
        visited_nodes = Vector([get_component(
                Bus, 
                sys, 
                setdiff(sys.bus_numbers, get_number.(visited_nodes)) |> first |> string
                    # a random node not visited yet starts a new island
            )])
        visited_branches, new_nodes = find_connected(branches, first(visited_nodes))
    end
    return islands
end

" Find nodes connected to a node "
function find_connected(branches::Vector{<: Branch}, node::Bus)
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

" An AffExpr nicely formatted to a string "
sprint_expr(expr::AffExpr) = join(Printf.@sprintf("%s%5.2f %s ", (x[2] > 0 ? "+" : "-"), abs(x[2]), x[1]) for x in expr.terms) * 
    (expr.constant > 0 ? " +" : " -") * string(abs(expr.constant)) * " <= 0"



# -------------------------------------- OLD --------------------------------------

# function run_benders(opfm::OPFmodel, voll, contingencies, prob;
#             time_limit_sec::Int64 = 600,
#             unit_commit::Bool = false,
#             max_shed::Float64 = 0.1,
#             max_curtail::Float64 = 1.0,
#             ratio::Float64= 0.5, 
#             circuit_breakers::Bool=false,
#             short_term_limit_multi::Float64 = 1.5,
#             ramp_minutes::Int64 = 10,
#             repair_time::Float64 = 1.0
#         )

#     # Set global variables
#     nodes = get_sorted_nodes(opfm.sys)
#     branches = get_sorted_branches(opfm.sys)
#     drop_idx = get_nodes_idx([x for x in nodes if x.bustype != BusTypes.REF])
#     drop_id = sort!(collect(keys(drop_idx)))
#     list_gen = make_list(opfm, get_ctrl_generation)
#     contanal = iterative_cont_anal(opfm.sys, nodes, branches, drop_idx, contingencies)

#     it = enumerate(contingencies)
#     next = iterate(it)
#     while next !== nothing
#         (c, cont), state = next
#         Pᵢ = get_Pᵢ(opfm)
#         ptdf = get_cont_ptdf(contanal, c)
#         overloads = get_overload(contanal, c, Pᵢ.data[drop_id])

#         for (i,ol) in enumerate(overloads)
#             abs(ol) < 0.001 && continue
#             expr = JuMP.AffExpr(-ol)

#             JuMP.add_to_expression!(expr, sum(
#                 -ptdf[i, drop_idx[n.number]] * 
#                 (sum((opfm.mod[:pg0][g.name] for g in list_gen[n.name]), init=0.0) - Pᵢ[n.name]) 
#                 for n in nodes if n.bustype != BusTypes.REF
#             ))
            
#             @info "Contingency on $(cont) resulted in overload on $(branches[i].name) of $(ol) \nCut added: $(sprint_expr(expr))\n"
#             set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
#             JuMP.@constraint(opfm.mod, expr .<= 0)
#             MOI.set(opfm.mod, MOI.Silent(), true) # supress output from the solver
#             opfm.mod = solve_model!(opfm.mod)
#             break
#         end
#         if abs(ol) < 0.001
#             next = iterate(it, state)
#         end
#     end
#     return opfm
# end