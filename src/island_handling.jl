# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022 

using PowerSystems
import SparseArrays


""" Find island(s) in the system returned in a nested Vector.
Each element of the Vector consists of two lists, nodes and branches, in that island. """
function get_islands(sys::System, branches::AbstractVector{<:Branch})
    islands = Vector{Tuple{Vector{Bus}, Vector{Branch}}}()
    visited_nodes = Vector{Bus}([branches[1].arc.from]) # start node on island 1 marked as visited
    visited_branches, new_nodes = find_connected(branches, first(visited_nodes))
        # all nodes connected are set as neighouring nodes not visited,
        # via visited branches

    bus_number_tot = length(sys.bus_numbers)
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
        # push!(islands, (sort_nodes!(visited_nodes), sort_branches!(visited_branches)))
        bus_number += length(visited_nodes)

        bus_number_tot == bus_number && break # all nodes are visited 
        if bus_number_tot < bus_number 
            for i in islands
                j = ""
                for k in sort!(i[1], by = x -> x.name)
                    if j == k
                        println(j.name)
                    else
                        j=k
                    end
                end
                #@show get_name.(i[2])
            end
            throw("More nodes counted, $(bus_number), than nodes in the system, $(bus_number_tot)!")
        end

        visited_nodes = Vector([setdiff(get_components(Bus, sys), visited_nodes) |> first])
                    # a random node not visited yet starts a new island
        visited_branches, new_nodes = find_connected(branches, first(visited_nodes))
    end
    return islands
end

function get_islands(nodes::AbstractVector{T}, branches::AbstractVector{T2}) where {T<:Integer, T2<:Tuple{Integer, Integer}}
    islands = Vector{Tuple{Vector{T}, Vector{T2}}}()
    visited_nodes = Vector{T}([branches[1][1]]) # start node on island 1 marked as visited
    visited_branches, new_nodes = find_connected(branches, first(visited_nodes))
        # all nodes connected are set as neighouring nodes not visited,
        # via visited branches

    bus_number_tot = length(nodes)
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

        println(visited_nodes, visited_branches)
        push!(islands, (visited_nodes, visited_branches))
        # push!(islands, (sort!(visited_nodes), sort!(visited_branches)))
        bus_number += length(visited_nodes)

        bus_number_tot == bus_number && break # all nodes are visited 
        bus_number_tot < bus_number && throw("More nodes counted, $(bus_number), than nodes in the system, $(bus_number_tot)!")

        visited_nodes = Vector([setdiff(nodes, visited_nodes) |> first])
                    # a random node not visited yet starts a new island
        visited_branches, new_nodes = find_connected(branches, first(visited_nodes))
    end
    return islands
end

" Find nodes connected to a node "
function find_connected(branches::AbstractVector{<:Branch}, node::Bus)
    new_branches = Vector{Branch}()
    new_nodes = Vector{Bus}()
    for branch in branches
        if branch.arc.from === node
            push!(new_branches, branch)
            push!(new_nodes, branch.arc.to)
        elseif branch.arc.to === node
            push!(new_branches, branch)
            push!(new_nodes, branch.arc.from)
        end
    end
    return new_branches, new_nodes
end

function find_connected(branches::AbstractVector{<:Tuple{Integer, Integer}}, node::Integer)
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

function find_connected(A::SparseArrays.SparseMatrixCSC{<:Integer, <:Integer}, branches::AbstractVector{<:Branch}, node::Integer)
    (n, b) = SparseArrays.findnz(A[:,node])
    new_branches = branches[n]
    return new_branches, union!(get_to.(get_arc.(new_branches)), get_from.(get_arc.(new_branches)))
end

#################################################################################################
# The code under this line is created by Sigurd Jakobsen (2023) in SINT_LF
#################################################################################################


"""
    create_connectivity_matrix(Y)
    Creates the connectivity matrix exaplined in [goderya_fast_1980](@cite).
    The code for creating the matrix is inspired by 
    https://github.com/NREL-SIIP/PowerSystems.jl
"""
function create_connectivity_matrix(Y::SparseArrays.SparseMatrixCSC)
    I, J, val = SparseArrays.findnz(Y)
    T = SparseArrays.dropzeros!(SparseArrays.sparse(I, J, val.!=0))
end

"""
    island_detection(T, isolated_buses)

    Uses Goderya's algorithm to detect islands [goderya_fast_1980}(@cite).
    The implmentation is adapted from [milano_power_2010](@cite).
    
    Inputs:
        T: Connectivity matrix.
        mul: Function used for the multiplication.

    Outputs:
        List of buses in the same island.
"""
function island_detection(T::SparseArrays.SparseMatrixCSC{<:Any, Ty}) where Ty <: Integer
    n_buses = size(T, 1)
    bus_islands = Vector{Vector{Ty}}()
    islands = 0
    conn = SparseArrays.spzeros(Bool, n_buses)
    idx = 1 
    n_isolated_buses = 0
    
    # Find the isolated buses
    for (i, bus) in enumerate(SparseArrays.diag(T))
        if bus == 0
            append!(bus_islands, i)
            n_isolated_buses += 1
            # In case the first buses are isolated
            if idx == i
                idx += 1
            end
        end
    end
    
    row = T[:, idx]
    nelm = SparseArrays.nnz(row) 
    
    while true
        while true
            row = T * row
            # row = mul(T, row)
            new_nelm = SparseArrays.nnz(row)
            if new_nelm == nelm
                break
            end
            nelm = new_nelm
        end
        append!(bus_islands, [row.nzind])
        conn += row 
        islands += 1
        if length(conn.nzind) >= (n_buses-n_isolated_buses)
            break
        end
        
        for element in (idx+1):(idx+length(conn.nzind))
            idx += 1
            if T[element, element] == false
                append!(bus_islands, [element])
                continue
            end
            if element ∉ conn.nzind
                break
            end
        end
        row = T[:, idx]
        nelm = SparseArrays.nnz(row) 
    end
    return bus_islands
end

function island_detection(T::SparseArrays.SparseMatrixCSC{Ty, <:Integer}, i::Integer, j::Integer) where Ty
    val = T[i, j]
    T[i, j] -= val
    T[j, i] -= val
    islands = island_detection(T)
    T[i, j] += val
    T[j, i] += val
    if sum(length(i) for i in islands) != size(T, 1) 
        @warn "Counted nodes $(sum(length(i) for i in islands)) != total nodes $size(T, 1) in contingency of line $i-$j !"
    end
    return islands
end
