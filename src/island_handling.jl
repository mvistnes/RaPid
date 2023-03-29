# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022 

using PowerSystems
import SparseArrays


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

function get_islands(nodes::Vector{T}, branches::Vector{T2}) where {T<:Integer, T2<:Tuple{Integer, Integer}}
    islands = Vector{Tuple{Vector{T}, Vector{T2}}}()
    visited_nodes = Vector{T}([branches[1][1]]) # start node on island 1 marked as visited
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

        push!(islands, (visited_nodes, visited_branches))
        # push!(islands, (sort!(visited_nodes), sort!(visited_branches)))
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

function find_connected(A::SparseArrays.SparseMatrixCSC{<:Integer, <:Integer}, branches::Vector{<:Branch}, node::Integer)
    (n, b) = SparseArrays.findnz(A[:,node])
    new_branches = branches[n]
    return new_branches, union!(get_to.(get_arc.(new_branches)), get_from.(get_arc.(new_branches)))
end
