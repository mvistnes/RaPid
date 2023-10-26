# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022 

# FIX: Sometimes breaks
""" Find island(s) in the system returned in a nested Vector.
Each element of the Vector consists of two lists, nodes and branches, in that island. """
function get_islands(sys::System, branches::AbstractVector{<:Branch})
    islands = Vector{Tuple{Vector{Bus},Vector{Branch}}}()
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
            # for i in islands
            #     j = ""
            #     for k in sort!(i[1], by = x -> x.name)
            #         if j == k
            #             println(j.name)
            #         else
            #             j=k
            #         end
            #     end
            #     #@show get_name.(i[2])
            # end
            throw("More nodes counted, $(bus_number), than nodes in the system, $(bus_number_tot)!")
        end

        visited_nodes = Vector([setdiff(get_components(Bus, sys), visited_nodes) |> first])
        # a random node not visited yet starts a new island
        visited_branches, new_nodes = find_connected(branches, first(visited_nodes))
    end
    return islands
end

# FIX: Sometimes breaks
function get_islands(nodes::AbstractVector{T}, branches::AbstractVector{T2}) where {T<:Integer,T2<:Tuple{Integer,Integer}}
    islands = Vector{Tuple{Vector{T},Vector{T2}}}()
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

        # println(visited_nodes, visited_branches)
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

function get_islands(T::SparseArrays.SparseMatrixCSC{<:Any,Ty}) where {Ty<:Integer}
    islands = Vector{Vector{Ty}}()
    visited_nodes = Vector{Ty}() # start node on island 1 marked as visited
    new_nodes = T[:, 1].nzind
    # all nodes connected are set as neighouring nodes not visited,
    # via visited branches

    bus_number_tot = size(T, 1)
    bus_number = 0
    while true

        # Visit new nodes until there are no neighouring nodes connected
        while !isempty(new_nodes)
            node = pop!(new_nodes)
            push!(visited_nodes, node)
            bn = T[:, node].nzind
            union!(new_nodes, setdiff(bn, visited_nodes))
        end

        # println(visited_nodes, visited_branches)
        push!(islands, visited_nodes)
        # push!(islands, (sort!(visited_nodes), sort!(visited_branches)))
        bus_number += length(visited_nodes)

        bus_number_tot == bus_number && break # all nodes are visited 
        bus_number_tot < bus_number && throw("More nodes counted, $(bus_number), than nodes in the system, $(bus_number_tot)!")

        visited_nodes = Vector([setdiff(nodes, visited_nodes) |> first])
        # a random node not visited yet starts a new island
        new_nodes = T[:, first(visited_nodes)].nzind
    end
    return islands
end

function is_islanded(
    DA::AbstractMatrix{T},
    B::AbstractMatrix{T},
    X::AbstractMatrix{T},
    cont::Tuple{Integer,Integer},
    branch::Integer;
    atol::Real=1e-5
) where {T<:Real}
    (fbus, tbus) = cont
    return isapprox(
        1 / B[fbus, tbus] + DA[branch, tbus] / B[fbus, tbus] *
                            ((X[fbus, fbus] - X[fbus, tbus]) - (X[tbus, fbus] - X[tbus, tbus])),
        zero(T), atol=atol)
end
is_islanded(pf::DCPowerFlow, cont::Tuple{Integer,Integer}, branch::Integer; atol::Real=1e-5) =
    is_islanded(pf.DA, pf.B, pf.X, cont, branch, atol=atol)

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

function find_connected(branches::AbstractVector{<:Tuple{Integer,Integer}}, node::Integer)
    new_branches = Vector{Tuple{Int,Int}}()
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

function find_connected(A::SparseArrays.SparseMatrixCSC{<:Integer,<:Integer}, branches::AbstractVector{<:Branch}, node::Integer)
    (n, b) = SparseArrays.findnz(A[:, node])
    new_branches = branches[n]
    return new_branches, union!(get_to.(get_arc.(new_branches)), get_from.(get_arc.(new_branches)))
end

"""
    Find islands after a contingency.
    Removes all connection between the two nodes i and j.
"""
function island_detection(T::SparseArrays.SparseMatrixCSC{Ty,<:Integer}, i::Integer, j::Integer; atol::Real=1e-5) where {Ty}
    val = T[i, j]
    SparseArrays.dropstored!(T, i, j)
    SparseArrays.dropstored!(T, j, i)
    if isapprox(T[i, i], val, atol=atol)
        SparseArrays.dropstored!(T, i, i)
        SparseArrays.dropstored!(T, j, j)
    else
        T[i, i] += val
        T[j, j] += val
    end
    # Values need to be exactly zero (or dropped) to be eliminated from the SparseMatrix in the algorithm
    islands = island_detection(create_connectivity_matrix(T))
    T[i, j] += val
    T[j, i] += val
    T[i, i] -= val
    T[j, j] -= val
    if sum(length(i) for i in islands) != size(T, 1)
        @warn "Counted nodes $(sum(length(i) for i in islands)) != total nodes $(size(T, 1)) in contingency of line $i-$j !"
    end
    return islands
end
function island_detection_thread_safe(Y::SparseArrays.SparseMatrixCSC{Ty,<:Integer}, i::Integer, j::Integer; atol::Real=1e-5) where {Ty}
    T = create_connectivity_matrix(Y)
    val = Y[i, j]
    SparseArrays.dropstored!(T, i, j)
    SparseArrays.dropstored!(T, j, i)
    isapprox(Y[i, i], val, atol=atol) && SparseArrays.dropstored!(T, i, i)
    isapprox(Y[j, j], val, atol=atol) && SparseArrays.dropstored!(T, j, j)
    # Values need to be exactly zero (or dropped) to be eliminated from the SparseMatrix in the algorithm
    islands = island_detection(T)
    if sum(length(i) for i in islands) != size(T, 1)
        @warn "Counted nodes $(sum(length(i) for i in islands)) != total nodes $(size(T, 1)) in contingency of line $i-$j !"
    end
    return islands
end
function island_detection_thread_safe(Y::SparseArrays.SparseMatrixCSC{Ty,<:Integer}, cont::AbstractVector{<:Tuple{Integer,Integer}}; atol::Real=1e-5) where {Ty}
    T = create_connectivity_matrix(Y)
    for (i,j) in cont
        val = Y[i, j]
        SparseArrays.dropstored!(T, i, j)
        SparseArrays.dropstored!(T, j, i)
        isapprox(Y[i, i], val, atol=atol) && SparseArrays.dropstored!(T, i, i)
        isapprox(Y[j, j], val, atol=atol) && SparseArrays.dropstored!(T, j, j)
    end
    # Values need to be exactly zero (or dropped) to be eliminated from the SparseMatrix in the algorithm
    islands = island_detection(T)
    if sum(length(i) for i in islands) != size(T, 1)
        @warn "Counted nodes $(sum(length(i) for i in islands)) != total nodes $(size(T, 1)) in contingency of line $i-$j !"
    end
    return islands
end

function find_ref_island(islands::Vector, slack::Integer)
    for (i, island) in enumerate(islands)
        if slack ∈ island
            return i
        end
    end
    return 0
end

function find_island_branches(island::Vector{<:Integer}, DA::SparseArrays.SparseMatrixCSC{<:Any,<:Integer}, c_branch::Integer)
    res = sort!(unique!(DA[:, island].rowval))
    if insorted(c_branch, res)
        deleteat!(res, searchsortedfirst(res, c_branch))
    end
    return res
end

function handle_islands(B::AbstractMatrix, DA::AbstractMatrix, contingency::Tuple{Integer,Integer}, branch::Integer, slack::Integer)
    islands = island_detection_thread_safe(B, contingency[1], contingency[2])
    island = find_ref_island(islands, slack)
    island_b = find_island_branches(islands[island], DA, branch)

    # Need at least one node to make a valid system
    if length(islands[island]) > 0
        return islands, island, island_b
    end
    return Vector{Vector{Int64}}(undef, 0), 0, Vector{Int64}(undef, 0)
end

function handle_islands(B::AbstractMatrix, contingency::Tuple{Integer,Integer}, slack::Integer)
    islands = island_detection_thread_safe(B, contingency[1], contingency[2])
    island = find_ref_island(islands, slack)
    return length(islands[island]) > 0 ? (islands, islands) : Vector{Vector{Int64}}(undef, 0), 0
    # Need at least one node to make a valid system
end

" Get all islands with the reference bus from all bx contingencies "
function get_all_islands(B::AbstractMatrix, bx::Vector{<:Tuple{Integer,Integer}}, slack::Integer)
    res = Vector{Vector{Int64}}(undef, length(bx))
    Threads.@threads for i in eachindex(bx)
        islands. island = handle_islands(B, bx[i], slack)
        res[i] = islands[island]
    end
    return res
end

function get_all_islands(opf::OPFsystem, slack::Integer)
    idx = get_nodes_idx(opf.nodes)
    bx = get_bus_idx.(opf.branches, [idx])
    cx = get_bus_idx.(opf.contingencies, [idx])
    A = calc_A(bx, length(opf.nodes))
    return get_all_islands(A' * A, cx, slack)
end

function find_islands(T::SparseArrays.SparseMatrixCSC{<:Any,<:Integer})
    n_buses = size(T, 1)
    bus_islands = []
    nodes = []
    node = 1
    n_isolated_buses = 0
    # Find the isolated buses
    for (i, bus) in enumerate(SparseArrays.diag(T))
        if bus == 0.0
            push!(bus_islands, [i])
            push!(nodes, i)
            n_isolated_buses += 1
            # In case the first buses are isolated
            if node == i
                node += 1
            end
        end
    end
    while true
        n_nodes = new_islands(T, node)
        push!(bus_islands, n_nodes)
        append!(nodes, n_nodes)
        if n_buses - n_isolated_buses > sum(length.(bus_islands), init=0)
            node = first(setdiff(1:n_buses, nodes))
        else
            break
        end
    end
    return bus_islands
end


function new_islands(T::SparseArrays.SparseMatrixCSC{<:Any,<:Integer}, node::Integer)
    island = T[:, node].nzind
    for i in island
        union!(island, T[:, i].nzind)
        # print(i, " ")
    end
    return island
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
    T = SparseArrays.dropzeros!(SparseArrays.sparse(I, J, val .!= 0))
end

"""
    island_detection(T)

    Uses Goderya's algorithm to detect islands [goderya_fast_1980}(@cite).
    The implmentation is adapted from [milano_power_2010](@cite).
    
    Inputs:
        T: Connectivity matrix.

    Outputs:
        List of buses in the same island.
"""
function island_detection(T::SparseArrays.SparseMatrixCSC{<:Any,Ty}) where {Ty<:Integer}
    n_buses = size(T, 1)
    bus_islands = Vector{Vector{Int64}}(undef, 0)
    islands = 0
    conn = SparseArrays.spzeros(Bool, n_buses)
    idx = 1
    n_isolated_buses = 0

    # Find the isolated buses
    for (i, bus) in enumerate(SparseArrays.diag(T))
        if bus == 0
            push!(bus_islands, [i])
            n_isolated_buses += 1
            # In case the first buses are isolated
            if idx == i
                idx += 1
            end
        end
    end

    row = T[idx, :]
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
        if length(conn.nzind) >= (n_buses - n_isolated_buses)
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
        row = T[idx, :]
        nelm = SparseArrays.nnz(row)
    end
    return bus_islands
end
