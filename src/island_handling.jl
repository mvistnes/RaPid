"""
    Fast check for islanding in DC power flow

    Arguments:
    - `DA::AbstractMatrix{T}`: The distribution factor matrix
    - `B::AbstractMatrix{T}`: The susceptance matrix
    - `X::AbstractMatrix{T}`: The inverse of the susceptance matrix
    - `cont::Tuple{Integer,Integer}`: The contingency as a tuple of (from, to) bus indices
    - `branch::Integer`: The index of the branch in the DA and B matrices
    - `atol::Real=1e-14`: Absolute tolerance for floating point comparison
"""
function is_islanded(
    DA::AbstractMatrix{T},
    B::AbstractMatrix{T},
    X::AbstractMatrix{T},
    cont::Tuple{Integer,Integer},
    branch::Integer;
    atol::Real=1e-14
) where {T<:Real}
    (fbus, tbus) = cont
    return isapprox(
        (1 + DA[branch, tbus] * ((X[fbus, fbus] - X[fbus, tbus]) - (X[tbus, fbus] - X[tbus, tbus]))) / B[fbus, tbus],
        zero(T), atol=atol)
end

"""
    Fast check for islanding in DC power flow.
    
    Arguments:
    - `pf::DCPowerFlow`: The DC power flow object
    - `cont::Tuple{Integer,Integer}`: The contingency as a tuple of (from, to) bus indices
    - `branch::Integer`: The index of the branch in the DA and B matrices
    - `atol::Real=1e-14`: Absolute tolerance for floating point comparison
"""
is_islanded(pf::DCPowerFlow, cont::Tuple{Integer,Integer}, branch::Integer; atol::Real=1e-14) =
    is_islanded(pf.DA, pf.B, pf.X, cont, branch, atol=atol)

"""
    Fast check for islanding in DC power flow.
    
    Always returns false as no islanding is possible with single bus contingencies.

    Arguments:
    - `pf::DCPowerFlow`: The DC power flow object
    - `cont::Integer`: The contingency as a bus index
    - `branch::Integer`: The index of the branch in the DA and B matrices
    - `atol::Real=1e-14`: Absolute tolerance for floating point comparison
"""
is_islanded(pf::DCPowerFlow, cont::Integer, branch::Integer; atol::Real=1e-14) = false

""" 
    Find nodes connected to a node 

    Arguments:
    - `branches::AbstractVector{<:Branch}`: List of branches in the system
    - `node::Bus`: The bus to find connections for
"""
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

""" 
    Find nodes connected to a node 

    Arguments:
    - `branches::AbstractVector{<:Branch}`: List of branches in the system
    - `node::Bus`: The bus to find connections for
"""
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

""" 
    Find nodes connected to a node 

    Arguments:
    - `A::SparseArrays.SparseMatrixCSC{<:Integer,<:Integer}`: The adjacency matrix of the system
    - `branches::AbstractVector{<:Branch}`: List of branches in the system
    - `node::Bus`: The bus to find connections for
"""
function find_connected(A::SparseArrays.SparseMatrixCSC{<:Integer,<:Integer}, branches::AbstractVector{<:Branch}, node::Integer)
    (n, b) = SparseArrays.findnz(A[:, node])
    new_branches = branches[n]
    return new_branches, union!(get_to.(get_arc.(new_branches)), get_from.(get_arc.(new_branches)))
end

"""
    Find nodes connected to a node 

    Arguments:
    - `system::System`: The power system
"""
function island_detection(system::System)
    branches = sort_components!(get_branches(system))
    nodes = sort_components!(get_nodes(system))
    A = calc_A(branches, length(nodes), get_nodes_idx(nodes))
    return island_detection(create_connectivity_matrix(A'A))
end

"""
    Find islands after a contingency.
    Removes all connection between the two nodes i and j.

    Arguments:
    - `T::SparseArrays.SparseMatrixCSC{Ty,<:Integer}`: The admittance or connectivity matrix of the system
    - `i::Integer`: The from bus index of the contingency
    - `j::Integer`: The to bus index of the contingency
    - `atol::Real=1e-14`: Absolute tolerance for floating point comparison
"""
function island_detection(T::SparseArrays.SparseMatrixCSC{Ty,<:Integer}, i::Integer, j::Integer; atol::Real=1e-14) where {Ty}
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

"""
    Find islands after a contingency, thread safe by not altering T.

    Arguments:
    - `T::SparseArrays.SparseMatrixCSC{Ty,<:Integer}`: The admittance or connectivity matrix of the system
    - `i::Integer`: The from bus index of the contingency
    - `j::Integer`: The to bus index of the contingency
    - `atol::Real=1e-14`: Absolute tolerance for floating point comparison
"""
function island_detection_thread_safe(Y::SparseArrays.SparseMatrixCSC{Ty,<:Integer}, i::Integer, j::Integer; atol::Real=1e-14) where {Ty}
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

"""
    Find islands after contingencies, thread safe by not altering T.

    Arguments:
    - `T::SparseArrays.SparseMatrixCSC{Ty,<:Integer}`: The admittance or connectivity matrix of the system
    - `cont::AbstractVector{<:Tuple{Integer,Integer}}`: The contingencies as a vector of (from, to) bus indices
    - `atol::Real=1e-14`: Absolute tolerance for floating point comparison
"""
function island_detection_thread_safe(Y::SparseArrays.SparseMatrixCSC{Ty,<:Integer}, cont::AbstractVector{<:Tuple{Integer,Integer}}; atol::Real=1e-14) where {Ty}
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

"""
    Find the island containing the slack bus

    Arguments:
    - `islands::Vector`: A vector of islands, each island being a vector of bus indices
    - `slack::Integer`: The index of the slack bus
"""
function find_ref_island(islands::Vector, slack::Integer)
    for (i, island) in enumerate(islands)
        if slack ∈ island
            return i
        end
    end
    return 0
end

"""
    Find branches in an island, excluding the outaged branch

    Arguments:
    - `island::Vector{<:Integer}`: A vector of bus indices in the island
    - `DA::SparseArrays.SparseMatrixCSC{<:Any,<:Integer}`: The branch-bus incidence matrix times diagonal matrix of branch susceptances
    - `c_branch::Integer`: The index of the outaged branch to exclude
"""
function find_island_branches(island::Vector{<:Integer}, DA::SparseArrays.SparseMatrixCSC{<:Any,<:Integer}, c_branch::Integer)
    res = sort!(unique!(DA[:, island].rowval))
    if insorted(c_branch, res)
        deleteat!(res, searchsortedfirst(res, c_branch))
    end
    return res
end

"""
    Find branches in an island, excluding the outaged branch

    Arguments:
    - `island::Vector{<:Integer}`: A vector of bus indices in the island
    - `DA::SparseArrays.SparseMatrixCSC{<:Any,<:Integer}`: The branch-bus incidence matrix times diagonal matrix of branch susceptances
    - `c_branchs::AbstractVector{<:Integer}`: The indices of outaged branches to exclude
"""
function find_island_branches(island::Vector{<:Integer}, DA::SparseArrays.SparseMatrixCSC{<:Any,<:Integer}, c_branches::AbstractVector{<:Integer})
    res = sort!(unique!(DA[:, island].rowval))
    for br in c_branches
        if insorted(br, res)
            deleteat!(res, searchsortedfirst(res, br))
        end
    end
    return res
end

"""
    Handle islands after a contingency.
    Returns all islands, the index of the island with the slack bus, and the branches in that island.

    Arguments:
    - `B::AbstractMatrix`: The susceptance matrix
    - `DA::AbstractMatrix`: The branch-bus incidence matrix times diagonal matrix of branch susceptances
    - `cont::Tuple{Integer,Integer}`: The contingency as a tuple of (from, to) bus indices
    - `branch::Integer`: The index of the branch in the DA and B matrices
    - `slack::Integer`: The index of the slack bus
"""
function handle_islands(B::AbstractMatrix, DA::AbstractMatrix, contingency::Tuple{Integer,Integer}, branch::Integer, slack::Integer)
    islands = island_detection(B, contingency[1], contingency[2])
    island = find_ref_island(islands, slack)
    island_b = find_island_branches(islands[island], DA, branch)

    # Need at least one node to make a valid system
    if length(islands[island]) > 0
        return islands, island, island_b
    end
    return Vector{Vector{Int64}}(undef, 0), 0, Vector{Int64}(undef, 0)
end

"""
    Handle islands after a contingency.
    Returns all islands, the index of the island with the slack bus, and the branches in that island.

    Arguments:
    - `B::AbstractMatrix`: The susceptance matrix
    - `DA::AbstractMatrix`: The branch-bus incidence matrix times diagonal matrix of branch susceptances
    - `contingencies::AbstractVector{<:Tuple{Integer,Integer}}`: The contingencies as a vector of tuples of (from, to) bus indices
    - `branch::Integer`: The index of the branch in the DA and B matrices
    - `slack::Integer`: The index of the slack bus
"""
function handle_islands(B::AbstractMatrix, DA::AbstractMatrix, contingencies::AbstractVector{<:Tuple{Integer,Integer}}, branch::Integer, slack::Integer)
    islands = island_detection(B, contingencies)
    island = find_ref_island(islands, slack)
    island_b = find_island_branches(islands[island], DA, branch)

    # Need at least one node to make a valid system
    if length(islands[island]) > 0
        return islands, island, island_b
    end
    return Vector{Vector{Int64}}(undef, 0), 0, Vector{Int64}(undef, 0)
end

"""
    Handle islands after a contingency.
    Returns all islands and the index of the island with the slack bus.

    Arguments:
    - `B::AbstractMatrix`: The susceptance matrix
    - `cont::Tuple{Integer,Integer}`: The contingency as a tuple of (from, to) bus indices
    - `slack::Integer`: The index of the slack bus
"""
function handle_islands(B::AbstractMatrix, contingency::Tuple{Integer,Integer}, slack::Integer)
    islands = island_detection(B, contingency[1], contingency[2])
    island = find_ref_island(islands, slack)
    return length(islands[island]) > 0 ? (islands, islands) : (Vector{Vector{Int64}}(undef, 0), 0)
    # Need at least one node to make a valid system
end

""" 
    Get all islands 

    Arguments:
    - `B::AbstractMatrix`: The susceptance matrix
    - `bx::AbstractVector{<:Tuple{Integer,Integer}}`: The contingencies as a vector of tuples of (from, to) bus indices
    - `slack::Integer`: The index of the slack bus
"""
function get_all_islands(B::AbstractMatrix, bx::Vector{<:Tuple{Integer,Integer}}, slack::Integer)
    res = Vector{Vector{Int64}}(undef, length(bx))
    Threads.@threads for i in eachindex(bx)
        islands. island = handle_islands(B, bx[i], slack)
        res[i] = islands[island]
    end
    return res
end

""" 
    Get all islands 

    Arguments:
    - `opf::OPFsystem`: The OPF system
    - `slack::Integer`: The index of the slack bus
"""
function get_all_islands(opf::OPFsystem, slack::Integer)
    idx = get_nodes_idx(opf.nodes)
    bx = get_bus_idx.(opf.branches, [idx])
    cx = get_bus_idx.(opf.contingencies, [idx])
    A = calc_A(bx, length(opf.nodes))
    return get_all_islands(A' * A, cx, slack)
end

"""
    Find the islands in the system

    Arguments:
    - `T::SparseArrays.SparseMatrixCSC{<:Any,<:Integer}`: The admittance or connectivity matrix of the system
"""
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

"""
    Find all nodes in the same island as the given node

    Arguments:
    - `T::SparseArrays.SparseMatrixCSC{<:Any,<:Integer}`: The admittance or connectivity matrix of the system
    - `node::Integer`: The bus index to find the island for
"""
function new_islands(T::SparseArrays.SparseMatrixCSC{<:Any,<:Integer}, node::Integer)
    island = T[:, node].nzind
    for i in island
        union!(island, T[:, i].nzind)
        # print(i, " ")
    end
    return island
end


# The code under this line is based on SINT_LF by Sigurd Jakobsen (2023)
######################################################

"""
    create_connectivity_matrix(Y)
    Creates the connectivity matrix exaplined in [goderya_fast_1980](@cite).
    The code for creating the matrix is inspired by 
    https://github.com/NREL-SIIP/PowerSystems.jl

    Arguments:
    - `Y::SparseArrays.SparseMatrixCSC`: The admittance matrix of the system
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
