# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
using PowerSystems

""" Return the net power injected at each node. """
function get_net_Pᵢ(opfm::OPFmodel, nodes::AbstractVector{Bus}, idx::Dict{<:Any, <:Int} = get_nodes_idx(nodes), Pᵢ = get_Pᵢ(opfm, nodes))
    p = JuMP.value.(opfm.mod[:ls0])
    for r in get_renewables(opfm.sys)
        Pᵢ[idx[r.bus.number]] += get_active_power(r) - p[get_name(r)]
    end
    for d in get_demands(opfm.sys)
        Pᵢ[idx[d.bus.number]] -= get_active_power(d) + p[get_name(d)]
    end
    # @assert abs(sum(Pᵢ)) < 0.001
    return Pᵢ
end

""" Return the power injected at each node. """
function get_Pᵢ(opfm::OPFmodel, nodes::AbstractVector{Bus}, idx::Dict{<:Any, <:Int} = get_nodes_idx(nodes))
    Pᵢ = zeros(length(nodes))
    p = JuMP.value.(opfm.mod[:pg0])
    for g in get_ctrl_generation(opfm.sys)
        Pᵢ[idx[g.bus.number]] += p[get_name(g)]
    end
    return Pᵢ
end

" Calculate the net power for each node from the voltage angles "
function calc_Pᵢ(branches::AbstractVector{<:Branch}, δ::AbstractVector{<:Real}, numnodes::Int64, idx::Dict{<:Any, <:Int}, outage::Tuple = (0,0))
    P = zeros(numnodes)
    for branch in branches
        (f, t) = get_bus_idx(branch, idx)
        if outage != (f, t)
            P[f] += (δ[f] - δ[t]) / branch.x
            P[t] -= (δ[f] - δ[t]) / branch.x
        end
    end
    return P
end

"""
Build the adjecency Matrix

Input:
    - branches: Tuples of the from- and to-node index for each branch
    - numnodes: The number of nodes in the system
"""
function build_adjacency(branches::AbstractVector{<:Tuple{Integer, Integer}}, numnodes::Integer)
    adj = SparseArrays.spzeros(Int8, numnodes, numnodes)
    for (f, t) in branches
        adj[f, t] = 1
        adj[t, f] = -1
        adj[f, f] = 1
        adj[t, t] = 1
    end
    return adj
end

build_adjacency(branches::AbstractVector{<:Branch}, numnodes::Integer, idx::Dict{<:Any, <:Integer}) =
    build_adjacency(get_bus_idx.(branches, [idx]), numnodes)

function connectivitymatrix(gens::AbstractVector{<:Integer}, numnodes::Integer) 
    mx = SparseArrays.spzeros(Int8, length(gens), numnodes)
    for (i, g) in enumerate(gens)
        mx[i, g] = 1
    end
    return mx
end

connectivitymatrix(system, numnodes::Integer, idx::Dict{<:Any, <:Integer}) =
    connectivitymatrix(getindex.([idx], get_number.(get_bus.(get_components(Generator, system)))), numnodes)

" Make the branch connectivity matrix "
function calc_A(branches::AbstractVector{<:Tuple{Integer, Integer}}, numnodes::Integer)
    A = SparseArrays.spzeros(Int8, length(branches), numnodes)
    for (i, (f, t)) in enumerate(branches)
        A[i,f] = one(Int8)
        A[i,t] = -one(Int8)
    end
    return A
end

calc_A(branches::AbstractVector{<:Branch}, numnodes::Integer, idx::Dict{<:Any, <:Integer}) =
    calc_A(get_bus_idx.(branches, [idx]), numnodes)

" Calculate the inverse of the adjecency matrix "
function calc_X(B::AbstractMatrix{T}, slack::Integer) where T<:Real
    X = Matrix(B)
    X[:,slack] .= zero(T)
    X[slack,:] .= zero(T)
    X[slack,slack] = one(T)
    X = inv(X)
    X[slack,slack] = zero(T)
    return X
end

calc_isf(A::AbstractMatrix{<:Real}, D::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real}) = D * A * X
calc_isf(DA::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real}) = DA * X

calc_B(A::AbstractMatrix{<:Real}, D::AbstractMatrix{<:Real}) = A' * D * A
fast_calc_B(A::AbstractMatrix{<:Real}, DA::AbstractMatrix{<:Real}) = A' * DA

" Builds an admittance matrix with the line series reactance of the lines. "
function calc_B(branches::AbstractVector{<:Tuple{Integer, Integer}}, X::AbstractVector{<:Real}, numnodes::Integer)
    B = SparseArrays.spzeros(numnodes, numnodes)
    for ((f, t), x) in zip(branches, X)
        B[f,f] += 1 / x
        B[t,t] += 1 / x
        B[f,t] -= 1 / x
        B[t,f] -= 1 / x
    end
    return B
end

calc_B(branches::AbstractVector{<:Branch}, numnodes::Integer, idx::Dict{<:Any, <:Integer}) =
    calc_B(get_bus_idx.(branches, [idx]), get_x.(branches), numnodes)

calc_D(x::AbstractVector{<:Real}) = LinearAlgebra.Diagonal(1 ./ x)
calc_D(branches::AbstractVector{<:Branch}) = calc_D(get_x.(branches))

function get_isf(branches::AbstractVector{<:Branch}, numnodes::Integer, idx::Dict{<:Any, <:Integer}, slack::Integer)
    A = calc_A(branches, numnodes, idx)
    D = calc_D(branches)
    return calc_isf(A, D, calc_X(calc_B(A, D), slack))
end

" Get the isf-matrix after a line outage "
function get_isf(A, D, B, from_bus_idx::Integer, to_bus_idx::Integer, i_branch::Integer, x::Real, slack::Integer)
    neutralize_line!(B, from_bus_idx, to_bus_idx, 1 / x)
    #d = LinearAlgebra.diag(D)
    isf = calc_isf(A, D, calc_X(B, slack))
    #         A[1:end .!= i_branch,:], 
    #         LinearAlgebra.Diagonal(d[1:end .!= i_branch]), 
    #         calc_X(B, slack)
    #     )
    isf[i_branch,:] .= 0
    neutralize_line!(B, from_bus_idx, to_bus_idx, -1 / x)
    return isf
end
function get_isf(A, D, B, idx::Dict{<:Any, <:Integer}, i_branch::Integer, branch::Branch, slack::Integer)
    (f, t) = get_bus_idx(branch, idx)
    return get_isf(A, D, B, f, t, i_branch, get_x(branch), slack)
end

function neutralize_line!(B::AbstractMatrix, i::Integer, j::Integer, val::Real)
    B[i,j] += val
    B[j,i] += val
    B[i,i] -= val
    B[j,j] -= val
end


" Make the PTDF matrix for using the input nodes and branches "
function get_ptdf(branches::AbstractVector{<:Branch}, numnodes::Integer, idx::Dict{<:Any, <:Integer}, slack::Integer)
    B = calc_B(branches, numnodes, idx)
    B[:,slack] .= zero(Float64)
    B[slack,:] .= zero(Float64)
    B[slack,slack] = one(Float64)
    return get_ptdf(LinearAlgebra.lu(B), branches, numnodes, idx, slack)
end
function get_ptdf(Bx, branches::AbstractVector{<:Branch}, numnodes::Integer, idx::Dict{<:Any, <:Integer}, slack::Integer)
    A = zeros(Float64, size(branches,1), numnodes) # Container for the distribution factors
    for (i, branch) in enumerate(branches)
        branch.x == 0 && continue
        ΔP = zeros(Float64, numnodes)
        (f, t) = get_bus_id(branch) # (f)rom and (t)o bus at this branch
        if f != slack
            ΔP[idx[f]] = 1 / branch.x
        end
        if t != slack # ∈ keys(idx)
            ΔP[idx[t]] = -1 / branch.x
        end
        A[i, :] = Bx \ ΔP # append factors to matrix
    end
    return A
end

" Return the overload of a line, else return 0.0 "
find_overload(flow::T, rate::Real) where {T<:Real} = abs(flow)-rate > 0.0 ? sign(flow)*(abs(flow)-rate) : zero(T)

" DC line flow calculation using Injection Shift Factors and Power Injection vector"
calculate_line_flows(isf::AbstractMatrix{<:Real}, Pᵢ::AbstractVector{<:Real}) = isf*Pᵢ

" Calculate the power flow on the lines from the voltage angles "
function calc_Pline(branches::AbstractVector{<:Branch}, δ::AbstractVector{<:Real}, idx::Dict{<:Any, <:Int})
    P = zeros(length(branches))
    for (i,branch) in enumerate(branches)
        (f, t) = get_bus_idx(branch, idx)
        P[i] = (δ[f] - δ[t]) / branch.x
    end
    return P
end

""" 
Calculate the power flow on the lines from the connectivity 
and the diagonal admittance matrices and the voltage angles 
"""
calc_Pline(A::AbstractMatrix{<:Real}, D::AbstractMatrix{<:Real}, δ::AbstractVector{<:Real}) = D * A * δ
calc_Pline(DA::AbstractMatrix{<:Real}, δ::AbstractVector{<:Real}) = DA * δ

""" 
Calculate the power flow on the lines from the connectivity 
and the diagonal admittance matrices and the voltage angles
in a contingency of the branch number.
"""
function calc_Pline(A::AbstractMatrix{<:Real}, D::AbstractMatrix{<:Real}, δ::AbstractVector{<:Real}, branch::Integer)
    P = calc_Pline(A, D, δ)
    P[branch] = 0.0
    return P
end
function calc_Pline(DA::AbstractMatrix{<:Real}, δ::AbstractVector{<:Real}, branch::Integer)
    P = calc_Pline(DA, δ)
    P[branch] = 0.0
    return P
end

" DC power flow calculation using the Admittance matrix and Power Injection vector returning the bus angles "
run_pf(B::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}) = B \ P

