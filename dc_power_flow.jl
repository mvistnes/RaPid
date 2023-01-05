# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
using PowerSystems

" Return the overload of a line, else return 0.0 "
find_overload(flow::T, rate::Real) where {T<:Real} = abs(flow)-rate > 0.0 ? sign(flow)*(abs(flow)-rate) : zero(T)

" DC line flow calculation using Injection Shift Factors and Power Injection vector"
calculate_line_flows(isf::AbstractMatrix{<:Real}, Pᵢ::AbstractVector{<:Real}) = isf*Pᵢ

" Calculate the power flow on the lines from the voltage angles "
function calc_Pline(branches::Vector{<:Branch}, δ::Vector{<:Real}, idx::Dict{<:Any, <:Int})
    P = zeros(length(branches))
    for (i,branch) in enumerate(branches)
        (f, t) = get_bus_idx(branch, idx)
        P[i] = (δ[f] - δ[t]) / branch.x
    end
    return P
end

calc_Pline2(branches::Vector{<:Branch}, δ::Vector{<:Real}, idx::Dict{<:Any, <:Int}) = 
        calc_pl.([δ], get_bus_idx.(branches, [idx]), get_x.(branches))

calc_Pline(δ::Vector{<:Real}, nodes::Tuple, x::Real) = (δ[nodes[1]] - δ[nodes[2]]) / x

""" 
Calculate the power flow on the lines from the connectivity 
and the diagonal admittance matrices and the voltage angles 
"""
calc_Pline(A::AbstractMatrix{<:Real}, D::AbstractMatrix{<:Real}, δ::Vector{<:Real}) = D * A * δ
calc_Pline(DA::AbstractMatrix{<:Real}, δ::Vector{<:Real}) = DA * δ

""" 
Calculate the power flow on the lines from the connectivity 
and the diagonal admittance matrices and the voltage angles
in a contingency of the branch number.
"""
function calc_Pline_contingency(A::AbstractMatrix{<:Real}, D::AbstractMatrix{<:Real}, δ::Vector{<:Real}, branch::Integer)
    P = calc_Pline(A, D, δ)
    P[branch] = 0.0
    return P
end
function calc_Pline_contingency(DA::AbstractMatrix{<:Real}, δ::Vector{<:Real}, branch::Integer)
    P = calc_Pline(DA, δ)
    P[branch] = 0.0
    return P
end

# Slower than the above by 5x
# calc_Pline_contingency(A::AbstractMatrix{<:Real}, D::AbstractMatrix{<:Real}, δ::Vector{<:Real}, branch::Integer) = 
#         D[:, 1:end .!= branch] * A[1:end .!= branch, :] * δ

" DC power flow calculation using the Admittance matrix and Power Injection vector returning the bus angles "
run_pf(B::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}) = B \ P

""" Return the net power injected at each node. """
function get_net_Pᵢ(opfm::OPFmodel, nodes::Vector{Bus}, idx::Dict{<:Any, <:Int} = get_nodes_idx(nodes), Pᵢ = get_Pᵢ(opfm, nodes))
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
function get_Pᵢ(opfm::OPFmodel, nodes::Vector{Bus}, idx::Dict{<:Any, <:Int} = get_nodes_idx(nodes))
    Pᵢ = zeros(length(nodes))
    p = JuMP.value.(opfm.mod[:pg0])
    for g in get_ctrl_generation(opfm.sys)
        Pᵢ[idx[g.bus.number]] += p[get_name(g)]
    end
    return Pᵢ
end

" Calculate the net power for each node from the voltage angles "
function calc_Pᵢ(branches::Vector{<:Branch}, δ::Vector{<:Real}, numnodes::Int64, idx::Dict{<:Any, <:Int}, outage::Tuple = (0,0))
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


" Make the PTDF matrix for using the input nodes and branches "
function get_ptdf(branches::Vector{<:Branch}, numnodes::Integer, idx::Dict{<:Any, <:Integer}, slack::Integer)
    B = build_B(branches, numnodes, idx)
    B[:,slack] .= zero(Float64)
    B[slack,:] .= zero(Float64)
    B[slack,slack] = one(Float64)
    return get_ptdf(LinearAlgebra.lu(B), branches, numnodes, idx, slack)
end
function get_ptdf(Bx, branches::Vector{<:Branch}, numnodes::Integer, idx::Dict{<:Any, <:Integer}, slack::Integer)
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


"""
Build the adjecency Matrix

Input:
    - branches: Tuples of the from- and to-node index for each branch
    - numnodes: The number of nodes in the system
"""
function build_adjacency(branches::Vector{<:Tuple{Integer, Integer}}, numnodes::Integer)
    adj = SparseArrays.spzeros(Int8, numnodes, numnodes)
    for (f, t) in branches
        adj[f, t] = 1
        adj[t, f] = -1
        adj[f, f] = 1
        adj[t, t] = 1
    end
    return adj
end

build_adjacency(branches::Vector{<:Branch}, numnodes::Integer, idx::Dict{<:Any, <:Integer}) =
    build_adjacency(get_bus_idx.(branches, [idx]), numnodes)


" Make the branch connectivity matrix "
function calc_A(branches::Vector{<:Tuple{Integer, Integer}}, numnodes::Integer)
    A = SparseArrays.spzeros(Int8, length(branches), numnodes)
    for (i, (f, t)) in enumerate(branches)
        A[i,f] = one(Float64)
        A[i,t] = -one(Float64)
    end
    return A
end

calc_A(branches::Vector{<:Branch}, numnodes::Integer, idx::Dict{<:Any, <:Integer}) =
    calc_A(get_bus_idx.(branches, [idx]), numnodes)

" Calculate the inverse of the adjecency matrix "
function calc_X(B::AbstractMatrix{T}, slack::Integer) where T<:Real
    B[:,slack] .= zero(T)
    B[slack,:] .= zero(T)
    B[slack,slack] = one(T)
    X = inv(Matrix(B))
    X[slack,slack] = zero(T)
    return X
end

calc_isf(A::AbstractMatrix{<:Real}, D::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real}) = D * A * X

calc_B(A::AbstractMatrix{<:Real}, D::AbstractMatrix{<:Real}) = A' * D * A

" Builds an admittance matrix with the line series reactance of the lines. "
function build_B(branches::Vector{<:Tuple{Integer, Integer}}, X::Vector{<:Real}, numnodes::Integer)
    B = SparseArrays.spzeros(numnodes, numnodes)
    for ((f, t), x) in zip(branches, X)
        B[f,f] += 1 / x
        B[t,t] += 1 / x
        B[f,t] -= 1 / x
        B[t,f] -= 1 / x
    end
    return B
end

build_B(branches::Vector{<:Branch}, numnodes::Integer, idx::Dict{<:Any, <:Integer}) =
    build_B(get_bus_idx.(branches, [idx]), get_x.(branches), numnodes)

calc_D(x::Vector{<:Real}) = LinearAlgebra.Diagonal(1 ./ x)
calc_D(branches::Vector{<:Branch}) = calc_D(get_x.(branches))

function get_isf(branches::Vector{<:Branch}, numnodes::Integer, idx::Dict{<:Any, <:Integer}, slack::Integer)
    A = calc_A(branches, numnodes, idx)
    D = calc_D(branches)
    return calc_isf(A, D, calc_X(calc_B(A, D), slack))
end

" Get the isf-matrix after a line outage "
function get_isf(A, D, B, from_bus_idx::Integer, to_bus_idx::Integer, i_branch::Integer, x::Real, slack::Integer)
    neutralize_line(B, t, f, 1 / x)
    #d = LinearAlgebra.diag(D)
    isf = calc_isf(A, D, calc_X(B, slack))
    #         A[1:end .!= i_branch,:], 
    #         LinearAlgebra.Diagonal(d[1:end .!= i_branch]), 
    #         calc_X(B, slack)
    #     )
    isf[i_branch,:] .= 0
    neutralize_line(B, t, f, -1 / x)
    return isf
end
function get_isf(A, D, B, idx::Dict{<:Any, <:Integer}, i_branch::Integer, branch::Branch, slack::Integer)
    (f, t) = get_bus_idx(branch, idx)
    return get_isf(A, D, B, f, t, i_branch, get_x(branch), slack)
end

function neutralize_line(B::AbstractMatrix, i::Integer, j::Integer, val::Real)
    B[i,j] += val
    B[j,i] += val
    B[i,i] -= val
    B[j,j] -= val
end

"""
    Calculation of voltage angles in a contingency case using IMML

    Input:
        - H_inv_from_bus: A column from the inverse admittance matrix
        - H_inv_to_bus: A column from the inverse admittance matrix
        - H: A value from the admittance matrix
        - δ₀: Inital voltage angles
        - from_bus: From bus index
        - to_bus: To bus index
        - slack: Slack bus index
        - change: The amount of reactance change of the line, <=1. 
          Default is 1 and this removes the line
"""
function get_changed_angles(
            H_inv_from_bus::AbstractVector{<:Real}, 
            H_inv_to_bus::AbstractVector{<:Real}, 
            H::Real, 
            δ₀::AbstractVector{<:Real}, 
            from_bus::Integer, 
            to_bus::Integer, 
            slack::Integer,
            change::Real = 1.0
        )

    x = zeros(length(δ₀))
    # x[1:end .!= slack] = change .* (H_inv_from_bus .- H_inv_to_bus)
    x[1:end .!= slack] = change .* (H_inv_from_bus[1:end .!= slack] .- H_inv_to_bus[1:end .!= slack])
    c_inv = 1/H + x[from_bus] - x[to_bus]
    delta = 1/c_inv * (δ₀[from_bus] - δ₀[to_bus])
    return δ₀ .- x .* delta
end

get_changed_angles(
        H_inv::AbstractMatrix{<:Real}, 
        H::AbstractMatrix{<:Real}, 
        δ₀::AbstractVector{<:Real}, 
        from_bus::Integer, 
        to_bus::Integer, 
        slack::Integer,
        change::Real = 1.0
    ) = get_changed_angles(
        H_inv[:,from_bus],
        H_inv[:,to_bus],
        H[from_bus,to_bus],
        δ₀,
        from_bus,
        to_bus,
        slack,
        change
    )