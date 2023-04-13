# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
using PowerSystems

abstract type PowerFlow end

mutable struct DCPowerFlow <: PowerFlow
    DA::AbstractMatrix{<:Real} # Diagonal admittance matrix times the connectivity matrix
    B::AbstractMatrix{<:Real} # The admittance matrix
    fact_B # A factorization of the admittance matrix
    X::AbstractMatrix{<:Real} # Inverse of the admittance matrix
    ϕ::AbstractMatrix{<:Real} # PTDF matrix
    Pᵢ::AbstractVector{<:Real} # Net power vector
    θ::AbstractVector{<:Real} # Bus voltage angles
    F::AbstractVector{<:Real} # Line power flow
    slack::Integer # Reference bus
end

function DCPowerFlow(nodes::AbstractVector{<:Bus}, branches::AbstractVector{<:Branch}, idx::Dict{<:Any, <:Int})
    slack = find_slack(nodes)[1]
    A = calc_A(branches, length(nodes), idx)
    D = calc_D(branches)
    DA = D * A
    B = fast_calc_B(A, DA)
    X = calc_X(B, slack)
    ϕ = calc_isf(DA, X)
    fact_B = LinearAlgebra.factorize(B)
    return DCPowerFlow(DA, B, fact_B, X, ϕ, Vector{typeof(first(B))}(), Vector{typeof(first(B))}(), Vector{typeof(first(B))}(), slack)
end

function DCPowerFlow(nodes::AbstractVector{<:Bus}, branches::AbstractVector{<:Branch}, idx::Dict{<:Any, <:Int}, Pᵢ::AbstractVector{<:Real})
    pf = DCPowerFlow(nodes, branches, idx)
    pf.Pᵢ = Pᵢ
    run_pf!(pf)
    calc_Pline!(pf)
    return pf
end

get_slack(pf::DCPowerFlow) = pf.slack
get_DA(pf::DCPowerFlow) = pf.DA
get_B(pf::DCPowerFlow) = pf.B
get_fact_B(pf::DCPowerFlow) = pf.fact_B
get_X(pf::DCPowerFlow) = pf.X
get_ϕ(pf::DCPowerFlow) = pf.ϕ
get_Pᵢ(pf::DCPowerFlow) = pf.Pᵢ
get_θ(pf::DCPowerFlow) = pf.θ

""" Calculate the net power for each node from the voltage angles """
function calc_Pᵢ(branches::AbstractVector{<:Branch}, θ::AbstractVector{<:Real}, numnodes::Int64, idx::Dict{<:Any, <:Int}, outage::Tuple = (0,0))
    P = zeros(numnodes)
    for branch in branches
        (f, t) = get_bus_idx(branch, idx)
        if outage != (f, t)
            P[f] += (θ[f] - θ[t]) / branch.x
            P[t] -= (θ[f] - θ[t]) / branch.x
        end
    end
    return P
end

""" Make the generator connectivity matrix """ # TODO: Fix
function connectivitymatrix(gens, numnodes::Integer, idx::Dict{<:Any, <:Integer})
    bus = getindex.([idx], get_number.(get_bus.(gens)))
    mx = SparseArrays.spzeros(String, length(gens), numnodes)
    for (i, (g, b)) in enumerate(zip(gens, bus))
        mx[i, b] = g.name
    end
    return mx
end

connectivitymatrix(system::System, numnodes::Integer, idx::Dict{<:Any, <:Integer}) =
    connectivitymatrix(get_components(Generator, system), numnodes, idx)

""" Make the branch connectivity matrix """
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

""" Make the diagonal suseptance matrix """
calc_D(x::AbstractVector{<:Real}) = LinearAlgebra.Diagonal(1 ./ x)
calc_D(branches::AbstractVector{<:Branch}) = calc_D(get_x.(branches))

calc_B(A::AbstractMatrix{<:Real}, D::AbstractMatrix{<:Real}) = A' * D * A
fast_calc_B(A::AbstractMatrix{<:Real}, DA::AbstractMatrix{<:Real}) = A' * DA

""" Builds an admittance matrix with the line series reactance of the lines. """
function calc_B(branches::AbstractVector{<:Tuple{Integer, Integer}}, x::AbstractVector{<:Real}, numnodes::Integer)
    B = SparseArrays.spzeros(numnodes, numnodes)
    for ((f, t), x) in zip(branches, x)
        B[f,f] += 1 / x
        B[t,t] += 1 / x
        B[f,t] -= 1 / x
        B[t,f] -= 1 / x
    end
    return B
end

calc_B(branches::AbstractVector{<:Branch}, numnodes::Integer, idx::Dict{<:Any, <:Integer}) =
    calc_B(get_bus_idx.(branches, [idx]), get_x.(branches), numnodes)

""" Calculate the inverse of the admittance matrix """
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

""" Make the isf-matrix """
function get_isf(branches::AbstractVector{<:Branch}, numnodes::Integer, idx::Dict{<:Any, <:Integer}, slack::Integer)
    A = calc_A(branches, numnodes, idx)
    D = calc_D(branches)
    return calc_isf(A, D, calc_X(calc_B(A, D), slack))
end

""" Make the isf-matrix after a line outage """
function get_isf(DA, B, from_bus_idx::Integer, to_bus_idx::Integer, i_branch::Integer, slack::Integer)
    x = B[from_bus_idx, to_bus_idx] * DA[i_branch, to_bus_idx] / B[from_bus_idx, to_bus_idx]
    neutralize_line!(B, from_bus_idx, to_bus_idx, -x)
    isf = calc_isf(DA, calc_X(B, slack))
    neutralize_line!(B, from_bus_idx, to_bus_idx, x)
    isf[i_branch,:] .= 0
    return isf
end
function get_isf(DA, B, idx::Dict{<:Any, <:Integer}, i_branch::Integer, branch::Branch, slack::Integer)
    (f, t) = get_bus_idx(branch, idx)
    return get_isf(DA, B, f, t, i_branch, slack)
end

function neutralize_line!(B::AbstractMatrix, i::Integer, j::Integer, val::Real)
    B[i,j] += val
    B[j,i] += val
    B[i,i] -= val
    B[j,j] -= val
end


""" Make the PTDF matrix for using the input nodes and branches """
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

""" Return the overload of a line, else return 0.0 """
find_overload(flow::T, rate::Real, atol::Real = 1e-6) where {T<:Real} = 
    abs(flow)-rate > atol ? sign(flow)*(abs(flow)-rate) : zero(T)

filter_overload(flow::AbstractVector{<:Real}, linerating::AbstractVector{<:Real}, atol::Real = 1e-6) = 
    [(i,ol) for (i,ol) in enumerate(find_overload.(flow, linerating)) if abs(ol) > atol]

filter_overload(Δflow::AbstractVector{<:Tuple}, linerating::AbstractVector{<:Real}, atol::Real = 1e-6) = 
[(i,find_overload(ol, linerating[i])) for (i,ol) in Δflow if abs(find_overload(ol, linerating[i])) > atol]

""" Calculate the power flow on the lines from the voltage angles """
function calc_Pline(branches::AbstractVector{<:Branch}, θ::AbstractVector{<:Real}, idx::Dict{<:Any, <:Int})
    P = zeros(length(branches))
    for (i,branch) in enumerate(branches)
        (f, t) = get_bus_idx(branch, idx)
        P[i] = (θ[f] - θ[t]) / branch.x
    end
    return P
end

""" 
Calculate the power flow on the lines from the connectivity 
and the diagonal admittance matrices and the voltage angles 
"""
calc_Pline(A::AbstractMatrix{<:Real}, D::AbstractMatrix{<:Real}, θ::AbstractVector{<:Real}) = D * A * θ
calc_Pline(DA::AbstractMatrix{<:Real}, θ::AbstractVector{<:Real}) = DA * θ
function calc_Pline!(pf::DCPowerFlow) 
    pf.F = pf.DA * pf.θ
end

""" DC line flow calculation using Injection Shift Factors and Power Injection vector"""
calculate_line_flows(isf::AbstractMatrix{<:Real}, Pᵢ::AbstractVector{<:Real}) = isf*Pᵢ
function calculate_line_flows!(pf::DCPowerFlow) 
    pf.F = pf.ϕ*pf.Pᵢ
end

""" DC power flow calculation using the Admittance matrix and Power Injection vector returning the bus angles """
run_pf(B::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}) = B \ P
function run_pf!(pf::DCPowerFlow) 
    pf.θ = pf.fact_B \ pf.Pᵢ
end

