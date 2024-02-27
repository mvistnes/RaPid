abstract type PowerFlow end

mutable struct DCPowerFlow{T1<:Real,T2<:Integer} <: PowerFlow
    DA::SparseArrays.SparseMatrixCSC{T1,T2} # Diagonal admittance matrix times the connectivity matrix
    B::SparseArrays.SparseMatrixCSC{T1,T2} # The admittance matrix
    K::KLU.KLUFactorization{T1, T2} # A factorization of the admittance matrix
    X::Matrix{T1} # Inverse of the admittance matrix
    ϕ::Matrix{T1} # PTDF matrix
    θ::Vector{T1} # Bus voltage angles
    F::Vector{T1} # Line power flow
    slack::T2 # Reference bus

    sp_tmp::SparseArrays.SparseMatrixCSC{T1,T2} # branch x bus sparse matrix
    K_tmp::KLU.KLUFactorization{T1, T2} # A factorization of the admittance matrix
    mnn_tmp::Matrix{T1} # bus x bus matrix
    mbn_tmp::Matrix{T1} # branch x bus matrix
    vn_tmp::Vector{T1} # bus vector
    vb_tmp::Vector{T1} # branch vector
end

function DCPowerFlow(branches::AbstractVector{<:Tuple{T2,T2}}, susceptance::AbstractVector{T1}, numnodes::Integer, slack::Integer
) where {T1<:Real,T2<:Integer}
    A = calc_A(branches, numnodes)
    DA = calc_DA(A, susceptance)
    B = calc_B(A, DA)
    K = get_klu(B, slack)
    ϕ = get_isf(K, DA, slack)
    set_tol_zero!(ϕ)
    X = calc_X(K, slack)
    return DCPowerFlow{T1,T2}(DA, B, K, X, ϕ, zeros(T1, numnodes), zeros(T1, length(branches)), slack,
        zeros(T1, size(DA)), get_klu(B, slack), zeros(T1, size(X)), zeros(T1, size(ϕ)), zeros(T1, numnodes), zeros(T1, length(branches)))
end
DCPowerFlow(nodes::AbstractVector{<:Bus}, branches::AbstractVector{<:Branch}, idx::Dict{<:Int,<:Int}) =
    DCPowerFlow(get_bus_idx.(branches, [idx]), PowerSystems.get_series_susceptance.(branches), length(nodes), find_slack(nodes)[1])

function DCPowerFlow(sys::System)
    nodes = sort_components!(get_nodes(sys))
    branches = sort_components!(get_branches(sys))
    idx = get_nodes_idx(nodes)
    return DCPowerFlow(nodes, branches, idx)
end

function DCPowerFlow(sys::System, Pᵢ::AbstractVector{<:Real})
    pf = DCPowerFlow(sys)
    calc_θ!(pf, Pᵢ)
    calc_Pline!(pf)
    return pf
end

function DCPowerFlow(nodes::AbstractVector{<:Bus}, branches::AbstractVector{<:Branch}, idx::Dict{<:Int,<:Int}, Pᵢ::AbstractVector{<:Real})
    pf = DCPowerFlow(nodes, branches, idx)
    calc_θ!(pf, Pᵢ)
    calc_Pline!(pf)
    return pf
end

get_slack(pf::DCPowerFlow) = pf.slack
get_DA(pf::DCPowerFlow) = pf.DA
get_B(pf::DCPowerFlow) = pf.B
get_K(pf::DCPowerFlow) = pf.K
get_X(pf::DCPowerFlow) = pf.X
get_ϕ(pf::DCPowerFlow) = pf.ϕ
get_θ(pf::DCPowerFlow) = pf.θ

set_θ!(pf::DCPowerFlow, model::Model) = copy!(pf.θ, get_sorted_angles(model))

""" Calculate the net power for each node from the voltage angles """
function calc_Pᵢ(branches::AbstractVector{<:Branch}, θ::AbstractVector{<:Real}, numnodes::Integer,
    idx::Dict{<:Int,<:Int}, outage::Tuple{<:Integer,<:Integer}=(0, 0)
)
    P = zeros(numnodes)
    for branch in branches
        (f, t) = get_bus_idx(branch, idx)
        val = (θ[f] - θ[t]) / branch.x
        if outage != (f, t)
            P[f] += val
            P[t] -= val
        end
    end
    return P
end

""" Calculate the net power for each node from the power flow """
function calc_Pᵢ_from_flow(branches::AbstractVector{<:Branch}, F::AbstractVector{<:Real}, numnodes::Integer, idx::Dict{<:Int,<:Int})
    P = zeros(numnodes)
    for (i, branch) in enumerate(branches)
        (f, t) = get_bus_idx(branch, idx)
        P[f] += F[i]
        P[t] -= F[i]
    end
    return P
end

""" Distributed slack """
function set_dist_slack!(ϕ::AbstractMatrix{<:Real}, mgx::AbstractMatrix{<:Real}, dist_slack::AbstractVector{<:Real})
    @assert !iszero(sum(dist_slack))
    slack_array = dist_slack / sum(dist_slack)
    ϕ = ϕ .- ((slack_array' * mgx) * ϕ')'
end

""" Make the component connectivity matrix """
function calc_connectivity(vals::AbstractVector{<:StaticInjection}, numnodes::Integer, idx::Dict{<:Int,<:Integer})
    num = length(vals)
    return SparseArrays.sparse(1:num, get_bus_idx.(vals, [idx]), one(Int8), num, numnodes)
end

""" Make the branch connectivity matrix """
function calc_A(branches::AbstractVector{<:Tuple{Integer,Integer}}, numnodes::Integer)
    num_b = length(branches)

    A_I = Vector{Int}(undef, 2 * num_b)
    @. A_I[1:num_b] = 1:num_b
    @. A_I[(num_b+1):end] = 1:num_b

    A_J = Vector{Int}(undef, 2 * num_b)
    @. A_J[1:num_b] = first(branches)
    @. A_J[(num_b+1):end] = last(branches)
    
    A_V = Vector{Int8}(undef, 2 * num_b)
    @. A_V[1:num_b] = one(Int8)
    @. A_V[(num_b+1):end] = -one(Int8)

    return SparseArrays.sparse(A_I, A_J, A_V, num_b, numnodes)
end
calc_A(branches::AbstractVector{<:Branch}, numnodes::Integer, idx::Dict{<:Int,<:Integer}) =
    calc_A(get_bus_idx.(branches, [idx]), numnodes)

""" Make the diagonal suseptance matrix """
calc_D(susceptance::AbstractVector{<:Real}) = LinearAlgebra.Diagonal(susceptance)
calc_D(branches::AbstractVector{<:Branch}) = calc_D(PowerSystems.get_series_susceptance.(branches))

calc_DA(A::AbstractMatrix{<:Real}, susceptance::AbstractVector{<:Real}) = calc_D(susceptance) * A
calc_DA(branches::AbstractVector{<:Branch}, numnodes::Integer, idx::Dict{<:Int,<:Integer}) = 
    calc_DA(calc_A(branches, idx), numnodes, PowerSystems.get_series_susceptance.(branches))

calc_B(A::AbstractMatrix{<:Real}, DA::AbstractMatrix{<:Real}) = A' * DA
function calc_B(branches::AbstractVector{<:Tuple{Integer,Integer}}, numnodes::Integer, susceptance::AbstractVector{<:Real})
    A = calc_A(branches, numnodes)
    D = calc_D(susceptance)
    return A' * (D * A)
end
calc_B(branches::AbstractVector{<:Branch}, numnodes::Integer, idx::Dict{<:Int,<:Integer}) =
    calc_B(get_bus_idx.(branches, [idx]), numnodes, PowerSystems.get_series_susceptance.(branches))

function _add_B!(B::SparseArrays.SparseMatrixCSC, br::ACBranch, idx::Dict{<:Int,<:Integer})
    (fbus, tbus) = get_bus_idx(br, idx)
    y = 1 / get_x(br)
    b = get_b(br)
    B[fbus,fbus] += (y - b.from)
    B[tbus,tbus] += (y - b.to)
    B[fbus,tbus] -= y
    B[tbus,fbus] -= y
    return
end

function _add_B!(B::SparseArrays.SparseMatrixCSC, br::Transformer2W, idx::Dict{<:Int,<:Integer})
    (fbus, tbus) = get_bus_idx(br, idx)
    y = 1 / get_x(br)
    b = get_primary_shunt(br)
    B[fbus,fbus] += (y + b)
    B[tbus,tbus] += y
    B[fbus,tbus] -= y
    B[tbus,fbus] -= y
    return
end

function _add_B!(B::SparseArrays.SparseMatrixCSC, br::TapTransformer, idx::Dict{<:Int,<:Integer})
    (fbus, tbus) = get_bus_idx(br, idx)
    y = 1 / get_x(br)
    tap = 1 / get_tap(br)
    b = get_primary_shunt(br)
    B[fbus,fbus] += ((y * tap^2) + b)
    B[tbus,tbus] += y
    B[fbus,tbus] -= (y * tap)
    B[tbus,fbus] -= (y * tap)
    return
end

function _add_B!(B::SparseArrays.SparseMatrixCSC, br::PhaseShiftingTransformer, idx::Dict{<:Int,<:Integer})
    (fbus, tbus) = get_bus_idx(br, idx)
    y = 1 / get_x(br)
    tap = get_tap(br) * exp(get_α(br))
    c_tap = get_tap(br) * exp(-1 * get_α(br))
    b = get_primary_shunt(br)
    B[fbus,fbus] += ((y * tap^2) + b)
    B[tbus,tbus] += y
    B[fbus,tbus] -= (y / c_tap)
    B[tbus,fbus] -= (y / tap)
    return
end

function _add_B!(B::SparseArrays.SparseMatrixCSC, br::FixedAdmittance, idx::Dict{<:Int,<:Integer})
    bus = idx[br.bus.number]
    B[bus,bus] -= imag(get_Y(br))
    return
end

# Slower than previous function
# """ Builds an admittance matrix with the line series reactance of the lines. """
function calc_B(branches::AbstractVector{<:Branch}, idx::Dict{<:Int,<:Integer}, fixed::AbstractVector{FixedAdmittance})
    numnodes = length(idx)
    B = SparseArrays.spzeros(numnodes, numnodes)
    for branch in branches
        _add_B!(B, branch, idx)
    end
    for fx in fixed
        _add_B!(B, fx, idx)
    end
    return B
end

function get_klu!(A::SparseArrays.SparseMatrixCSC{T1, T2}, slack::Integer) where {T1<:Real,T2<:Integer}
    A[:, slack] .= zero(T1)
    A[slack, :] .= zero(T1)
    A[slack, slack] = one(T1)
    return KLU.klu(A)
end
get_klu(A::SparseArrays.SparseMatrixCSC, slack::Integer) =
    get_klu!(copy(A), slack)

""" Calculate the inverse of the admittance matrix.
    X must be filled with the values of B """
function calc_X!(X::AbstractMatrix{T}, K::KLU.KLUFactorization{T,<:Integer}, slack::Integer
) where {T<:Real}
    fill!(X, zero(T))
    for i in axes(X, 1)
        X[i,i] = one(T)
    end
    KLU.solve!(K, X)
    X[slack, slack] = zero(T)
    return X
end

""" Calculate the inverse of the admittance matrix """
function calc_X!(X::Matrix{<:Real}, B::AbstractMatrix{T}, slack::Integer) where {T<:Real}
    copy!(X, B)
    calc_X!(X, get_klu(B, slack), slack)
    return X
end
calc_X(B::AbstractMatrix{<:Real}, slack::Integer) = calc_X!(Matrix(B), B, slack)
calc_X(K::KLU.KLUFactorization{T,<:Integer}, slack::Integer) where {T<:Real} = calc_X!(Matrix{T}(undef, size(K)), K, slack)

function calc_X_vec!(x::Vector{T}, K::KLU.KLUFactorization{T,<:Integer}, bus::Integer, slack::Integer
) where {T<:Real}
    x .= zero(T)
    x[bus] = one(T)
    KLU.solve!(K, x)
    x[slack] = zero(T)
    return x
end
calc_X_vec(pf::DCPowerFlow, bus::Integer) = calc_X_vec!(pf.vn_tmp, pf.K, bus, pf.slack)
    
calc_isf(DA::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real}) = DA * X
calc_isf!(ϕ::AbstractMatrix{<:Real}, DA::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real}) =
    LinearAlgebra.mul!(ϕ, DA, X)

""" Make the isf-matrix """
function get_isf(K::KLU.KLUFactorization{T,<:Integer}, DA::AbstractMatrix{T}, slack::Integer
) where {T<:Real}
    ϕ = Matrix(DA')
    KLU.solve!(K, ϕ)
    ϕ[slack, :] .= zero(T)
    return ϕ'
end
function get_isf!(B::SparseArrays.SparseMatrixCSC{T,<:Integer}, DA::AbstractMatrix{T}, slack::Integer
) where {T<:Real}
    return get_isf(get_klu(B, slack), DA, slack)
end
function get_isf(branches::AbstractVector{<:Branch}, nodes::AbstractVector{<:Bus},
    idx::Dict{<:Int,<:Integer}=get_nodes_idx(nodes), slack::Integer=find_slack(nodes)[1])
    A = calc_A(branches, length(nodes), idx)
    DA = calc_DA(A, PowerSystems.get_series_susceptance.(branches))
    B = calc_B(A, DA)
    return get_isf!(B, DA, slack)
end

function calc_isf_vec!(ϕ_vec::Vector{T}, K::KLU.KLUFactorization{T,<:Integer}, DA::AbstractMatrix{T}, branch::Integer
) where {T<:Real}
    copyto!(ϕ_vec, Vector(DA[branch,:]))
    KLU.solve!(K, ϕ_vec)
    return ϕ_vec
end
calc_isf_vec(pf::DCPowerFlow, branch::Integer) = calc_isf_vec!(pf.vn_tmp, pf.K, pf.DA, branch)

""" Find voltage angles from the factorization of the B-matrix and injected power. Change both θ and K """
function _calc_θ!(θ::AbstractVector{T}, K::KLU.KLUFactorization{T,<:Integer}, P::AbstractVector{T}, slack::Integer
) where {T<:Real}
    copyto!(θ, P)
    KLU.solve!(K, θ)
    θ[slack] = zero(T)
    return θ
end
""" Find voltage angles from the B-matrix and injected power. Change θ """
_calc_θ!(θ::AbstractVector{T}, B::SparseArrays.SparseMatrixCSC{T,<:Integer}, P::AbstractVector{T}, slack::Integer
) where {T<:Real} =
    _calc_θ!(θ, get_klu(B, slack), P, slack)
calc_θ!(B::AbstractMatrix{T}, P::AbstractVector{T}, slack::Integer) where {T<:Real} =
    _calc_θ!(Vector{T}(undef, length(P)), B, P, slack)

""" Calculate the power flow on the lines from the voltage angles """
function calc_Pline(branches::AbstractVector{<:Branch}, θ::AbstractVector{<:Real}, idx::Dict{<:Int,<:Int})
    P = zeros(length(branches))
    for (i, branch) in enumerate(branches)
        (f, t) = get_bus_idx(branch, idx)
        P[i] = (θ[f] - θ[t]) / branch.x
    end
    return P
end

""" 
Calculate the power flow on the lines from the connectivity 
and the diagonal admittance matrices and the voltage angles 
"""
calc_Pline(DA::AbstractMatrix{<:Real}, θ::AbstractVector{<:Real}) = DA * θ
calc_Pline!(pf::DCPowerFlow) = LinearAlgebra.mul!(pf.F, pf.DA, pf.θ)

""" DC line flow calculation using Injection Shift Factors and Power Injection vector"""
calculate_line_flows(isf::AbstractMatrix{<:Real}, Pᵢ::AbstractVector{<:Real}) = isf * Pᵢ
calculate_line_flows!(pf::DCPowerFlow, Pᵢ::AbstractVector{<:Real}) = LinearAlgebra.mul!(pf.F, pf.ϕ, Pᵢ)

""" DC power flow calculation using the Admittance matrix factorization and Power Injection vector returning the bus angles """
function run_pf!(θ::AbstractVector{<:Real}, K::KLU.KLUFactorization, P::AbstractVector{<:Real}, slack::Integer)
    copy!(θ, P)
    KLU.solve!(K, θ)
    θ[slack] = 0.0
    return θ
end
run_pf(K::KLU.KLUFactorization, P::AbstractVector{<:Real}, slack::Integer) = run_pf!(similar(P), K, P, slack)
function run_pf!(pf::DCPowerFlow, Pᵢ::AbstractVector{<:Real})
    copy!(pf.θ, Pᵢ)
    KLU.solve!(pf.K, pf.θ)
    pf.θ[pf.slack] = 0.0
    return pf.θ
end

""" DC power flow calculation using the inverse Admittance matrix and Power Injection vector returning the bus angles """
calc_θ(X::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}) = X * P
calc_θ!(pf::DCPowerFlow, Pᵢ::AbstractVector{<:Real}) = LinearAlgebra.mul!(pf.θ, pf.X, Pᵢ)

""" Active Power Injection vector calculation using the Admittance matrix and the bus angles """
calc_Pᵢ(B::AbstractMatrix{<:Real}, θ::AbstractVector{<:Real}) = B * θ
calc_Pᵢ!(Pᵢ::AbstractVector{<:Real}, pf::DCPowerFlow) = LinearAlgebra.mul!(Pᵢ, pf.B, pf.θ)
calc_Pᵢ(pf::DCPowerFlow) = calc_Pᵢ!(Vector{eltype(pf.θ)}(undef, size(pf.θ)), pf)
