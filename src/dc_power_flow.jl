abstract type PowerFlow end

mutable struct DCPowerFlow{TR<:Real,TI<:Integer} <: PowerFlow
    DA::SparseArrays.SparseMatrixCSC{TR,TI} # Diagonal admittance matrix times the connectivity matrix
    B::SparseArrays.SparseMatrixCSC{TR,TI} # The admittance matrix
    K::KLU.KLUFactorization{TR, TI} # A factorization of the admittance matrix
    # X::Matrix{TR} # Inverse of the admittance matrix
    # ϕ::Matrix{TR} # PTDF matrix
    θ::Vector{TR} # Bus voltage angles
    F::Vector{TR} # Line power flow
    dist_slack::Vector{TR}
    slack::TI # Reference bus

    # sp_tmp::SparseArrays.SparseMatrixCSC{TR,TI} # branch x bus sparse matrix
    # K_tmp::KLU.KLUFactorization{TR, TI} # A factorization of the admittance matrix
    # mnn_tmp::Matrix{TR} # bus x bus matrix
    # mbn_tmp::Matrix{TR} # branch x bus matrix
    vn_tmp::Vector{TR} # bus vector
    vb_tmp::Vector{TR} # branch vector
end

function DCPowerFlow(branches::AbstractVector{<:Tuple{TI,TI}}, susceptance::AbstractVector{TR}, numnodes::Integer, dist_slack::AbstractVector{TR}, slack::Integer
) where {TR<:Real,TI<:Integer}
    A = calc_A(branches, numnodes)
    DA = calc_DA(A, susceptance)
    B = calc_B(A, DA)
    K = get_klu(B, slack)
    # ϕ = get_isf(K, DA, slack)
    # set_tol_zero!(ϕ)
    # X = calc_X(K, slack)
    dist = dist_slack / sum(dist_slack)
    return DCPowerFlow{TR,TI}(DA, B, K, zeros(TR, numnodes), zeros(TR, length(branches)), dist, slack,
        zeros(TR, numnodes), zeros(TR, length(branches)))
end
DCPowerFlow(nodes::AbstractVector{<:Bus}, branches::AbstractVector{<:Branch}, idx::Dict{<:Int,<:Int}) =
    DCPowerFlow(get_bus_idx.(branches, [idx]), PowerSystems.get_series_susceptance.(branches), length(nodes), find_slack(nodes)[1])

function DCPowerFlow(sys::System)
    nodes = sort_components!(get_nodes(sys))
    branches = sort_components!(get_branches(sys))
    idx = get_nodes_idx(nodes)
    return DCPowerFlow(nodes, branches, idx)
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
    dist_slack /= sum(dist_slack)
    ϕ = ϕ .- ((dist_slack' * mgx) * ϕ')'
end
function set_dist_slack!(ϕ::AbstractMatrix{<:Real}, opf::OPFsystem, dist_slack::AbstractVector{<:Real} = Float64[])
    if isempty(dist_slack)
        dist_slack = getproperty.(get_active_power_limits.(opf.ctrl_generation), [:max])
    end
    mgx = calc_connectivity(opf.ctrl_generation, length(opf.nodes), opf.idx)
    set_dist_slack!(ϕ, mgx, dist_slack)
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

function get_klu!(A::SparseArrays.SparseMatrixCSC{TR, TI}, slack::Integer) where {TR<:Real,TI<:Integer}
    A[:, slack] .= zero(TR)
    A[slack, :] .= zero(TR)
    A[slack, slack] = one(TR)
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

function calc_X_vec!(X::Vector{T}, K::KLU.KLUFactorization{T,<:Integer}, branch::Integer, slack::Integer
) where {T<:Real}
    X .= zero(T)
    X[branch] = one(T)
    KLU.solve!(K, X)
    X[slack] = zero(T)
    return X
end
calc_X_vec(pf::DCPowerFlow, branch::Integer, slack::Integer) = calc_X_vec!(pf.vn_tmp, pf.K, branch, slack)
    
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
    copyto!(ϕ_vec, DA[branch,:])
    KLU.solve!(K, ϕ_vec)
    return ϕ_vec
end

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
calc_Pline!(F::AbstractMatrix{<:Real}, DA::AbstractMatrix{<:Real}, θ::AbstractVector{<:Real}) = LinearAlgebra.mul!(F, DA, θ)
calc_Pline(DA::AbstractMatrix{T}, θ::AbstractVector{T}) where {T<:Real} = calc_Pline!(Vector{T}(undef, size(DA,1)), DA, θ)
calc_Pline!(pf::DCPowerFlow) = LinearAlgebra.mul!(pf.F, pf.DA, pf.θ)

""" DC line flow calculation using Injection Shift Factors and Power Injection vector"""
calculate_line_flows(isf::AbstractMatrix{<:Real}, Pᵢ::AbstractVector{<:Real}) = isf * Pᵢ
function calculate_line_flows!(pf::DCPowerFlow, Pᵢ::AbstractVector{<:Real})
    run_pf!(pf, Pᵢ)
    return calc_Pline!(pf)
end

function calculate_line_flows!(F::AbstractMatrix{<:Real}, θ::AbstractVector{<:Real}, K::KLU.KLUFactorization, 
    DA::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, slack::Integer
)
    run_pf!(θ, K, P, slack)
    calc_Pline!(F, DA, θ)
    return F
end

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
