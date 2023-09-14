# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

abstract type PowerFlow end

mutable struct DCPowerFlow{T1<:Real,T2<:Integer} <: PowerFlow
    DA::SparseArrays.SparseMatrixCSC{T1,T2} # Diagonal admittance matrix times the connectivity matrix
    B::SparseArrays.SparseMatrixCSC{T1,T2} # The admittance matrix
    fact_B::LinearAlgebra.Factorization{T1} # A factorization of the admittance matrix
    X::Matrix{T1} # Inverse of the admittance matrix
    ϕ::Matrix{T1} # PTDF matrix
    θ::Vector{T1} # Bus voltage angles
    F::Vector{T1} # Line power flow
    slack::T2 # Reference bus
end

function DCPowerFlow(branches::AbstractVector{<:Tuple{T2,T2}}, x::AbstractVector{T1}, numnodes::Integer, slack::Integer) where {T1<:Real,T2<:Integer}
    A = calc_A(branches)
    DA = calc_DA(A, x)
    B = calc_B(A, DA)
    X = calc_X(B, slack)
    ϕ = Matrix{T1}(undef, size(A))
    calc_isf!(ϕ, DA, X)
    fact_B = LinearAlgebra.factorize(Matrix(B))
    return DCPowerFlow{T1,T2}(DA, B, fact_B, X, ϕ, zeros(T1, numnodes), zeros(T1, length(branches)), slack)
end
DCPowerFlow(nodes::AbstractVector{<:Bus}, branches::AbstractVector{<:Branch}, idx::Dict{<:Int,<:Int}) =
    DCPowerFlow(get_bus_idx.(branches, [idx]), get_x.(branches), length(nodes), find_slack(nodes)[1])

function DCPowerFlow(sys::System)
    nodes = sort_components!(get_nodes(sys))
    branches = sort_components!(get_branches(sys))
    idx = get_nodes_idx(nodes)
    return DCPowerFlow(nodes, branches, idx)
end

function DCPowerFlow(model::Model, nodes::AbstractVector{<:Bus}, branches::AbstractVector{<:Branch}, idx::Dict{<:Int,<:Int}=get_nodes_idx(nodes))
    pf = DCPowerFlow(nodes, branches, idx)
    set_θ!(pf, model)
    calc_Pline!(pf)
    return pf
end

function DCPowerFlow(nodes::AbstractVector{<:Bus}, branches::AbstractVector{<:Branch}, Pᵢ::AbstractVector{<:Real}, idx::Dict{<:Int,<:Int}=get_nodes_idx(nodes))
    pf = DCPowerFlow(nodes, branches, idx)
    calc_θ!(pf, Pᵢ)
    calc_Pline!(pf)
    return pf
end
get_slack(pf::DCPowerFlow) = pf.slack
get_DA(pf::DCPowerFlow) = pf.DA
get_B(pf::DCPowerFlow) = pf.B
get_fact_B(pf::DCPowerFlow) = pf.fact_B
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

""" Make the generator connectivity matrix """ # TODO: Fix
function connectivitymatrix(gens, numnodes::Integer, idx::Dict{<:Int,<:Integer})
    bus = getindex.([idx], get_number.(get_bus.(gens)))
    mx = SparseArrays.spzeros(String, length(gens), numnodes)
    for (i, (g, b)) in enumerate(zip(gens, bus))
        mx[i, b] = g.name
    end
    return mx
end

connectivitymatrix(system::System, numnodes::Integer, idx::Dict{<:Int,<:Integer}) =
    connectivitymatrix(get_components(Generator, system), numnodes, idx)

""" Make the branch connectivity matrix """
# function calc_A(branches::AbstractVector{<:Tuple{Integer, Integer}}, numnodes::Integer)
#     A = SparseArrays.spzeros(Int8, length(branches), numnodes)
#     for (i, (f, t)) in enumerate(branches)
#         A[i,f] = one(Int8)
#         A[i,t] = -one(Int8)
#     end
#     return A
# end
# calc_A(branches::AbstractVector{<:Branch}, numnodes::Integer, idx::Dict{<:Int, <:Integer}) =
#     calc_A(get_bus_idx.(branches, [idx]), numnodes)

function calc_A(branches::AbstractVector{<:Tuple{Integer,Integer}})
    num_b = length(branches)

    A_I = Vector{Int}(undef, 2 * num_b)
    A_J = Vector{Int}(undef, 2 * num_b)
    A_V = Vector{Int8}(undef, 2 * num_b)

    # build incidence matrix A (lines x buses)
    for (ix, b) in pairs(branches)
        (fbus, tbus) = b

        A_I[ix] = ix
        A_J[ix] = fbus
        A_V[ix] = one(Int8)

        A_I[num_b + ix] = ix
        A_J[num_b + ix] = tbus
        A_V[num_b + ix] = -one(Int8)
    end

    return SparseArrays.sparse(A_I, A_J, A_V)
end
calc_A(branches::AbstractVector{<:Branch}, idx::Dict{<:Int,<:Integer}) =
    calc_A(get_bus_idx.(branches, [idx]))

""" Make the diagonal suseptance matrix """
calc_D(x::AbstractVector{<:Real}) = LinearAlgebra.Diagonal(1.0 ./ x)
calc_D(branches::AbstractVector{<:Branch}) = calc_D(get_x.(branches))

calc_DA(A::AbstractMatrix{<:Real}, x::AbstractVector{<:Real}) = calc_D(x) * A
calc_DA(branches::AbstractVector{<:Branch}, idx::Dict{<:Int,<:Integer}) = calc_DA(calc_A(branches, idx), get_x.(branches))

calc_B(A::AbstractMatrix{<:Real}, DA::AbstractMatrix{<:Real}) = A' * DA
function calc_B(branches::AbstractVector{<:Tuple{Integer,Integer}}, x::AbstractVector{<:Real})
    A = calc_A(branches)
    D = calc_D(x)
        # TODO: add support for TapTransformer and PhaseShiftingTransformer
    return A' * (D * A)
end
calc_B(branches::AbstractVector{<:Branch}, idx::Dict{<:Int,<:Integer}) =
    calc_B(get_bus_idx.(branches, [idx]), get_x.(branches))

# Slower than previous function
# """ Builds an admittance matrix with the line series reactance of the lines. """
# function calc_B(branches::AbstractVector{<:Tuple{Integer, Integer}}, x::AbstractVector{<:Real}, numnodes::Integer)
#     B = SparseArrays.spzeros(numnodes, numnodes)
#     for ((f, t), x) in zip(branches, x)
#         B[f,f] += 1 / x
#         B[t,t] += 1 / x
#         B[f,t] -= 1 / x
#         B[t,f] -= 1 / x
#     end
#     return B
# end

""" Calculate the inverse of the admittance matrix.
    X must be filled with the values of B """
function _calc_X!(X::Matrix{T}, slack::Integer) where {T<:Real}
    X[:, slack] .= zero(T)
    X[slack, :] .= zero(T)
    X[slack, slack] = one(T)
    LinearAlgebra.inv!(LinearAlgebra.lu!(X))
    # LinearAlgebra.inv!(LinearAlgebra.cholesky!(X))
    # X should always be positive definite and thus Cholesky can be used (OR IS THIS TRUE??)
    X[slack, slack] = zero(T)
    return X
end

""" Calculate the inverse of the admittance matrix """
function calc_X!(X::Matrix{<:Real}, B::AbstractMatrix{T}, slack::Integer) where {T<:Real}
    copy!(X, B)
    _calc_X!(X, slack)
    return X
end
calc_X(B::AbstractMatrix{<:Real}, slack::Integer) = calc_X!(Matrix(B), B, slack)

""" 
    Calculate the inverse of the admittance matrix after a line outage. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function calc_X!(X::Matrix{<:Real}, DA::AbstractMatrix{<:Real}, B::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer}, cont_branch::Integer, slack::Integer)
    copy!(X, B)
    neutralize_line!(X, cont[1], cont[2], DA[cont_branch, cont[1]])
    _calc_X!(X, slack)
    return X
end
calc_X(DA::AbstractMatrix{<:Real}, B::AbstractMatrix{T}, cont::Tuple{Integer,Integer}, cont_branch::Integer, slack::Integer) where {T<:Real} =
    calc_X!(Matrix{T}(undef, size(B)), DA, B, cont, cont_branch, slack)

""" 
    Calculate the inverse of the admittance matrix after a line outage which splits the system. 
    cont[1] (from_bus), cont[2] (to_bus), cont_branch branch number, and island is sorted index numbers 
"""
function calc_X(DA::AbstractMatrix{<:Real}, B::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer, island::AbstractVector{<:Integer})
    (fbus, tbus) = cont
    X = Matrix(B[island, island])
    c = ifelse(insorted(fbus, island), fbus, tbus)
    i = searchsortedfirst(island, c)
    X[i, i] += DA[cont_branch, tbus]
    _calc_X!(X, searchsortedfirst(island, slack))
    return X
end

function neutralize_line!(B::AbstractMatrix, i::Integer, j::Integer, val::Real)
    B[i, j] += val
    B[j, i] += val
    B[i, i] -= val
    B[j, j] -= val
end

calc_isf(D::AbstractMatrix{<:Real}, A::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real}) = D * A * X
calc_isf(DA::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real}) = DA * X
calc_isf!(ϕ::AbstractMatrix{<:Real}, DA::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real}) =
    LinearAlgebra.mul!(ϕ, DA, X)

""" Make the isf-matrix """
function get_isf!(B::SparseArrays.SparseMatrixCSC{T,<:Integer}, DA::AbstractMatrix{T}, slack::Integer) where {T<:Real}
    B[:, slack] .= zero(T)
    B[slack, :] .= zero(T)
    B[slack, slack] = one(T)
    K = KLU.klu(B) # TODO: Option to use old KLU containers or symbolics
    ϕ = Matrix(DA')
    KLU.solve!(K, ϕ)
    ϕ[slack, :] .= zero(T)
    return ϕ'
end
function get_isf(branches::AbstractVector{<:Branch}, nodes::AbstractVector{<:Bus},
    idx::Dict{<:Int,<:Integer}=get_nodes_idx(nodes), slack::Integer=find_slack(nodes)[1])
    A = calc_A(branches, idx)
    DA = calc_DA(A, get_x.(branches))
    B = calc_B(A, DA)
    return get_isf!(B, DA, slack)
end

""" Find voltage angles from the B-matrix and injected power """
function _calc_θ!(θ::AbstractVector{T}, B::SparseArrays.SparseMatrixCSC{T,<:Integer}, P::AbstractVector{T}, slack::Integer) where {T<:Real}
    B[:, slack] .= zero(T)
    B[slack, :] .= zero(T)
    B[slack, slack] = one(T)
    K = KLU.klu(B) # TODO: Option to use old KLU containers or symbolics
    copyto!(θ, P)
    KLU.solve!(K, θ)
    θ[slack] = zero(T)
    return θ
end
calc_θ!(B::AbstractMatrix{T}, P::AbstractVector{T}, slack::Integer) where {T<:Real} =
    _calc_θ!(Vector{T}(undef, length(P)), B, P, slack)

""" 
    Make the isf-matrix after a line outage using base case D*A and B. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function get_isf!(ϕ::AbstractMatrix{<:Real}, DA::AbstractMatrix{<:Real}, B::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer
)
    X = calc_X(DA, B, cont, cont_branch, slack)
    calc_isf!(ϕ, DA, X)
    ϕ[cont_branch, :] .= 0
    return ϕ
end
get_isf(DA::AbstractMatrix{T}, B::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer) where {T<:Real} =
    get_isf!(Matrix{T}(undef, size(DA)), DA, B, cont, cont_branch, slack)

""" 
    Find voltage angles after a line outage using base case D*A and B. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function get_θ!(θ::AbstractVector{<:Real}, B::AbstractMatrix{<:Real}, DA::AbstractMatrix{<:Real}, 
    B0::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, 
    cont::Tuple{Integer,Integer}, cont_branch::Integer, slack::Integer
)
    copyto!(B, B0)
    neutralize_line!(B, cont[1], cont[2], DA[cont_branch, cont[1]])
    _calc_θ!(θ, B, P, slack)
    return θ
end
get_θ(DA::AbstractMatrix{T}, B0::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, 
    cont::Tuple{Integer,Integer}, cont_branch::Integer, slack::Integer
) where {T<:Real} = get_θ!(Vector{T}(undef, size(DA, 2)), similar(B0), DA, B0, P, cont, cont_branch, slack)

function calculate_line_flows!(F::AbstractVector{T}, θ::AbstractVector{<:Real}, B::AbstractMatrix{<:Real}, 
    DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, 
    cont::Tuple{Integer,Integer}, cont_branch::Integer, slack::Integer
) where {T<:Real}
    get_θ!(θ, B, DA, B0, P, cont, cont_branch, slack)
    LinearAlgebra.mul!(F, DA, θ)
    F[cont_branch] = zero(T)
    return F
end

""" 
    Make the isf-matrix after a line outage, which splits the system, using base case D*A and B. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function get_isf!(ϕ::AbstractMatrix{T}, DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer, nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}
) where {T<:Real}
    B = get_cont_B(DA, B0, cont, cont_branch, nodes)
    # ϕ[sorted_missing(branches, size(ϕ,1)), sorted_missing(nodes, size(ϕ,2))] .= zero(T)
    fill!(ϕ, zero(T))
    ϕ[branches, nodes] = get_isf!(B, view(DA, branches, nodes), searchsortedfirst(nodes, slack))
    return ϕ
end
get_isf(DA::AbstractMatrix{T}, B0::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer, nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}
) where {T<:Real} = get_isf!(Matrix{T}(undef, size(DA)), DA, B0, cont, cont_branch, slack, nodes, branches)

""" 
    Make the B-matrix after a line outage, which splits the system, using base case D*A and B. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function get_cont_B(DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, nodes::AbstractVector{<:Integer}
)
    (fbus, tbus) = cont
    B = copy(view(B0, nodes, nodes))
    c = ifelse(insorted(fbus, nodes), fbus, tbus)
    i = searchsortedfirst(nodes, c)
    B[i, i] += DA[cont_branch, tbus]
    return B
end

""" 
    Find voltage angles after a line outage, which splits the system, using base case B. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function get_θ(DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, 
    cont::Tuple{Integer,Integer}, cont_branch::Integer, slack::Integer, nodes::AbstractVector{<:Integer}
)
    B = get_cont_B(DA, B0, cont, cont_branch, nodes)
    θ = calc_θ!(B, view(P, nodes), searchsortedfirst(nodes, slack))
    return θ
end

function calculate_line_flows!(F::AbstractVector{T}, θ::AbstractVector{<:Real}, DA::AbstractMatrix{<:Real}, 
    B::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, cont::Tuple{Integer,Integer}, cont_branch::Integer, 
    slack::Integer, nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}
) where {T<:Real}
    θ[nodes] = get_θ(DA, B, P, cont, cont_branch, slack, nodes)
    LinearAlgebra.mul!(F, DA, θ)
    F[setdiff(1:size(DA,1), branches)] .= zero(T)
    return F
end

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
calc_Pline(A::AbstractMatrix{<:Real}, D::AbstractMatrix{<:Real}, θ::AbstractVector{<:Real}) = D * A * θ
calc_Pline(DA::AbstractMatrix{<:Real}, θ::AbstractVector{<:Real}) = DA * θ
calc_Pline!(pf::DCPowerFlow) = LinearAlgebra.mul!(pf.F, pf.DA, pf.θ)

""" DC line flow calculation using Injection Shift Factors and Power Injection vector"""
calculate_line_flows(isf::AbstractMatrix{<:Real}, Pᵢ::AbstractVector{<:Real}) = isf * Pᵢ
calculate_line_flows!(pf::DCPowerFlow, Pᵢ::AbstractVector{<:Real}) = LinearAlgebra.mul!(pf.F, pf.ϕ, Pᵢ)

""" DC power flow calculation using the Admittance matrix and Power Injection vector returning the bus angles """
run_pf(B::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}) = B \ P
run_pf!(pf::DCPowerFlow, Pᵢ::AbstractVector{<:Real}) = LinearAlgebra.ldiv!(pf.θ, pf.fact_B, Pᵢ)

""" DC power flow calculation using the inverse Admittance matrix and Power Injection vector returning the bus angles """
calc_θ(X::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}) = X * P
calc_θ!(pf::DCPowerFlow, Pᵢ::AbstractVector{<:Real}) = LinearAlgebra.mul!(pf.θ, pf.X, Pᵢ)

""" Active Power Injection vector calculation using the Admittance matrix and the bus angles """
calc_Pᵢ(B::AbstractMatrix{<:Real}, θ::AbstractVector{<:Real}) = B * θ
calc_Pᵢ!(Pᵢ::AbstractVector{<:Real}, pf::DCPowerFlow) = LinearAlgebra.mul!(Pᵢ, pf.B, pf.θ)
calc_Pᵢ(pf::DCPowerFlow) = calc_Pᵢ!(Vector{eltype(pf.θ)}(undef, size(pf.θ)), pf)

function calculate_island_line_flows(pf::DCPowerFlow, cont::Tuple{Integer,Integer}, cont_branch::Integer, Pᵢ::AbstractVector{<:Real})
    island, island_b = handle_islands(pf.B, pf.DA, cont, cont_branch, pf.slack)
    ptdf = get_isf(pf.DA, pf.B, cont, cont_branch, pf.slack, island, island_b)
    return ptdf * Pᵢ[island]
end

function get_contingency_ptdf(opf::OPFsystem, pf::DCPowerFlow)
    idx = get_nodes_idx(opf.nodes)
    contids = [(x, get_bus_idx(opf.branches[x], idx)) for x in indexin(opf.contingencies, opf.branches)]
    ptdf = Array{Float64}(undef, size(pf.ϕ, 1), size(pf.ϕ, 2), length(opf.contingencies))
    # for i in eachindex(contids)
    Threads.@threads for i in eachindex(contids)
        (cont_branch, cont) = contids[i]
        if !is_islanded(pf, cont, cont_branch)
            get_isf!(ptdf[:, :, i], pf.X, pf.B, pf.DA, cont, cont_branch)
        else
            islands = island_detection_thread_safe(pf.B, cont[1], cont[2])
            island = find_ref_island(islands, pf.slack)
            island_b = find_island_branches(islands[island], pf.DA, cont_branch)
            fill!(ptdf[:, :, i], zero(Float64))
            get_isf!(ptdf[island_b, islands[island], i], pf.DA, pf.B, cont, cont_branch, pf.slack)
        end
    end
    return ptdf
end