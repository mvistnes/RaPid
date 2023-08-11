# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

abstract type PowerFlow end

mutable struct DCPowerFlow{T1<:Real, T2<:Integer} <: PowerFlow 
    DA::SparseArrays.SparseMatrixCSC{T1, T2} # Diagonal admittance matrix times the connectivity matrix
    B::SparseArrays.SparseMatrixCSC{T1, T2} # The admittance matrix
    fact_B::LinearAlgebra.Factorization{T1} # A factorization of the admittance matrix
    X::Matrix{T1} # Inverse of the admittance matrix
    ϕ::Matrix{T1} # PTDF matrix
    θ::Vector{T1} # Bus voltage angles
    F::Vector{T1} # Line power flow
    slack::T2 # Reference bus
end

function DCPowerFlow(branches::AbstractVector{<:Tuple{T2, T2}}, x::AbstractVector{T1}, numnodes::Integer, slack::Integer) where {T1<:Real, T2<:Integer}
    A = calc_A(branches, numnodes)
    D = calc_D(x)
    DA = D * A
    B = calc_B(A, DA)
    X = calc_X(B, slack)
    ϕ = zeros(eltype(X), size(A))
    calc_isf!(ϕ, DA, X)
    fact_B = LinearAlgebra.factorize(Matrix(B))
    return DCPowerFlow{T1, T2}(DA, B, fact_B, X, ϕ, zeros(T1, numnodes), zeros(T1, length(branches)), slack)
end
DCPowerFlow(nodes::AbstractVector{<:Bus}, branches::AbstractVector{<:Branch}, idx::Dict{<:Int, <:Int} = get_nodes_idx(nodes)) =
    DCPowerFlow(get_bus_idx.(branches, [idx]), get_x.(branches), length(nodes), find_slack(nodes)[1])

function DCPowerFlow(sys::System)
    nodes = sort_components!(get_nodes(sys))
    branches = sort_components!(get_branches(sys))
    return DCPowerFlow(nodes, branches)
end

function DCPowerFlow(model::Model, nodes::AbstractVector{<:Bus}, branches::AbstractVector{<:Branch}, idx::Dict{<:Int, <:Int} = get_nodes_idx(nodes))
    pf = DCPowerFlow(nodes, branches, idx)
    set_θ!(pf, model)
    calc_Pline!(pf)
    return pf
end

function DCPowerFlow(nodes::AbstractVector{<:Bus}, branches::AbstractVector{<:Branch}, Pᵢ::AbstractVector{<:Real}, idx::Dict{<:Int, <:Int} = get_nodes_idx(nodes))
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
        idx::Dict{<:Int, <:Int}, outage::Tuple{<:Integer, <:Integer} = (0,0)
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
function calc_Pᵢ_from_flow(branches::AbstractVector{<:Branch}, F::AbstractVector{<:Real}, numnodes::Integer, idx::Dict{<:Int, <:Int})
    P = zeros(numnodes)
    for (i,branch) in enumerate(branches)
        (f, t) = get_bus_idx(branch, idx)
        P[f] += F[i]
        P[t] -= F[i]
    end
    return P
end

""" Make the generator connectivity matrix """ # TODO: Fix
function connectivitymatrix(gens, numnodes::Integer, idx::Dict{<:Int, <:Integer})
    bus = getindex.([idx], get_number.(get_bus.(gens)))
    mx = SparseArrays.spzeros(String, length(gens), numnodes)
    for (i, (g, b)) in enumerate(zip(gens, bus))
        mx[i, b] = g.name
    end
    return mx
end

connectivitymatrix(system::System, numnodes::Integer, idx::Dict{<:Int, <:Integer}) =
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

calc_A(branches::AbstractVector{<:Branch}, numnodes::Integer, idx::Dict{<:Int, <:Integer}) =
    calc_A(get_bus_idx.(branches, [idx]), numnodes)

""" Make the diagonal suseptance matrix """
calc_D(x::AbstractVector{<:Real}) = LinearAlgebra.Diagonal(1. ./ x)
calc_D(branches::AbstractVector{<:Branch}) = calc_D(get_x.(branches))

# calc_B(A::AbstractMatrix{<:Real}, D::AbstractMatrix{<:Real}) = A' * D * A
calc_B(A::AbstractMatrix{<:Real}, DA::AbstractMatrix{<:Real}) = A' * DA

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

calc_B(branches::AbstractVector{<:Branch}, numnodes::Integer, idx::Dict{<:Int, <:Integer}) =
    calc_B(get_bus_idx.(branches, [idx]), get_x.(branches), numnodes)

""" Calculate the inverse of the admittance matrix.
    X must be filled with the values of B """
function _calc_X!(X::Matrix{T}, slack::Integer) where T<:Real
    X[:,slack] .= zero(T)
    X[slack,:] .= zero(T)
    X[slack,slack] = one(T)
    LinearAlgebra.inv!(LinearAlgebra.cholesky!(X))
        # X should always be positive definite and thus Cholesky can be used
    X[slack,slack] = zero(T)
    return X
end

""" Calculate the inverse of the admittance matrix """
function calc_X!(X::Matrix{<:Real}, B::AbstractMatrix{T}, slack::Integer) where T<:Real 
    copy!(X, B)
    _calc_X!(X, slack)
    return X
end
calc_X(B::AbstractMatrix{<:Real}, slack::Integer) = calc_X!(Matrix(B), B, slack)

""" 
    Calculate the inverse of the admittance matrix after a line outage. 
    cont[1] (from_bus), cont[2] (to_bus), and i_branch are index numbers 
"""
function calc_X!(X::Matrix{<:Real}, DA::AbstractMatrix{<:Real}, B::AbstractMatrix{<:Real}, cont::Tuple{Integer, Integer}, i_branch::Integer, slack::Integer)
    x = B[cont[1], cont[2]] * DA[i_branch, cont[2]] / B[cont[1], cont[2]]
    copy!(X, B)
    neutralize_line!(X, cont[1], cont[2], -x)
    _calc_X!(X, slack)
    return X
end
calc_X(DA::AbstractMatrix{<:Real}, B::AbstractMatrix{T}, cont::Tuple{Integer, Integer}, i_branch::Integer, slack::Integer) where T<:Real =
    calc_X!(Matrix{T}(undef, size(B)), DA, B, cont, i_branch, slack)

""" 
    Calculate the inverse of the admittance matrix after a line outage which splits the system. 
    cont[1] (from_bus), cont[2] (to_bus), i_branch branch number, and island is index numbers 
"""
function calc_X(DA::AbstractMatrix{<:Real}, B::AbstractMatrix{<:Real}, cont::Tuple{Integer, Integer}, 
        i_branch::Integer, slack::Integer, island::AbstractVector{<:Integer})
    (fbus, tbus) = cont
    x = B[fbus, tbus] * DA[i_branch, tbus] / B[fbus, tbus]
    X = Matrix(B[island, island])
    c = ifelse(insorted(fbus, island), fbus, tbus)
    i = searchsortedfirst(island, c)
    X[i, i] += x
    _calc_X!(X, searchsortedfirst(island, slack))
    return X
end

function neutralize_line!(B::AbstractMatrix, i::Integer, j::Integer, val::Real)
    B[i,j] += val
    B[j,i] += val
    B[i,i] -= val
    B[j,j] -= val
end

calc_isf(A::AbstractMatrix{<:Real}, D::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real}) = D * A * X
calc_isf(DA::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real}) = DA * X
calc_isf!(ϕ::AbstractMatrix{<:Real}, DA::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real}) = 
    LinearAlgebra.mul!(ϕ, DA, X)

""" Make the isf-matrix """
function get_isf(branches::AbstractVector{<:Branch}, nodes::AbstractVector{<:Bus}, 
        idx::Dict{<:Int, <:Integer} = get_nodes_idx(nodes), slack::Integer = find_slack(nodes)[1])
    A = calc_A(branches, length(nodes), idx)
    D = calc_D(branches)
    DA = D * A
    B = calc_B(A, DA)
    X = calc_X(B, slack)
    ϕ = zeros(eltype(X), size(A))
    calc_isf!(ϕ, DA, X)
    return ϕ
end

""" 
    Make the isf-matrix after a line outage using base case D*A and B. 
    cont[1] (from_bus), cont[2] (to_bus), and i_branch are index numbers 
"""
function get_isf!(ϕ::AbstractMatrix{<:Real}, DA::AbstractMatrix{<:Real}, B::AbstractMatrix{<:Real}, cont::Tuple{Integer, Integer}, 
        i_branch::Integer, slack::Integer)
    X = calc_X(DA, B, cont, i_branch, slack)
    calc_isf!(ϕ, DA, X)
    ϕ[i_branch,:] .= 0
    return ϕ
end
function get_isf(DA::AbstractMatrix{<:Real}, B::AbstractMatrix{<:Real}, cont::Tuple{Integer, Integer}, 
        i_branch::Integer, slack::Integer)
    ϕ = similar(size(DA))
    get_isf!(ϕ, DA, B, cont, i_branch, slack)
    return ϕ
end

""" 
    Make the isf-matrix after a line outage, which splits the system, using base case D*A and B. 
    cont[1] (from_bus), cont[2] (to_bus), and i_branch are index numbers 
"""
function get_isf(DA::AbstractMatrix{<:Real}, B::AbstractMatrix{<:Real}, cont::Tuple{Integer, Integer}, 
        c::Integer, slack::Integer, nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer})
    X = calc_X(DA, B, cont, c, slack, nodes)
    return calc_isf(DA[branches, nodes], X)
end

""" Return the overload of a line, else return 0.0 """
find_overload(flow::T, rate::Real, atol::Real = 1e-6) where {T<:Real} = 
    abs(flow)-rate > atol ? sign(flow)*(abs(flow)-rate) : zero(T)

filter_overload(flow::AbstractVector{<:Real}, linerating::AbstractVector{<:Real}, atol::Real = 1e-6) = 
    [(i,ol) for (i,ol) in enumerate(find_overload.(flow, linerating)) if abs(ol) > atol]

filter_overload(Δflow::AbstractVector{<:Tuple}, linerating::AbstractVector{<:Real}, atol::Real = 1e-6) = 
    [(i,find_overload(ol, linerating[i])) for (i,ol) in Δflow if abs(find_overload(ol, linerating[i])) > atol]

""" Calculate the power flow on the lines from the voltage angles """
function calc_Pline(branches::AbstractVector{<:Branch}, θ::AbstractVector{<:Real}, idx::Dict{<:Int, <:Int})
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
calc_Pline!(pf::DCPowerFlow) = LinearAlgebra.mul!(pf.F, pf.DA, pf.θ)

""" DC line flow calculation using Injection Shift Factors and Power Injection vector"""
calculate_line_flows(isf::AbstractMatrix{<:Real}, Pᵢ::AbstractVector{<:Real}) = isf*Pᵢ
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

function calculate_island_line_flows(pf::DCPowerFlow, cont::Tuple{Integer, Integer}, c::Integer, Pᵢ::AbstractVector{<:Real})
    island, island_b = handle_islands(pf.B, pf.DA, cont, c, pf.slack)
    ptdf = get_isf(pf.DA, pf.B, cont, c, pf.slack, island, island_b)
    return ptdf * Pᵢ[island]
end

function get_contingency_ptdf(opfm::OPFmodel, pf::DCPowerFlow)
    idx = get_nodes_idx(opfm.nodes)
    contids = [(x, get_bus_idx(opfm.branches[x], idx)) for x in indexin(opfm.contingencies, opfm.branches)]
    ptdf = Array{Float64}(undef, size(pf.ϕ, 1), size(pf.ϕ, 2), length(opfm.contingencies))
    # for i in eachindex(contids)
    Threads.@threads for i in eachindex(contids)
        (c, cont) = contids[i]
        if !is_islanded(pf, cont, c)
            get_isf!(ptdf[:,:,i], pf.X, pf.B, pf.DA, cont, c)
        else
            islands = island_detection_thread_safe(pf.B, cont[1], cont[2])
            island = find_ref_island(islands, pf.slack)
            island_b = find_island_branches(islands[island], pf.DA, c)
            fill!(ptdf[:,:,i], zero(Float64))
            get_isf!(ptdf[island_b, islands[island],i], pf.DA, pf.B, cont, c, pf.slack) 
        end
    end
    return ptdf
end