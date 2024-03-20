function neutralize_line!(B::AbstractMatrix, i::Integer, j::Integer, val::Real)
    B[i, j] += val
    B[j, i] += val
    B[i, i] -= val
    B[j, j] -= val
end

# TODO
function calc_dist_slack!(ϕ::AbstractMatrix{<:Real}, ϕ₀::AbstractMatrix{<:Real}, mgx::AbstractMatrix{<:Real}, dist_slack::AbstractVector{<:Real}, c::Integer)
    @error "not implemented"
    @assert !iszero(sum(dist_slack))
    slack_array = dist_slack / sum(dist_slack)
    c_val = slack_array[c]
    slack_array *= (c_val / (1 - c_val))
    slack_array[c] /= (c_val / (1 - c_val))
    ϕ = ϕ₀ .+ ((slack_array' * mgx) * ϕ₀')'
    return ϕ
end

""" 
    Make the B-matrix after a line outage, using base case D*A and B. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function calc_cont_B(DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{T}, cont::AbstractVector{<:Tuple{Integer,Integer}}, 
    cont_branch::AbstractVector{<:Integer}, slack::Integer
) where {T<:Real}
    B = copy(B0)
    B[:, slack] .= zero(T)
    B[slack, :] .= zero(T)
    B[slack, slack] = one(T)
    neutralize_line!.([B], first.(cont), last.(cont), getindex.([DA], cont_branch, first.(cont)))
    return B
end

""" 
    Make the B-matrix after a line outage, which splits the system, using base case D*A and B. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function calc_cont_B(DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, cont::AbstractVector{<:Tuple{Integer,Integer}}, 
    cont_branch::AbstractVector{<:Integer}, nodes::AbstractVector{<:Integer}
)
    B = copy(view(B0, nodes, nodes))
    for ((fbus, tbus), branch) in zip(cont, cont_branch)
        c = ifelse(insorted(fbus, nodes), fbus, tbus)
        i = searchsortedfirst(nodes, c)
        B[i, i] -= DA[branch, fbus] # Only value of the contingency left inside nodes
    end
    return B
end

""" 
    Calculate the inverse of the admittance matrix after a line outage which splits the system. 
    cont[1] (from_bus), cont[2] (to_bus), cont_branch branch number, and island is sorted index numbers 
"""
function calc_X(DA::AbstractMatrix{<:Real}, B::AbstractMatrix{<:Real}, cont::AbstractVector{<:Tuple{Integer,Integer}}, 
    cont_branch::AbstractVector{<:Integer}, slack::Integer, island::AbstractVector{<:Integer}
)
    X = Matrix(B[island, island])
    for ((fbus, tbus), branch) in zip(cont, cont_branch)
        c = ifelse(insorted(fbus, island), fbus, tbus)
        i = searchsortedfirst(island, c)
        X[i, i] -= DA[branch, fbus]
    end
    _calc_X!(X, searchsortedfirst(island, slack))
    return X
end
function calc_X!(X::AbstractMatrix{<:Real}, X₀::AbstractMatrix{<:Real}, cont::AbstractVector{<:Tuple{Integer,Integer}},
    island::AbstractVector{<:Integer}
)
    copy!(X, X₀)
    for (fbus, tbus) in cont
        c = ifelse(!insorted(fbus, island), fbus, tbus)
        bus = searchsortedfirst(island, c)
        X[bus, :] .= 0.
        X[:, bus] .= 0.
    end
    return X
end

""" 
    Find voltage angles after a line outage using base case D*A and B. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function calc_θ!(θ::AbstractVector{<:Real}, B::AbstractMatrix{<:Real}, DA::AbstractMatrix{<:Real}, 
    B0::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, 
    cont::AbstractVector{<:Tuple{Integer,Integer}}, cont_branch::AbstractVector{<:Integer}, slack::Integer
)
    copyto!(B, B0)
    neutralize_line!.([B], first.(cont), last.(cont), getindex.([DA], cont_branch, first.(cont)))
    _calc_θ!(θ, B, P, slack)
    return θ
end
calc_θ(DA::AbstractMatrix{T}, B0::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, 
    cont::AbstractVector{<:Tuple{Integer,Integer}}, cont_branch::AbstractVector{<:Integer}, slack::Integer
) where {T<:Real} = calc_θ!(Vector{T}(undef, size(DA, 2)), similar(B0), DA, B0, P, cont, cont_branch, slack)

function calculate_line_flows!(F::AbstractVector{T}, θ::AbstractVector{<:Real}, B::AbstractMatrix{<:Real}, 
    DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, 
    cont::AbstractVector{<:Tuple{Integer,Integer}}, cont_branch::AbstractVector{<:Integer}, slack::Integer
) where {T<:Real}
    calc_θ!(θ, B, DA, B0, P, cont, cont_branch, slack)
    LinearAlgebra.mul!(F, DA, θ)
    F[cont_branch] .= zero(T)
    return F
end

""" 
    Find voltage angles after a line outage, which splits the system, using base case B. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function calc_θ(DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, 
    cont::AbstractVector{<:Tuple{Integer,Integer}}, cont_branch::AbstractVector{<:Integer}, slack::Integer, 
    nodes::AbstractVector{<:Integer}
)
    B = calc_cont_B(DA, B0, cont, cont_branch, nodes)
    θ = calc_θ!(B, view(P, nodes), searchsortedfirst(nodes, slack))
    return θ
end

function calculate_line_flows!(F::AbstractVector{<:Real}, θ::AbstractVector{<:Real}, DA::AbstractMatrix{<:Real}, 
    B::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, cont::AbstractVector{<:Tuple{Integer,Integer}}, 
    cont_branch::AbstractVector{<:Integer}, slack::Integer, nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}
)
    θ[nodes] = calc_θ(DA, B, P, cont, cont_branch, slack, nodes)
    LinearAlgebra.mul!(F, DA, θ)
    zero_not_in_array!(F, branches)
    return F
end

""" 
    Make the isf-matrix after a line outage using base case D*A and B. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function calc_isf!(ϕ::AbstractMatrix{<:Real}, DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, 
    cont::AbstractVector{<:Tuple{Integer,Integer}}, cont_branch::AbstractVector{<:Integer}, slack::Integer
)
    B = calc_cont_B(DA, B0, cont, cont_branch, slack)
    K = calc_klu!(B, slack)
    copy!(ϕ, calc_isf(K, DA, slack))
    ϕ[cont_branch, :] .= 0.0
    return ϕ
end
calc_isf(DA::AbstractMatrix{T}, B::AbstractMatrix{<:Real}, cont::AbstractVector{<:Tuple{Integer,Integer}}, 
    cont_branch::AbstractVector{<:Integer}, slack::Integer
) where {T<:Real} =
    calc_isf!(Matrix{T}(undef, size(DA)), DA, B, cont, cont_branch, slack)

""" 
    Make the isf-matrix after a line outage, which splits the system, using base case D*A and B. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function calc_isf!(ϕ::AbstractMatrix{T}, DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, cont::AbstractVector{<:Tuple{Integer,Integer}}, 
    cont_branch::AbstractVector{<:Integer}, slack::Integer, nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}
) where {T<:Real}
    B = calc_cont_B(DA, B0, cont, cont_branch, nodes)
    # ϕ[sorted_missing(branches, size(ϕ,1)), sorted_missing(nodes, size(ϕ,2))] .= zero(T)
    fill!(ϕ, zero(T))
    ϕ[branches, nodes] = calc_isf!(B, view(DA, branches, nodes), searchsortedfirst(nodes, slack))
    return ϕ
end
calc_isf(DA::AbstractMatrix{T}, B0::AbstractMatrix{<:Real}, cont::AbstractVector{<:Tuple{Integer,Integer}}, 
    cont_branch::AbstractVector{<:Integer}, slack::Integer, nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}
) where {T<:Real} = calc_isf!(Matrix{T}(undef, size(DA)), DA, B0, cont, cont_branch, slack, nodes, branches)

""" Only for single branch contingnecies on a radial """
function calc_isf!(ϕ::AbstractMatrix{<:Real}, ϕ₀::AbstractMatrix{<:Real},
    nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}
)
    copy!(ϕ, ϕ₀)
    zero_not_in_array!(ϕ, nodes, Val(2))
    zero_not_in_array!(ϕ, branches, Val(1))
    return ϕ
end

""" Only for single branch contingnecies on a radial """
function calculate_line_flows!(F::AbstractVector{<:Real}, ϕ::AbstractMatrix{<:Real}, ϕ₀::AbstractMatrix{<:Real}, Pᵢ::AbstractVector{<:Real}, 
    nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}
)
    calc_isf!(ϕ, ϕ₀, nodes, branches)
    LinearAlgebra.mul!(F, ϕ, Pᵢ)
    zero_not_in_array!(F, branches)
    return F
end
