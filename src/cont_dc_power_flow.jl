function neutralize_line!(B::AbstractMatrix, i::Integer, j::Integer, val::Real)
    B[i, j] += val
    B[j, i] += val
    B[i, i] -= val
    B[j, j] -= val
end

""" 
    Make the B-matrix after a line outage, which splits the system, using base case D*A and B. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function calc_cont_B(DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
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
    Calculate the inverse of the admittance matrix after a line outage which splits the system. 
    cont[1] (from_bus), cont[2] (to_bus), cont_branch branch number, and island is sorted index numbers 
"""
function calc_X(DA::AbstractMatrix{<:Real}, B::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer, island::AbstractVector{<:Integer}
)
    (fbus, tbus) = cont
    X = Matrix(B[island, island])
    c = ifelse(insorted(fbus, island), fbus, tbus)
    i = searchsortedfirst(island, c)
    X[i, i] += DA[cont_branch, tbus]
    _calc_X!(X, searchsortedfirst(island, slack))
    return X
end
function calc_X!(X::AbstractMatrix{<:Real}, X₀::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    island::AbstractVector{<:Integer}
)
    copy!(X, X₀)
    (fbus, tbus) = cont
    c = ifelse(!insorted(fbus, island), fbus, tbus)
    bus = searchsortedfirst(island, c)
    X[bus, :] .= 0.
    X[:, bus] .= 0.
    return X
end

""" 
    Find voltage angles after a line outage using base case D*A and B. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function calc_θ!(θ::AbstractVector{<:Real}, B::AbstractMatrix{<:Real}, DA::AbstractMatrix{<:Real}, 
    B0::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, 
    cont::Tuple{Integer,Integer}, cont_branch::Integer, slack::Integer
)
    copyto!(B, B0)
    neutralize_line!(B, cont[1], cont[2], DA[cont_branch, cont[1]])
    _calc_θ!(θ, B, P, slack)
    return θ
end
calc_θ(DA::AbstractMatrix{T}, B0::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, 
    cont::Tuple{Integer,Integer}, cont_branch::Integer, slack::Integer
) where {T<:Real} = calc_θ!(Vector{T}(undef, size(DA, 2)), similar(B0), DA, B0, P, cont, cont_branch, slack)

function calculate_line_flows!(F::AbstractVector{T}, θ::AbstractVector{<:Real}, B::AbstractMatrix{<:Real}, 
    DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, 
    cont::Tuple{Integer,Integer}, cont_branch::Integer, slack::Integer
) where {T<:Real}
    calc_θ!(θ, B, DA, B0, P, cont, cont_branch, slack)
    LinearAlgebra.mul!(F, DA, θ)
    F[cont_branch] = zero(T)
    return F
end

""" 
    Find voltage angles after a line outage, which splits the system, using base case B. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function calc_θ(DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, 
    cont::Tuple{Integer,Integer}, cont_branch::Integer, slack::Integer, nodes::AbstractVector{<:Integer}
)
    B = calc_cont_B(DA, B0, cont, cont_branch, nodes)
    θ = calc_θ!(B, view(P, nodes), searchsortedfirst(nodes, slack))
    return θ
end

function calculate_line_flows!(F::AbstractVector{T}, θ::AbstractVector{<:Real}, DA::AbstractMatrix{<:Real}, 
    B::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, cont::Tuple{Integer,Integer}, cont_branch::Integer, 
    slack::Integer, nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}
) where {T<:Real}
    θ[nodes] = calc_θ(DA, B, P, cont, cont_branch, slack, nodes)
    LinearAlgebra.mul!(F, DA, θ)
    zero_not_in_array!(F, branches)
    return F
end
""" 
    Make the isf-matrix after a line outage using base case D*A and B. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function calc_isf!(ϕ::AbstractMatrix{<:Real}, K::KLU.KLUFactorization{T,<:Integer}, DA::AbstractMatrix{<:Real}, 
    B0::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer}, cont_branch::Integer, slack::Integer
) where {T<:Real}
    B = calc_cont_B(DA, B0, cont, cont_branch, slack)
    KLU.klu!(K, B)
    copy!(ϕ, calc_isf(K, DA, slack))
    ϕ[cont_branch, :] .= 0.0
    return ϕ
end
calc_isf(DA::AbstractMatrix{T}, B::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer) where {T<:Real} =
    calc_isf!(Matrix{T}(undef, size(DA)), calc_klu!(B, slack), DA, B, cont, cont_branch, slack)

function calc_isf_vec!(ϕ::AbstractVector{<:Real}, K::KLU.KLUFactorization{T,<:Integer}, DA::AbstractMatrix{<:Real}, 
    B0::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer}, cont_branch::Integer, slack::Integer, branch::Integer
) where {T<:Real}
    B = calc_cont_B(DA, B0, cont, cont_branch, slack)
    KLU.klu!(K, B)
    calc_isf_vec!(ϕ, K, DA, slack, branch)
    return ϕ
end
calc_isf_vec(DA::AbstractMatrix{T}, B::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer, branch::Integer) where {T<:Real} =
    calc_isf_vec!(Vector{T}(undef, size(DA,2)), calc_klu!(B, slack), DA, B, cont, cont_branch, slack, branch)
    
function calc_cont_B(DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{T}, cont::Tuple{Integer,Integer}, cont_branch::Integer,
    slack::Integer
) where {T<:Real}
    B = copy(B0)
    B[:, slack] .= zero(T)
    B[slack, :] .= zero(T)
    B[slack, slack] = one(T)
    neutralize_line!(B, cont[1], cont[2], DA[cont_branch, cont[1]])
    return B
end

""" 
    Make the isf-matrix after a line outage, which splits the system, using base case D*A and B. 
    cont[1] (from_bus), cont[2] (to_bus), and cont_branch are index numbers 
"""
function calc_isf!(ϕ::AbstractMatrix{T}, DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer, nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}
) where {T<:Real}
    B = calc_cont_B(DA, B0, cont, cont_branch, nodes)
    # ϕ[sorted_missing(branches, size(ϕ,1)), sorted_missing(nodes, size(ϕ,2))] .= zero(T)
    fill!(ϕ, zero(T))
    ϕ[branches, nodes] = calc_isf!(B, view(DA, branches, nodes), searchsortedfirst(nodes, slack))
    return ϕ
end
calc_isf(DA::AbstractMatrix{T}, B0::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer, nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}
) where {T<:Real} = calc_isf!(Matrix{T}(undef, size(DA)), DA, B0, cont, cont_branch, slack, nodes, branches)

function calc_isf_vec!(ϕ::AbstractVector{T}, DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer, nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}, branch::Integer
) where {T<:Real}
    B = calc_cont_B(DA, B0, cont, cont_branch, nodes)
    # ϕ[sorted_missing(branches, size(ϕ,1)), sorted_missing(nodes, size(ϕ,2))] .= zero(T)
    fill!(ϕ, zero(T))
    ϕ[nodes] = calc_isf_vec!(B, view(DA, branches, nodes), searchsortedfirst(nodes, slack), searchsortedfirst(branches, branch))
    return ϕ
end
calc_isf_vec(DA::AbstractMatrix{T}, B0::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer, nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}, branch::Integer
) where {T<:Real} = calc_isf_vec!(Vector{T}(undef, size(DA,2)), DA, B0, cont, cont_branch, slack, nodes, branches, branch)

function calc_isf!(ϕ::AbstractMatrix{<:Real}, ϕ₀::AbstractMatrix{<:Real},
    nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}
)
    copy!(ϕ, ϕ₀)
    zero_not_in_array!(ϕ, nodes, Val(2))
    zero_not_in_array!(ϕ, branches, Val(1))
    return ϕ
end

function calc_isf(pf::DCPowerFlow, cont::Tuple{Real,Real}, c::Integer, islands::Vector, 
    island::Integer, island_b::Vector{<:Integer}
)
    return calc_isf!(similar(pf.ϕ), pf.ϕ, islands[island], island_b)
end

function calculate_line_flows!(F::AbstractVector{T}, ϕ::AbstractMatrix{<:Real}, ϕ₀::AbstractMatrix{<:Real}, Pᵢ::AbstractVector{<:Real}, 
    nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}
) where {T<:Real}
    ix = ones(size(ϕ,2))
    zero_not_in_array!(ix, nodes)
    LinearAlgebra.mul!(ϕ, ϕ₀, LinearAlgebra.Diagonal(ix))
    LinearAlgebra.mul!(F, ϕ, Pᵢ)
    zero_not_in_array!(F, branches)
    return F
end
