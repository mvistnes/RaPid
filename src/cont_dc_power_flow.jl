"""
    Remove a line from the B-matrix by neutralizing its effect.

    Parameters:
    - `B`: The B-matrix to modify.
    - `i`: From bus index.
    - `j`: To bus index.
    - `val`: The value of the line to remove.
"""
function neutralize_line!(B::AbstractMatrix, i::Integer, j::Integer, val::Real)
    B[i, j] += val
    B[j, i] += val
    B[i, i] -= val
    B[j, j] -= val
end

""" 
    Make the B-matrix after a line outage, which splits the system, using base case D*A and B. 

    Parameters:
    - `DA`: The D*A matrix.
    - `B0`: The base B-matrix.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `cont_branch`: The index of the outaged branch.
    - `nodes`: Sorted indices of the nodes in the system (island).
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

    Parameters:
    - `DA`: The D*A matrix.
    - `B`: The base B-matrix.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `cont_branch`: The index of the outaged branch.
    - `slack`: The index of the slack bus.
    - `island`: Sorted indices of the nodes in the system (island).
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

""" 
    Calculate the inverse of the admittance matrix after a line outage which splits the system. 

    Parameters:
    - `X`: Matrix to store the result.
    - `X₀`: Base inverse admittance matrix.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `island`: Sorted indices of the nodes in the system (island).
"""
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
    
    Parameters:
    - `θ`: Vector to store the voltage angles.
    - `B`: The B-matrix.
    - `DA`: The D*A matrix.
    - `B0`: The base B-matrix.
    - `P`: The active power injection vector.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `cont_branch`: The index of the outaged branch.
    - `slack`: The index of the slack bus.
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

""" 
    Find voltage angles after a line outage using base case D*A and B. 
    
    Parameters:
    - `DA`: The D*A matrix.
    - `B0`: The base B-matrix.
    - `P`: The active power injection vector.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `cont_branch`: The index of the outaged branch.
    - `slack`: The index of the slack bus.
"""
calc_θ(DA::AbstractMatrix{T}, B0::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, 
    cont::Tuple{Integer,Integer}, cont_branch::Integer, slack::Integer
) where {T<:Real} = calc_θ!(Vector{T}(undef, size(DA, 2)), similar(B0), DA, B0, P, cont, cont_branch, slack)

""" 
    Calculate line flows after a line outage using base case D*A and B. 
    
    Parameters:
    - `F`: Vector to store the line flows.
    - `θ`: Voltage angles.
    - `B`: The B-matrix.
    - `DA`: The D*A matrix.
    - `B0`: The base B-matrix.
    - `P`: The active power injection vector.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `cont_branch`: The index of the outaged branch.
    - `slack`: The index of the slack bus.
"""
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
    
    Parameters:
    - `DA`: The D*A matrix.
    - `B0`: The base B-matrix.
    - `P`: The active power injection vector.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `cont_branch`: The index of the outaged branch.
    - `slack`: The index of the slack bus.
    - `nodes`: Sorted indices of the nodes in the system (island).
"""
function calc_θ(DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, P::AbstractVector{<:Real}, 
    cont::Tuple{Integer,Integer}, cont_branch::Integer, slack::Integer, nodes::AbstractVector{<:Integer}
)
    B = calc_cont_B(DA, B0, cont, cont_branch, nodes)
    θ = calc_θ!(B, view(P, nodes), searchsortedfirst(nodes, slack))
    return θ
end

""" 
    Calculate line flows after a line outage, which splits the system, using base case D*A and B. 
    
    Parameters:
    - `F`: Vector to store the line flows.
    - `θ`: Voltage angles.
    - `DA`: The D*A matrix.
    - `B`: The B-matrix.
    - `P`: The active power injection vector.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `cont_branch`: The index of the outaged branch.
    - `slack`: The index of the slack bus.
    - `nodes`: Sorted indices of the nodes in the system (island).
    - `branches`: Sorted indices of the branches in the system (island).
"""
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
    Make the ptdf-matrix after a line outage using base case D*A and B. 
    
    Parameters:
    - `ϕ`: Matrix to store the ptdf-matrix.
    - `K`: KLU factorization of the B-matrix.
    - `DA`: The D*A matrix.
    - `B0`: The base B-matrix.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `cont_branch`: The index of the outaged branch.
    - `slack`: The index of the slack bus.
"""
function calc_ptdf!(ϕ::AbstractMatrix{<:Real}, K::KLU.KLUFactorization{T,<:Integer}, DA::AbstractMatrix{<:Real}, 
    B0::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer}, cont_branch::Integer, slack::Integer
) where {T<:Real}
    B = calc_cont_B(DA, B0, cont, cont_branch, slack)
    KLU.klu!(K, B)
    copy!(ϕ, calc_ptdf(K, DA, slack))
    ϕ[cont_branch, :] .= 0.0
    return ϕ
end

""" 
    Make the ptdf-matrix after a line outage using base case D*A and B. 
    
    Parameters:
    - `DA`: The D*A matrix.
    - `B`: The base B-matrix.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `cont_branch`: The index of the outaged branch.
    - `slack`: The index of the slack bus.
"""
calc_ptdf(DA::AbstractMatrix{T}, B::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer) where {T<:Real} =
    calc_ptdf!(Matrix{T}(undef, size(DA)), calc_klu!(B, slack), DA, B, cont, cont_branch, slack)

"""
    Make a ptdf-vector after a line outage using base case D*A and B. 
    
    Parameters:
    - `ϕ`: Vector to store the ptdf-vector.
    - `K`: KLU factorization of the B-matrix.
    - `DA`: The D*A matrix.
    - `B0`: The base B-matrix.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `cont_branch`: The index of the outaged branch.
    - `slack`: The index of the slack bus.
    - `branch`: The index of the branch for which to calculate the ptdf-vector.
"""
function calc_ptdf_vec!(ϕ::AbstractVector{<:Real}, K::KLU.KLUFactorization{T,<:Integer}, DA::AbstractMatrix{<:Real}, 
    B0::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer}, cont_branch::Integer, slack::Integer, branch::Integer
) where {T<:Real}
    B = calc_cont_B(DA, B0, cont, cont_branch, slack)
    KLU.klu!(K, B)
    calc_ptdf_vec!(ϕ, K, DA, slack, branch)
    return ϕ
end

"""
    Make a ptdf-vector after a line outage using base case D*A and B. 
    
    Parameters:
    - `DA`: The D*A matrix.
    - `B`: The base B-matrix.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `cont_branch`: The index of the outaged branch.
    - `slack`: The index of the slack bus.
    - `branch`: The index of the branch for which to calculate the ptdf-vector.
"""
calc_ptdf_vec(DA::AbstractMatrix{T}, B::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer, branch::Integer) where {T<:Real} =
    calc_ptdf_vec!(Vector{T}(undef, size(DA,2)), calc_klu!(B, slack), DA, B, cont, cont_branch, slack, branch)

"""
    Make the B-matrix after a line outage, which splits the system, using base case D*A and B. 
    
    Parameters:
    - `DA`: The D*A matrix.
    - `B0`: The base B-matrix.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `cont_branch`: The index of the outaged branch.
    - `slack`: The index of the slack bus.
"""
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
    Make the ptdf-matrix after a line outage, which splits the system, using base case D*A and B. 
    
    Parameters:
    - `ϕ`: Matrix to store the ptdf-matrix.
    - `DA`: The D*A matrix.
    - `B0`: The base B-matrix.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `cont_branch`: The index of the outaged branch.
    - `slack`: The index of the slack bus.
    - `nodes`: Sorted indices of the nodes in the system (island).
    - `branches`: Sorted indices of the branches in the system (island). 
"""
function calc_ptdf!(ϕ::AbstractMatrix{T}, DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer, nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}
) where {T<:Real}
    B = calc_cont_B(DA, B0, cont, cont_branch, nodes)
    # ϕ[sorted_missing(branches, size(ϕ,1)), sorted_missing(nodes, size(ϕ,2))] .= zero(T)
    fill!(ϕ, zero(T))
    ϕ[branches, nodes] = calc_ptdf!(B, view(DA, branches, nodes), searchsortedfirst(nodes, slack))
    return ϕ
end

""" 
    Make the ptdf-matrix after a line outage, which splits the system, using base case D*A and B. 
    
    Parameters:
    - `DA`: The D*A matrix.
    - `B0`: The base B-matrix.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `cont_branch`: The index of the outaged branch.
    - `slack`: The index of the slack bus.
    - `nodes`: Sorted indices of the nodes in the system (island).
    - `branches`: Sorted indices of the branches in the system (island). 
"""
calc_ptdf(DA::AbstractMatrix{T}, B0::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer, nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}
) where {T<:Real} = calc_ptdf!(Matrix{T}(undef, size(DA)), DA, B0, cont, cont_branch, slack, nodes, branches)

""" 
    Make a ptdf-vector after a line outage, which splits the system, using base case D*A and B. 
    
    Parameters:
    - `ϕ`: Vector to store the ptdf-vector.
    - `DA`: The D*A matrix.
    - `B0`: The base B-matrix.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `cont_branch`: The index of the outaged branch.
    - `slack`: The index of the slack bus.
    - `nodes`: Sorted indices of the nodes in the system (island).
    - `branches`: Sorted indices of the branches in the system (island).
"""
function calc_ptdf_vec!(ϕ::AbstractVector{T}, DA::AbstractMatrix{<:Real}, B0::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer, nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}, branch::Integer
) where {T<:Real}
    B = calc_cont_B(DA, B0, cont, cont_branch, nodes)
    # ϕ[sorted_missing(branches, size(ϕ,1)), sorted_missing(nodes, size(ϕ,2))] .= zero(T)
    fill!(ϕ, zero(T))
    ϕ[nodes] = calc_ptdf_vec!(B, view(DA, branches, nodes), searchsortedfirst(nodes, slack), searchsortedfirst(branches, branch))
    return ϕ
end

""" 
    Make a ptdf-vector after a line outage, which splits the system, using base case D*A and B. 
    
    Parameters:
    - `DA`: The D*A matrix.
    - `B0`: The base B-matrix.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `cont_branch`: The index of the outaged branch.
    - `slack`: The index of the slack bus.
    - `nodes`: Sorted indices of the nodes in the system (island).
    - `branches`: Sorted indices of the branches in the system (island).
    - `branch`: The index of the branch for which to calculate the ptdf-vector.
"""
calc_ptdf_vec(DA::AbstractMatrix{T}, B0::AbstractMatrix{<:Real}, cont::Tuple{Integer,Integer},
    cont_branch::Integer, slack::Integer, nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}, branch::Integer
) where {T<:Real} = calc_ptdf_vec!(Vector{T}(undef, size(DA,2)), DA, B0, cont, cont_branch, slack, nodes, branches, branch)

"""
    Calculate the ptdf-matrix after a line outage, which splits the system, using base case D*A and B. 

    Parameters:
    - `ϕ`: Matrix to store the ptdf-matrix.
    - `ϕ₀`: Base ptdf-matrix.
    - `nodes`: Sorted indices of the nodes in the system (island).
    - `branches`: Sorted indices of the branches in the system (island).
"""
function calc_ptdf!(ϕ::AbstractMatrix{<:Real}, ϕ₀::AbstractMatrix{<:Real},
    nodes::AbstractVector{<:Integer}, branches::AbstractVector{<:Integer}
)
    copy!(ϕ, ϕ₀)
    zero_not_in_array!(ϕ, nodes, Val(2))
    zero_not_in_array!(ϕ, branches, Val(1))
    return ϕ
end

"""
    Calculate the ptdf-matrix after a line outage, which splits the system, using base case D*A and B. 

    Parameters:
    - `pf`: A `DCPowerFlow` object containing the power flow data.
    - `cont`: A tuple containing the from and to bus indices of the line outage.
    - `c`: The index of the contingency.
    - `islands`: Vector of island indices.
    - `island`: The index of the island.
    - `island_b`: Vector of branch indices in the island.
"""
function calc_ptdf(pf::DCPowerFlow, cont::Tuple{Real,Real}, c::Integer, islands::Vector, 
    island::Integer, island_b::Vector{<:Integer}
)
    return calc_ptdf!(similar(pf.ϕ), pf.ϕ, islands[island], island_b)
end

"""
    Calculate line flows after a line outage, which splits the system, using base case D*A and B. 

    Parameters:
    - `F`: Vector to store the line flows.
    - `ϕ`: Voltage angles.
    - `ϕ₀`: Base ptdf-matrix.
    - `Pᵢ`: Active power injection vector.
    - `nodes`: Sorted indices of the nodes in the system (island).
    - `branches`: Sorted indices of the branches in the system (island).
"""
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
