# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
# using ReusePatterns

# mutable struct IMML
#     pf::DCPowerFlow
#     contingencies::Vector
#     linerating::Vector{<:Real}
# end

# function IMML(nodes::AbstractVector{<:Bus}, branches::AbstractVector{<:Branch}, 
#         idx::Dict{<:Any, <:Int}, Pᵢ::AbstractVector{<:Real}, contingencies
#     )
#     return IMML(DCPowerFlow(nodes, branches, idx, Pᵢ), contingencies, get_rate.(branches))
# end
# @forward((IMML, :pf), DCPowerFlow)
# contingencies(imml::IMML) = imml.contingencies
# linerating(imml::IMML) = imml.linerating

"""
Calculation of voltage angles in a contingency case using IMML

Input:
    - X: The inverse admittance matrix
    - DA: Diagonal suseptance matrix times the connectivity matrix
    - θ₀: Inital voltage angles
    - from_bus: From bus index
    - to_bus: To bus index
    - slack: Slack bus index
    - change: The amount of reactance change between the buses, <=1. 
        Default is 1 and this removes all lines
"""
@views function get_changed_angles(
        X::AbstractMatrix{<:Real}, 
        B::AbstractMatrix{<:Real}, 
        DA::AbstractMatrix{<:Real},
        θ₀::AbstractVector{<:Real}, 
        from_bus::Integer, 
        to_bus::Integer, 
        branch::Integer; 
        atol=1e-5
    )
    change = DA[branch, to_bus] / B[from_bus, to_bus]
    x = change .* (X[:,from_bus] .- X[:,to_bus])
    # x[slack] = 0.0
    c⁻¹ = 1/B[from_bus,to_bus] + x[from_bus] - x[to_bus] 
    delta = 1/c⁻¹ * (θ₀[from_bus] - θ₀[to_bus])
    if isapprox(delta, 0.0; atol=atol)
        throw(DivideError())
    end
    return θ₀ .- x .* delta
end


"""
Calculation of the inverse admittance matrix in a contingency case using IMML

Input:
    - X: The inverse admittance matrix
    - DA: Diagonal suseptance matrix times the connectivity matrix
    - from_bus: From bus index
    - to_bus: To bus index
    - slack: Slack bus index
    - change: The amount of reactance change between the buses, <=1. 
        Default is 1 and this removes all lines
"""
@views function get_changed_X(
        X::AbstractMatrix{<:Real}, 
        B::AbstractMatrix{<:Real}, 
        DA::AbstractMatrix{<:Real},
        from_bus::Integer, 
        to_bus::Integer, 
        branch::Integer; 
        atol=1e-5
    )

    # X_new = X - (X*A[from_bus,:]*DA[from_bus,to_bus]*x)/(1+DA[from_bus,to_bus]*x*A[from_bus,:])
    change = DA[branch, to_bus] / B[from_bus, to_bus]
    x = X[:,from_bus] .- X[:,to_bus]
    c⁻¹ = 1/B[from_bus,to_bus] + change * (x[from_bus] - x[to_bus])
    if isapprox(c⁻¹, 0.0; atol=atol)
        throw(DivideError())
    end
    delta = 1/c⁻¹ .* x
    return X - change .* x .* delta'
end

" Get the isf-matrix after a line outage using IMML "
function get_isf(pf::DCPowerFlow, from_bus_idx::Integer, to_bus_idx::Integer, i_branch::Integer)
    isf = calc_isf(pf.DA, get_changed_X(pf.X, pf.B, pf.DA, from_bus_idx, to_bus_idx, i_branch))
    isf[i_branch,:] .= 0
    return isf
end

"""
Calculation of line flow in a contingency case using IMML.
If DivideError, there are island(s) in the system.

Input:
    - Pl0: Initial line power flow
    - ptdf: Power Transfer Distribution Factor matrix
    - B: Suseptance matrix
    - DA: Diagonal suseptance matrix times the connectivity matrix
    - θ₀: Inital voltage angles
    - from_bus: From bus index
    - to_bus: To bus index
    - slack: Slack bus index
    - change: The amount of reactance change between the buses, <=1. 
        Default is 1 and this removes all lines
"""
@views function calculate_line_flows(
        Pl0::AbstractVector{<:Real},
        ptdf::AbstractMatrix{<:Real},
        B::AbstractMatrix{<:Real}, 
        DA::AbstractMatrix{<:Real},
        X::AbstractMatrix{<:Real}, 
        θ₀::AbstractVector{<:Real}, 
        from_bus::Integer, 
        to_bus::Integer, 
        branch::Integer; 
        atol=1e-5
    )
    change = DA[branch, to_bus] / B[from_bus, to_bus]
    x = change .* (X[:,from_bus] .- X[:,to_bus])
    c⁻¹ = 1/B[from_bus,to_bus] + x[from_bus] - x[to_bus]
    delta = 1/c⁻¹ * (θ₀[from_bus] - θ₀[to_bus])
    if isapprox(delta, 0.0; atol=atol)
        throw(DivideError())
    end
    Pl = Pl0 .- (ptdf[:, from_bus] .- ptdf[:, to_bus]) .* change .* delta
    Pl[branch] = 0.0
    return Pl
end

function get_overload(
        Pl0::AbstractVector{<:Real}, 
        ptdf::AbstractMatrix{<:Real}, 
        B::AbstractMatrix{<:Real},
        DA::AbstractMatrix{<:Real}, 
        X::AbstractMatrix{<:Real}, 
        θ::AbstractVector{<:Real}, 
        branch::Integer,
        cont::Tuple{Integer, Integer},
        linerating::AbstractVector{<:Real}
    )
    find_overload.(
                calculate_line_flows(Pl0, ptdf, B, DA, X, θ, cont[1], cont[2], branch), 
                linerating
            )
end

function get_overload(
        Pl0::AbstractVector{<:Real}, 
        pf::DCPowerFlow,
        branch::Integer,
        cont::Tuple{Integer, Integer},
        linerating::AbstractVector{<:Real}
    )
    find_overload.(
                calculate_line_flows(Pl0, pf.ptdf, pf.B, pf.DA, pf.X, pf.θ, cont[1], cont[2], branch), 
                linerating
            )
end

""" 
LODF value for a contingency at line l_mn change in line k_pq 
    From the book Optimization of power system operation 
"""
@views get_lodf(x_l::Real, m::Integer, n::Integer, x_k::Real, p::Integer, q::Integer, X::AbstractMatrix) = 
    (x_l / x_k) * (X[p,m] - X[q,m] - X[p,n] + X[q,n]) / 
    (x_l - (X[m,m] + X[n,n] - 2 * X[m,n]))
@views get_lodf(x_l::Real, m::Integer, n::Integer, x_k::AbstractVector{<:Real}, A::AbstractMatrix, X::AbstractMatrix) = 
    (x_l ./ x_k) .* A * (X[:,m] - X[:,n]) ./ 
    (x_l - (X[m,m] + X[n,n] - 2 * X[m,n]))

function get_lodf(branch_l::Branch, branch_k::Branch, X::AbstractMatrix, idx::Dict{<:Any, <:Int}) 
    (m, n) = get_bus_idx(branch_l, idx)
    (p, q) = get_bus_idx(branch_k, idx)
    return get_lodf(get_x(branch_l), m, n, get_x(branch_k), p, q, X)
end
function get_lodf(branch_l::Branch, branches::AbstractVector{<:Branch}, A::AbstractMatrix, X::AbstractMatrix, idx::Dict{<:Any, <:Int}) 
    (m, n) = get_bus_idx(branch_l, idx)
    return get_lodf(get_x(branch_l), m, n, get_x.(branches), A, X)
end
function get_lodf(from_bus, to_bus, x::AbstractVector{<:Real}, A::AbstractMatrix, X::AbstractMatrix)
    mx = reshape(reduce(vcat, get_lodf.(x, from_bus, to_bus, [x], [A], [X])), (length(x),length(x)))
    return mx - LinearAlgebra.Diagonal(mx) - LinearAlgebra.I
end
