# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

mutable struct IMML
    A::AbstractMatrix{<:Real}
    D::AbstractMatrix{<:Real}
    DA::AbstractMatrix{<:Real}
    B::AbstractMatrix{<:Real}
    X::AbstractMatrix{<:Real}
    ptdf::AbstractMatrix{<:Real}
    contingencies::Vector
    linerating::Vector{<:Real}
end

function IMML(nodes, branches, idx, contingencies)
    numnodes = length(nodes)
    A = calc_A(branches, numnodes, idx)
    D = calc_D(branches)
    DA = D*A
    B = fast_calc_B(A, DA)
    slack = find_slack(nodes)[1]
    X = calc_X(B, slack)
    ptdf = calc_isf(DA, X)
    return IMML(A, D, DA, B, X, ptdf, contingencies, get_rate.(branches))
end

"""
    Calculation of voltage angles in a contingency case using IMML

    Input:
        - X: The inverse admittance matrix
        - B: The admittance matrix
        - δ₀: Inital voltage angles
        - from_bus: From bus index
        - to_bus: To bus index
        - slack: Slack bus index
        - change: The amount of reactance change between the buses, <=1. 
          Default is 1 and this removes all lines
"""
function get_changed_angles(
        X::AbstractMatrix{<:Real}, 
        B::AbstractMatrix{<:Real}, 
        δ₀::AbstractVector{<:Real}, 
        from_bus::Integer, 
        to_bus::Integer, 
        slack::Integer,
        change::Real = 1.0; 
        atol=1e-5
    )

    x = change .* (X[:,from_bus] .- X[:,to_bus])
    x[slack] = 0.0
    c⁻¹ = 1/B[from_bus,to_bus] + x[from_bus] - x[to_bus] # Sjekk om dette blir ~inf ved øyer
    if isapprox(c⁻¹, inf; atol=atol)
        throw(DivideError())
    end
    delta = 1/c⁻¹ * (δ₀[from_bus] - δ₀[to_bus])
    # return .- x .* delta
    return δ₀ .- x .* delta
end

"""
    Calculation of line flow in a contingency case using IMML.
    If the output is approx. zero there are islands in the system.

    Input:
        - Pl0: Initial line power flow
        - ptdf: Power Transfer Distribution Factor matrix
        - B: The admittance matrix
        - δ₀: Inital voltage angles
        - from_bus: From bus index
        - to_bus: To bus index
        - slack: Slack bus index
        - change: The amount of reactance change between the buses, <=1. 
          Default is 1 and this removes all lines
"""
@views function calculate_delta_line_flows(
        ptdf::AbstractMatrix{<:Real},
        X::AbstractMatrix{<:Real}, 
        B::AbstractMatrix{<:Real}, # could be changed to D-matrix, need only the negative admittance of the line
        δ₀::AbstractVector{<:Real}, 
        from_bus::Integer, 
        to_bus::Integer, 
        branch::Integer,
        change::Real = 1.0
    )
    x = change .* (X[:,from_bus] .- X[:,to_bus])
    c⁻¹ = 1/B[from_bus,to_bus] + x[from_bus] - x[to_bus]
    delta = 1/c⁻¹ * (δ₀[from_bus] - δ₀[to_bus])
    Pl = (ptdf[:, from_bus] .- ptdf[:, to_bus]) .* change .* delta
    # Pl[branch] = 0.0 # OBS: double check!!
    return Pl
end

function calculate_line_flows(
        Pl0::AbstractVector{<:Real},
        ptdf::AbstractMatrix{<:Real},
        X::AbstractMatrix{<:Real}, 
        B::AbstractMatrix{<:Real}, # could be changed to D-matrix, need only the negative admittance of the line
        δ₀::AbstractVector{<:Real}, 
        from_bus, to_bus, branch
    )
    @assert length(from_bus) == length(to_bus) == length(branch)
    Pl = Pl0 .- sum(calculate_delta_line_flows(ptdf, X, B, δ₀, f, t, b, 1.0) 
                    for (f,t,b) in zip(from_bus, to_bus, branch))
    return Pl
end

get_overload(
        Pl0::AbstractVector{<:Real}, 
        ptdf::AbstractMatrix{<:Real}, 
        B::AbstractMatrix{<:Real}, 
        X::AbstractMatrix{<:Real}, 
        δ::AbstractVector{<:Real}, 
        branch::Integer,
        cont::Tuple{Integer, Integer},
        linerating::AbstractVector{<:Real}
    ) = 
        find_overload.(
                calculate_line_flows(Pl0, ptdf, X, B, δ, cont[1], cont[2], branch), 
                linerating
            )

get_overload(
        Pl0::AbstractVector{<:Real}, 
        δ::AbstractVector{<:Real}, 
        imml::IMML,
        branch::Integer,
        cont::Tuple{Integer, Integer}
    ) = 
        find_overload.(
                calculate_line_flows(Pl0, imml.ptdf, imml.X, imml.B, δ, cont[1], cont[2], branch), 
                imml.linerating
            )

""" LODF value for a contingency at line l_mn change in line k_pq 
    From the book Optimization of power system operation """
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
