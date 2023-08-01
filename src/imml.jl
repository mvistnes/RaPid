# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

""" Primitive IMML """
get_power_flow_change(F::AbstractVector{<:Real}, ϕ::AbstractMatrix{<:Real}, A::AbstractMatrix{<:Integer}, branch) = 
    F .+ get_change(ϕ, A, branch) * F[branch]

function get_change(ϕ::AbstractMatrix{<:Real}, A::AbstractMatrix{<:Integer}, branch; atol::Real = 1e-5) 
    x = LinearAlgebra.I - ϕ[branch,:]'*A[branch,:]
    if isapprox(x, zero(typeof(x)); atol=atol)
        return zero(typeof(x))
    end
    return ϕ*A[branch,:]*inv(x)
end

"""
Calculation of voltage angles in a contingency case using IMML

Input:
    - X: The inverse admittance matrix
    - B: Suseptance matrix
    - DA: Diagonal suseptance matrix times the connectivity matrix
    - θ₀: Inital voltage angles
    - from_bus: From bus index
    - to_bus: To bus index
    - branch: Branch index
    - change: The amount of reactance change between the buses, <=1. 
        Default is 1 and this removes all lines
"""
@views function get_changed_angles(
        X::AbstractMatrix{<:Real}, 
        b::Real, 
        change::Real,
        θ₀::AbstractVector{<:Real}, 
        from_bus::Integer, 
        to_bus::Integer; 
        atol::Real = 1e-5
    )
    # change = DA[branch, to_bus] / B[from_bus, to_bus]
    x = change * (X[:,from_bus] - X[:,to_bus])
    c⁻¹ = 1/b + x[from_bus] - x[to_bus] 
    if isapprox(c⁻¹, zero(typeof(c⁻¹)); atol=atol)
        return throw(DivideError())
    end
    delta = 1/c⁻¹ * (θ₀[from_bus] - θ₀[to_bus])
    return θ₀ - x * delta
end

""" 
Calculate the power flow on the lines from the connectivity 
and the diagonal admittance matrices and the voltage angles
in a contingency of the branch number.
"""
function calc_Pline(pf::DCPowerFlow, cont::Tuple{Integer, Integer}, branch::Integer)
    θ = get_changed_angles(pf.X, pf.B[cont[1], cont[2]], pf.DA[branch, cont[2]] / pf.B[cont[1], cont[2]], 
        pf.θ, cont[1], cont[2])
    if θ == zero(typeof(pf.X))
        return θ
    else
        P = calc_Pline(pf.DA, θ)
        P[branch] = 0.0
        return P
    end
end

"""
Calculation of the inverse admittance matrix in a contingency case using IMML

Input:
    - X: The inverse admittance matrix
    - B: Suseptance matrix
    - DA: Diagonal suseptance matrix times the connectivity matrix
    - from_bus: From bus index
    - to_bus: To bus index
    - branch: Branch index
    - change: The amount of reactance change between the buses, <=1. 
        Default is 1 and this removes all lines
"""
@views function get_changed_X!(
        X::AbstractMatrix{<:Real}, 
        X0::AbstractMatrix{<:Real}, 
        B::AbstractMatrix{<:Real}, 
        DA::AbstractMatrix{<:Real}, 
        from_bus::Integer, 
        to_bus::Integer, 
        branch::Integer; 
        atol::Real = 1e-5
    )

    # X_new = X - (X*A[from_bus,:]*DA[from_bus,to_bus]*x)/(1+DA[from_bus,to_bus]*x*A[from_bus,:])
    change = DA[branch, to_bus] / B[from_bus, to_bus]
    x = X0[:,from_bus] - X0[:,to_bus]
    c⁻¹ = 1/B[to_bus, from_bus] + change * (x[from_bus] - x[to_bus])
    if isapprox(c⁻¹, zero(typeof(c⁻¹)); atol=atol)
        throw(DivideError())
    end
    # delta = 1/c⁻¹ * x
    copy!(X, X0)
    LinearAlgebra.mul!(X, x, x', -change * 1/c⁻¹, true) # mul!(C, A, B, α, β) -> C == $A B α + C β$
end

" Get the isf-matrix after a line outage using IMML "
function get_isf!(
        isf::AbstractMatrix{<:Real}, 
        X::AbstractMatrix{<:Real}, 
        X0::AbstractMatrix{<:Real}, 
        B::AbstractMatrix{<:Real}, 
        DA::AbstractMatrix{<:Real}, 
        cont::Tuple{Integer, Integer}, 
        branch::Integer
    )
    get_changed_X!(X, X0, B, DA, cont[1], cont[2], branch)
    calc_isf!(isf, DA, X)
    isf[branch,:] .= 0
end

function get_isf(
        pf::DCPowerFlow, 
        cont::Tuple{Integer, Integer}, 
        branch::Integer
    )
    isf = similar(pf.ϕ)
    X = similar(pf.X)
    get_isf!(isf, X, pf.X, pf.B, pf.DA, cont, branch)
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
@views function calculate_line_flows!(
        Pl::AbstractVector{<:Real},
        Pl0::AbstractVector{<:Real},
        ptdf::AbstractMatrix{<:Real},
        B::AbstractMatrix{<:Real},
        DA::AbstractMatrix{<:Real},
        X::AbstractMatrix{<:Real}, 
        θ₀::AbstractVector{<:Real}, 
        from_bus::Integer, 
        to_bus::Integer, 
        branch::Integer; 
        atol::Real = 1e-5
    )
    change = DA[branch, to_bus] / B[from_bus, to_bus]
    # x = change * (X[:,from_bus] - X[:,to_bus])
    c⁻¹ = 1/B[from_bus, to_bus] + change * (X[from_bus,from_bus] - X[from_bus,to_bus] - X[to_bus,from_bus] + X[to_bus,to_bus])
    if isapprox(c⁻¹, zero(typeof(c⁻¹)); atol=atol)
        return throw(DivideError())
    end
    delta = 1/c⁻¹ * (θ₀[from_bus] - θ₀[to_bus])
    @. Pl = Pl0 - (ptdf[:, from_bus] - ptdf[:, to_bus]) * change * delta
    Pl[branch] = 0.0
end

function calculate_line_flows!(
        Pl::AbstractVector{<:Real},
        pf::DCPowerFlow,
        cont::Tuple{Integer, Integer}, 
        branch::Integer
    )
    return calculate_line_flows!(Pl, pf.F, pf.ϕ, pf.B, pf.DA, pf.X, pf.θ, cont[1], cont[2], branch)
end
function calculate_line_flows(
        pf::DCPowerFlow,
        cont::Tuple{Integer, Integer}, 
        branch::Integer
    )
    Pl = similar(pf.F)
    calculate_line_flows!(Pl, pf, cont, branch)
    return Pl
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
