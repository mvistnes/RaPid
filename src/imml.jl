""" Primitive IMML """
get_power_flow_change(F::AbstractVector{<:Real}, ϕ::AbstractMatrix{<:Real}, A::AbstractMatrix{<:Integer}, branch) =
    F .+ get_change(ϕ, A, branch) * F[branch]

function get_change(ϕ::AbstractMatrix{<:Real}, A::AbstractMatrix{<:Integer}, branch; atol::Real=1e-10)
    x = LinearAlgebra.I - ϕ[branch, :]' * A[branch, :]
    if isapprox(x, zero(typeof(x)); atol=atol)
        return zeros(typeof(x), size(x))
    end
    return ϕ * A[branch, :] * inv(x)
end

""" Multi contingency Woodbury """
@views function get_changed_X!(
    X::AbstractMatrix{T},
    X₀::AbstractMatrix{T},
    B::AbstractMatrix{T},
    DA::AbstractMatrix{T},
    bx::AbstractVector{<:Tuple{Integer, Integer}},
    branches::AbstractVector{<:Integer};
    atol::Real=1e-10
) where {T<:Real}
    iE = X₀
    F = X₀[:, last.(bx)]
    iG = inv(B[first.(bx), last.(bx)])
    H = F'
    X = iE - iE*F*inv(iG + H*iE*F)*H*iE
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
"""
@views function get_changed_angles!(
    θ::AbstractVector{T},
    X::AbstractMatrix{T},
    B::AbstractMatrix{T},
    DA::AbstractMatrix{T},
    θ₀::AbstractVector{T},
    from_bus::Integer,
    to_bus::Integer,
    branch::Integer;
    atol::Real=1e-10
) where {T<:Real}
    change = DA[branch, to_bus] / B[from_bus, to_bus]
    # x = change * (X[:, from_bus] - X[:, to_bus])
    c⁻¹ = inv(B[from_bus, to_bus]) + change * (X[from_bus, from_bus] - X[from_bus, to_bus] - X[to_bus, from_bus] + X[to_bus, to_bus])
    if isapprox(c⁻¹, zero(T); atol=atol)
        if size(SparseArrays.getindex(B,from_bus,:).nzind, 1) <= 2 # Assummes only one islanded bus
            c⁻¹ = inv(B[from_bus, to_bus]) + change * (- X[from_bus, to_bus] + X[to_bus, to_bus])
            delta = inv(c⁻¹) * (θ₀[from_bus] - θ₀[to_bus])
            @. θ = θ₀ - (- X[:, to_bus]) * delta * change
        elseif size(SparseArrays.getindex(B,to_bus,:).nzind, 1) <= 2
            c⁻¹ = inv(B[from_bus, to_bus]) + change * (X[from_bus, from_bus] - X[to_bus, from_bus])
            delta = inv(c⁻¹) * (θ₀[from_bus] - θ₀[to_bus])
            @. θ = θ₀ - (X[:, from_bus]) * delta * change
        else
            return throw(DivideError())
        end
    else
        delta = inv(c⁻¹) * (θ₀[from_bus] - θ₀[to_bus])
        @. θ = θ₀ - (X[:, from_bus] - X[:, to_bus]) * delta * change
    end
    return θ
end

""" 
Calculate the power flow on the lines from the connectivity 
and the diagonal admittance matrices and the voltage angles
in a contingency of the branch number.
"""
function calc_Pline!(
    F::AbstractVector{<:Real},
    θ::AbstractVector{<:Real},
    X::AbstractMatrix{<:Real},
    B::AbstractMatrix{<:Real},
    DA::AbstractMatrix{<:Real},
    θ₀::AbstractVector{<:Real},
    cont::Tuple{Integer,Integer},
    branch::Integer
)
    get_changed_angles!(θ, X, B, DA, θ₀, cont[1], cont[2], branch)
    LinearAlgebra.mul!(F, DA, θ)
    F[branch] = 0.0
    return F
end
function calc_Pline(pf::DCPowerFlow, cont::Tuple{Integer,Integer}, branch::Integer)
    θ = similar(pf.θ)
    F = similar(pf.F)
    calc_Pline!(F, θ, pf.X, pf.B, pf.DA, pf.θ, cont, branch)
    return F
end

function calc_Pline!(
    F::AbstractVector{<:Real},
    θ::AbstractVector{<:Real},
    X::AbstractMatrix{<:Real},
    B::AbstractMatrix{<:Real},
    DA::AbstractMatrix{<:Real},
    θ₀::AbstractVector{<:Real},
    Pᵢ::AbstractVector{<:Real},
    cont::Tuple{Integer,Integer},
    branch::Integer
)
    θ₂ = calc_θ(X, Pᵢ)
    return calc_Pline!(F, θ, X, B, DA, θ₂, cont, branch)
end

"""
Calculation of the inverse admittance matrix in a contingency case using IMML

Input:
    - X: Container for the output matrix
    - X0: The inverse admittance matrix
    - B: Suseptance matrix
    - DA: Diagonal suseptance matrix times the connectivity matrix
    - from_bus: From bus index
    - to_bus: To bus index
    - branch: Branch index
"""
@views function get_changed_X!(
    X::AbstractMatrix{T},
    X₀::AbstractMatrix{T},
    B::AbstractMatrix{T},
    DA::AbstractMatrix{T},
    from_bus::Integer,
    to_bus::Integer,
    branch::Integer;
    atol::Real=1e-10
) where {T<:Real}

    # X_new = X - (X*A[from_bus,:]*DA[from_bus,to_bus]*x)/(1+DA[from_bus,to_bus]*x*A[from_bus,:])
    change = DA[branch, to_bus] / B[from_bus, to_bus]
    x = X₀[:, from_bus] - X₀[:, to_bus]
    c⁻¹ = inv(B[to_bus, from_bus]) + change * (x[from_bus] - x[to_bus])
    if isapprox(c⁻¹, zero(T); atol=atol)
        if size(SparseArrays.getindex(B,from_bus,:).nzind, 1) <= 2 # Assummes only one islanded bus
            c⁻¹ = inv(B[to_bus, from_bus]) + change * (- x[to_bus])
            copy!(X, X₀)
            X[:,from_bus] .= zero(T)
            X[from_bus,:] .= zero(T)
            X[from_bus,from_bus] = one(T)
        elseif size(SparseArrays.getindex(B,to_bus,:).nzind, 1) <= 2
            c⁻¹ = inv(B[to_bus, from_bus]) + change * (x[from_bus])
            copy!(X, X₀)
            X[:,to_bus] .= zero(T)
            X[to_bus,:] .= zero(T)
            X[to_bus,to_bus] = one(T)
        else
            throw(DivideError())
        end
    else
        # delta = 1/c⁻¹ * x
        copy!(X, X₀)
    end
    LinearAlgebra.mul!(X, x, x', -change * inv(c⁻¹), true) # mul!(C, A, B, α, β) -> C == $A B α + C β$
end

""" Get the isf-matrix after a line outage using IMML. 
    isf and X are containers for output and calculation and will be overwritten """
function get_isf!(
    isf::AbstractMatrix{<:Real},
    X::AbstractMatrix{<:Real},
    X₀::AbstractMatrix{<:Real},
    B::AbstractMatrix{<:Real},
    DA::AbstractMatrix{<:Real},
    cont::Tuple{Integer,Integer},
    branch::Integer
)
    get_changed_X!(X, X₀, B, DA, cont[1], cont[2], branch)
    calc_isf!(isf, DA, X)
    isf[branch, :] .= 0.0
    nothing
end

""" Get the isf-matrix after a line outage using IMML. 
     """
function get_isf(
    pf::DCPowerFlow,
    cont::Tuple{Integer,Integer},
    branch::Integer
)
    get_isf!(pf.mbn_tmp, pf.mnn_tmp, pf.X, pf.B, pf.DA, cont, branch)
    return pf.mbn_tmp
end

function get_isf(pf::DCPowerFlow, cont::Real, c::Integer)
    return pf.ϕ
end

"""
Calculation of line flow in a contingency case using IMML.
If DivideError, there are island(s) in the system.

Input:
    - Pl: Container for output
    - Pl0: Initial line power flow
    - ptdf: Power Transfer Distribution Factor matrix
    - B: Suseptance matrix
    - DA: Diagonal suseptance matrix times the connectivity matrix
    - X: The inverse admittance matrix
    - θ: Voltage angles
    - from_bus: From bus index
    - to_bus: To bus index
    - branch: Contingency branch index
"""
@views function calculate_line_flows!(
    Pl::AbstractVector{T},
    Pl0::AbstractVector{T},
    ptdf::AbstractMatrix{T},
    B::AbstractMatrix{T},
    DA::AbstractMatrix{T},
    X::AbstractMatrix{T},
    θ::AbstractVector{T},
    from_bus::Integer,
    to_bus::Integer,
    branch::Integer;
    atol::Real=1e-10
) where {T<:Real}
    change = DA[branch, to_bus] / B[from_bus, to_bus]
    # x = change * (X[:,from_bus] - X[:,to_bus])
    c⁻¹ = inv(B[from_bus, to_bus]) + change * (X[from_bus, from_bus] - X[from_bus, to_bus] - X[to_bus, from_bus] + X[to_bus, to_bus])
    if isapprox(c⁻¹, zero(T); atol=atol)
        return throw(DivideError())
    end
    delta = inv(c⁻¹) * (θ[from_bus] - θ[to_bus])
    @. Pl = Pl0 - (ptdf[:, from_bus] - ptdf[:, to_bus]) * change * delta
    Pl[branch] = 0.0
    return nothing
end
@views function calculate_line_flows!(
    Pl::AbstractVector{T},
    Pl0::AbstractVector{T},
    ptdf::AbstractMatrix{T},
    B::AbstractMatrix{T},
    DA::AbstractMatrix{T},
    X::AbstractMatrix{T},
    θ::AbstractVector{T},
    from_bus::Integer,
    to_bus::Integer,
    branch::Integer,
    ::Val{1}; # from_bus is not connected to the system
    atol::Real=1e-10
) where {T<:Real}
    change = DA[branch, to_bus] / B[from_bus, to_bus]
    c⁻¹ = inv(B[from_bus, to_bus]) + change * (- X[from_bus, to_bus] + X[to_bus, to_bus])
    if isapprox(c⁻¹, zero(T); atol=atol)
        return throw(DivideError())
    end
    delta = inv(c⁻¹) * (θ[from_bus] - θ[to_bus])
    @. Pl = Pl0 - (- ptdf[:, to_bus]) * change * delta
    Pl[branch] = 0.0
    return nothing
end
@views function calculate_line_flows!(
    Pl::AbstractVector{T},
    Pl0::AbstractVector{T},
    ptdf::AbstractMatrix{T},
    B::AbstractMatrix{T},
    DA::AbstractMatrix{T},
    X::AbstractMatrix{T},
    θ::AbstractVector{T},
    from_bus::Integer,
    to_bus::Integer,
    branch::Integer,
    ::Val{2}; # to_bus is not connected to the system
    atol::Real=1e-10
) where {T<:Real}
    change = DA[branch, to_bus] / B[from_bus, to_bus]
    c⁻¹ = inv(B[from_bus, to_bus]) + change * (X[from_bus, from_bus] - X[to_bus, from_bus])
    if isapprox(c⁻¹, zero(T); atol=atol)
        return throw(DivideError())
    end
    delta = inv(c⁻¹) * (θ[from_bus] - θ[to_bus])
    @. Pl = Pl0 - (ptdf[:, from_bus]) * change * delta
    Pl[branch] = 0.0
    return nothing
end

function calculate_line_flows!(
    Pl::AbstractVector{<:Real},
    pf::DCPowerFlow,
    cont::Tuple{Integer,Integer},
    branch::Integer;
    Pᵢ::AbstractVector{<:Real} = Float64[],
    nodes::AbstractVector{<:Integer} = Int[],
    branches::AbstractVector{<:Integer} = Int[]
)
    if !isempty(Pᵢ)
        θ = run_pf!(pf.vn_tmp, pf.K, Pᵢ, pf.slack)
        F = LinearAlgebra.mul!(pf.vb_tmp, pf.DA, pf.vn_tmp)
    else
        θ = pf.θ
        F = pf.F
    end
    if !isempty(nodes)
        island = not_insorted_nodes(cont[1], cont[2], nodes)
        calculate_line_flows!(Pl, F, pf.ϕ, pf.B, pf.DA, pf.X, θ, cont[1], cont[2], branch, Val(island))
        zero_not_in_array!(Pl, branches)
    else
        calculate_line_flows!(Pl, F, pf.ϕ, pf.B, pf.DA, pf.X, θ, cont[1], cont[2], branch)
    end
    return nothing
end
function calculate_line_flows(
    pf::DCPowerFlow,
    cont::Tuple{Integer,Integer},
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
    (x_l / x_k) * (X[p, m] - X[q, m] - X[p, n] + X[q, n]) /
    (x_l - (X[m, m] + X[n, n] - 2 * X[m, n]))
@views get_lodf(x_l::Real, m::Integer, n::Integer, x_k::AbstractVector{<:Real}, A::AbstractMatrix, X::AbstractMatrix) =
    (x_l ./ x_k) .* A * (X[:, m] - X[:, n]) ./
    (x_l - (X[m, m] + X[n, n] - 2 * X[m, n]))

function get_lodf(branch_l::Branch, branch_k::Branch, X::AbstractMatrix, idx::Dict{<:Any,<:Int})
    (m, n) = get_bus_idx(branch_l, idx)
    (p, q) = get_bus_idx(branch_k, idx)
    return get_lodf(get_x(branch_l), m, n, get_x(branch_k), p, q, X)
end
function get_lodf(branch_l::Branch, branches::AbstractVector{<:Branch}, A::AbstractMatrix, X::AbstractMatrix, idx::Dict{<:Any,<:Int})
    (m, n) = get_bus_idx(branch_l, idx)
    return get_lodf(get_x(branch_l), m, n, get_x.(branches), A, X)
end
function get_lodf(from_bus, to_bus, x::AbstractVector{<:Real}, A::AbstractMatrix, X::AbstractMatrix)
    mx = reshape(reduce(vcat, get_lodf.(x, from_bus, to_bus, [x], [A], [X])), (length(x), length(x)))
    return mx - LinearAlgebra.Diagonal(mx) - LinearAlgebra.I
end
