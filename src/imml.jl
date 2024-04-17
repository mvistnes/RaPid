""" Sherman Morrison """
function calc_X(X::AbstractMatrix{<:Real}, A::AbstractMatrix{<:Integer}, br_x::Real, branch::Integer; atol::Real=1e-10)
    u = A[branch,:]
    v = -u
    mx = X - (br_x * X * u * v' * X) / (1 + br_x * v' * X * u)
    set_tol_zero!(mx, atol)
    return mx
end

""" Sherman Morrison Woodbury """
function calc_X(X::AbstractMatrix{T}, A::AbstractMatrix{<:Real}, D::AbstractMatrix{<:Real}, 
    branches::AbstractVector{<:Integer}; atol::Real=1e-10
) where {T<:Real}
    H = A[branches,:]
    G = D[branches, branches]
    XF = X * H'
    cinv = inv(G) - H * XF
    if isapprox(cinv, zero(T); atol=atol)
        return throw(DivideError())
    end
    mx = X + XF * inv(cinv) * XF'
    set_tol_zero!(mx, atol)
    return mx
end

""" Get the isf-matrix after a line outage using Sherman Morrison Woodbury. """
function calc_isf(
    pf::DCPowerFlow,
    cont::AbstractVector{<:Tuple{Integer,Integer}},
    branch::AbstractVector{<:Integer}
)
    ϕ = similar(pf.ϕ)
    if length(cont) < 2
        calc_isf!(ϕ, pf.mnn_tmp, pf.X, pf.B, pf.DA, cont[1], branch[1])
    else
        X = calc_X(pf.X, pf.A, pf.D, branch)
        calc_isf!(ϕ, pf.DA, X)
        ϕ[branch, :] .= 0.0
        set_tol_zero!(ϕ)
    end
    return ϕ
end

function calc_ptdf_vec(Xf::AbstractVector{T}, Xt::AbstractVector{T}, Xk::AbstractVector{T}, Xl::AbstractVector{T}, 
    x_c::Real, x_l::Real, fbus::Integer, tbus::Integer, kbus::Integer, lbus::Integer; atol::Real=1e-10
) where {T<:Real}
    cinv = 1 + x_c * (Xf[tbus] - Xf[fbus] + Xt[fbus] - Xt[tbus])
    if isapprox(cinv, zero(T); atol=atol)
        return throw(DivideError())
    end
    xft = -Xf[kbus] + Xt[kbus] + Xf[lbus] - Xt[lbus]
    vec = x_l * ((Xk .- Xl) .- x_c * xft / cinv * (Xf .- Xt))
    # set_tol_zero!(vec)
    # xftl = Xl .- x_c * (Xt .- Xf) * (Xf[lbus] - Xt[lbus]) / cinv
    # xftk = Xk .- x_c * (Xt .- Xf) * (Xf[kbus] - Xt[kbus]) / cinv
    # vec = x_l * (xftk .- xftl)
    return vec
end

function calc_ptdf_vec(B::AbstractMatrix{<:Real}, K::KLU.KLUFactorization{T,<:Integer}, slack::Integer, 
    fbus::Integer, tbus::Integer, kbus::Integer, lbus::Integer; atol::Real=1e-10
) where {T<:Real}
    n = size(B, 1)
    Xf = calc_X_vec!(Vector{T}(undef, n), K, fbus, slack)
    Xt = calc_X_vec!(Vector{T}(undef, n), K, tbus, slack)
    Xk = calc_X_vec!(Vector{T}(undef, n), K, kbus, slack)
    Xl = calc_X_vec!(Vector{T}(undef, n), K, lbus, slack)
    return calc_ptdf_vec(Xf, Xt, Xk, Xl, -B[fbus, tbus], -B[kbus, lbus], fbus, tbus, kbus, lbus, atol=atol)
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
@views function calc_changed_X!(
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
        throw(DivideError())
    end
    # delta = 1/c⁻¹ * x
    copy!(X, X₀)
    LinearAlgebra.mul!(X, x, x', -change * inv(c⁻¹), true) # mul!(C, A, B, α, β) -> C == $A B α + C β$
end

""" Get the isf-matrix after a line outage using IMML. 
    isf and X are containers for output and calculation and will be overwritten """
function calc_isf!(
    isf::AbstractMatrix{<:Real},
    X::AbstractMatrix{<:Real},
    X₀::AbstractMatrix{<:Real},
    B::AbstractMatrix{<:Real},
    DA::AbstractMatrix{<:Real},
    cont::Tuple{Integer,Integer},
    branch::Integer
)
    calc_changed_X!(X, X₀, B, DA, cont[1], cont[2], branch)
    calc_isf!(isf, DA, X)
    isf[branch, :] .= 0.0
    set_tol_zero!(isf)
    return isf
end

""" Get the isf-matrix after a line outage using IMML. 
     """
function calc_isf(
    pf::DCPowerFlow,
    cont::Tuple{Integer,Integer},
    branch::Integer
)
    calc_isf!(pf.mbn_tmp, pf.mnn_tmp, pf.X, pf.B, pf.DA, cont, branch)
    return pf.mbn_tmp
end

function calc_isf(pf::DCPowerFlow, cont::Real, c::Integer)
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
    Xf::AbstractVector{T},
    Xt::AbstractVector{T},
    θ::AbstractVector{T},
    from_bus::Integer,
    to_bus::Integer,
    branch::Integer;
    atol::Real=1e-10
) where {T<:Real}
    change = DA[branch, to_bus] / B[from_bus, to_bus]
    # x = change * (X[:,from_bus] - X[:,to_bus])
    c⁻¹ = inv(B[from_bus, to_bus]) + change * (Xf[from_bus] - Xf[to_bus] - Xt[from_bus] + Xt[to_bus])
    if isapprox(c⁻¹, zero(T); atol=atol)
        return throw(DivideError())
    end
    delta = inv(c⁻¹) * (θ[from_bus] - θ[to_bus])
    @. Pl = Pl0 - (ptdf[:, from_bus] - ptdf[:, to_bus]) * change * delta
    Pl[branch] = zero(T)
    return Pl
end
function calculate_line_flows!(
    Pl::AbstractVector{<:Real},
    pf::DCPowerFlow,
    cont::Tuple{Integer,Integer},
    branch::Integer,
    nodes::AbstractVector{<:Integer} = Int64[],
    branches::AbstractVector{<:Integer} = Int64[]
)
    if isempty(nodes) || size(pf.B, 1) == length(nodes)
        calculate_line_flows!(Pl, pf.F, pf.ϕ, pf.B, pf.DA, view(pf.X,:,cont[1]), view(pf.X,:,cont[2]), pf.θ, cont[1], cont[2], branch)
    elseif insorted(cont[1], nodes)
        calculate_line_flows!(Pl, pf.F, pf.ϕ, pf.B, pf.DA, view(pf.X,:,cont[1]), zeros(size(pf.X,1)), pf.θ, cont[1], cont[2], branch)
        zero_not_in_array!(Pl, branches)
    elseif insorted(cont[2], nodes)
        calculate_line_flows!(Pl, pf.F, pf.ϕ, pf.B, pf.DA, zeros(size(pf.X,1)), view(pf.X,:,cont[2]), pf.θ, cont[1], cont[2], branch)
        zero_not_in_array!(Pl, branches)
    else
        throw(DivideError())
    end
    return Pl
end

@views function calculate_line_flows!(
    Pl::AbstractVector{T},
    ϕ::AbstractMatrix{T},
    B::AbstractMatrix{T},
    DA::AbstractMatrix{T},
    Xf::AbstractVector{T},
    Xt::AbstractVector{T},
    Pᵢ::AbstractVector{T},
    from_bus::Integer,
    to_bus::Integer,
    branch::Integer;
    atol::Real=1e-10
) where {T<:Real}
    change = DA[branch, to_bus] / B[from_bus, to_bus]
    x = Xf - Xt
    c⁻¹ = inv(B[from_bus, to_bus]) + change * (x[from_bus] - x[to_bus])
    if isapprox(c⁻¹, zero(T); atol=atol)
        return throw(DivideError())
    end
    delta = change * inv(c⁻¹) * LinearAlgebra.dot(x, Pᵢ)
    LinearAlgebra.mul!(Pl, DA, x) # mul!(C, A, B) -> C == $A B$
    LinearAlgebra.mul!(Pl, ϕ, Pᵢ, true, -delta) # mul!(C, A, B, α, β) -> C == $A B α + C β$
    Pl[branch] = zero(T)
    return Pl
end
function calculate_line_flows!(
    Pl::AbstractVector{<:Real},
    pf::DCPowerFlow,
    cont::Tuple{Integer,Integer},
    branch::Integer,
    Pᵢ::AbstractVector{<:Real}
)
    calculate_line_flows!(Pl, pf.ϕ, pf.B, pf.DA, view(pf.X,:,cont[1]), view(pf.X,:,cont[2]), Pᵢ, cont[1], cont[2], branch) 
    return Pl
end

function calculate_line_flows!(
    Pl::AbstractVector{T},
    ϕ::AbstractMatrix{T},
    B::AbstractMatrix{T},
    A::AbstractMatrix{T},
    DA::AbstractMatrix{T},
    X::AbstractMatrix{T},
    Pᵢ::AbstractVector{T},
    fnodes::AbstractVector{<:Integer},
    tnodes::AbstractVector{<:Integer},
    branches::AbstractVector{<:Integer};
    atol::Real=1e-10
) where {T<:Real}
    change = LinearAlgebra.Diagonal(getindex.([DA],branches, tnodes) .* inv.(getindex.([B], fnodes, tnodes)))
    Ab = A[branches,:]
    XF = X * Ab' * change
    c⁻¹ = LinearAlgebra.Diagonal(inv.(getindex.([DA], branches, fnodes))) - Ab * XF
    if isapprox(LinearAlgebra.det(c⁻¹), zero(T); atol=atol)
        return throw(DivideError())
    end
    LinearAlgebra.mul!(Pl, ϕ, (Pᵢ .+ Ab'inv(c⁻¹) * XF'Pᵢ))
    Pl[branches] .= zero(T)
    return Pl
end
function calculate_line_flows!(
    Pl::AbstractVector{<:Real},
    pf::DCPowerFlow,
    conts::AbstractVector{<:Tuple{Integer,Integer}},
    branches::AbstractVector{<:Integer},
    Pᵢ::AbstractVector{<:Real}
)
    calculate_line_flows!(Pl, pf.ϕ, pf.B, pf.A, pf.DA, pf.X, Pᵢ, first.(conts), last.(conts), branches) 
    return Pl
end

function calculate_line_flows(
    pf::DCPowerFlow,
    cont,
    branch,
    Pᵢ::AbstractVector{<:Real},
)
    Pl = similar(pf.F)
    calculate_line_flows!(Pl, pf, cont, branch, Pᵢ)
    return Pl
end