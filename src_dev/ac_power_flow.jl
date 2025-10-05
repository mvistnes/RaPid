mutable struct ACPowerFlow{T1<:Real,T2<:Integer} <: PowerFlow
    DA::SparseArrays.SparseMatrixCSC{T1,T2} # Diagonal admittance matrix times the connectivity matrix
    B::SparseArrays.SparseMatrixCSC{T1,T2} # The admittance matrix
    K::KLU.KLUFactorization{T1, T2} # A factorization of the admittance matrix
    X::Matrix{T1} # Inverse of the admittance matrix
    ϕ::Matrix{T1} # PTDF matrix
    θ::Vector{T1} # Bus voltage angles
    F::Vector{T1} # Line power flow
    slack::T2 # Reference bus

    sp_tmp::SparseArrays.SparseMatrixCSC{T1,T2} # branch x bus sparse matrix
    K_tmp::KLU.KLUFactorization{T1, T2} # A factorization of the admittance matrix
    m_tmp::Matrix{T1} # bus x bus matrix
    vn_tmp::Vector{T1} # bus vector
    vb_tmp::Vector{T1} # branch vector
end

function ACPowerFlow(branches::AbstractVector{<:Tuple{T2,T2}}, susceptance::AbstractVector{T1}, numnodes::Integer, slack::Integer) where {T1<:Real,T2<:Integer}
    A = calc_A(branches)
    DA = calc_DA(A, susceptance)
    B = calc_B(A, DA)
    K = get_klu(B, slack)
    ϕ = get_ptdf(K, DA, slack)
    X = calc_X(K, slack)
    return ACPowerFlow{T1,T2}(DA, B, K, X, ϕ, zeros(T1, numnodes), zeros(T1, length(branches)), slack,
        zeros(T1, size(DA)), get_klu(B, slack), zeros(T1, size(X)), zeros(T1, numnodes), zeros(T1, length(branches)))
end
ACPowerFlow(nodes::AbstractVector{<:Bus}, branches::AbstractVector{<:Branch}, idx::Dict{<:Int,<:Int}) =
    ACPowerFlow(get_bus_idx.(branches, [idx]), PowerSystems.get_series_susceptance.(branches), length(nodes), find_slack(nodes)[1])

function ACPowerFlow(sys::System)
    nodes = sort_components!(get_nodes(sys))
    branches = sort_components!(get_branches(sys))
    idx = get_nodes_idx(nodes)
    return ACPowerFlow(nodes, branches, idx)
end

function ACPowerFlow(nodes::AbstractVector{<:Bus}, branches::AbstractVector{<:Branch}, idx::Dict{<:Int,<:Int}, Pᵢ::AbstractVector{<:Real})
    pf = ACPowerFlow(nodes, branches, idx)
    calc_θ!(pf, Pᵢ)
    calc_Pline!(pf)
    return pf
end

get_slack(pf::ACPowerFlow) = pf.slack
get_DA(pf::ACPowerFlow) = pf.DA
get_B(pf::ACPowerFlow) = pf.B
get_K(pf::ACPowerFlow) = pf.K
get_X(pf::ACPowerFlow) = pf.X
get_ϕ(pf::ACPowerFlow) = pf.ϕ
get_θ(pf::ACPowerFlow) = pf.θ

set_θ!(pf::ACPowerFlow, model::Model) = copy!(pf.θ, get_sorted_angles(model))

function _add_Ybus!(ybus::SparseArrays.SparseMatrixCSC, br::ACBranch, idx::Dict{<:Int,<:Integer})
    (fbus, tbus) = get_bus_idx(br, idx)
    y = 1 / (get_r(br) + get_x(br) * im)
    b = get_b(br)
    ybus[fbus,fbus] += (y + b.from * im)
    ybus[tbus,tbus] += (y + b.to * im)
    ybus[fbus,tbus] -= y
    ybus[tbus,fbus] -= y
    return
end

function _add_Ybus!(ybus::SparseArrays.SparseMatrixCSC, br::Transformer2W, idx::Dict{<:Int,<:Integer})
    (fbus, tbus) = get_bus_idx(br, idx)
    y = 1 / (get_r(br) + get_x(br) * im)
    b = get_primary_shunt(br)
    ybus[fbus,fbus] += (y + b * im)
    ybus[tbus,tbus] += y
    ybus[fbus,tbus] -= y
    ybus[tbus,fbus] -= y
    return
end

function _add_Ybus!(ybus::SparseArrays.SparseMatrixCSC, br::TapTransformer, idx::Dict{<:Int,<:Integer})
    (fbus, tbus) = get_bus_idx(br, idx)
    y = 1 / (get_r(br) + get_x(br) * im)
    tap = 1 / get_tap(br)
    b = get_primary_shunt(br)
    ybus[fbus,fbus] += ((y * tap^2) + b * im)
    ybus[tbus,tbus] += y
    ybus[fbus,tbus] -= (y * tap)
    ybus[tbus,fbus] -= (y * tap)
    return
end

function _add_Ybus!(ybus::SparseArrays.SparseMatrixCSC, br::PhaseShiftingTransformer, idx::Dict{<:Int,<:Integer})
    (fbus, tbus) = get_bus_idx(br, idx)
    y = 1 / (get_r(br) + get_x(br) * im)
    tap = get_tap(br) * exp(get_α(br) * im)
    c_tap = get_tap(br) * exp(get_α(br) * -im)
    b = get_primary_shunt(br)
    ybus[fbus,fbus] += ((y * tap^2) + b * im)
    ybus[tbus,tbus] += y
    ybus[fbus,tbus] -= (y / c_tap)
    ybus[tbus,fbus] -= (y / tap)
    return
end

function _add_Ybus!(ybus::SparseArrays.SparseMatrixCSC, br::FixedAdmittance, idx::Dict{<:Int,<:Integer})
    bus = idx[br.bus.number]
    ybus[bus,bus] += get_Y(br)
    return
end

# """ Builds an admittance matrix with the line series reactance of the lines. """
function calc_Ybus(branches::AbstractVector{<:Branch}, idx::Dict{<:Int,<:Integer}, fixed::AbstractVector{FixedAdmittance})
    numnodes = length(idx)
    ybus = SparseArrays.spzeros(ComplexF64, numnodes, numnodes)
    for branch in branches
        _add_Ybus!(ybus, branch, idx)
    end
    for fx in fixed
        _add_Ybus!(ybus, fx, idx)
    end
    return ybus
end

function calc_B_mark(branches::AbstractVector{<:Branch}, idx::Dict{<:Int,<:Integer}, fixed::AbstractVector{FixedAdmittance})
    return calc_B(branches, idx, fixed)
end

function calc_B_doublemark(ybus::SparseArrays.SparseMatrixCSC{ComplexF64, Int})
    return -imag(ybus)
end
function calc_B_doublemark(branches::AbstractVector{<:Branch}, idx::Dict{<:Int,<:Integer}, fixed::AbstractVector{FixedAdmittance})
    return calc_B_doublemark(calc_Ybus(branches, idx, fixed))
end

#Intermediate calculations
Tij(Gij::Real, Bij::Real, thetaij::Real) = Gij*cos(thetaij) + Bij*sin(thetaij)
Uij(Gij::Real, Bij::Real, thetaij::Real) = Gij*sin(thetaij) - Bij*cos(thetaij)

"""Power equation at a bus"""
function injected_power(type::Function, ybus::SparseArrays.SparseMatrixCSC{ComplexF64, Int}, vang::AbstractVector{<:Real}, 
    vmag::AbstractVector{<:Real}, bus::Integer
)
    val = 0.0
    ycol = ybus[bus,:]
    for (i, y) in zip(ycol.nzind, ycol.nzval)
        val += vmag[i] * type(real(y), imag(y), vang[bus] - vang[i])
    end
    return val*vmag[bus]
end
injected_active_power(ybus::SparseArrays.SparseMatrixCSC{ComplexF64, Int}, vang::AbstractVector{<:Real}, 
    vmag::AbstractVector{<:Real}, bus::Integer
) = injected_power(Tij, ybus, vang, vmag, bus)
injected_reactive_power(ybus::SparseArrays.SparseMatrixCSC{ComplexF64, Int}, vang::AbstractVector{<:Real}, 
    vmag::AbstractVector{<:Real}, bus::Integer
) = injected_power(Uij, ybus, vang, vmag, bus)


function injected_power(ybus::SparseArrays.SparseMatrixCSC{ComplexF64, Int}, vang::AbstractVector{<:Real}, 
    vmag::AbstractVector{<:Real}
)
    return injected_active_power.([ybus], [vang], [vmag], axes(ybus, 1)), 
        injected_reactive_power.([ybus], [vang], [vmag], axes(ybus, 1))
end

function injected_power!(p::AbstractVector{<:Real}, q::AbstractVector{<:Real}, 
    p0::AbstractVector{<:Real}, q0::AbstractVector{<:Real}, 
    ybus::SparseArrays.SparseMatrixCSC{ComplexF64, Int}, vang::AbstractVector{<:Real}, 
    vmag::AbstractVector{<:Real}
)
    p .= p0 .- injected_active_power.([ybus], [vang], [vmag], axes(ybus, 1))
    q .= q0 .- injected_reactive_power.([ybus], [vang], [vmag], axes(ybus, 1))
    return p, q
end