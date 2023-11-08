# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

"""
    Container for OPF state types.
    Base: base case.
    P: Only preventive actions.
    C1: Preventive and corrective actions short term. 
    C2: Preventive and corrective actions long term.
    C2F: Preventive and corrective actions failures long term.

    P cannot be combined with C1.
    C2F cannot be used without C2.
    Use ex.:
        - OPF(true, false, true, true, true)
        - OPF(true, true, false, true, false)
"""
struct OPF 
    Base::Bool
    P::Bool
    C1::Bool
    C2::Bool
    C2F::Bool
end

"OPF with only case"
const Base_SCOPF = OPF(true, false, false, false, false)

"OPF with only preventive actions"
const P_SCOPF = OPF(true, true, false, false, false)

"OPF with preventive for one post-contingency state and corrective actions for one post-contingency state"
const PC_SCOPF = OPF(true, true, false, true, false)

"OPF with preventive for one post-contingency state and corrective actions for one post-contingency state with corrective failures"
const PCF_SCOPF = OPF(true, true, false, true, true)

"OPF with preventive and corrective actions for two post-contingency state"
const PC2_SCOPF = OPF(true, false, true, true, false)

"OPF with preventive and corrective actions for two post-contingency state with corrective failures"
const PC2F_SCOPF = OPF(true, false, true, true, true)

function assert(opf::OPF)
    @assert opf.P ⊼ opf.C1
    @assert !opf.C2F | (opf.C2 & opf.C2F)
    return
end

function +(x::OPF, y::OPF)
    return OPF(x.Base | y.Base, x.P | y.P, x.C1 | y.C1, x.C2 | y.C2, x.C2F | y.C2F)
end

function -(x::OPF, y::OPF)
    return OPF(x.Base > y.Base, x.P > y.P, x.C1 > y.C1, x.C2 > y.C2, x.C2F > y.C2F)
end


""" Holds the short term variables for contingencies """
mutable struct ExprC{T<:Vector{JuMP.VariableRef}}
    pc::T
    pgc::T
    prc::T
    lsc::T
end

""" Holds the long term variables for contingencies """
mutable struct ExprCC{T<:Vector{JuMP.VariableRef}}
    pcc::T
    pgu::T
    pgd::T
    pfdccc::T
    prcc::T
    lscc::T
end

""" Holds the long term variables for contingencies, no ramp up allowed """
mutable struct ExprCCX{T<:Vector{JuMP.VariableRef}}
    pccx::T
    pgdx::T
    prccx::T
    lsccx::T
end


""" OPF system type """
mutable struct OPFsystem{T<:Real}
    cost_ctrl_gen::Vector{<:Tuple{T, T}}
    cost_renewables::Vector{T}
    voll::Vector{T}
    prob::Vector{T}

    ctrl_generation::Vector{Generator}
    branches::Vector{ACBranch}
    dc_branches::Vector{TwoTerminalHVDCLine}
    nodes::Vector{ACBus}
    demands::Vector{StaticLoad}
    renewables::Vector{RenewableGen}
    contingencies::Vector{Component}
end

""" Constructor for OPFsystem """
function opfsystem(sys::System, voll::Vector{T}, contingencies::Vector{<:Component}=Component[], prob::Vector{T}=[];
    check=false
) where {T<:Real}

    ctrl_generation = sort_components!(get_ctrl_generation(sys))
    branches = sort_components!(get_branches(sys))
    dc_branches = sort_components!(get_dc_branches(sys))
    nodes = sort_components!(get_nodes(sys))
    demands = sort_components!(get_demands(sys))
    renewables = sort_components!(get_renewables(sys))

    cost_ctrl_gen = get_generator_cost.(ctrl_generation)::Vector{Tuple{Float64, Float64}}
    cost_renewables = Vector{Float64}() # Vector{Float64}([get_generator_cost(g)[2] for g in renewables])

    if check
        check_values(getindex.(cost_ctrl_gen, 1), ctrl_generation, "cost_1")
        check_values(getindex.(cost_ctrl_gen, 2), ctrl_generation, "cost_2")
        check_values(voll, demands, "voll")
        check_values.(ctrl_generation)
        check_values.(branches)
        check_values.(dc_branches)
    end

    A = calc_B(branches, nodes)
    islands = island_detection(A'A)
    if length(islands) > 1
        @warn "The system is separated into islands" islands
    end

    return OPFsystem{T}(cost_ctrl_gen, cost_renewables, voll, prob,
        ctrl_generation, branches, dc_branches, nodes, demands, renewables, contingencies)
end

""" Automatic constructor for OPFsystem where voll, prob, and contingencies are automatically computed. """
function opfsystem(sys::System)
    voll, prob, contingencies = setup(sys, 1, 4)
    fix_generation_cost!(sys)
    return opfsystem(sys, voll, contingencies, prob)
end

get_cost_ctrl_gen(opf::OPFsystem) = opf.cost_ctrl_gen
get_cost_renewables(opf::OPFsystem) = opf.cost_renewables
get_voll(opf::OPFsystem) = opf.voll
get_prob(opf::OPFsystem) = opf.prob

get_ctrl_generation(opf::OPFsystem) = opf.ctrl_generation
get_branches(opf::OPFsystem) = opf.branches
get_dc_branches(opf::OPFsystem) = opf.dc_branches
get_nodes(opf::OPFsystem) = opf.nodes
get_demands(opf::OPFsystem) = opf.demands
get_renewables(opf::OPFsystem) = opf.renewables
get_contingencies(opf::OPFsystem) = opf.contingencies
