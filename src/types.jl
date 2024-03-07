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

abstract type ContExpr end

""" Holds the short term variables for contingencies """
mutable struct ExprC{T<:Vector{JuMP.VariableRef}} <: ContExpr
    pc::T
    pgu::T
    pgd::T
    lsc::T
end
ExprC() = ExprC(JuMP.VariableRef[], JuMP.VariableRef[], JuMP.VariableRef[], JuMP.VariableRef[])
get_name(::ExprC) = "Corr1"

""" Holds the long term variables for contingencies """
mutable struct ExprCC{T<:Vector{JuMP.VariableRef}} <: ContExpr
    pcc::T
    pgu::T
    pgd::T
    pfdccc::T
    lscc::T
end
ExprCC() = ExprCC(JuMP.VariableRef[], JuMP.VariableRef[], JuMP.VariableRef[], JuMP.VariableRef[], JuMP.VariableRef[])
get_name(::ExprCC) = "Corr2"

""" Holds the long term variables for contingencies, no ramp up allowed """
mutable struct ExprCCX{T<:Vector{JuMP.VariableRef}} <: ContExpr
    pccx::T
    pgdx::T
    lsccx::T
end
ExprCCX() = ExprCCX(JuMP.VariableRef[], JuMP.VariableRef[], JuMP.VariableRef[])
get_name(::ExprCCX) = "Corr2x"
