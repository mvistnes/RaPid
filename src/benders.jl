""" Benders type """
mutable struct Benders{TR<:Real}
    obj::JuMP.AffExpr
    Pᵢ::Vector{TR}
    Pg::Vector{TR}
    Pd::Vector{TR}
end

""" Constructor for Benders type """
function benders(opf::OPFsystem, mod::Model)
    obj = objective_function(mod)
    # Pᵢ = get_net_Pᵢ(opf)
    Pᵢ = get_value(mod, :p0)
    Pg = get_controllable(opf, mod)
    Pd = get_Pd(opf, mod) # Fixed injection
    @assert isapprox(Pg, (Pᵢ - Pd); atol=1e-14) string(Pg - (Pᵢ - Pd))
    return Benders{Float64}(obj, Pᵢ, Pg, Pd)
end

""" 
Solve the optimization model using Benders decomposition.

Function creates extra production constraints on generators in the system
based on the power transfer distribution factors of the generators in the
system. 
"""
function run_benders(
    type::OPF, 
    system::System, 
    optimizer, 
    voll::Vector{<:Number}, 
    prob::Vector{Float64}, 
    contingencies::Vector{<:Component};
    dist_slack=Float64[],
    time_limit_sec::Int64=600, 
    ramp_minutes=10, 
    ramp_mult=10, 
    max_shed=1.0, 
    max_curtail=1.0, 
    atol=1e-6,
    short_term_multi::Real=1.5, 
    long_term_multi::Real=1.2, 
    max_itr=10,
    p_failure=0.0,
    silent=true,
    debug=false
)
    case = Case(opf_base(OPF(true, false, false, false, false), system, optimizer, voll=voll, contingencies=contingencies, prob=prob,
        dist_slack=dist_slack, time_limit_sec=time_limit_sec, ramp_minutes=ramp_minutes, ramp_mult=ramp_mult, max_shed=max_shed, max_curtail=max_curtail,
        short_term_multi=short_term_multi, long_term_multi=long_term_multi, p_failure=p_failure, silent=silent, debug=debug)...)

    total_solve_time = constrain_branches!(case.model, case.pf, case.oplim, 0.0)
    if !type.P & !type.C1 & !type.C2 & !type.C2F
        return case, total_solve_time
    end
    return run_benders!(type, case, atol, max_itr)
end

run_benders!(type::OPF, case::Case, atol=1e-6, max_itr=max(length(case.opf.contingencies), 5)) =
    run_benders!(type, case.model, case.opf, case.pf, case.oplim, case.Pc, case.Pcc, case.Pccx, atol, max_itr)

"""
Solve the optimization model using Benders decomposition.
"""
function run_benders!(
    type::OPF, 
    mod::Model, 
    opf::OPFsystem, 
    pf::DCPowerFlow, 
    oplim::Oplimits, 
    Pc::Dict{<:Integer, ExprC}, 
    Pcc::Dict{<:Integer, ExprCC}, 
    Pccx::Dict{<:Integer, ExprCCX},
    atol=1e-6,
    max_itr=10
)
    assert(type)
    @assert !type.C1 || isempty(Pc)
    @assert !type.C2 || isempty(Pcc)
    @assert !type.C2F || isempty(Pccx)
    
    total_solve_time = solve_time(mod)
    !has_values(mod) && return Case(mod, opf, pf, oplim, Pc, Pcc, Pccx), total_solve_time
    @debug "lower_bound = $(objective_value(mod))"

    # Set variables
    !isempty(oplim.dist_slack) && set_dist_slack!(pf.ϕ, opf.mgx, oplim.dist_slack)
    bd = benders(opf, mod)

    overloads = zeros(length(opf.contingencies))
    # PI = zeros(length(opf.contingencies))
    pg = get_value(mod, :pg0)
    for (i, c_obj) in enumerate(opf.contingencies)
        cont = typesort_component(c_obj, opf)
        # cont  = get_bus_idx(opf.contingencies[c], opf.idx)
        if is_islanded(pf, cont[2], cont[1])
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont[2], cont[1], pf.slack)
            type.C1 && init_P!(Pc, opf, oplim, mod, bd.obj, islands, island, i)
            type.C2 && init_P!(Pcc, opf, oplim, mod, bd.obj, islands, island, i)
            type.C2F && init_P!(Pccx, opf, oplim, mod, bd.obj, islands, island, i)
            @debug "Island: Contingency $(string(typeof(c_obj))) $(get_name(c_obj))"
            if isempty(islands)
                fill!(flow, 0.0)
            else
                calculate_line_flows!(pf.vb_tmp, pf, cont[2], cont[1], nodes=islands[island], branches=island_b)
                flow = pf.vb_tmp
            end
        else
            if typeof(c_obj) <: ACBranch || (typeof(c_obj) <: Pair && first(c_obj) == "branch")
                flow = calculate_contingency_line_flows(mod, pf, bd.Pᵢ, cont, i, c_obj, Int[], Int[])
            else
                flow = calculate_contingency_line_flows(mod, pf, bd.Pᵢ, cont, i, c_obj, Int[], Int[], pg[cont[1]])
            end
        end

        # PI[i] = maximum((flow ./ (oplim.branch_rating * oplim.short_term_multi)).^2)

        # Calculate the power flow with the new outage and find if there are any overloads
        overload = filter_overload(flow, oplim.branch_rating * oplim.short_term_multi)

        if !isempty(overload)
            overloads[i] = maximum(x -> abs(x[2]), overload)
        end
    end

    total_solve_time = update_model!(mod, pf, oplim, bd, total_solve_time)
    !has_values(mod) && return Case(mod, opf, pf, oplim, Pc, Pcc, Pccx), total_solve_time

    ΔPc = zeros(length(bd.Pg))
    ΔPcc = zeros(length(bd.Pg))
    ΔPccx = zeros(length(bd.Pg))
    pre = 0
    corr1 = 0
    corr2 = 0
    corr2f = 0

    brst = oplim.branch_rating * oplim.short_term_multi
    brlt = oplim.branch_rating * oplim.long_term_multi

    olc = Vector{Tuple{Int,Float64}}()
    olcc = Vector{Tuple{Int,Float64}}()
    olccx = Vector{Tuple{Int,Float64}}()
    islands = Vector{Vector{Int}}()
    island = 0
    island_b = Int[]

    permutation = sortperm(overloads, rev=true)
    # permutation = sortperm(PI, rev=true)
    cut_added = 1
    for iterations in 1:max_itr
        if cut_added == 0 # loops until no new cuts are added for the contingencies
            # @printf "\nEND: Total solve time %.4f.\n" total_solve_time
            print_cuts(type, pre, corr1, corr2, corr2f)
            return Case(mod, opf, pf, oplim, Pc, Pcc, Pccx), total_solve_time
        end
        cut_added = 0
        @info "\n------------------\nIteration $iterations"

        for i in permutation
            c_obj = opf.contingencies[i]
            cont = typesort_component(c_obj, opf)
            if is_islanded(pf, cont[2], cont[1])
                islands, island, island_b = handle_islands(pf.B, pf.DA, cont[2], cont[1], pf.slack)
                inodes = islands[island]
            else
                empty!(islands)
                island = 0
                empty!(island_b)
                inodes = Int[]
            end

            if type.P
                olc = calculate_contingency_overload(brst, mod, pf, bd, cont, i, c_obj, inodes, island_b, atol)
            end
            if type.C1
                olc = calculate_contingency_overload!(ΔPc, brst, Pc, opf, mod, pf, bd, cont, i, c_obj, inodes, island_b, atol)
            end
            if type.C2
                olcc = calculate_contingency_overload!(ΔPcc, brlt, Pcc, opf, mod, pf, bd, cont, i, c_obj, inodes, island_b, atol)
            end
            if type.C2F
                olccx = calculate_contingency_overload!(ΔPccx, brlt, Pccx, opf, mod, pf, bd, cont, i, c_obj, inodes, island_b, atol)
            end
            if !isempty(olc) || !isempty(olcc) || !isempty(olccx) # ptdf calculation is more computational expensive than line flow
                if is_islanded(pf, cont[2], cont[1])
                    ptdf = get_isf(pf, cont[2], cont[1], islands, island, island_b)
                else
                    ptdf = get_isf(pf, cont[2], cont[1])
                    set_tol_zero!(ptdf)
                end

                type.C1 && get(Pc, i, 0) == 0 && fill!(ΔPc, 0.0)
                type.C2 && get(Pcc, i, 0) == 0 && fill!(ΔPcc, 0.0)
                type.C2F && get(Pccx, i, 0) == 0 && fill!(ΔPccx, 0.0)
            end

            # Cannot change the model before all data is exctracted!
            if !isempty(olc)
                if type.P
                    cut_added, pre = add_cut(opf, mod, bd, ptdf, olc, c_obj, cut_added, atol, pre)
                elseif type.C1
                    cut_added, corr1 = add_cut(Pc, opf, oplim, mod, bd, ΔPc, ptdf, olc, islands, island, c_obj, i, cut_added, atol, corr1)
                end
            end
            if !isempty(olcc)
                cut_added, corr2 = add_cut(Pcc, opf, oplim, mod, bd, ΔPcc, ptdf, olcc, islands, island, c_obj, i, cut_added, atol, corr2)
            end
            if !isempty(olccx)
                cut_added, corr2f = add_cut(Pccx, opf, oplim, mod, bd, ΔPccx, ptdf, olccx, islands, island, c_obj, i, cut_added, atol, corr2f)
            end
            if cut_added > 1
                total_solve_time = update_model!(mod, pf, oplim, bd, total_solve_time)
                cut_added = 1
            end
            !has_values(mod) && return Case(mod, opf, pf, oplim, Pc, Pcc, Pccx), total_solve_time
        end

    end
    @warn "Reached $(max_itr) iterations without a stable solution."
    return Case(mod, opf, pf, oplim, Pc, Pcc, Pccx), total_solve_time
end

""" Solve model and update the power flow object """
function update_model!(mod::Model, pf::DCPowerFlow, oplim::Oplimits, bd::Benders, total_solve_time::Real)
    # set_warm_start!(mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
    set_objective_function(mod, bd.obj)
    total_solve_time = constrain_branches!(mod, pf, oplim, total_solve_time)
    bd.Pᵢ = get_value(mod, :p0)
    @. bd.Pg = bd.Pᵢ - bd.Pd
    return total_solve_time
end

"""
    Calculate the contingency power flow with corrective actions.
    Assummes that island-contingencies have active variables from the pre-procedure.
"""
function calculate_contingency_line_flows!(ΔP::Vector{<:Real}, P::Dict{<:Integer, T}, opf::OPFsystem, 
    mod::Model, pf::DCPowerFlow, Pᵢ::AbstractVector{<:Real}, cont::Tuple{Real, Tuple{Real,Real}}, 
    c::Integer, c_obj::Component, island::Vector, island_b::Vector{<:Integer}
) where {T}
    if get(P, c, 0) == 0
        return calculate_contingency_line_flows(mod, pf, Pᵢ, cont, c, c_obj, island, island_b)
    end
    ΔP = get_ΔP!(ΔP, mod, opf, P, c)
    if iszero(ΔP)
        calculate_line_flows!(pf.vb_tmp, pf, cont[2], cont[1], nodes=island, branches=island_b)
    else
        calculate_line_flows!(pf.vb_tmp, pf, cont[2], cont[1], Pᵢ=(Pᵢ .+ ΔP), nodes=island, branches=island_b)
    end
    return pf.vb_tmp
end
function calculate_contingency_line_flows!(ΔP::Vector{<:Real}, P::Dict{<:Integer, T}, opf::OPFsystem, 
    mod::Model, pf::DCPowerFlow, Pᵢ::AbstractVector{<:Real}, cont::Tuple{Real,Real},
    c::Integer, c_obj::Component, island::Vector, island_b::Vector{<:Integer}
) where {T}
    if get(P, c, 0) == 0
        return calculate_contingency_line_flows(mod, pf, Pᵢ, cont, c, c_obj, island, island_b)
    end
    ΔP = get_ΔP!(ΔP, mod, opf, P, c)
    ΔP[cont[2]] -= get_ΔP(mod, cont, c_obj, P, c)
    LinearAlgebra.mul!(pf.vb_tmp, pf.ϕ, (Pᵢ .+ ΔP))
    return pf.vb_tmp
end

"""
    Calculate the contingency power flow without corrective actions.
"""
function calculate_contingency_line_flows(
    mod::Model, pf::DCPowerFlow, P::Vector{<:Real}, cont::Tuple{Real, Tuple{Real,Real}}, c::Integer, c_obj::Component,
    island::Vector, island_b::Vector{<:Integer}, val=0.0
)
    calculate_line_flows!(pf.vb_tmp, pf, cont[2], cont[1], nodes=island, branches=island_b)
    return pf.vb_tmp
end
function calculate_contingency_line_flows(
    mod::Model, pf::DCPowerFlow, P::Vector{<:Real}, cont::Tuple{Real,Real}, c::Integer, c_obj::StaticInjection,
    island::Vector, island_b::Vector{<:Integer}, val=get_ΔP(mod, cont, c_obj)
)
    P[cont[2]] -= val
    LinearAlgebra.mul!(pf.vb_tmp, pf.ϕ, P)
    P[cont[2]] += val
    return pf.vb_tmp
end

"""
    Calculate the contingency overflow with corrective actions.
"""
function calculate_contingency_overload!(ΔP::Vector{<:Real}, branch_rating::Vector{<:Real},
    Pc::Dict, opf::OPFsystem, mod::Model, pf::DCPowerFlow, bd::Benders, 
    cont::Union{Tuple{Real, Tuple{Real,Real}}, Tuple{Real,Real}}, c::Integer, c_obj::Component, 
    island::Vector, island_b::Vector{<:Integer}, atol::Real=1e-6
)
    flow = calculate_contingency_line_flows!(ΔP, Pc, opf, mod, pf, bd.Pᵢ, cont, c, c_obj, island, island_b)
    return filter_overload(flow, branch_rating, atol)
end

"""
    Calculate the contingency overflow without corrective actions.
"""
function calculate_contingency_overload(branch_rating::Vector{<:Real}, mod::Model, 
    pf::DCPowerFlow, bd::Benders, cont::Union{Tuple{Real, Tuple{Real,Real}}, Tuple{Real,Real}}, c::Integer, c_obj::Component,
    island::Vector, island_b::Vector{<:Integer}, atol::Real=1e-6
)
    flow = calculate_contingency_line_flows(mod, pf, bd.Pᵢ, cont, c, c_obj, island, island_b)
    return filter_overload(flow, branch_rating, atol)
end

get_ΔP(mod::Model, cont::Tuple{Real,Real}, ::Generator) = 
    value(mod[:pg0][cont[1]])
get_ΔP(mod::Model, cont::Tuple{Real,Real}, ::Generator, Pc::Dict{<:Integer, ExprC}, c::Integer) = 
    value(mod[:pg0][cont[1]]) + value(Pc[c].pgu[cont[1]]) - value(Pc[c].pgd[cont[1]])
get_ΔP(mod::Model, cont::Tuple{Real,Real}, ::Generator, Pcc::Dict{<:Integer, ExprCC}, c::Integer) = 
    value(mod[:pg0][cont[1]]) + value(Pcc[c].pgu[cont[1]]) - value(Pcc[c].pgd[cont[1]])
get_ΔP(mod::Model, cont::Tuple{Real,Real}, ::Generator, Pccx::Dict{<:Integer, ExprCCX}, c::Integer) = 
    value(mod[:pg0][cont[1]]) - value(Pccx[c].pgdx[cont[1]])

get_ΔP(mod::Model, cont::Tuple{Real,Real}, ::StaticLoad) = 
    value(mod[:pd0][cont[1]]) - value(mod[:ls0][cont[1]])
get_ΔP(mod::Model, cont::Tuple{Real,Real}, ::StaticLoad, Pc::Dict{<:Integer, ExprC}, c::Integer) = 
    value(mod[:pd0][cont[1]]) - value(Pc[c].lsc[cont[1]])
get_ΔP(mod::Model, cont::Tuple{Real,Real}, ::StaticLoad, Pcc::Dict{<:Integer, ExprCC}, c::Integer) = 
    value(mod[:pd0][cont[1]]) - value(Pcc[c].lscc[cont[1]])
get_ΔP(mod::Model, cont::Tuple{Real,Real}, ::StaticLoad, Pccx::Dict{<:Integer, ExprCCX}, c::Integer) = 
    value(mod[:pd0][cont[1]]) - value(Pccx[c].lsccx[cont[1]])

""" Return the contingency power injection change at each node changing ΔP. """
function get_ΔP!(ΔP::Vector{T}, mod::Model, opf::OPFsystem, 
    P::Dict{<:Integer, <:ContExpr}, c::Integer
) where {T<:Real}
    x = get(P, c, 0)
    if x != 0
        copy!(ΔP, get_ΔP(x, opf, mod))
    else
        fill!(ΔP, zero(T))
    end
    return ΔP
end

""" Return the short term post-contingency power injection change at each node. """
@inline function get_ΔP(pc::ExprC, opf::OPFsystem, mod::Model)
    return opf.mgx' * (get_value(mod, pc.pgu) - get_value(mod, pc.pgd)) - 
        opf.mrx' * get_value(mod, pc.prc) + opf.mdx' * get_value(mod, pc.lsc)
end

""" Return the long term post-contingency power injection change at each node. """
@inline function get_ΔP(pcc::ExprCC, opf::OPFsystem, mod::Model)
    return opf.mgx' * (get_value(mod, pcc.pgu) - get_value(mod, pcc.pgd)) - 
        opf.mrx' * get_value(mod, pcc.prcc) + opf.mdx' * get_value(mod, pcc.lscc) +
        opf.mdcx' * get_value(mod, pcc.pfdccc)
end

""" Return the long term post-contingency with failure power injection change at each node. """
@inline function get_ΔP(pccx::ExprCCX, opf::OPFsystem, mod::Model)
    return -opf.mgx' * get_value(mod, pccx.pgdx) + opf.mdx' * get_value(mod, pccx.lsccx)
end

function add_overload_expr!(mod::Model, expr::JuMP.AbstractJuMPScalar, ol::Real, type::String, id::Integer,
    c_obj::Component, opf::OPFsystem, i::Integer, atol::Real
)
    @info @sprintf "%s %d: Contingency %s %s; overload on %s of %.6f" type id string(typeof(c_obj)) c_obj.name opf.branches[i].name ol
    @debug "Cut added: $(sprint_expr(expr,atol))\n"
    if ol < 0
        JuMP.@constraint(mod, expr <= ol)
    else
        JuMP.@constraint(mod, expr >= ol)
    end
end

function add_cut(opf::OPFsystem, mod::Model, bd::Benders, ptdf::AbstractMatrix{<:Real}, overloads::Vector{<:Tuple{Integer,Real}},
    c_obj::Component, cut_added::Integer, atol::Real, id::Integer
)
    for (i, ol) in overloads
        expr = AffExpr() + sum(ptdf[i, :] .* bd.Pᵢ)
        for j in axes(ptdf,2)
            add_to_expression!(expr, -ptdf[i, j], mod[:inj_p0][j])
        end

        id += 1
        add_overload_expr!(mod, expr, ol, "Pre", id, c_obj, opf, i, atol)
        cut_added = 2
    end
    return cut_added, id
end

""" Finding and adding Benders cuts to mitigate overloads from one contingency """
function add_cut(P::Dict{<:Integer, <:ContExpr}, opf::OPFsystem, oplim::Oplimits, mod::Model, bd::Benders, ΔP::Vector{<:Real},
    ptdf::AbstractMatrix{<:Real}, overloads::Vector{<:Tuple{Integer,Real}}, islands::Vector, island::Integer,
    c_obj::Component, c::Integer, cut_added::Integer, atol::Real, id::Integer
)
    p = get(P, c, 0)
    if p == 0  # If the contingency is not run before, a set of corrective variables is added
        init_P!(P, opf, oplim, mod, bd.obj, islands, island, c)
        p = get(P, c, 0)
    end
    for (i, ol) in overloads
        expr = make_cut(p, opf, mod)
        add_to_expression!.(expr, bd.Pᵢ + ΔP)
        expr2 = JuMP.@expression(mod, sum((ptdf[i, :] .* expr)))

        id += 1
        add_overload_expr!(mod, expr2, ol, get_name(p), id, c_obj, opf, i, atol)
        cut_added = 2
    end
    return cut_added, id
end

function make_cut(pc::ExprC, opf::OPFsystem, mod::Model)
    expr = JuMP.@expression(mod, -opf.mgx' * (mod[:pg0] + pc.pgu - pc.pgd))
    add_to_expression!.(expr, -opf.mrx' * (mod[:pr] - pc.prc))
    add_to_expression!.(expr, -opf.mdcx' * mod[:pfdc0])
    add_to_expression!.(expr, -opf.mdx' * (pc.lsc - mod[:pd]))
    return expr
end

function make_cut(pcc::ExprCC, opf::OPFsystem, mod::Model)
    expr = JuMP.@expression(mod, -opf.mgx' * (mod[:pg0] + pcc.pgu - pcc.pgd))
    add_to_expression!.(expr, -opf.mrx' * (mod[:pr] - pcc.prcc))
    add_to_expression!.(expr, -opf.mdcx' * pcc.pfdccc)
    add_to_expression!.(expr, -opf.mdx' * (pcc.lscc - mod[:pd]))
    return expr
end

function make_cut(pccx::ExprCCX, opf::OPFsystem, mod::Model)
    expr = JuMP.@expression(mod, -opf.mgx' * (mod[:pg0] - pccx.pgdx))
    add_to_expression!.(expr, -opf.mrx' * (mod[:pr] - pccx.prccx))
    add_to_expression!.(expr, -opf.mdcx' * mod[:pfdc0])
    add_to_expression!.(expr, -opf.mdx' * (pccx.lsccx - mod[:pd]))
    return expr
end