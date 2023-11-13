""" Benders type """
mutable struct Benders{TR<:Real,TI<:Integer}
    idx::Dict{TI,TI}
    list::Vector{CTypes{TI}}

    obj::AbstractJuMPScalar

    Pᵢ::Vector{TR}
    Pg::Vector{TR}
    Pd::Vector{TR}
end

""" Constructor for Benders type """
function benders(opf::OPFsystem, mod::Model)
    idx = get_nodes_idx(opf.nodes)
    list = make_list(opf, idx, opf.nodes)
    # cgen = connectivitymatrix(system, length(nodes), idx)

    obj = objective_function(mod)
    # Pᵢ = get_net_Pᵢ(opf, idx)
    Pᵢ = get_value(mod, :p0)
    Pg = get_controllable(opf, mod, idx)
    Pd = get_Pd(opf, idx) # Fixed injection
    @assert isapprox(Pg, (Pᵢ - Pd); atol=1e-8) string(Pg - (Pᵢ - Pd))
    return Benders(idx, list, obj, Pᵢ, Pg, Pd)
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
    time_limit_sec::Int64=600, 
    ramp_minutes=10, 
    ramp_mult=10, 
    max_shed=1.0, 
    max_curtail=1.0, 
    lim=1e-6,
    short_term_multi::Real=1.5, 
    long_term_multi::Real=1.2, 
    max_itr=length(contingencies),
    p_failure=0.0,
    silent=true,
    debug=false
)
    mod, opf, pf, oplim, Pc, Pcc, Pccx = opf_base(OPF(true, false, false, false, false), system, optimizer, voll=voll, contingencies=contingencies, prob=prob,
        time_limit_sec=time_limit_sec, ramp_minutes=ramp_minutes, ramp_mult=ramp_mult, max_shed=max_shed, max_curtail=max_curtail,
        short_term_multi=short_term_multi, long_term_multi=long_term_multi, p_failure=p_failure, silent=silent, debug=debug)

    solve_model!(mod)
    if !type.P & !type.C1 & !type.C2 & !type.C2F
        return mod, opf, pf, oplim, Pc, Pcc, Pccx, solve_time(mod)
    end
    return run_benders!(type, mod, opf, pf, oplim, Pc, Pcc, Pccx, lim, max_itr)
end

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
    lim=1e-6,
    max_itr=length(opf.contingencies)
)
    assert(type)
    @assert !type.C1 || isempty(Pc)
    @assert !type.C2 || isempty(Pcc)
    @assert !type.C2F || isempty(Pccx)
    
    total_solve_time = solve_time(mod)
    !has_values(mod) && return mod, opf, pf, oplim, Pc, Pcc, Pccx, total_solve_time
    @debug "lower_bound = $(objective_value(mod))"

    # Set variables
    bd = benders(opf, mod)
    calc_θ!(pf, bd.Pᵢ)
    calc_Pline!(pf)
    overloads = zeros(length(opf.contingencies))

    pg = get_value(mod, :pg0)
    for c_obj in opf.contingencies
        (typelist, c, cont) = typesort_component(c_obj, opf, bd.idx)
        # cont  = get_bus_idx(opf.contingencies[c], bd.idx)
        if is_islanded(pf, cont, c)
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont, c, pf.slack)
            type.C1 && init_P!(Pc, opf, oplim, mod, bd.obj, bd.list, islands, island, c)
            type.C2 && init_P!(Pcc, opf, oplim, mod, bd.obj, bd.list, islands, island, c)
            type.C2F && init_P!(Pccx, opf, oplim, mod, bd.obj, bd.list, islands, island, c)
            @debug "Island: Contingency $(string(typeof(c_obj))) $(get_name(c_obj))"
            if isempty(islands)
                fill!(flow, 0.0)
            else
                calculate_line_flows!(pf.vb_tmp, pf, cont, c, nodes=islands[island], branches=island_b)
                flow = pf.vb_tmp
            end
        else
            if typeof(c_obj) <: ACBranch
                flow = calculate_contingency_line_flows(mod, pf, bd.Pᵢ, cont, c, c_obj, Int[], Int[])
            else
                flow = calculate_contingency_line_flows(mod, pf, bd.Pᵢ, cont, c, c_obj, Int[], Int[], pg[c])
            end
        end

        # Calculate the power flow with the new outage and find if there are any overloads
        overload = filter_overload(flow, oplim.branch_rating * oplim.short_term_multi)

        if !isempty(overload)
            overloads[c] = maximum(x -> abs(x[2]), overload)
        end
    end

    total_solve_time = update_model!(mod, pf, bd, total_solve_time)
    !has_values(mod) && return mod, opf, pf, oplim, Pc, Pcc, Pccx, total_solve_time

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
    cut_added = 1
    for iterations in 1:max_itr
        if cut_added == 0 # loops until no new cuts are added for the contingencies
            # @printf "\nEND: Total solve time %.4f.\n" total_solve_time
            print_cuts(type, pre, corr1, corr2, corr2f)
            return mod, opf, pf, oplim, Pc, Pcc, Pccx, total_solve_time
        end
        cut_added = 0
        @info "Iteration $iterations"

        for c_obj in opf.contingencies[permutation]
            (typelist, c, cont) = typesort_component(c_obj, opf, bd.idx)
            if is_islanded(pf, cont, c)
                islands, island, island_b = handle_islands(pf.B, pf.DA, cont, c, pf.slack)
                inodes = islands[island]
            else
                empty!(islands)
                island = 0
                empty!(island_b)
                inodes = Int[]
            end

            if type.P
                olc = calculate_contingency_overload(brst, mod, pf, bd, cont, c, c_obj, inodes, island_b)
            end
            if type.C1
                olc = calculate_contingency_overload!(ΔPc, brst, Pc, opf, mod, pf, bd, cont, c, c_obj, inodes, island_b)
            end
            if type.C2
                olcc = calculate_contingency_overload!(ΔPcc, brlt, Pcc, opf, mod, pf, bd, cont, c, c_obj, inodes, island_b)
            end
            if type.C2F
                olccx = calculate_contingency_overload!(ΔPccx, brlt, Pccx, opf, mod, pf, bd, cont, c, c_obj, inodes, island_b)
            end
            if !isempty(olc) || !isempty(olcc) || !isempty(olccx) # ptdf calculation is more computational expensive than line flow
                if is_islanded(pf, cont, c)
                    ptdf = get_isf(pf, cont, c, islands, island, island_b)
                else
                    ptdf = get_isf(pf, cont, c)
                end
                set_tol_zero!(ptdf)

                type.C1 && get(Pc, c, 0) == 0 && fill!(ΔPc, 0.0)
                type.C2 && get(Pcc, c, 0) == 0 && fill!(ΔPcc, 0.0)
                type.C2F && get(Pccx, c, 0) == 0 && fill!(ΔPccx, 0.0)

                # # Calculate the power flow with the new outage and find if there are any overloads
                # overloads_c, overloads_cc, overloads_ccx = find_overloads(Val(type), flow, ptdf, bd.Pᵢ, oplim, ΔPc, ΔPcc, ΔPccx)
                # if !isempty(olc) || !isempty(overloads_c) 
                #     if isempty(olc)
                #         println([o[2] for o in overloads_c])
                #     elseif isempty(overloads_c) 
                #         println([o for o in olc])
                #     else
                #         println([(o1[2], o2[2], o1[2] - o2[2]) for (o1, o2) in zip(overloads_c, olc)])
                #     end
                #     # error(olc)
                # end
                # if !isempty(olcc) || !isempty(overloads_cc) 
                #     if isempty(olcc)
                #         println([o[2] for o in overloads_cc])
                #     elseif isempty(overloads_cc) 
                #         println([o for o in olcc])
                #     else
                #         println([(o1[2], o2[2], o1[2] - o2[2]) for (o1, o2) in zip(overloads_cc, olcc)])
                #     end
                #     # error(olcc)
                # end
            end

            # Cannot change the model before all data is exctracted!
            if !isempty(olc)
                if type.P
                    cut_added, pre = add_cut(opf, mod, bd, ptdf, olc, c_obj, cut_added, lim, pre)
                elseif type.C1
                    cut_added, corr1 = add_cut(Pc, opf, oplim, mod, bd, ΔPc, ptdf, olc, islands, island, c_obj, c, cut_added, lim, corr1)
                end
            end
            if !isempty(olcc)
                cut_added, corr2 = add_cut(Pcc, opf, oplim, mod, bd, ΔPcc, ptdf, olcc, islands, island, c_obj, c, cut_added, lim, corr2)
            end
            if !isempty(olccx)
                cut_added, corr2f = add_cut(Pccx, opf, oplim, mod, bd, ΔPccx, ptdf, olccx, islands, island, c_obj, c, cut_added, lim, corr2f)
            end
            if cut_added > 1
                total_solve_time = update_model!(mod, pf, bd, total_solve_time)
                cut_added = 1
            end
            !has_values(mod) && return mod, opf, pf, oplim, Pc, Pcc, Pccx, total_solve_time
        end

    end
    @warn "Reached $(length(opf.contingencies)) iterations without a stable solution."
    return mod, opf, pf, oplim, Pc, Pcc, Pccx, total_solve_time
end

""" Solve model and update the power flow object """
function update_model!(mod::Model, pf::DCPowerFlow, bd::Benders, total_solve_time::Real)
    # set_warm_start!(mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
    set_objective_function(mod, bd.obj)
    solve_model!(mod)
    total_solve_time += solve_time(mod)

    bd.Pᵢ = get_value(mod, :p0)
    @. bd.Pg = bd.Pᵢ - bd.Pd
    calc_θ!(pf, bd.Pᵢ)
    calc_Pline!(pf)
    return total_solve_time
end

"""
    Calculate the contingency power flow with corrective actions.
    Assummes that island-contingencies have active variables from the pre-procedure.
"""
function calculate_contingency_line_flows!(ΔP::Vector{<:Real}, P::Dict{<:Integer, T}, opf::OPFsystem, 
    mod::Model, pf::DCPowerFlow, Pᵢ::AbstractVector{<:Real}, list::Vector{<:CTypes{Int}}, cont::Tuple{Real,Real}, 
    c::Integer, c_obj::Component, island::Vector, island_b::Vector{<:Integer}
) where {T}
    if get(P, c, 0) == 0
        return calculate_contingency_line_flows(mod, pf, Pᵢ, cont, c, c_obj, island, island_b)
    end
    ΔP = get_ΔP!(ΔP, mod, opf, list, P, c)
    if iszero(ΔP)
        calculate_line_flows!(pf.vb_tmp, pf, cont, c, nodes=island, branches=island_b)
    else
        calculate_line_flows!(pf.vb_tmp, pf, cont, c, Pᵢ=(Pᵢ .+ ΔP), nodes=island, branches=island_b)
    end
    return pf.vb_tmp
end
function calculate_contingency_line_flows!(ΔP::Vector{<:Real}, P::Dict{<:Integer, T}, opf::OPFsystem, 
    mod::Model, pf::DCPowerFlow, Pᵢ::AbstractVector{<:Real}, list::Vector{<:CTypes{Int}}, cont::Real,
    c::Integer, c_obj::Component, island::Vector, island_b::Vector{<:Integer}
) where {T}
    if get(P, c, 0) == 0
        return calculate_contingency_line_flows(mod, pf, Pᵢ, cont, c, c_obj, island, island_b)
    end
    ΔP = get_ΔP!(ΔP, mod, opf, list, P, c)
    ΔP[cont] -= get_ΔP(mod, c, c_obj, P)
    LinearAlgebra.mul!(pf.vb_tmp, pf.ϕ, (Pᵢ .+ ΔP))
    return pf.vb_tmp
end

"""
    Calculate the contingency power flow without corrective actions.
"""
function calculate_contingency_line_flows(
    mod::Model, pf::DCPowerFlow, P::Vector{<:Real}, cont::Tuple{Real,Real}, c::Integer, c_obj::Component,
    island::Vector, island_b::Vector{<:Integer}, val=0.0
)
    calculate_line_flows!(pf.vb_tmp, pf, cont, c, nodes=island, branches=island_b)
    return pf.vb_tmp
end
function calculate_contingency_line_flows(
    mod::Model, pf::DCPowerFlow, P::Vector{<:Real}, cont::Real, c::Integer, c_obj::StaticInjection,
    island::Vector, island_b::Vector{<:Integer}, val=get_ΔP(mod, c, c_obj)
)
    P[cont] -= val
    LinearAlgebra.mul!(pf.vb_tmp, pf.ϕ, P)
    P[cont] += val
    return pf.vb_tmp
end

"""
    Calculate the contingency overflow with corrective actions.
"""
function calculate_contingency_overload!(ΔP::Vector{<:Real}, branch_rating::Vector{<:Real},
    Pc::Dict, opf::OPFsystem, mod::Model, pf::DCPowerFlow, bd::Benders, 
    cont::Union{Real, Tuple{Real,Real}}, c::Integer, c_obj::Component, 
    island::Vector, island_b::Vector{<:Integer}
)
    flow = calculate_contingency_line_flows!(ΔP, Pc, opf, mod, pf, bd.Pᵢ, bd.list, cont, c, c_obj, island, island_b)
    return filter_overload(flow, branch_rating)
end

"""
    Calculate the contingency overflow without corrective actions.
"""
function calculate_contingency_overload(branch_rating::Vector{<:Real}, mod::Model, 
    pf::DCPowerFlow, bd::Benders, cont::Union{Real, Tuple{Real,Real}}, c::Integer, c_obj::Component,
    island::Vector, island_b::Vector{<:Integer}
)
    flow = calculate_contingency_line_flows(mod, pf, bd.Pᵢ, cont, c, c_obj, island, island_b)
    return filter_overload(flow, branch_rating)
end

get_ΔP(mod::Model, c::Integer, ::Generator) = 
    value(mod[:pg0][c])
get_ΔP(mod::Model, c::Integer, ::Generator, Pc::Dict{<:Integer, ExprC}) = 
    value(mod[:pg0][c]) - value(Pc[c].pgc[c])
get_ΔP(mod::Model, c::Integer, ::Generator, Pcc::Dict{<:Integer, ExprCC}) = 
    value(mod[:pg0][c]) + value(Pcc[c].pgu[c]) - value(Pcc[c].pgd[c])
get_ΔP(mod::Model, c::Integer, ::Generator, Pccx::Dict{<:Integer, ExprCCX}) = 
    value(mod[:pg0][c]) - value(Pccx[c].pgdx[c])

get_ΔP(mod::Model, c::Integer, ::StaticLoad) = 
    value(mod[:pd0][c]) - value(mod[:ls0][c])
get_ΔP(mod::Model, c::Integer, ::StaticLoad, Pc::Dict{<:Integer, ExprC}) = 
    value(mod[:pd0][c]) - value(Pc[c].lsc[c])
get_ΔP(mod::Model, c::Integer, ::StaticLoad, Pcc::Dict{<:Integer, ExprCC}) = 
    value(mod[:pd0][c]) - value(Pcc[c].lscc[c])
get_ΔP(mod::Model, c::Integer, ::StaticLoad, Pccx::Dict{<:Integer, ExprCCX}) = 
    value(mod[:pd0][c]) - value(Pccx[c].lsccx[c])

""" Return the short term post-contingency power injection change at each node. """
function get_ΔP!(ΔP::Vector{T}, mod::Model, opf::OPFsystem, list::Vector{<:CTypes{Int}}, 
    Pc::Dict{<:Integer, ExprC}, c::Integer
) where {T<:Real}
    fill!(ΔP, zero(T))
    x = get(Pc, c, 0)
    if x != 0
        pgc = get_value(mod, x.pgc)
        prc = get_value(mod, x.prc)
        lsc = get_value(mod, x.lsc)
        for (i, n) in enumerate(list)
            for g in n.ctrl_generation
                ΔP[i] -= pgc[g]
            end
            for g in n.renewables
                ΔP[i] -= prc[g]
            end
            for g in n.demands
                ΔP[i] += lsc[g]
            end
        end
    end
    return ΔP
end

""" Return the long term post-contingency power injection change at each node. """
function get_ΔP!(ΔP::Vector{T}, mod::Model, opf::OPFsystem, list::Vector{<:CTypes{Int}}, 
    Pcc::Dict{<:Integer, ExprCC}, c::Integer
) where {T<:Real}
    fill!(ΔP, zero(eltype(ΔP)))
    x = get(Pcc, c, 0)
    if x != 0
        pgu = get_value(mod, x.pgu)
        pgd = get_value(mod, x.pgd)
        pfdccc = get_value(mod, x.pfdccc)
        prcc = get_value(mod, x.prcc)
        lscc = get_value(mod, x.lscc)
        for (i, n) in enumerate(list)
            for g in n.ctrl_generation
                ΔP[i] += (pgu[g] - pgd[g])
            end
            for g in n.dc_branches
                ΔP[i] += beta(opf.nodes[i], opf.dc_branches[g]) * pfdccc[g]
            end
            for g in n.renewables
                ΔP[i] -= prcc[g]
            end
            for g in n.demands
                ΔP[i] += lscc[g]
            end
        end
    end
    return ΔP
end

""" Return the long term post-contingency with failure power injection change at each node. """
function get_ΔP!(ΔP::Vector{T}, mod::Model, opf::OPFsystem, list::Vector{<:CTypes{Int}}, 
    Pccx::Dict{<:Integer, ExprCCX}, c::Integer
) where {T<:Real}
    fill!(ΔP, zero(eltype(ΔP)))
    x = get(Pccx, c, 0)
    if x != 0
        pgdx = get_value(mod, x.pgdx)
        lsccx = get_value(mod, x.lsccx)
        for (i, n) in enumerate(list)
            for g in n.ctrl_generation
                ΔP[i] -= pgdx[g]
            end
            for g in n.demands
                ΔP[i] += lsccx[g]
            end
        end
    end
    return ΔP
end

function get_P(P::Dict, opf::OPFsystem, oplim::Oplimits, mod::Model, obj::JuMP.AbstractJuMPScalar, list::Vector{<:CTypes{Int}},
    islands::Vector, island::Integer, c::Integer
)
    if get(P, c, 0) == 0  # If the contingency is not run before, a set of corrective variables is added
        init_P!(P, opf, oplim, mod, obj, list, islands, island, c)
    end
    return get(P, c, 0)
end

function add_overload_expr!(mod::Model, expr::JuMP.AbstractJuMPScalar, ol::Real, type::String, id::Integer,
    c_obj::Component, opf::OPFsystem, i::Integer, lim::Real
)
    @info @sprintf "%s %d: Contingency %s %s; overload on %s of %.4f" type id string(typeof(c_obj)) get_name(c_obj) opf.branches[i].name ol
    @debug "Cut added: $(sprint_expr(expr,lim))\n"
    if ol < 0
        JuMP.@constraint(mod, expr <= ol)
    else
        JuMP.@constraint(mod, expr >= ol)
    end
end

function add_cut(opf::OPFsystem, mod::Model, bd::Benders, ptdf::AbstractMatrix{<:Real}, overloads::Vector{<:Tuple{Integer,Real}},
    c_obj::Component, cut_added::Integer, lim::Real, id::Integer
)
    for (i, ol) in overloads
        expr = JuMP.@expression(mod, sum((ptdf[i, inode] * (
                bd.Pg[inode] -
                sum((mod[:pg0][ctrl] for ctrl in sublist.ctrl_generation), init=0.0) -
                sum((beta(sublist.node, opf.dc_branches[d]) * mod[:pfdc0][d] for d in sublist.dc_branches), init=0.0) +
                sum((mod[:pr0][r] for r in sublist.renewables), init=0.0) -
                sum((mod[:ls0][d] for d in sublist.demands), init=0.0)
            ) for (inode, sublist) in enumerate(bd.list)), init=0.0
        ))

        id += 1
        add_overload_expr!(mod, expr, ol, "Pre", id, c_obj, opf, i, lim)
        cut_added = 2
    end
    return cut_added, id
end

function add_cut(Pc::Dict{<:Integer, ExprC}, opf::OPFsystem, oplim::Oplimits, mod::Model, bd::Benders, ΔPc::Vector{<:Real},
    ptdf::AbstractMatrix{<:Real}, overloads::Vector{<:Tuple{Integer,Real}}, islands::Vector, island::Integer,
    c_obj::Component, c::Integer, cut_added::Integer, lim::Real, id::Integer
)
    pc = get_P(Pc, opf, oplim, mod, bd.obj, bd.list, islands, island, c)
    for (i, ol) in overloads
        expr = JuMP.@expression(mod, sum((ptdf[i, inode] * (
                bd.Pg[inode] + ΔPc[inode] -
                sum((mod[:pg0][ctrl] - pc.pgc[ctrl] for ctrl in sublist.ctrl_generation), init=0.0) -
                sum((beta(sublist.node, opf.dc_branches[d]) * mod[:pfdc0][d] for d in sublist.dc_branches), init=0.0) +
                sum((pc.prc[r] for r in sublist.renewables), init=0.0) -
                sum((pc.lsc[d] for d in sublist.demands), init=0.0)
            ) for (inode, sublist) in enumerate(bd.list)), init=0.0
        ))

        id += 1
        add_overload_expr!(mod, expr, ol, "Corr1", id, c_obj, opf, i, lim)
        cut_added = 2
    end
    return cut_added, id
end

function add_cut(Pcc::Dict{<:Integer, ExprCC}, opf::OPFsystem, oplim::Oplimits, mod::Model, bd::Benders, ΔPcc::Vector{<:Real},
    ptdf::AbstractMatrix{<:Real}, overloads::Vector{<:Tuple{Integer,Real}}, islands::Vector, island::Integer,
    c_obj::Component, c::Integer, cut_added::Integer, lim::Real, id::Integer
)
    pcc = get_P(Pcc, opf, oplim, mod, bd.obj, bd.list, islands, island, c)
    for (i, ol) in overloads
        # Finding and adding the Benders cut
        expr = JuMP.@expression(mod, sum((ptdf[i, inode] * (
                bd.Pg[inode] + ΔPcc[inode] -
                sum((mod[:pg0][ctrl] + pcc.pgu[ctrl] - pcc.pgd[ctrl] for ctrl in sublist.ctrl_generation), init=0.0) -
                sum((beta(sublist.node, opf.dc_branches[d]) * pcc.pfdccc[d] for d in sublist.dc_branches), init=0.0) +
                sum((pcc.prcc[r] for r in sublist.renewables), init=0.0) -
                sum((pcc.lscc[d] for d in sublist.demands), init=0.0)
            ) for (inode, sublist) in enumerate(bd.list)), init=0.0
        ))

        id += 1
        add_overload_expr!(mod, expr, ol, "Corr2", id, c_obj, opf, i, lim)
        cut_added = 3
    end
    return cut_added, id
end

function add_cut(Pccx::Dict{<:Integer, ExprCCX}, opf::OPFsystem, oplim::Oplimits, mod::Model, bd::Benders, ΔPccx::Vector{<:Real},
    ptdf::AbstractMatrix{<:Real}, overloads::Vector{<:Tuple{Integer,Real}}, islands::Vector, island::Integer,
    c_obj::Component, c::Integer, cut_added::Integer, lim::Real, id::Integer
)
    pccx = get_P(Pccx, opf, oplim, mod, bd.obj, bd.list, islands, island, c)
    # sort!(overloads, rev = true, by = x -> abs(x[2]))
    for (i, ol) in overloads
        # Finding and adding the Benders cut
        expr = JuMP.@expression(mod, sum((ptdf[i, inode] * (
                bd.Pg[inode] + ΔPccx[inode] -
                sum((mod[:pg0][ctrl] - pccx.pgdx[ctrl] for ctrl in sublist.ctrl_generation), init=0.0) +
                sum((pccx.prccx[r] for r in sublist.renewables), init=0.0) -
                sum((pccx.lsccx[d] for d in sublist.demands), init=0.0)
            ) for (inode, sublist) in enumerate(bd.list)), init=0.0
        ))

        id += 1
        add_overload_expr!(mod, expr, ol, "Corr2x", id, c_obj, opf, i, lim)
        cut_added = 3
    end
    return cut_added, id
end

" An AbstractJuMPScalar nicely formatted to a string "
sprint_expr(expr::AbstractJuMPScalar, lim=1e-6) =
    join(Printf.@sprintf("%s%5.2f %s ", (x[2] > 0 ? "+" : "-"), abs(x[2]), x[1])
         for x in expr.terms if abs(x[2]) > lim) *
    Printf.@sprintf("<= %s%.2f", (expr.constant > 0 ? "-" : " "), abs(expr.constant)
    )

function print_benders_results(opf::OPFsystem, mod::Model, Pc::Dict=Dict(), Pcc::Dict=Dict(), Pccx::Dict=Dict(), lim::Real=1e-6)
    function print_c(itr, symb::String, i_g::Int, lim::Real)
        for i in 1:length(opf.contingencies)
            c = get(itr, i, 0)
            if c != 0 && JuMP.value(getfield(c, Symbol(symb))[i_g]) > lim
                @printf("          c %12s: %s: %.3f\n", opf.contingencies[i].name, symb, JuMP.value(getfield(c, Symbol(symb))[i_g]))
            end
        end
    end
    for (i_g, g) in enumerate(opf.ctrl_generation)
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(mod[:pg0][i_g]), get_active_power_limits(g).max)
        print_c(Pc, "pgc", i_g, lim)
        print_c(Pcc, "pgu", i_g, lim)
        print_c(Pcc, "pgd", i_g, lim)
        print_c(Pccx, "pgdx", i_g, lim)
    end
    for (i_g, g) in enumerate(opf.dc_branches)
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(mod[:pfdc0][i_g]), get_active_power_limits(g).max)
        print_c(Pcc, "pfdccc", i_g, lim)
    end
    for (i_g, g) in enumerate(opf.renewables)
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(mod[:pr0][i_g]), get_active_power_limits(g).max)
        print_c(Pc, "prc", i_g, lim)
        print_c(Pcc, "prcc", i_g, lim)
        print_c(Pccx, "prccx", i_g, lim)
    end
    for (i_g, g) in enumerate(opf.demands)
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(mod[:ls0][i_g]), get_active_power(g))
        print_c(Pc, "lsc", i_g, lim)
        print_c(Pcc, "lscc", i_g, lim)
        print_c(Pccx, "lsccx", i_g, lim)
    end
end

function print_cuts(type::OPF, pre, corr1, corr2, corr2f)
    print("Cuts added:")
    type.P && print(" pre=", pre)
    type.C1 && print(" corr1=", corr1)
    type.C2 && print(" corr2=", corr2)
    type.C2F && print(" corr2f=", corr2f)
    println("")
end