""" ModelState type """
mutable struct ModelState{TR<:Real}
    Pᵢ::Vector{TR}
    Pg::Vector{TR}
    Pd::Vector{TR}
end

""" Constructor for ModelState type """
function decomposition(opf::OPFsystem, m::Model)
    # Pᵢ = get_net_Pᵢ(opf)
    Pᵢ = get_value(m, :p0)
    Pg = get_controllable(opf, m)
    Pd = get_Pd(opf, m) # Fixed injection
    @assert isapprox(Pg, (Pᵢ - Pd); atol=1e-10) string(Pg - (Pᵢ - Pd))
    return ModelState{Float64}(Pᵢ, Pg, Pd)
end

""" 
Solve the optimization model using decomposition.

Function creates extra production constraints on generators in the system
based on the power transfer distribution factors of the generators in the
system. 
"""
function run_decomposition(
    type::OPF, 
    system::System, 
    optimizer, 
    voll::Vector{<:Number}, 
    prob::Vector{<:Number}, 
    contingencies::Vector{Pair{String, Vector{Int64}}};
    dist_slack::Vector{<:Real}=Float64[],
    time_limit_sec::Number=600, 
    ramp_minutes::Number=10, 
    ramp_mult::Number=10, 
    renew_cost::Number=0.0,
    renew_ramp::Number=1000.0,
    max_shed::Number=1.0, 
    max_curtail::Number=1.0, 
    reserve=0.02,
    atol::Number=1e-6,
    short_term_multi::Real=1.5, 
    long_term_multi::Real=1.2, 
    max_itr::Integer=10,
    p_failure::Number=0.0,
    zero_c_cost::Bool = false,
    silent::Bool=true,
    debug::Bool=false
)
    case = Case(opf_base(OPF(true, false, false, false, false), system, optimizer, voll=voll, contingencies=contingencies, prob=prob,
        dist_slack=dist_slack, time_limit_sec=time_limit_sec, ramp_minutes=ramp_minutes, ramp_mult=ramp_mult, renew_cost=renew_cost, 
        renew_ramp=renew_ramp, max_shed=max_shed, max_curtail=max_curtail, reserve=reserve, short_term_multi=short_term_multi, long_term_multi=long_term_multi, 
        p_failure=p_failure, zero_c_cost=zero_c_cost, silent=silent, debug=debug)...)

    total_solve_time = constrain_branches!(case, 0.0)
    if !type.P & !type.C1 & !type.C2 & !type.C2F
        return case, total_solve_time
    end
    return run_decomposition!(type, case, atol, max_itr)
end

run_decomposition!(type::OPF, case::Case, atol=1e-6, max_itr=max(length(case.opf.contingencies), 5)) =
    run_decomposition!(type, case.model, case.opf, case.pf, case.oplim, case.brc_up, case.brc_down, case.Pc, case.Pcc, case.Pccx, atol, max_itr)

"""
Solve the optimization model using decomposition.
"""
function run_decomposition!(
    type::OPF, 
    m::Model, 
    opf::OPFsystem, 
    pf::DCPowerFlow, 
    oplim::Oplimits,
    brc_up::Dict{<:Integer, ConstraintRef}, 
    brc_down::Dict{<:Integer, ConstraintRef},
    Pc::Dict{<:Integer, ExprC}, 
    Pcc::Dict{<:Integer, ExprCC}, 
    Pccx::Dict{<:Integer, ExprCCX},
    atol::Number=1e-6,
    max_itr::Integer=10
)
    # assert(type)
    # @assert !type.C1 || isempty(Pc)
    # @assert !type.C2 || isempty(Pcc)
    # @assert !type.C2F || isempty(Pccx)
    total_solve_time = solve_time(m)
    !is_solved_and_feasible(m) && return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time
    @debug "lower_bound = $(objective_value(m))"

    # Set variables
    !isempty(oplim.dist_slack) && set_dist_slack!(pf.ϕ, opf.mgx, oplim.dist_slack)
    state = decomposition(opf, m)
    if type.P
        @constraint(m, sum(m[:pgru]) >= oplim.reserve)
        @constraint(m, sum(m[:pgrd]) >= oplim.reserve)
    end

    ΔPc = zeros(length(state.Pg))
    ΔPcc = zeros(length(state.Pg))
    ΔPccx = zeros(length(state.Pg))
    pre = 0
    corr1 = 0
    corr2 = 0
    corr2f = 0

    brst = oplim.branch_rating * oplim.short_term_multi
    brlt = oplim.branch_rating * oplim.long_term_multi

    olp = Vector{Tuple{Int64,Float64}}()
    olc = Vector{Tuple{Int64,Float64}}()
    olcc = Vector{Tuple{Int64,Float64}}()
    olccx = Vector{Tuple{Int64,Float64}}()
    islands = [Vector{Int64}()]
    islands_b = [Vector{Int64}()]
    ptdf = similar(pf.vn_tmp)

    overloads = zeros(length(opf.contingencies))
    # PI = zeros(length(opf.contingencies))
    for (c, c_obj) in enumerate(opf.contingencies)
        c_i = c_obj.second
        length(c_i) == 0 && continue
        if c_obj.first == "branch"
            c_n = [get_bus_idx(opf.mbx[i,:]) for i in c_i]
            if is_islanded(pf, c_n, c_i)
                islands, islands_b = handle_islands(pf.B, pf.DA, c_n, c_i)
                type.C1 && init_P!(Pc, opf, oplim, m, islands, c)
                type.C2 && init_P!(Pcc, opf, oplim, m, islands, c)
                type.C2F && init_P!(Pccx, opf, oplim, m, islands, c)
                # @debug "Island: Contingency $(string(typeof(c_obj))) $(get_name(c_obj))"
                overload = []
                if !isempty(islands)
                    for (island, island_b) in zip(islands, islands_b)
                        s = find_slack(island, pf.slack, oplim.pg_lim_max, opf.mgx)
                        s < 1 && continue
                        calculate_line_flows!(pf.vb_tmp, pf.vn_tmp, pf.DA, pf.B, state.Pᵢ, c_n, c_i, s, island, island_b)
                        push!(overload, filter_overload(pf.vb_tmp, brst))
                    end
                    overload = reduce(vcat, overload)
                end
            else
                calculate_line_flows!(pf.vb_tmp, pf, c_n, c_i, state.Pᵢ)
                overload = filter_overload(pf.vb_tmp, brst)
            end
            # PI[i] = maximum((flow ./ (brst)).^2)
        elseif c_obj.first == "gen"
            c_n = [opf.mgx[i,:].nzind[1] for i in c_i]
            overload = calculate_contingency_overload!(brst, pf, c_n, c_i, islands[1], islands_b[1], pf.slack, atol)
        end

        if !isempty(overload)
            overloads[c] = maximum(x -> abs(x[2]), overload)
        end
    end

    # total_solve_time = update_model!(m, pf, oplim, brc_up, brc_down, state, total_solve_time)
    # !is_solved_and_feasible(m) && return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time

    permutation = sortperm(overloads, rev=true)
    # permutation = sortperm(PI, rev=true)
    obj_val = Inf
    cut_added = []
    for iterations in 1:max_itr
        total_solve_time = update_model!(m, pf, oplim, brc_up, brc_down, state, total_solve_time)
        !is_solved_and_feasible(m) && return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time
        if obj_val == JuMP.objective_value(m) # loops until no new cuts are added for the contingencies
            # @printf "\nEND: Total solve time %.4f.\n" total_solve_time
            print_cuts(iterations, type, pre, corr1, corr2, corr2f)
            return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time
        end
        obj_val = JuMP.objective_value(m)
        @info "\n------------------\nIteration $iterations"

        for c in permutation
            c_obj = opf.contingencies[c]
            c_i = c_obj.second
            length(c_i) == 0 && continue
            c_n = c_obj.first == "branch" ? [get_bus_idx(opf.mbx[i,:]) for i in c_i] : [opf.mgx[i,:].nzind[1] for i in c_i]
            if is_islanded(pf, c_n, c_i)
                islands, islands_b = handle_islands(pf.B, pf.DA, c_n, c_i)
                slack = find_slack.(islands, [pf.slack], [oplim.pg_lim_max], [opf.mgx])
            else
                islands = [Vector{Int64}()]
                islands_b = [Vector{Int64}()]
                slack = [pf.slack]
            end

            if type.P
                olp = [s > 0 ? calculate_contingency_overload!(brst, pf, state.Pᵢ, c_n, c_i, inodes, ibranches, s, atol) : Vector{Tuple{Int64, Float64}}() for (inodes, ibranches, s) in zip(islands, islands_b, slack) if isempty(inodes)]
            end
            if type.C1
                olc = [s > 0 ? calculate_contingency_overload!(ΔPc, brst, Pc, opf, m, pf, state, c_n, c_i, c, inodes, ibranches, s, atol) : Vector{Tuple{Int64, Float64}}() for (inodes, ibranches, s) in zip(islands, islands_b, slack)]
            end
            if type.C2
                olcc = [s > 0 ? calculate_contingency_overload!(ΔPcc, brlt, Pcc, opf, m, pf, state, c_n, c_i, c, inodes, ibranches, s, atol) : Vector{Tuple{Int64, Float64}}() for (inodes, ibranches, s) in zip(islands, islands_b, slack)]
            end
            if type.C2F
                olccx = [s > 0 ? calculate_contingency_overload!(ΔPccx, brlt, Pccx, opf, m, pf, state, c_n, c_i, c, inodes, ibranches, s, atol) : Vector{Tuple{Int64, Float64}}() for (inodes, ibranches, s) in zip(islands, islands_b, slack)]
            end

            # Cannot change the model before all data is exctracted!
            if sum(.!isempty.(olp)) > 0
                for (j,ol) in enumerate(olp)
                    pre = add_cut(pf, m, brst, ptdf, ol, c_obj, cut_added, atol, pre, islands, j, islands_b, c_n, c_i)
                end
            end
            if sum(.!isempty.(olc)) > 0
                get(Pc, c, 0) == 0 && fill!(ΔPc, 0.0)
                for (j,ol) in enumerate(olc)
                    corr1 = add_cut(Pc, pf, opf, oplim, m, state, brst, ptdf, ol, islands, j, islands_b, c_n, c_i, c_obj, c, cut_added, atol, corr1)
                end
            end
            if sum(.!isempty.(olcc)) > 0
                type.C2 && get(Pcc, c, 0) == 0 && fill!(ΔPcc, 0.0)
                for (j,ol) in enumerate(olcc)
                    corr2 = add_cut(Pcc, pf, opf, oplim, m, state, brlt, ptdf, ol, islands, j, islands_b, c_n, c_i, c_obj, c, cut_added, atol, corr2)
                end
            end
            if sum(.!isempty.(olccx)) > 0
                type.C2F && get(Pccx, c, 0) == 0 && fill!(ΔPccx, 0.0)
                for (j,ol) in enumerate(olccx)
                    corr2f = add_cut(Pccx, pf, opf, oplim, m, state, brlt, ptdf, ol, islands, j, islands_b, c_n, c_i, c_obj, c, cut_added, atol, corr2f)
                end
            end
            # if !isempty(cut_added)
            #     total_solve_time = update_model!(m, pf, oplim, brc_up, brc_down, state, total_solve_time)
            #     if !is_solved_and_feasible(m)
            #         @warn "Infeasible solution for Contingency $(c_obj.first) $(c_obj.second). Deleting latest cuts!"
            #         delete.([m], cut_added)
            #         total_solve_time = update_model!(m, pf, oplim, brc_up, brc_down, state, total_solve_time)
            #         !is_solved_and_feasible(m) && return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time
            #     end
            #     cut_added = []
            # end
            # !is_solved_and_feasible(m) && return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time
        end

    end
    @warn "Reached $(max_itr) iterations without a stable solution."
    return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time
end

""" Solve model and update the power flow object """
function update_model!(m::Model, pf::DCPowerFlow, oplim::Oplimits, brc_up::Dict{<:Integer, ConstraintRef}, 
    brc_down::Dict{<:Integer, ConstraintRef}, state::ModelState, total_solve_time::Real
)
    # set_warm_start!(m, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
    set_objective_function(m, m[:obj])
    total_solve_time = constrain_branches!(m, pf, oplim, brc_up, brc_down, total_solve_time)
    state.Pᵢ = get_value(m, :p0)
    @. state.Pg = state.Pᵢ - state.Pd
    return total_solve_time
end

function calculate_contingency_line_flows!(pf::DCPowerFlow, c_n::Vector{<:Tuple{Real,Real}}, c_i::Vector{<:Real},
    island::Vector, island_b::Vector{<:Integer}, slack::Integer, Pᵢ::AbstractVector{<:Real}
)
    if length(c_i) == 1
        calculate_line_flows!(pf.vb_tmp, pf, c_n[1], c_i[1], island, island_b)
    else
        if isempty(island)
            calculate_line_flows!(pf.vb_tmp, pf, c_n, c_i, Pᵢ)
        else
            calculate_line_flows!(pf.vb_tmp, pf.vn_tmp, pf.DA, pf.B, Pᵢ, c_n, c_i, slack, island, island_b)
        end
    end
    return pf.vb_tmp
end
function calculate_contingency_line_flows!(
    m::Model, pf::DCPowerFlow, Pᵢ::AbstractVector{<:Real}, c_n::Vector{<:Integer}, c_i::Vector{<:Integer}
) 
    val = get_gen_ΔP.([m], c_i)
    @. Pᵢ[c_n] -= val
    LinearAlgebra.mul!(pf.vb_tmp, pf.ϕ, Pᵢ)
    @. Pᵢ[c_n] += val
    return pf.vb_tmp
end


"""
    Calculate the contingency power flow accounting for corrective actions.
    Assummes that island-contingencies have active variables from the pre-procedure.
"""
function calculate_contingency_line_flows!(ΔP::Vector{<:Real}, P::Dict{<:Integer, T}, opf::OPFsystem, 
    m::Model, pf::DCPowerFlow, Pᵢ::AbstractVector{<:Real}, c_n::Vector{<:Tuple{Integer,Integer}}, c_i::Vector{<:Integer},
    c::Integer, island::Vector, island_b::Vector{<:Integer}, slack::Integer
) where {T}
    if isempty(island)
        if length(c_i) == 1
            if get(P, c, 0) == 0
                calculate_line_flows!(pf.vb_tmp, pf, c_n[1], c_i[1])
            else
                ΔP = get_ΔP!(ΔP, m, opf, P, c)
                calculate_line_flows!(pf.vb_tmp, pf, c_n[1], c_i[1], (Pᵢ .+ ΔP))
            end
        else
            if get(P, c, 0) == 0
                calculate_line_flows!(pf.vb_tmp, pf.vn_tmp, similar(pf.B), pf.DA, pf.B, Pᵢ, c_n, c_i, slack)
            else
                ΔP = get_ΔP!(ΔP, m, opf, P, c)
                calculate_line_flows!(pf.vb_tmp, pf.vn_tmp, similar(pf.B), pf.DA, pf.B, (Pᵢ .+ ΔP), c_n, c_i, slack)
            end
        end
    else
        if get(P, c, 0) == 0
            calculate_line_flows!(pf.vb_tmp, pf.vn_tmp, pf.DA, pf.B, Pᵢ, c_n, c_i, slack, island, island_b)
        else
            ΔP = get_ΔP!(ΔP, m, opf, P, c)
            calculate_line_flows!(pf.vb_tmp, pf.vn_tmp, pf.DA, pf.B, (Pᵢ .+ ΔP), c_n, c_i, slack, island, island_b)
        end
    end
    return pf.vb_tmp
end
function calculate_contingency_line_flows!(ΔP::Vector{<:Real}, P::Dict{<:Integer, T}, opf::OPFsystem, 
    m::Model, pf::DCPowerFlow, Pᵢ::AbstractVector{<:Real}, c_n::Vector{<:Integer}, c_i::Vector{<:Integer},
    c::Integer, island::Vector, island_b::Vector{<:Integer}, slack::Integer
) where {T}
    if get(P, c, 0) == 0
        return calculate_contingency_line_flows!(m, pf, Pᵢ, c_n, c_i)
    end
    ΔP = get_ΔP!(ΔP, m, opf, P, c)
    ΔP[c_n] .-= get_gen_ΔP.([m], c_i, [P], [c])
    LinearAlgebra.mul!(pf.vb_tmp, pf.ϕ, (Pᵢ .+ ΔP))
    return pf.vb_tmp
end

"""
    Calculate the contingency overflow.
    c_n = affected node(s)
    c_i = contingecy index in it's component vector
"""
function calculate_contingency_overload!(branch_rating::Vector{<:Real}, pf::DCPowerFlow, 
    Pᵢ::AbstractVector{<:Real}, c_n, c_i, 
    island::Vector, island_b::Vector{<:Integer}, slack::Integer, atol::Real=1e-6
)
    flow = calculate_contingency_line_flows!(pf, c_n, c_i, island, island_b, slack, Pᵢ)
    flow = zero_outage!(flow, island_b)
    return filter_overload(flow, branch_rating, atol)
end

"""
    Calculate the contingency overflow accounting for corrective actions.
    c_n = affected node(s)
    c_i = contingecy index in it's component vector
    c = contingecy index
"""
function calculate_contingency_overload!(ΔP::Vector{<:Real}, branch_rating::Vector{<:Real},
    Pc::Dict, opf::OPFsystem, m::Model, pf::DCPowerFlow, state::ModelState, 
    c_n, c_i, c::Integer, 
    island::Vector, island_b::Vector{<:Integer}, slack::Integer, atol::Real=1e-6
)
    flow = calculate_contingency_line_flows!(ΔP, Pc, opf, m, pf, state.Pᵢ, c_n, c_i, c, island, island_b, slack)
    flow = zero_outage!(flow, island_b)
    return filter_overload(flow, branch_rating, atol)
end

function zero_outage!(flow::Vector, branches::Vector)
    if !isempty(branches)
        zero_not_in_array!(flow, branches)
    end
    return flow
end

get_gen_ΔP(m::Model, cont::Real) = 
    value(m[:pg0][cont])
get_gen_ΔP(m::Model, cont::Real, Pc::Dict{<:Integer, ExprC}, c::Integer) = 
    value(m[:pg0][cont]) + value(Pc[c].pgu[cont]) - value(Pc[c].pgd[cont])
get_gen_ΔP(m::Model, cont::Real, Pcc::Dict{<:Integer, ExprCC}, c::Integer) = 
    value(m[:pg0][cont]) + value(Pcc[c].pgu[cont]) - value(Pcc[c].pgd[cont])
get_gen_ΔP(m::Model, cont::Real, Pccx::Dict{<:Integer, ExprCCX}, c::Integer) = 
    value(m[:pg0][cont]) - value(Pccx[c].pgdx[cont])

get_load_ΔP(m::Model, cont::Real) = 
    value(m[:pd0][cont]) - value(m[:ls0][cont])
get_load_ΔP(m::Model, cont::Real, Pc::Dict{<:Integer, ExprC}, c::Integer) = 
    value(m[:pd0][cont]) - value(Pc[c].lsc[cont])
get_load_ΔP(m::Model, cont::Real, Pcc::Dict{<:Integer, ExprCC}, c::Integer) = 
    value(m[:pd0][cont]) - value(Pcc[c].lscc[cont])
get_load_ΔP(m::Model, cont::Real, Pccx::Dict{<:Integer, ExprCCX}, c::Integer) = 
    value(m[:pd0][cont]) - value(Pccx[c].lsccx[cont])

""" Return the contingency power injection change at each node changing ΔP. """
function get_ΔP!(ΔP::Vector{T}, m::Model, opf::OPFsystem, 
    P::Dict{<:Integer, <:ContExpr}, c::Integer
) where {T<:Real}
    x = get(P, c, 0)
    if x != 0
        copy!(ΔP, get_ΔP(x, opf, m))
        sumP = sum(ΔP)
        if sumP > 1e-6
            ΔP .+= (get_value(m, m[:pgru])' * opf.mgx) .* sumP
        elseif sumP < -1e-6
            ΔP .+= (get_value(m, m[:pgrd])' * opf.mgx) .* sumP
        end
    else
        fill!(ΔP, zero(T))
    end
    return ΔP
end

""" Return the short term post-contingency power injection change at each node. """
@inline function get_ΔP(pc::ExprC, opf::OPFsystem, m::Model)
    return opf.mgx' * (get_value(m, pc.pgu) - get_value(m, pc.pgd)) + opf.mdx' * get_value(m, pc.lsc)
end

""" Return the long term post-contingency power injection change at each node. """
@inline function get_ΔP(pcc::ExprCC, opf::OPFsystem, m::Model)
    return opf.mgx' * (get_value(m, pcc.pgu) - get_value(m, pcc.pgd)) + opf.mdx' * get_value(m, pcc.lscc) +
        opf.mdcx' * get_value(m, pcc.pfdccc)
end

""" Return the long term post-contingency with failure power injection change at each node. """
@inline function get_ΔP(pccx::ExprCCX, opf::OPFsystem, m::Model)
    return -opf.mgx' * get_value(m, pccx.pgdx) + opf.mdx' * get_value(m, pccx.lsccx)
end

function add_overload_expr!(m::Model, expr::JuMP.AbstractJuMPScalar, ol::Real, rate::Real, type::String, id::Integer,
    c_obj::Pair{String, Vector{Int64}}, i::Integer, atol::Real
)
    @info "$type $id: C $(c_obj.first) $(c_obj.second); ol br $i of $(abs(ol))"
    @debug "Cut added: $(sprint_expr(expr,atol))\n"
    if ol < 0
        return JuMP.@constraint(m, expr >= -rate)
    else
        return JuMP.@constraint(m, expr <= rate)
    end
end

function add_cut(pf::DCPowerFlow, m::Model, rate::Vector{<:Real}, ptdf::AbstractVector{<:Real}, overloads::Vector{<:Tuple{Integer,Real}},
    c_obj::Pair{String, Vector{Int64}}, cut_added::Vector, atol::Real, id::Integer, islands::Vector, island::Integer, islands_b::Vector, c_n, c_i
)
    for (i, ol) in overloads
        if is_islanded(pf, c_n, c_i)
            calc_isf_vec!(ptdf, pf.DA, pf.B, c_n, c_i, i, pf.slack, islands[island], islands_b[island])
        else
            calc_isf_vec!(ptdf, pf, c_n, c_i, Tuple(pf.DA[i,:].nzind), i)
        end
        expr = AffExpr() #+ sum(ptdf[i, :] .* state.Pᵢ)
        add_to_expression!.(expr, ptdf, m[:inj_p0])

        id += 1
        push!(cut_added, add_overload_expr!(m, expr, ol, rate[i], "Pre", id, c_obj, i, atol))
    end
    return id
end

""" Finding and adding cuts to mitigate overloads from one contingency """
function add_cut(P::Dict{<:Integer, <:ContExpr}, pf::DCPowerFlow, opf::OPFsystem, oplim::Oplimits, m::Model, state::ModelState, rate::Vector{<:Real},
    ptdf::AbstractVector{<:Real}, overloads::Vector{<:Tuple{Integer,Real}}, islands::Vector, island::Integer, islands_b::Vector, c_n, c_i,
    c_obj::Pair{String, Vector{Int64}}, c::Integer, cut_added::Vector, atol::Real, id::Integer
)
    p = get(P, c, 0)
    if p == 0  # If the contingency is not run before, a set of corrective variables is added
        init_P!(P, opf, oplim, m, islands, c)
        p = get(P, c, 0)
    end
    for (i, ol) in overloads
        if is_islanded(pf, c_n, c_i)
            calc_isf_vec!(ptdf, pf.DA, pf.B, c_n, c_i, i, pf.slack, islands[island], islands_b[island])
        else
            calc_isf_vec!(ptdf, pf, c_n, c_i, Tuple(pf.DA[i,:].nzind), i)
        end
        expr = make_cut(p, opf, m)
        # add_to_expression!.(expr, state.Pᵢ + ΔP)
        expr2 = JuMP.@expression(m, sum((ptdf .* expr)))

        id += 1
        push!(cut_added, add_overload_expr!(m, expr2, ol, rate[i], get_name(p), id, c_obj, i, atol))
    end
    return id
end

function make_cut(pc::ExprC, opf::OPFsystem, m::Model)
    expr = JuMP.@expression(m, opf.mgx' * (m[:pg0] + pc.pgu - pc.pgd))
    add_to_expression!.(expr, opf.mdcx' * m[:pfdc0])
    add_to_expression!.(expr, opf.mdx' * (pc.lsc - m[:pd]))
    return expr
end

function make_cut(pcc::ExprCC, opf::OPFsystem, m::Model)
    expr = JuMP.@expression(m, opf.mgx' * (m[:pg0] + pcc.pgu - pcc.pgd))
    add_to_expression!.(expr, opf.mdcx' * pcc.pfdccc)
    add_to_expression!.(expr, opf.mdx' * (pcc.lscc - m[:pd]))
    return expr
end

function make_cut(pccx::ExprCCX, opf::OPFsystem, m::Model)
    expr = JuMP.@expression(m, opf.mgx' * (m[:pg0] - pccx.pgdx))
    add_to_expression!.(expr, opf.mdcx' * m[:pfdc0])
    add_to_expression!.(expr, opf.mdx' * (pccx.lsccx - m[:pd]))
    return expr
end

function initialize_islands(opf::OPFsystem, islands::Vector, island::Int)
    for sep in islands[1:end .!= island]
        if length(sep) > 1
            ix = 1:length(sep)
            bx = [Tuple(opf.mbx[i,:].nzind) for i in axes(opf.mbx,1)]
            ibx = first.(bx) .∈ [ix] .&& last.(bx) .∈ [ix]
            return SCOPF.DCPowerFlow(opf.nodes[ix], opf.branches[ibx], SCOPF.get_nodes_idx(opf.nodes[ix]))
        end
    end
    return nothing
end
