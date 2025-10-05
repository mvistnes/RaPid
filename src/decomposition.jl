""" 
    ModelState type

    Stores the objective function and variables needed for decomposition.
    
    Fields:
    - `obj`: Objective function of the optimization model.
    - `Pᵢ`: Vector of the initial power injections.
    - `Pg`: Vector of the controllable generation.
    - `Pd`: Vector of the fixed power injections.
"""
mutable struct ModelState{TR<:Real}
    obj::JuMP.AffExpr
    Pᵢ::Vector{TR}
    Pg::Vector{TR}
    Pd::Vector{TR}
end

""" Constructor for ModelState type """
function model_state(opf::OPFsystem, m::Model)
    obj = objective_function(m)
    # Pᵢ = get_net_Pᵢ(opf)
    Pᵢ = get_value(m, :p0)
    Pg = get_controllable(opf, m)
    Pd = get_Pd(opf, m) # Fixed injection
    @assert isapprox(Pg, (Pᵢ - Pd); atol=1e-14) string(Pg - (Pᵢ - Pd))
    return ModelState{Float64}(obj, Pᵢ, Pg, Pd)
end

""" 
    Solve the optimization model using decomposition.

    Function creates extra production constraints on generators in the system
    based on the power transfer distribution factors of the generators in the
    system. 

    Arguments:
    - `type::OPF`: Type of SCOPF to solve.
    - `system::System`: Power system to model the SCOPF on.
    - `optimizer`: Optimizer to use for the optimization model.
    - `voll::Vector{<:Number}`: Value of lost load for each demand in the system.
    - `prob::Vector{Float64}`: Probability of failure for each contingency.
    - `contingencies::Vector{<:Component}`: Vector of contingencies to consider.
    - `dist_slack=Float64[]`: Vector of distributed slack factors.
    - `time_limit_sec::Int64=600`: Time limit for the optimization model in seconds.
    - `ramp_minutes=10`: Ramp time for corrective actions in minutes for generators.
    - `ramp_mult=10`: Ramp (post-contingency) cost multiplier for generators.
    - `max_shed=1.0`: Maximum load shedding allowed in MW.
    - `max_curtail=1.0`: Maximum renewable curtailment allowed in MW.
    - `atol=1e-6`: Absolute tolerance for convergence of the decomposition algorithm.
    - `short_term_multi::Real=1.5`: Short term rating multiplier for branches.
    - `long_term_multi::Real=1.2`: Long term rating multiplier for branches.
    - `max_itr=10`: Maximum number of iterations for the decomposition algorithm.
    - `all_post_c=true`: Whether to initialize all post-contingency variables for islands before starting the decomposition.
    - `p_failure=0.0`: Probability of failure for corrective actions.
    - `silent=true`: Whether to suppress output from the optimizer.
    - `debug=false`: Whether to enable debug output.
"""
function run_decomposed_optimization(
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
    all_post_c=true,
    p_failure=0.0,
    silent=true,
    debug=false
)
    case = @timeit timeo "base" Case(opf_base(OPF(true, false, false, false, false), system, optimizer, voll=voll, contingencies=contingencies, prob=prob,
        dist_slack=dist_slack, time_limit_sec=time_limit_sec, ramp_minutes=ramp_minutes, ramp_mult=ramp_mult, max_shed=max_shed, max_curtail=max_curtail,
        short_term_multi=short_term_multi, long_term_multi=long_term_multi, p_failure=p_failure, silent=silent, debug=debug)...)

    total_solve_time = @timeit timeo "base" constrain_branches!(case, 0.0)
    if !type.P & !type.C1 & !type.C2 & !type.C2F
        return case, total_solve_time
    end
    res = @timeit timeo "contingencies" run_decomposed_optimization!(type, case, total_solve_time, atol, max_itr, all_post_c)
    # show(timeo)
    return res
end

run_decomposed_optimization!(type::OPF, case::Case, total_solve_time::Float64, atol=1e-6, max_itr=max(length(case.opf.contingencies), 5), all_post_c=true) =
    run_decomposed_optimization!(type, case.model, case.opf, case.pf, case.oplim, case.brc_up, case.brc_down, case.Pc, case.Pcc, case.Pccx, total_solve_time, atol, max_itr, all_post_c)

"""
    Solve the optimization model using decomposition.

    Arguments:
    - `type::OPF`: Type of SCOPF to solve.
    - `m::Model`: JuMP model of the optimization problem.
    - `opf::OPFsystem`: OPF system data structure.
    - `pf::DCPowerFlow`: DC power flow data structure.
    - `oplim::Oplimits`: Operating limits data structure.
    - `brc_up::Dict{<:Integer, ConstraintRef}`: Dictionary of branch flow upper limit constraints.
    - `brc_down::Dict{<:Integer, ConstraintRef}`: Dictionary of branch flow lower limit constraints.
    - `Pc::Dict{<:Integer, ExprC}`: Dictionary of short-term corrective action expressions.
    - `Pcc::Dict{<:Integer, ExprCC}`: Dictionary of long-term corrective action expressions.
    - `Pccx::Dict{<:Integer, ExprCCX}`: Dictionary of long-term corrective action after corrective failure expressions.
    - `total_solve_time=0.0`: Total time spent solving the optimization model.
    - `atol=1e-6`: Absolute tolerance for convergence of the decomposition algorithm.
    - `max_itr=10`: Maximum number of iterations for the decomposition algorithm.
    - `all_post_c=true`: Whether to initialize all post-contingency variables for islands before starting the decomposition.
"""
function run_decomposed_optimization!(
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
    total_solve_time=0.0,
    atol=1e-6,
    max_itr=10,
    all_post_c=true
)
    assert(type)
    @assert !type.C1 || isempty(Pc)
    @assert !type.C2 || isempty(Pcc)
    @assert !type.C2F || isempty(Pccx)
    
    !is_solved_and_feasible(m) && return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time
    @debug "lower_bound = $(objective_value(m))"

    # Set variables
    state = model_state(opf, m)

    brst = oplim.branch_rating * oplim.short_term_multi
    brlt = oplim.branch_rating * oplim.long_term_multi
    bx = get_bus_idx.(opf.branches, [opf.idx])
    flow = zeros(length(pf.vb_tmp))

    # Initialize post-contingency variables for islands if required
    if all_post_c
        @timeit timeo "pre" begin
        for (i, c_obj) in enumerate(opf.contingencies)
            cont = typesort_component(c_obj, opf)
            # cont  = get_bus_idx(opf.contingencies[c], opf.idx)
            if is_islanded(pf, cont[2], cont[1])
                islands, island, island_b = handle_islands(pf.B, pf.DA, cont[2], cont[1], pf.slack)
                length(islands[island]) < 2 && continue
                type.C1 && init_P!(Pc, opf, oplim, m, state.obj, islands, island, i)
                type.C2 && init_P!(Pcc, opf, oplim, m, state.obj, islands, island, i)
                type.C2F && init_P!(Pccx, opf, oplim, m, state.obj, islands, island, i)
                @debug "Island: Contingency $(string(typeof(c_obj))) $(get_name(c_obj))"
            end
        end
        end

        total_solve_time = update_model!(m, pf, oplim, brc_up, brc_down, state, total_solve_time)
        !is_solved_and_feasible(m) && return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time
    end
    old_obj = -1.0

    # Initialize variables
    ΔPc = zeros(length(state.Pᵢ))
    ΔPcc = zeros(length(state.Pᵢ))
    ΔPccx = zeros(length(state.Pᵢ))
    pre = 0
    corr1 = 0
    corr2 = 0
    corr2f = 0

    olc = Vector{Tuple{Int,Float64}}()
    olcc = Vector{Tuple{Int,Float64}}()
    olccx = Vector{Tuple{Int,Float64}}()
    islands = Vector{Vector{Int}}()
    island = 0
    island_b = Int[]
    ptdf = similar(pf.vn_tmp)

    # Main iteration loop
    cut_added = 1
    for iterations in 1:max_itr
        if cut_added == 0 # loops until no new cuts are added for the contingencies
            # @printf "\nEND: Total solve time %.4f.\n" total_solve_time
            print_cuts(type, pre, corr1, corr2, corr2f)
            return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time
        elseif old_obj == objective_value(m)
            @warn "Reached $(iterations) iterations without change to the solution."
            return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time
        end
        cut_added = 0
        old_obj = objective_value(m)
        @info "\n------------------\nIteration $iterations"

        # Iterate over all contingencies adding restrictions as needed
        for (i, c_obj) in enumerate(opf.contingencies)
            cont = typesort_component(c_obj, opf)

            # Find if the system is islanded by the contingency
            @timeit timeo "split" begin
            if is_islanded(pf, cont[2], cont[1])
                islands, island, island_b = handle_islands(pf.B, pf.DA, cont[2], cont[1], pf.slack)
                inodes = islands[island]
                length(inodes) < 2 && continue
                # initialize_islands(opf, islands, island)
            else
                empty!(islands)
                island = 0
                empty!(island_b)
                inodes = Int[]
            end
            end

            # Calculate the line flows for the contingency, then find overloaded lines
            @timeit timeo "power flow" begin
            if type.P
                olc = calculate_contingency_overload!(flow, brlt, m, pf, state, cont, i, c_obj, inodes, island_b, atol)
            end
            if type.C1
                olc = calculate_contingency_overload!(flow, ΔPc, brst, Pc, opf, m, pf, state, cont, i, c_obj, inodes, island_b, atol)
            end
            if type.C2
                olcc = calculate_contingency_overload!(flow, ΔPcc, brlt, Pcc, opf, m, pf, state, cont, i, c_obj, inodes, island_b, atol)
            end
            if type.C2F
                olccx = calculate_contingency_overload!(flow, ΔPccx, brlt, Pccx, opf, m, pf, state, cont, i, c_obj, inodes, island_b, atol)
            end
            end

            # If any corrective actions are initialized, use them, else fill in zeros
            if !isempty(olc) || !isempty(olcc) || !isempty(olccx) 
                type.C1 && get(Pc, i, 0) == 0 && fill!(ΔPc, 0.0)
                type.C2 && get(Pcc, i, 0) == 0 && fill!(ΔPcc, 0.0)
                type.C2F && get(Pccx, i, 0) == 0 && fill!(ΔPccx, 0.0)
            end

            # NB: Cannot change the model before all data is exctracted!

            # Add restrictions for the contingency if there are overloads
            @timeit timeo "cut" begin
            if !isempty(olc)
                if type.P
                    cut_added, pre = add_cut(ptdf, opf, pf, bx, m, state, brlt, olc, islands, island, island_b, c_obj, cont, cut_added, atol, pre)
                elseif type.C1
                    cut_added, corr1 = add_cut(ptdf, Pc, opf, pf, bx, oplim, m, state, brst, ΔPc, olc, islands, island, island_b, c_obj, cont, i, cut_added, atol, corr1)
                end
            end
            if !isempty(olcc)
                cut_added, corr2 = add_cut(ptdf, Pcc, opf, pf, bx, oplim, m, state, brlt, ΔPcc, olcc, islands, island, island_b, c_obj, cont, i, cut_added, atol, corr2)
            end
            if !isempty(olccx)
                cut_added, corr2f = add_cut(ptdf, Pccx, opf, pf, bx, oplim, m, state, brlt, ΔPccx, olccx, islands, island, island_b, c_obj, cont, i, cut_added, atol, corr2f)
            end
            end
        end

        # Update the model if restrictions were added
        if cut_added > 1
            total_solve_time = update_model!(m, pf, oplim, brc_up, brc_down, state, total_solve_time)
            cut_added = 1
        end
        if !is_solved_and_feasible(m) 
            return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time
        end
    end
    @warn "Reached $(max_itr) iterations without a stable solution."
    return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time
end

""" Solve model and update the power flow object """
@timeit timeo "update model" function update_model!(m::Model, pf::DCPowerFlow, oplim::Oplimits, brc_up::Dict{<:Integer, ConstraintRef}, 
    brc_down::Dict{<:Integer, ConstraintRef}, state::ModelState, total_solve_time::Real
)
    # set_warm_start!(m, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
    set_objective_function(m, state.obj)
    total_solve_time = constrain_branches!(m, pf, oplim, brc_up, brc_down, total_solve_time)
    state.Pᵢ .= get_value(m, :p0)
    @. state.Pg = state.Pᵢ - state.Pd
    return total_solve_time
end

"""
    Calculate the contingency power flow with corrective actions.
    Assummes that island-contingencies have active variables from the pre-procedure.
"""
function calculate_contingency_line_flows!(flow::Vector{<:Real}, ΔP::Vector{<:Real}, P::Dict{<:Integer, T}, opf::OPFsystem, 
    m::Model, pf::DCPowerFlow, Pᵢ::AbstractVector{<:Real}, cont::Tuple{Real, Tuple{Real,Real}}, 
    c::Integer, c_obj::Component, island::Vector, island_b::Vector{<:Integer}; to=0
) where {T}
    if get(P, c, 0) == 0
        return calculate_contingency_line_flows!(flow, m, pf, Pᵢ, cont, c, c_obj, island, island_b)
    end
    @timeit timeo "ΔP" get_ΔP!(ΔP, m, opf, P, c)
    if iszero(ΔP)
        # calculate_line_flows!(flow, pf, cont[2], cont[1], nodes=island, branches=island_b)
        @timeit timeo "Calculate flow" calc_Pline!(flow, pf.vn_tmp, pf, cont[2], cont[1], nodes=island, branches=island_b)
    else
        # calculate_line_flows!(flow, pf, cont[2], cont[1], Pᵢ=(Pᵢ .+ ΔP), nodes=island, branches=island_b)
        @timeit timeo "Calculate flow" calc_Pline!(flow, pf.vn_tmp, pf, cont[2], cont[1], Pᵢ=(Pᵢ .+ ΔP), nodes=island, branches=island_b)
    end
    return flow
end
function calculate_contingency_line_flows!(flow::Vector{<:Real}, ΔP::Vector{<:Real}, P::Dict{<:Integer, T}, opf::OPFsystem, 
    m::Model, pf::DCPowerFlow, Pᵢ::AbstractVector{<:Real}, cont::Tuple{Real,Real},
    c::Integer, c_obj::Component, island::Vector, island_b::Vector{<:Integer}
) where {T}
    if get(P, c, 0) == 0
        return calculate_contingency_line_flows!(flow, m, pf, Pᵢ, cont, c, c_obj, island, island_b)
    end
    @timeit timeo "ΔP" get_ΔP!(ΔP, m, opf, P, c)
    ΔP[cont[2]] -= get_ΔP(m, cont, c_obj, P, c)
    @timeit timeo "Calculate flow" LinearAlgebra.mul!(flow, pf.ϕ, (Pᵢ .+ ΔP))
    return flow
end

"""
    Calculate the contingency power flow without corrective actions.
"""
function calculate_contingency_line_flows!(flow::Vector{<:Real}, 
    m::Model, pf::DCPowerFlow, P::Vector{<:Real}, cont::Tuple{Real, Tuple{Real,Real}}, c::Integer, c_obj::Component,
    island::Vector, island_b::Vector{<:Integer}, val=0.0
)
    # return calculate_line_flows!(flow, pf, cont[2], cont[1], nodes=island, branches=island_b)
    return @timeit timeo "Calculate flow" calc_Pline!(flow, pf.vn_tmp, pf, cont[2], cont[1], nodes=island, branches=island_b)
end
function calculate_contingency_line_flows!(flow::Vector{<:Real}, 
    m::Model, pf::DCPowerFlow, P::Vector{<:Real}, cont::Tuple{Real,Real}, c::Integer, c_obj::StaticInjection,
    island::Vector, island_b::Vector{<:Integer}, val=get_ΔP(m, cont, c_obj)
)
    P[cont[2]] -= val
    @timeit timeo "Calculate flow" LinearAlgebra.mul!(flow, pf.ϕ, P)
    P[cont[2]] += val
    return flow
end

"""
    Calculate the contingency overflow with corrective actions.
"""
function calculate_contingency_overload!(flow::Vector{<:Real}, ΔP::Vector{<:Real}, branch_rating::Vector{<:Real},
    Pc::Dict, opf::OPFsystem, m::Model, pf::DCPowerFlow, state::ModelState, 
    cont::Union{Tuple{Real, Tuple{Real,Real}}, Tuple{Real,Real}}, c::Integer, c_obj::Component, 
    island::Vector, island_b::Vector{<:Integer}, atol::Real=1e-6
)
    calculate_contingency_line_flows!(flow, ΔP, Pc, opf, m, pf, state.Pᵢ, cont, c, c_obj, island, island_b)
    return filter_overload(flow, branch_rating, atol)
end

"""
    Calculate the contingency overflow without corrective actions.
"""
function calculate_contingency_overload!(flow::Vector{<:Real}, branch_rating::Vector{<:Real}, m::Model, 
    pf::DCPowerFlow, state::ModelState, cont::Union{Tuple{Real, Tuple{Real,Real}}, Tuple{Real,Real}}, c::Integer, c_obj::Component,
    island::Vector, island_b::Vector{<:Integer}, atol::Real=1e-6
)
    calculate_contingency_line_flows!(flow, m, pf, state.Pᵢ, cont, c, c_obj, island, island_b)
    return filter_overload(flow, branch_rating, atol)
end

get_ΔP(m::Model, cont::Tuple{Real,Real}, ::Generator) = 
    value(m[:pg0][cont[1]])
get_ΔP(m::Model, cont::Tuple{Real,Real}, ::Generator, Pc::Dict{<:Integer, ExprC}, c::Integer) = 
    value(m[:pg0][cont[1]]) + value(Pc[c].pgu[cont[1]]) - value(Pc[c].pgd[cont[1]])
get_ΔP(m::Model, cont::Tuple{Real,Real}, ::Generator, Pcc::Dict{<:Integer, ExprCC}, c::Integer) = 
    value(m[:pg0][cont[1]]) + value(Pcc[c].pgu[cont[1]]) - value(Pcc[c].pgd[cont[1]])
get_ΔP(m::Model, cont::Tuple{Real,Real}, ::Generator, Pccx::Dict{<:Integer, ExprCCX}, c::Integer) = 
    value(m[:pg0][cont[1]]) - value(Pccx[c].pgdx[cont[1]])

get_ΔP(m::Model, cont::Tuple{Real,Real}, ::StaticLoad) = 
    value(m[:pd0][cont[1]]) - value(m[:ls0][cont[1]])
get_ΔP(m::Model, cont::Tuple{Real,Real}, ::StaticLoad, Pc::Dict{<:Integer, ExprC}, c::Integer) = 
    value(m[:pd0][cont[1]]) - value(Pc[c].lsc[cont[1]])
get_ΔP(m::Model, cont::Tuple{Real,Real}, ::StaticLoad, Pcc::Dict{<:Integer, ExprCC}, c::Integer) = 
    value(m[:pd0][cont[1]]) - value(Pcc[c].lscc[cont[1]])
get_ΔP(m::Model, cont::Tuple{Real,Real}, ::StaticLoad, Pccx::Dict{<:Integer, ExprCCX}, c::Integer) = 
    value(m[:pd0][cont[1]]) - value(Pccx[c].lsccx[cont[1]])

""" Return the contingency power injection change at each node changing ΔP. """
function get_ΔP!(ΔP::Vector{T}, m::Model, opf::OPFsystem, 
    P::Dict{<:Integer, <:ContExpr}, c::Integer
) where {T<:Real}
    x = get(P, c, 0)
    if x != 0
        copy!(ΔP, get_ΔP(x, opf, m))
    else
        fill!(ΔP, zero(T))
    end
    return ΔP
end

""" Return the short term post-contingency power injection change at each node. """
@inline function get_ΔP(pc::ExprC, opf::OPFsystem, m::Model)
    return opf.mgx' * (get_value(m, pc.pgu) - get_value(m, pc.pgd)) - 
        opf.mrx' * get_value(m, pc.prc) + opf.mdx' * get_value(m, pc.lsc)
end

""" Return the long term post-contingency power injection change at each node. """
@inline function get_ΔP(pcc::ExprCC, opf::OPFsystem, m::Model)
    return opf.mgx' * (get_value(m, pcc.pgu) - get_value(m, pcc.pgd)) - 
        opf.mrx' * get_value(m, pcc.prcc) + opf.mdx' * get_value(m, pcc.lscc) +
        opf.mdcx' * get_value(m, pcc.pfdccc)
end

""" Return the long term post-contingency with failure power injection change at each node. """
@inline function get_ΔP(pccx::ExprCCX, opf::OPFsystem, m::Model)
    return -opf.mgx' * get_value(m, pccx.pgdx) + opf.mdx' * get_value(m, pccx.lsccx)
end

function add_overload_expr!(m::Model, expr::JuMP.AbstractJuMPScalar, ol::Real, rate::Real, type::String, id::Integer,
    c_obj::Component, opf::OPFsystem, i::Integer, atol::Real
)
    @info @sprintf "%s %d: Contingency %s %s; overload on %s of %.6f" type id string(typeof(c_obj)) c_obj.name opf.branches[i].name ol
    @debug "Cut added: $(sprint_expr(expr,atol))\n"
    if ol < 0
        return JuMP.@constraint(m, expr >= -rate)
    else
        return JuMP.@constraint(m, expr <= rate)
    end
end

""" Finding and adding preventive cuts to mitigate overloads from one contingency """
function add_cut(ptdf::AbstractVector, opf::OPFsystem, pf::DCPowerFlow, bx::Vector, m::Model, state::ModelState, rate::Vector{<:Real}, overloads::Vector{<:Tuple{Integer,Real}},
    islands::Vector, island::Integer, island_b::Vector, c_obj::Component, cont, cut_added::Integer, atol::Real, id::Integer
)
    for (i, ol) in overloads
        if is_islanded(pf, cont[2], cont[1])
            calc_ptdf_vec!(ptdf, pf.DA, pf.B, cont[2], cont[1], pf.slack, islands[island], island_b, 1)
            # ptdf = calc_ptdf(pf, cont[2], cont[1], islands, island, island_b)
        else
            # ptdf = calc_ptdf(pf, cont[2], cont[1])
            calc_ptdf_vec!(ptdf, pf, cont[1], cont[2][1], cont[2][2], i, bx[i][1], bx[i][2])
        end
        expr = AffExpr() #+ sum(ptdf .* state.Pᵢ)
        for j in axes(ptdf,1)
            add_to_expression!(expr, ptdf[j], m[:inj_p0][j])
        end

        id += 1
        add_overload_expr!(m, expr, ol, rate[i], "Pre", id, c_obj, opf, i, atol)
        cut_added = 2
    end
    return cut_added, id
end

""" Finding and adding corrective cuts to mitigate overloads from one contingency """
function add_cut(ptdf::AbstractVector, P::Dict{<:Integer, <:ContExpr}, opf::OPFsystem, pf::DCPowerFlow, bx::Vector, oplim::Oplimits, m::Model, state::ModelState, rate::Vector{<:Real}, ΔP::Vector{<:Real},
    overloads::Vector{<:Tuple{Integer,Real}}, islands::Vector, island::Integer, island_b::Vector, 
    c_obj::Component, cont, c::Integer, cut_added::Integer, atol::Real, id::Integer
)
    p = get(P, c, 0)
    if p == 0  # If the contingency is not run before, a set of corrective variables is added
        init_P!(P, opf, oplim, m, state.obj, islands, island, c)
        p = get(P, c, 0)
    end
    for (i, ol) in overloads
        if is_islanded(pf, cont[2], cont[1])
            calc_ptdf_vec!(ptdf, pf.DA, pf.B, cont[2], cont[1], pf.slack, islands[island], island_b, 1)
            # ptdf = calc_ptdf(pf, cont[2], cont[1], islands, island, island_b)
        else
            # ptdf = calc_ptdf(pf, cont[2], cont[1])
            calc_ptdf_vec!(ptdf, pf, cont[1], cont[2][1], cont[2][2], i, bx[i][1], bx[i][2])
        end
        expr = make_cut(p, opf, m)
        # add_to_expression!.(expr, state.Pᵢ + ΔP)
        expr2 = JuMP.@expression(m, sum((ptdf .* expr)))
        # println("c_line", cont[1], " (", cont[2][1], ", ", cont[2][2], ") o_line", i, " (", bx[i][1], ", ", bx[i][2], ")")
        id += 1
        add_overload_expr!(m, expr2, ol, rate[i], get_name(p), id, c_obj, opf, i, atol)
        cut_added = 2
    end
    return cut_added, id
end

function make_cut(pc::ExprC, opf::OPFsystem, m::Model)
    expr = JuMP.@expression(m, opf.mgx' * (m[:pg0] + pc.pgu - pc.pgd))
    add_to_expression!.(expr, opf.mrx' * (m[:pr] - pc.prc))
    add_to_expression!.(expr, opf.mdcx' * m[:pfdc0])
    add_to_expression!.(expr, opf.mdx' * (pc.lsc - m[:pd]))
    return expr
end

function make_cut(pcc::ExprCC, opf::OPFsystem, m::Model)
    expr = JuMP.@expression(m, opf.mgx' * (m[:pg0] + pcc.pgu - pcc.pgd))
    add_to_expression!.(expr, opf.mrx' * (m[:pr] - pcc.prcc))
    add_to_expression!.(expr, opf.mdcx' * pcc.pfdccc)
    add_to_expression!.(expr, opf.mdx' * (pcc.lscc - m[:pd]))
    return expr
end

function make_cut(pccx::ExprCCX, opf::OPFsystem, m::Model)
    expr = JuMP.@expression(m, opf.mgx' * (m[:pg0] - pccx.pgdx))
    add_to_expression!.(expr, opf.mrx' * (m[:pr] - pccx.prccx))
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