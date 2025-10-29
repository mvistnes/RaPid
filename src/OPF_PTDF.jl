""" 
    Operational limits type 
        
    Parameters:
    - `ramp_mult`: Multiplier for ramping limits.
    - `ramp_minutes`: Ramping time in minutes.
    - `p_failure`: Probability of failure of corrective actions.
    - `branch_rating`: Thermal rating of branches.
    - `short_term_multi`: Short-term multiplier for branch ratings.
    - `long_term_multi`: Long-term multiplier for branch ratings.
    - `pg_lim_min`: Minimum active power limits for generators.
    - `pg_lim_max`: Maximum active power limits for generators.
    - `rampup`: Ramp-up limits for generators.
    - `rampdown`: Ramp-down limits for generators.
    - `pr_lim`: Active power limits for renewables.
    - `max_curtail`: Maximum curtailment limits for renewables.
    - `dc_lim_min`: Minimum (negative) active power limits for DC branches.
    - `dc_lim_max`: Maximum active power limits for DC branches.
    - `pd_lim`: Active power limits for demands.
    - `max_shed`: Maximum load shedding limits.
"""
struct Oplimits{TR<:Real}
    ramp_mult::TR
    ramp_minutes::TR
    p_failure::TR

    branch_rating::Vector{TR}
    short_term_multi::Union{TR,Vector{TR}}
    long_term_multi::Union{TR,Vector{TR}}
    pg_lim_min::Vector{TR}
    pg_lim_max::Vector{TR}
    rampup::Vector{TR}
    rampdown::Vector{TR}
    pr_lim::Vector{TR}
    max_curtail::Union{TR,Vector{TR}}
    dc_lim_min::Vector{TR}
    dc_lim_max::Vector{TR}
    pd_lim::Vector{TR}
    max_shed::Union{TR,Vector{TR}}
end

""" 
    Constructor for Oplimits 
    
    Parameters:
    - `opf`: An instance of `OPFsystem` containing the power system data.
    - `max_shed`: Maximum load shedding limits, can be a single value or a vector.
    - `max_curtail`: Maximum curtailment limits for renewables, can be a single value or a vector.
    - `ramp_mult`: Multiplier for ramping limits.
    - `ramp_minutes`: Ramping time in minutes.
    - `p_failure`: Probability of failure of corrective actions.
    - `short_term_multi`: Short-term multiplier for branch ratings, can be a single value or a vector.
    - `long_term_multi`: Long-term multiplier for branch ratings, can be a single value or a vector.
"""
function oplimits(
    opf::OPFsystem,
    max_shed::Union{TR,Vector{TR}},
    max_curtail::Union{TR,Vector{TR}},
    ramp_mult::TR,
    ramp_minutes::TR,
    p_failure::TR,
    short_term_multi::Union{TR,Vector{TR}},
    long_term_multi::Union{TR,Vector{TR}}
) where {TR<:Real}
    branch_rating::Vector{TR} = get_rating.(opf.branches)
    (pg_lim_min::Vector{TR}, pg_lim_max::Vector{TR}) = split_pair(get_active_power_limits.(opf.ctrl_generation))
    (rampup::Vector{TR}, rampdown::Vector{TR}) = split_pair(get_ramp_limits.(opf.ctrl_generation))
    pr_lim::Vector{TR} = get_active_power.(opf.renewables)
    (dc_lim_min::Vector{TR}, dc_lim_max::Vector{TR}) = split_pair(get_active_power_limits_from.(opf.dc_branches))
    pd_lim::Vector{TR} = get_active_power.(opf.demands)
    if any(pg_lim_min .< 0.0)
        @error "Negative pg_lim_min, enable min limits on generators!"
    end
    for v in (branch_rating, pg_lim_min, pg_lim_max, rampup, rampdown, pr_lim, dc_lim_min, dc_lim_max, pd_lim)
        for i in eachindex(v)
            if isnan(v[i])
                v[i] = zero(eltype(v))
            end
        end
    end
    return Oplimits{TR}(ramp_mult, ramp_minutes, p_failure, branch_rating, short_term_multi, long_term_multi, 
        zeros(length(pg_lim_max)), pg_lim_max, rampup, rampdown, pr_lim, max_curtail, dc_lim_min, dc_lim_max, pd_lim, max_shed)
end

"""
    Case{TR,TI} <: AbstractCase

    A mutable struct representing a case for OPF analysis.
    
    Parameters:
    - `model`: The optimization model.
    - `opf`: The power system data.
    - `pf`: The DC power flow object.
    - `oplim`: Operational limits for the case.
    - `brc_up`: Dictionary of upper branch constraints.
    - `brc_down`: Dictionary of lower branch constraints.
    - `Pc`: Dictionary of short-term corrective actions.
    - `Pcc`: Dictionary of long-term corrective actions.
    - `Pccx`: Dictionary of long-term corrective actions after failure of actions.
"""
mutable struct Case{TR<:Real,TI<:Integer}
    model::Model
    opf::OPFsystem{TR}
    pf::DCPowerFlow{TI,TR}
    oplim::Oplimits{TR}
    brc_up::Dict{TI,ConstraintRef}
    brc_down::Dict{TI,ConstraintRef}
    Pc::Dict{TI,ExprC}
    Pcc::Dict{TI,ExprCC}
    Pccx::Dict{TI,ExprCCX}
end

""" 
    Initialize an OPF of a power system 

    This function sets up the optimization model for the OPF problem, including variables, constraints, and objectives.

    Parameters:
    - `type`: The type of OPF to be solved.
    - `system`: Power system data.
    - `optimizer`: The optimizer to be used for solving the OPF.
    - `voll`: Value of lost load.
    - `contingencies`: List of contingencies to be considered.
    - `prob`: Probability of each contingency.
    - `dist_slack`: Distribution of generator slack.
    - `time_limit_sec`: Time limit for the optimization.
    - `ramp_minutes`: Ramping time in minutes.
    - `ramp_mult`: Multiplier for ramping limits.
    - `max_shed`: Maximum load shedding limits.
    - `max_curtail`: Maximum curtailment limits for renewables.
    - `short_term_multi`: Short-term multiplier for branch ratings.
    - `long_term_multi`: Long-term multiplier for branch ratings.
    - `p_failure`: Probability of failure of corrective actions.
    - `silent`: Suppresses output messages.
    - `debug`: Enables debug mode.
    """
function opf_base(type::OPF, system::System, optimizer;
    voll=Float64[],
    contingencies=Component[],
    prob=Float64[],
    dist_slack=Float64[],
    time_limit_sec::Int64=600,
    ramp_minutes::Real=10.0,
    ramp_mult::Real=10.0,
    max_shed=1.0,
    max_curtail=1.0,
    short_term_multi=1.5,
    long_term_multi=1.2,
    p_failure=0.0,
    silent=true,
    debug=false
)
    @assert type.Base
    assert(type)
    @assert !type.C2F | !iszero(p_failure)
    m = create_model(optimizer, time_limit_sec=time_limit_sec, silent=silent, debug=debug)
    opf = isempty(voll) ? opfsystem(system) : opfsystem(system, voll, contingencies, prob, ramp_mult)
    pf = isempty(dist_slack) ? DCPowerFlow(opf.nodes, opf.branches, opf.idx) : DCPowerFlow(opf.nodes, opf.branches, opf.idx, slack_array=dist_slack / sum(dist_slack))
    
    brc_up = Dict{Int64,ConstraintRef}()
    brc_down = Dict{Int64,ConstraintRef}()
    Pc = Dict{Int64,ExprC}()
    Pcc = Dict{Int64,ExprCC}()
    Pccx = Dict{Int64,ExprCCX}()

    oplim = oplimits(opf, max_shed, max_curtail, ramp_mult, ramp_minutes, p_failure, short_term_multi, long_term_multi)

    @variables(m, begin
        p0[n in 1:length(opf.nodes)]
        # active power injection on each node
        oplim.pg_lim_min[g] <= pg0[g in 1:length(opf.ctrl_generation)] <= oplim.pg_lim_max[g]
        # active power variables for the generators
        oplim.dc_lim_min[l] <= pfdc0[l in 1:length(opf.dc_branches)] <= oplim.dc_lim_max[l]
        # power flow on DC branches
        pd[d in 1:length(opf.demands)]
        0.0 <= ls0[d in 1:length(opf.demands)] <= oplim.pd_lim[d] * (typeof(oplim.max_shed) <: Real ? 1.0 : oplim.max_shed[d])
        # demand curtailment variables
        pr[d in 1:length(opf.renewables)]
        0.0 <= pr0[d in 1:length(opf.renewables)] <= oplim.pr_lim[d] * (typeof(oplim.max_curtail) <: Real ? 1.0 : oplim.max_curtail[d])
        # renewable curtailment variables
    end)

    @objective(m, Min, sum(c.var * g for (c, g) in zip(opf.cost_ctrl_gen, pg0)) + sum(opf.voll' * ls0))
    # @objective(m, Min, sum(c[1] * g^2 + c[2] * g for (c, g) in zip(opf.cost_ctrl_gen, pg0)) + sum(opf.voll' * ls0) + sum(pr0 * 1))

    JuMP.fix.(pd, oplim.pd_lim)
    JuMP.fix.(pr, oplim.pr_lim)
    if typeof(oplim.max_shed) <: Real
        @constraint(m, sum_max_shed, sum(ls0) <= oplim.max_shed)
    end

    @expression(m, inj_p0, opf.mgx' * pg0)
    add_to_expression!.(inj_p0, opf.mrx' * (pr - pr0))
    add_to_expression!.(inj_p0, opf.mdcx' * pfdc0)
    add_to_expression!.(inj_p0, opf.mdx' * (ls0 - pd))
    @constraint(m, injected_power, inj_p0 .== p0)

    # add_branch_constraints!(m, pf.ϕ, p0, oplim.branch_rating)
    @expression(m, balance, sum(pg0, init=0.0))
    add_to_expression!.(balance, pr)
    add_to_expression!.(balance, -1, pr0)
    add_to_expression!.(balance, -1, pd)
    add_to_expression!.(balance, ls0)
    @constraint(m, power_balance, balance == 0.0)

    if any([type.P, type.C1, type.C2, type.C2F])
        return add_all_contingencies!(type, opf, oplim, m, pf, brc_up, brc_down, Pc, Pcc, Pccx)
    end
    return m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx
end

""" 
    Add all branch constraints to the model.
    
    Parameters:
    - `m`: Optimization model.
    - `ptdf`: PTDF matrix.
    - `p`: Node active power injection.
    - `brc_up`: Upper branch constraints.
    - `brc_down`: Lower branch constraints.
    - `rating`: Branch thermal ratings.
"""
function add_branch_constraints!(m::Model, ptdf::AbstractMatrix{<:Real}, p::AbstractVector{VariableRef}, 
    brc_up::Dict{<:Integer, ConstraintRef}, brc_down::Dict{<:Integer, ConstraintRef}, rating::AbstractVector{<:Real}
)
    for i = axes(ptdf,1)
        # ptdf0 = GenericAffExpr(0.0, Pair.(p, ptdf[i,:])) 
        ptdf0 = @expression(m, AffExpr())
        for j in axes(ptdf,2)
            add_to_expression!(ptdf0, ptdf[i,j], p[j])
        end
        brc_down[i] = @constraint(m, ptdf0 + rating[i] >= 0.0)
        brc_up[i] = @constraint(m, ptdf0 - rating[i] <= 0.0)
    end
    return mod
end

""" 
    Add a branch constraint to the model.
    
    Parameters:
    - `m`: Optimization model.
    - `pf`: DC power flow object.
    - `p`: Node active power injection.
    - `brc_up`: Upper branch constraints.
    - `brc_down`: Lower branch constraints.
    - `branch`: Index of the branch constraint to add to the model.
    - `rating`: Branch thermal ratings.
"""
function add_branch_constraint!(m::Model, pf::DCPowerFlow, p::AbstractVector{VariableRef}, 
    brc_up::Dict{<:Integer, ConstraintRef}, brc_down::Dict{<:Integer, ConstraintRef}, branch::Integer, rating::Real
)
    # ptdf = calc_ptdf_vec(pf, branch)
    # ptdf0 = GenericAffExpr(0.0, Pair.(p, ptdf[i,:])) 
    ptdf0 = @expression(m, AffExpr())
    add_to_expression!.(ptdf0, pf.ϕ[branch, :], p)
    brc_down[branch] = @constraint(m, ptdf0 + rating >= 0.0)
    brc_up[branch] = @constraint(m, ptdf0 - rating <= 0.0)
    return mod
end

""" 
    Iteratively add branch constraints to the model until all overloaded branches are constrained.
    
    Parameters:
    - `m`: Optimization model.
    - `pf`: DC power flow object.
    - `oplim`: Operational limits for the case.
    - `brc_up`: Upper branch constraints.
    - `brc_down`: Lower branch constraints.
    - `total_solve_time`: Total solve time for the model.
    - `atol`: Absolute tolerance for branch overload detection.
"""
function constrain_branches!(m::Model, pf::DCPowerFlow, oplim::Oplimits, brc_up::Dict{<:Integer, ConstraintRef}, brc_down::Dict{<:Integer, ConstraintRef}, 
    total_solve_time::Real, atol::Real=1e-6
)
    # if !is_solved_and_feasible(m)
        # Note: While using a direct_model, this check fails after the model is modified for some solvers
        total_solve_time = update_model!(m, pf, total_solve_time)
    # end
    while true
        ol_br = find_overloaded_branches(pf.F, oplim.branch_rating, atol)
        isempty(ol_br) && break 
        !is_solved_and_feasible(m) && break
        for br in ol_br
            add_branch_constraint!(m, pf, m[:p0], brc_up, brc_down, br, oplim.branch_rating[br])
            @info "Branch $br added"
        end
        total_solve_time = update_model!(m, pf, total_solve_time)
    end
    return total_solve_time
end
""" 
    Iteratively add branch constraints to the model until all overloaded branches are constrained.
    
    Parameters:
    - `case`: A `Case` object.
    - `total_solve_time`: Total solve time for the model.
    - `atol`: Absolute tolerance for branch overload detection.
"""
constrain_branches!(case::Case, total_solve_time::Real, atol::Real=1e-6) = 
    constrain_branches!(case.model, case.pf, case.oplim, case.brc_up, case.brc_down, total_solve_time, atol)

""" 
    Solve model and update the power flow object 

    Parameters:
    - `m`: Optimization model.
    - `pf`: DC power flow object.
    - `total_solve_time`: Total solve time for the model.
"""
function update_model!(m::Model, pf::DCPowerFlow, total_solve_time::Real)
    solve_model!(m)
    total_solve_time += solve_time(m)
    Pᵢ = get_value(m, :p0)
    run_pf!(pf, Pᵢ)
    calc_Pline!(pf)
    return total_solve_time
end

"""
    Add all contingencies to the model and update the power flow object.

    Parameters:
    - `type`: The type of OPF to be solved.
    - `opf`: The power system data.
    - `oplim`: Operational limits for the case.
    - `m`: Optimization model.
    - `pf`: DC power flow object.
    - `brc_up`: Upper branch constraints.
    - `brc_down`: Lower branch constraints.
    - `Pc`: Short-term corrective actions.
    - `Pcc`: Long-term corrective actions.
    - `Pccx`: Long-term corrective actions after failure of actions.
"""
function add_all_contingencies!(type::OPF, opf::OPFsystem, oplim::Oplimits, m::Model,
    pf::DCPowerFlow, brc_up::Dict{<:Integer, ConstraintRef}, brc_down::Dict{<:Integer, ConstraintRef}, 
    Pc::Dict{<:Integer,ExprC}, Pcc::Dict{<:Integer,ExprCC}, Pccx::Dict{<:Integer,ExprCCX}
)
    obj = objective_function(m)
    for (i, c_obj) in enumerate(opf.contingencies)
        cont = typesort_component(c_obj, opf)
        if is_islanded(pf, cont[2], cont[1])
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont[2], cont[1], pf.slack)
            length(islands[island]) < 2 && continue
            ptdf = calc_ptdf(pf, cont[2], cont[1], islands, island, island_b)
            set_tol_zero!(ptdf)
        else
            ptdf = calc_ptdf(pf, cont[2], cont[1])
            islands = Vector{Vector{Int64}}[]
            island = 0
        end
        type.P && add_contingency!(opf, pf, oplim, m, brc_up, brc_down, ptdf, i)
        type.C1 && add_contingency!(Pc, opf, pf, oplim, m, brc_up, brc_down, obj, islands, island, ptdf, i)
        type.C2 && add_contingency!(Pcc, opf, pf, oplim, m, brc_up, brc_down, obj, islands, island, ptdf, i)
        type.C2F && add_contingency!(Pccx, opf, pf, oplim, m, brc_up, brc_down, obj, islands, island, ptdf, i)

        @debug "Contingency $(get_name(c_obj)) is added"
    end
    set_objective_function(m, obj)
    return m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx
end

"""
    Add all contingencies to the model and update the power flow object.

    Parameters:
    - `type`: The type of OPF to be solved.
    - `case`: A `Case` object containing the power system data and model.
"""
add_all_contingencies!(type::OPF, case::Case) = 
    add_all_contingencies!(type, case.opf, case.oplim, case.model, case.pf, case.brc_up, case.brc_down, case.Pc, case.Pcc, case.Pccx)


"""
    Add a contingency to the model for the base case.

    Parameters:
    - `opf`: The power system data.
    - `pf`: DC power flow object.
    - `oplim`: Operational limits for the case.
    - `m`: Optimization model.
    - `brc_up`: Upper branch constraints.
    - `brc_down`: Lower branch constraints.
    - `ptdf`: PTDF matrix.
    - `c`: Contingency index.
"""
function add_contingency!(opf::OPFsystem, pf::DCPowerFlow, oplim::Oplimits, m::Model, 
    brc_up::Dict{<:Integer, ConstraintRef}, brc_down::Dict{<:Integer, ConstraintRef}, 
    ptdf::AbstractMatrix{<:Real}, c::Integer
)
    p = @variable(m, [n in 1:length(opf.nodes)], base_name = "p"*string(c))
    @constraint(m, m[:inj_p0] .== p)
    add_branch_constraints!(m, ptdf, p, brc_up, brc_down, oplim.branch_rating * oplim.long_term_multi)
    # add_branch_constraints!(m, opf, pf, p, oplim.branch_rating * oplim.long_term_multi, c)
end

""" 
    Add a contingency to the model for short-term corrective actions.

    Parameters:
    - `Pc`: Dictionary of short-term corrective actions.
    - `opf`: The power system data.
    - `pf`: DC power flow object.
    - `oplim`: Operational limits for the case.
    - `m`: Optimization model.
    - `brc_up`: Upper branch constraints.
    - `brc_down`: Lower branch constraints.
    - `obj`: Objective function.
    - `islands`: List of islands in the power system.
    - `island`: Index of the island to add.
    - `ptdf`: PTDF matrix.
    - `c`: Contingency index.
"""
function add_contingency!(Pc::Dict{<:Integer,ExprC}, opf::OPFsystem, pf::DCPowerFlow, oplim::Oplimits, m::Model, 
    brc_up::Dict{<:Integer, ConstraintRef}, brc_down::Dict{<:Integer, ConstraintRef}, obj::AbstractJuMPScalar, islands::Vector,
    island::Integer, ptdf::AbstractMatrix{<:Real}, c::Integer
)
    pc = @variable(m, [n in 1:length(opf.nodes)], base_name = "pc"*string(c))
    pgu, pgd, prc, lsc = init_P!(Pc, opf, oplim, m, obj, islands, island, c)
    Pc[c].pc = pc

    inj_pc = @expression(m, opf.mgx' * (m[:pg0] + pgu - pgd))
    add_to_expression!.(inj_pc, opf.mrx' * (m[:pr] - prc))
    add_to_expression!.(inj_pc, opf.mdcx' * m[:pfdc0])
    add_to_expression!.(inj_pc, opf.mdx' * (lsc - m[:pd]))
    @constraint(m, inj_pc .== pc)

    add_branch_constraints!(m, ptdf, pc, brc_up, brc_down, oplim.branch_rating * oplim.short_term_multi)
    # add_branch_constraints!(m, opf, pf, pc, oplim.branch_rating * oplim.short_term_multi, c)
end

""" 
    Add a contingency to the model for long-term corrective actions.

    Parameters:
    - `Pcc`: Dictionary of long-term corrective actions.
    - `opf`: The power system data.
    - `pf`: DC power flow object.
    - `oplim`: Operational limits for the case.
    - `m`: Optimization model.
    - `brc_up`: Upper branch constraints.
    - `brc_down`: Lower branch constraints.
    - `obj`: Objective function.
    - `islands`: List of islands in the power system.
    - `island`: Index of the island to add.
    - `ptdf`: PTDF matrix.
    - `c`: Contingency index.
"""
function add_contingency!(Pcc::Dict{<:Integer,ExprCC}, opf::OPFsystem, pf::DCPowerFlow, oplim::Oplimits, m::Model, 
    brc_up::Dict{<:Integer, ConstraintRef}, brc_down::Dict{<:Integer, ConstraintRef}, obj::AbstractJuMPScalar, islands::Vector,
    island::Integer, ptdf::AbstractMatrix{<:Real}, c::Integer
)
    pcc = @variable(m, [n in 1:length(opf.nodes)], base_name = "pcc"*string(c))
    pgu, pgd, pfdccc, prcc, lscc = init_P!(Pcc, opf, oplim, m, obj, islands, island, c)
    Pcc[c].pcc = pcc

    inj_pcc = @expression(m, opf.mgx' * (m[:pg0] + pgu - pgd))
    add_to_expression!.(inj_pcc, opf.mrx' * (m[:pr] - prcc))
    add_to_expression!.(inj_pcc, opf.mdcx' * pfdccc)
    add_to_expression!.(inj_pcc, opf.mdx' * (lscc - m[:pd]))
    @constraint(m, inj_pcc .== pcc)

    add_branch_constraints!(m, ptdf, pcc, brc_up, brc_down, oplim.branch_rating * oplim.long_term_multi)
    # add_branch_constraints!(m, opf, pf, pcc, oplim.branch_rating * oplim.long_term_multi, c)
end

""" 
    Add a contingency to the model for long-term corrective actions after failure of actions.

    Parameters:
    - `Pccx`: Dictionary of long-term corrective actions after failure of actions.
    - `opf`: The power system data.
    - `pf`: DC power flow object.
    - `oplim`: Operational limits for the case.
    - `m`: Optimization model.
    - `brc_up`: Upper branch constraints.
    - `brc_down`: Lower branch constraints.
    - `obj`: Objective function.
    - `islands`: List of islands in the power system.
    - `island`: Index of the island to add.
    - `ptdf`: PTDF matrix.
    - `c`: Contingency index.
"""
function add_contingency!(Pccx::Dict{<:Integer,ExprCCX}, opf::OPFsystem, pf::DCPowerFlow, oplim::Oplimits, m::Model, 
    brc_up::Dict{<:Integer, ConstraintRef}, brc_down::Dict{<:Integer, ConstraintRef}, obj::AbstractJuMPScalar, islands::Vector,
    island::Integer, ptdf::AbstractMatrix{<:Real}, c::Integer
)
    pccx = @variable(m, [n in 1:length(opf.nodes)], base_name = "pccx"*string(c))
    pgd, prcc, lscc = init_P!(Pccx, opf, oplim, m, obj, islands, island, c)
    Pccx[c].pccx = pccx

    inj_pccx = @expression(m, opf.mgx' * (m[:pg0] - pgd))
    add_to_expression!.(inj_pccx, opf.mrx' * (m[:pr] - prcc))
    add_to_expression!.(inj_pccx, opf.mdcx' * m[:pfdc0])
    add_to_expression!.(inj_pccx, opf.mdx' * (lscc - m[:pd]))
    @constraint(m, inj_pccx .== pccx)

    add_branch_constraints!(m, ptdf, pccx, brc_up, brc_down, oplim.branch_rating * oplim.long_term_multi)
    # add_branch_constraints!(m, opf, pf, pccx, oplim.branch_rating * oplim.long_term_multi, c)
end

""" 
    Initialize variables for the short-term corrective state.

    Parameters:
    - `Pc`: Dictionary of short-term corrective actions.
    - `opf`: The power system data.
    - `oplim`: Operational limits for the case.
    - `m`: Optimization model.
    - `obj`: Objective function.
    - `islands`: List of islands in the power system.
    - `island`: Index of the island to add.
    - `c`: Contingency index.
"""
function init_P!(Pc::Dict{<:Integer,ExprC}, opf::OPFsystem, oplim::Oplimits, m::Model, obj::AbstractJuMPScalar, 
    islands::Vector, island::Integer, c::Integer
)
    pgu = @variable(m, [g in 1:length(opf.ctrl_generation)], base_name = "pguc"*string(c),
        lower_bound = 0.0, upper_bound = oplim.rampup[g] * 0.0)
    # active power variables for the generators in contingencies ramp up 
    pgd = @variable(m, [g in 1:length(opf.ctrl_generation)], base_name = "pgdc"*string(c),
        lower_bound = 0.0) #, upper_bound = rampdown[g] * 0.0)
    # and ramp down
    prc = @variable(m, [d in 1:length(opf.renewables)], base_name = "prc"*string(c),
        lower_bound = 0.0, upper_bound = oplim.pr_lim[d] * (typeof(oplim.max_curtail) <: Real ? 1.0 : oplim.max_curtail[d]))
    lsc = @variable(m, [d in 1:length(opf.demands)], base_name = "lsc"*string(c),
        lower_bound = 0.0, upper_bound = oplim.pd_lim[d] * (typeof(oplim.max_shed) <: Real ? 1.0 : oplim.max_shed[d]))
    Pc[c] = ExprC(VariableRef[], pgu, pgd, prc, lsc)

    p_survive = 1.0 - oplim.p_failure
    add_to_expression!(obj, opf.prob[c], sum(opf.voll' * lsc))
    add_to_expression!(obj, opf.prob[c], sum(oplim.ramp_mult * 1 * prc))
    for (cost, gu, gd) in zip(opf.cost_ctrl_gen, pgu, pgd)
        add_to_expression!(obj, opf.prob[c] * p_survive * cost.ramp, gu)
        add_to_expression!(obj, opf.prob[c] * p_survive * cost.ramp, gd)
        # add_to_expression!.(obj, opf.prob[c], p_survive * oplim.ramp_mult * (cost[1] * g^2 + cost[2] * g))
    end

    # Add new constraints that limit the corrective variables within operating limits
    balance_pc = @expression(m, sum(lsc))
    add_to_expression!.(balance_pc, -1, prc)
    add_to_expression!.(balance_pc, pgu)
    add_to_expression!.(balance_pc, -1, pgd)
    @constraint(m, balance_pc == 0.0)
    if isempty(islands)
        @constraint(m, m[:pg0] .+ pgu .- pgd .>= oplim.pg_lim_min)
        # @constraint(m, m[:pg0] .+ pgu .- pgd .<= oplim.pg_lim_max)
        if typeof(oplim.max_shed) <: Real
            @constraint(m, sum(lsc) <= oplim.max_shed)
        end
    else
        itr = length(islands[island]) < 2 ? Int64[] : islands[island]
        for n in itr
            for g = opf.mgx[:,n].nzind
                @constraints(m, begin
                    m[:pg0][g] + pgu[g] - pgd[g] >= oplim.pg_lim_min[g]
                    # m[:pg0][g] + pgu[g] - pgd[g] <= oplim.pg_lim_max[g]
                end)
            end
        end
        if typeof(oplim.max_shed) <: Real
            expr = AffExpr()
            for n in itr, d in opf.mdx[:,n].nzind
                add_to_expression!(expr, lsc[d])
            end
            @constraint(m, expr <= oplim.max_shed)
        end
        itr = isempty(itr) ? (1:length(opf.nodes)) : islands[1:end.!=island]
        for in_vec in itr
            for n in in_vec
                for g in opf.mgx[:,n].nzind
                    set_upper_bound(pgu[g], oplim.pg_lim_max[g] - oplim.pg_lim_min[g])
                    @constraint(m, m[:pg0][g] + pgu[g] - pgd[g] == 0.0)
                end
                fix!(prc, oplim.pr_lim, opf.mrx[:,n].nzind)
                fix!(lsc, oplim.pd_lim, opf.mdx[:,n].nzind)
            end
        end
    end
    return pgu, pgd, prc, lsc
end

"""
    Initialize variables for the long-term corrective state.

    Parameters:
    - `Pcc`: Dictionary of long-term corrective actions.
    - `opf`: The power system data.
    - `oplim`: Operational limits for the case.
    - `m`: Optimization model.
    - `obj`: Objective function.
    - `islands`: List of islands in the power system.
    - `island`: Index of the island to add.
    - `c`: Contingency index.
"""
function init_P!(Pcc::Dict{<:Integer,ExprCC}, opf::OPFsystem, oplim::Oplimits, m::Model, obj::AbstractJuMPScalar, 
    islands::Vector, island::Integer, c::Integer
)
    # Add corrective variables
    pgu = @variable(m, [g in 1:length(opf.ctrl_generation)], base_name = @sprintf("pgucc%s", c),
        lower_bound = 0.0, upper_bound = oplim.rampup[g] * oplim.ramp_minutes)
    # active power variables for the generators in contingencies ramp up 
    pgd = @variable(m, [g in 1:length(opf.ctrl_generation)], base_name = @sprintf("pgdcc%s", c),
        lower_bound = 0.0) #, upper_bound = rampdown[g] * ramp_minutes)
    # and ramp down
    pfdccc = @variable(m, [d in 1:length(opf.dc_branches)], base_name = @sprintf("pfdccc%s", c),
        lower_bound = oplim.dc_lim_min[d], upper_bound = oplim.dc_lim_max[d])
    prcc = @variable(m, [r in 1:length(opf.renewables)], base_name = @sprintf("prcc%s", c),
        lower_bound = 0.0, upper_bound = oplim.pr_lim[r] * (typeof(oplim.max_curtail) <: Real ? 1.0 : oplim.max_curtail[d]))
    lscc = @variable(m, [d in 1:length(opf.demands)], base_name = @sprintf("lscc%s", c),
        lower_bound = 0.0, upper_bound = oplim.pd_lim[d] * (typeof(oplim.max_shed) <: Real ? 1.0 : oplim.max_shed[d]))
    # load curtailment variables in in contingencies
    Pcc[c] = ExprCC(VariableRef[], pgu, pgd, pfdccc, prcc, lscc)

    p_survive = 1.0 - oplim.p_failure
    # Extend the objective with the corrective variables
    add_to_expression!.(obj, opf.prob[c],
        # (1.0 - p_failure) * (sum(opf.voll' * lscc) + sum(60 * (pgu + pgd)) # + # TODO: remove 60 and uncomment next lines for non-4-area analysis!!!!
        p_survive * (sum(opf.voll' * lscc) # + sum(opf.cost_ctrl_gen' * ramp_mult * (pgu + pgd))
        # (sum(opf.voll' * lscc) +
        # sum(opf.cost_ctrl_gen' * ramp_mult * pgu) # +
        # sum(opf.cost_ctrl_gen' * pgd)
        ))
    add_to_expression!(obj, opf.prob[c], sum(p_survive * oplim.ramp_mult * 1 * prcc))
    for (cost, gu, gd) in zip(opf.cost_ctrl_gen, pgu, pgd)
        add_to_expression!(obj, opf.prob[c] * p_survive * cost.ramp, gu)
        add_to_expression!(obj, opf.prob[c] * p_survive * cost.ramp, gd)
        # add_to_expression!.(obj, opf.prob[c], p_survive * oplim.ramp_mult * (cost[1] * gu^2 + cost[2] * gu))
        # add_to_expression!.(obj, opf.prob[c], p_survive * oplim.ramp_mult * (cost[1] * gd^2 + cost[2] * gd))
    end

    # Add new constraints that limit the corrective variables within operating limits
    balance_pcc = @expression(m, sum(pgu))
    add_to_expression!.(balance_pcc, lscc)
    add_to_expression!.(balance_pcc, -1, pgd)
    add_to_expression!.(balance_pcc, -1, prcc)
    @constraint(m, balance_pcc == 0.0)
    if isempty(islands)
        @constraints(m, begin
            m[:pg0] .+ pgu .- pgd .>= oplim.pg_lim_min
            m[:pg0] .+ pgu .- pgd .<= oplim.pg_lim_max
        end)
        if typeof(oplim.max_shed) <: Real
            @constraint(m, sum(lscc) <= oplim.max_shed)
        end
    else
        itr = length(islands[island]) < 2 ? Int64[] : islands[island]
        for n in itr
            for g in opf.mgx[:,n].nzind
                expr = @expression(m, m[:pg0][g] + pgu[g] - pgd[g])
                @constraint(m, expr >= oplim.pg_lim_min[g])
                @constraint(m, expr <= oplim.pg_lim_max[g])
            end
        end
        if typeof(oplim.max_shed) <: Real
            expr = AffExpr()
            for n in itr, d in opf.mdx[:,n].nzind
                add_to_expression!(expr, lscc[d])
            end
            @constraint(m, expr <= oplim.max_shed)
        end
        itr = isempty(itr) ? (1:length(opf.nodes)) : islands[1:end.!=island]
        for in_vec in itr
            for n in in_vec
                for g in opf.mgx[:,n].nzind
                    set_upper_bound(pgu[g], oplim.pg_lim_max[g] - oplim.pg_lim_min[g])
                    @constraint(m, m[:pg0][g] + pgu[g] - pgd[g] == 0.0)
                end
                fix!(pfdccc, opf.mdcx[:,n].nzind)
                fix!(prcc, oplim.pr_lim, opf.mrx[:,n].nzind)
                fix!(lscc, oplim.pd_lim, opf.mdx[:,n].nzind)
            end
        end
    end
    return pgu, pgd, pfdccc, prcc, lscc
end

"""
    Initialize variables for the long-term corrective state with contingencies.

    Parameters:
    - `Pcc`: Dictionary of long-term corrective actions.
    - `opf`: The power system data.
    - `oplim`: Operational limits for the case.
    - `m`: Optimization model.
    - `obj`: Objective function.
    - `islands`: List of islands in the power system.
    - `island`: Index of the island to add.
    - `c`: Contingency index.
"""
function init_P!(Pccx::Dict{<:Integer,ExprCCX}, opf::OPFsystem, oplim::Oplimits, m::Model, obj::AbstractJuMPScalar, 
    islands::Vector, island::Integer, c::Integer
)
    # Add corrective variables
    pgd = @variable(m, [g in 1:length(opf.ctrl_generation)], base_name = @sprintf("pgdx%s", c),
        lower_bound = 0.0)#, upper_bound = oplim.rampdown[g] * oplim.ramp_minutes)
    # and ramp down
    prcc = @variable(m, [r in 1:length(opf.renewables)], base_name = @sprintf("prccx%s", c),
        lower_bound = 0.0, upper_bound = oplim.pr_lim[d] * (typeof(oplim.max_curtail) <: Real ? 1.0 : oplim.max_curtail[d]))
    lscc = @variable(m, [d in 1:length(opf.demands)], base_name = @sprintf("lsccx%s", c),
        lower_bound = 0.0, upper_bound = oplim.pd_lim[d] * (typeof(oplim.max_shed) <: Real ? 1.0 : oplim.max_shed[d]))
    # load curtailment variables in in contingencies
    Pccx[c] = ExprCCX(VariableRef[], pgd, prcc, lscc)

    # Extend the objective with the corrective variables
    add_to_expression!(obj, opf.prob[c] * oplim.p_failure, sum(opf.voll' * lscc))
    add_to_expression!(obj, opf.prob[c], sum(oplim.p_failure .* 1 .* prcc))
    add_to_expression!(obj, opf.prob[c] * oplim.p_failure, sum(60 * pgd))

    # Add new constraints that limit the corrective variables within operating limits
    balance_pccx = @expression(m, sum(lscc))
    add_to_expression!.(balance_pccx, -pgd)
    @constraint(m, balance_pccx == 0.0)
    if isempty(islands)
        @constraint(m, m[:pg0] .- pgd .>= 0.0)
    else
        itr = length(islands[island]) < 2 ? Int64[] : islands[island]
        for n in itr
            for g = opf.mgx[:,n].nzind
                @constraint(m, m[:pg0][g] - pgd[g] >= 0.0)
            end
        end
        if typeof(oplim.max_shed) <: Real
            expr = AffExpr()
            for n in itr, d in opf.mdx[:,n].nzind
                add_to_expression!(expr, lscc[d])
            end
            @constraint(m, expr <= oplim.max_shed)
        end
        itr = isempty(itr) ? (1:length(opf.nodes)) : islands[1:end.!=island]
        for in_vec in itr
            for n in in_vec
                for g = opf.mgx[:,n].nzind
                    @constraint(m, m[:pg0][g] - pgd[g] == 0.0)
                end
                fix!(prcc, oplim.pr_lim, opf.mrx[:,n].nzind)
                fix!(lscc, oplim.pd_lim, opf.mdx[:,n].nzind)
            end
        end
    end
    return pgd, prcc, lscc
end

""" 
    Force variables to equal a value 

    Parameters:
    - `var`: Variables to be fixed.
    - `val`: Values to which the variables should be fixed.
    - `vec`: Indices indicating which variables to fix in the vectors.
"""
fix!(var::AbstractVector{VariableRef}, val::AbstractVector{<:Real}, vec::AbstractVector) =
    JuMP.fix.(var[vec], val[vec]; force=true)

""" 
    Force variables to equal zero 

    Parameters:
    - `var`: Variables to be fixed.
    - `vec`: Indices indicating which variables to fix in the vectors.
"""
fix!(var::AbstractVector{VariableRef}, vec::AbstractVector) =
    JuMP.fix.(var[vec], 0.0; force=true)

""" Force variables to equal its current value in the model """
fix_values!(m::Model, symb::Symbol) = JuMP.fix.(m[symb], JuMP.value.(m[symb]), force=true)
fix_values!(m::Model, var::AbstractVector{VariableRef}) = JuMP.fix.(var, JuMP.value.(m[var]), force=true)

""" Fix all base case varibles to its current values in the model """
function fix_base_case!(m::Model)
    fix_values!(m, :pg0)
    fix_values!(m, :pfdc0)
    fix_values!(m, :ls0)
    fix_values!(m, :pr0)
end

""" Fix all contingency variables to its current values in the model """
function fix_contingencies!(m::Model, P::Dict{<:Integer,<:ContExpr})
    for (_, c) in P
        for symb in propertynames(c)
            fix_values!(m, getproperty(c, symb))
        end
    end
end

""" Calculate the cost of VOLL for the provided variables """
calc_cens(m::Model, opf::OPFsystem, var::AbstractVector{VariableRef}) = sum(opf.voll .* get_value(m, var))

""" Calculate the cost of generation for the provided variables """
calc_ctrl_cost(m::Model, opf::OPFsystem, var::AbstractVector{VariableRef}, symb=:var) =
    sum(getproperty(c, symb) * g for (c, g) in zip(opf.cost_ctrl_gen, get_value(m, var)))
# sum(c[1] * g^2 + c[2] * g for (c, g) in zip(opf.cost_ctrl_gen, get_value(m, var)))

""" Calculate the cost of the modelled slack for all contingencies. """
function calc_slack_cost(case::Case)
    m = case.model
    opf = case.opf
    pf = case.pf
    costs = zeros(length(opf.contingencies))
    pg0 = get_value(m, :pg0)
    if isempty(case.pf.ϕ.slack_array)
        slack_cost = sum(last.(get_generator_cost.(opf.ctrl_generation[opf.mgx[:,pf.slack].nzind])))
    else
        slack_cost = opf.mgx' * (getproperty.(opf.cost_ctrl_gen, :ramp) .* case.pf.ϕ.slack_array .* pg0)
    end
    pg0 = opf.mgx' * pg0
    pr = opf.mrx' * (get_value(m, :pr) - get_value(m, :pr0))
    # pr_cost = pr
    pd = opf.mdx' * (get_value(m, :pd) - get_value(m, :ls0))
    ls_cost = opf.mdx' * (opf.voll' * (get_value(m, :pd) - get_value(m, :ls0)))
    for (i, c_obj) in enumerate(opf.contingencies)
        cont = typesort_component(c_obj, opf)
        if is_islanded(pf, cont[2], cont[1])
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont[2], cont[1], pf.slack)
            nodes = reduce(vcat, islands[1:end.!=island], init=[])
            ΔP = sum(pg0[nodes]) + sum(pr[nodes]) - sum(pd[nodes])
            costs[i] = sum(slack_cost * abs(ΔP)) + sum(ls_cost[nodes]) #+ sum(pr_cost[nodes])
        elseif typeof(c_obj) <: StaticInjection
            costs[i] = sum(slack_cost * value(m[:pg0][cont[1]]))
        end
    end
    return costs
end

""" Calculate the base case cost """
calc_objective(m::Model, opf::OPFsystem) =
    calc_ctrl_cost(m, opf, m[:pg0]) + calc_cens(m, opf, m[:ls0])

""" Calculate the short term state cost """
calc_objective(m::Model, opf::OPFsystem, expr::ExprC) =
    calc_ctrl_cost(m, opf, expr.pgu, :ramp) + calc_ctrl_cost(m, opf, expr.pgd, :ramp) + calc_cens(m, opf, expr.lsc)

""" Calculate the long term state cost """
calc_objective(m::Model, opf::OPFsystem, expr::ExprCC) =
    calc_ctrl_cost(m, opf, expr.pgu, :ramp) + calc_ctrl_cost(m, opf, expr.pgd, :ramp) + calc_cens(m, opf, expr.lscc)

""" Calculate the long term with corrective failure state cost """
calc_objective(m::Model, opf::OPFsystem, expr::ExprCCX) =
    60 * calc_ctrl_cost(m, opf, expr.pgdx, :ramp) + calc_cens(m, opf, expr.lsccx)

""" Calculate a state cost """
calc_objective(m::Model, opf::OPFsystem, P::Dict) =
    sum((opf.prob[i] * calc_objective(m, opf, c) for (i, c) in P), init=0.0)

""" Calculate the total cost """
calc_objective(m::Model, opf::OPFsystem, Pc::Dict{<:Integer,ExprC},
    Pcc::Dict{<:Integer,ExprCC}, Pccx::Dict{<:Integer,ExprCCX}
) = calc_objective(m, opf) + calc_objective(m, opf, Pc) +
    calc_objective(m, opf, Pcc) + calc_objective(m, opf, Pccx)

"""
    Print the costs of the base case and slack variable.

    Parameters:
    - `case`: A `Case` object containing the power system data and model.

    Returns:
    - `base_cost`: The cost of the base case.
    - `slack_cost`: The cost of the slack variable.
"""
function print_costs(case::Case)
    base_cost = calc_objective(case.model, case.opf)
    slack_cost = sum(case.opf.prob .* calc_slack_cost(case))
    @printf "Base cost %.2f; " base_cost
    @printf "Slack cost %.2f; " slack_cost
    @printf "Total cost %.2f \n" base_cost + slack_cost
    return base_cost, slack_cost
end

"""
    Get the costs for the base and all contingencies in the case.

    Parameters:
    - `case`: A `Case` object containing the power system data and model.

    Returns:
    - `costs`: A sparse matrix containing the costs for each contingency.
    - `labels`: A vector of labels for each column in the costs matrix.
"""
function get_costs(case::Case)
    costs = SparseArrays.spzeros(length(case.opf.contingencies), 5)
    costs[:,1] .= calc_objective(case.model, case.opf)
    for c in keys(case.Pc)
        costs[c,2] = calc_objective(case.model, case.opf, case.Pc[c])
    end
    for c in keys(case.Pcc)
        costs[c,3] = calc_objective(case.model, case.opf, case.Pcc[c])
    end
    for c in keys(case.Pccx)
        costs[c,4] = calc_objective(case.model, case.opf, case.Pccx[c])
    end
    costs[:,5] = calc_slack_cost(case)
    return costs, [:Base, :Pc, :Pcc, :Pccx, :distslack]
end