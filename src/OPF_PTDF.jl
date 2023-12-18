""" Operational limits type """
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
    dist_slack::Vector{TR}
    pr_lim::Vector{TR}
    max_curtail::Union{TR,Vector{TR}}
    dc_lim_min::Vector{TR}
    dc_lim_max::Vector{TR}
    pd_lim::Vector{TR}
    max_shed::Union{TR,Vector{TR}}
end

""" Constructor for Oplimits """
function oplimits(
    opf::OPFsystem,
    max_shed::Union{TR,Vector{TR}},
    max_curtail::Union{TR,Vector{TR}},
    dist_slack::Vector{TR},
    ramp_mult::TR,
    ramp_minutes::TR,
    p_failure::TR,
    short_term_multi::Union{TR,Vector{TR}},
    long_term_multi::Union{TR,Vector{TR}}
) where {TR<:Real}
    branch_rating::Vector{TR} = get_rate.(opf.branches)
    (pg_lim_min::Vector{TR}, pg_lim_max::Vector{TR}) = split_pair(get_active_power_limits.(opf.ctrl_generation))
    (rampup::Vector{TR}, rampdown::Vector{TR}) = split_pair(get_ramp_limits.(opf.ctrl_generation))
    pr_lim::Vector{TR} = get_active_power.(opf.renewables)
    (dc_lim_min::Vector{TR}, dc_lim_max::Vector{TR}) = split_pair(get_active_power_limits_from.(opf.dc_branches))
    pd_lim::Vector{TR} = get_active_power.(opf.demands)
    if !isempty(dist_slack)
        dist_slack = dist_slack / sum(dist_slack)
    end
    return Oplimits{TR}(ramp_mult, ramp_minutes, p_failure, branch_rating, short_term_multi, long_term_multi, 
        pg_lim_min, pg_lim_max, rampup, rampdown, dist_slack, pr_lim, max_curtail, dc_lim_min, dc_lim_max, pd_lim, max_shed)
end

mutable struct Case{TR<:Real,TI<:Integer}
    model::Model
    opf::OPFsystem{TR}
    pf::DCPowerFlow{TR,TI}
    oplim::Oplimits{TR}
    Pc::Dict{Int,ExprC}
    Pcc::Dict{Int,ExprCC}
    Pccx::Dict{Int,ExprCCX}
end

""" Run an OPF of a power system """
function opf_base(type::OPF, system::System, optimizer;
    voll=Float64[],
    contingencies=Component[],
    prob=Float64[],
    dist_slack=Float64[],
    time_limit_sec::Int64=600,
    unit_commit::Bool=false,
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
    mod = create_model(optimizer, time_limit_sec=time_limit_sec, silent=silent, debug=debug)
    opf = isempty(voll) ? opfsystem(system) : opfsystem(system, voll, contingencies, prob, ramp_mult)
    pf = DCPowerFlow(opf.nodes, opf.branches, opf.idx)
    opf.dc_branches = DCBranch[]
    
    Pc = Dict{Int,ExprC}()
    Pcc = Dict{Int,ExprCC}()
    Pccx = Dict{Int,ExprCCX}()

    oplim = oplimits(opf, max_shed, max_curtail, dist_slack, ramp_mult, ramp_minutes, p_failure, short_term_multi, long_term_multi)

    @variables(mod, begin
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

    @objective(mod, Min, sum(c.var * g for (c, g) in zip(opf.cost_ctrl_gen, pg0)) + sum(opf.voll' * ls0) + sum(pr0 * 1))
    # @objective(mod, Min, sum(c[1] * g^2 + c[2] * g for (c, g) in zip(opf.cost_ctrl_gen, pg0)) + sum(opf.voll' * ls0) + sum(pr0 * 1))

    JuMP.fix.(pd, oplim.pd_lim)
    JuMP.fix.(pr, oplim.pr_lim)
    if typeof(oplim.max_shed) <: Real
        @constraint(mod, sum(ls0) <= oplim.max_shed)
    end

    @expression(mod, inj_p0, opf.mgx' * pg0)
    add_to_expression!.(inj_p0, opf.mrx' * (pr - pr0))
    add_to_expression!.(inj_p0, opf.mdcx' * pfdc0)
    add_to_expression!.(inj_p0, opf.mdx' * (ls0 - pd))
    @constraint(mod, inj_p0 .== p0)

    @expression(mod, ptdf0, pf.ϕ * p0)
    @constraint(mod, branch_lim_n, -oplim.branch_rating .<= ptdf0)
    @constraint(mod, branch_lim_p, ptdf0 .<= oplim.branch_rating)
    @expression(mod, balance, sum(pg0, init=0.0))
    add_to_expression!.(balance, pr)
    add_to_expression!.(balance, -pr0)
    add_to_expression!.(balance, -pd)
    add_to_expression!.(balance, ls0)
    @constraint(mod, power_balance, balance == 0.0)

    if unit_commit
        add_unit_commit!(opf)
    end

    case = Case(mod, opf, pf, oplim, Pc, Pcc, Pccx)
    if any([type.P, type.C1, type.C2, type.C2F])
        return add_all_contingencies!(type, opf, oplim, mod, pf, Pc, Pcc, Pccx)
    end
    return mod, opf, pf, oplim, Pc, Pcc, Pccx
end

function add_all_contingencies!(type::OPF, opf::OPFsystem, oplim::Oplimits, mod::Model,
    pf::DCPowerFlow, Pc::Dict{<:Integer,ExprC}, Pcc::Dict{<:Integer,ExprCC}, Pccx::Dict{<:Integer,ExprCCX}
)
    obj = objective_function(mod)
    set_dist_slack!(pf, opf, oplim.dist_slack)
    for (i, c_obj) in enumerate(opf.contingencies)
        cont = typesort_component(c_obj, opf)
        if is_islanded(pf, cont[2], cont[1])
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont[2], cont[1], pf.slack)
            ptdf = get_isf(pf, cont[2], cont[1], islands, island, island_b)
            set_tol_zero!(ptdf)
        else
            ptdf = get_isf(pf, cont[2], cont[1])
            islands = Vector{Vector{Int}}[]
            island = 0
        end
        type.P && add_contingencies!(opf, oplim, mod, ptdf, i)
        type.C1 && add_contingencies!(Pc, opf, oplim, mod, obj, islands, island, ptdf, i)
        type.C2 && add_contingencies!(Pcc, opf, oplim, mod, obj, islands, island, ptdf, i)
        type.C2F && add_contingencies!(Pccx, opf, oplim, mod, obj, islands, island, ptdf, i)

        @debug "Contingency $(get_name(c_obj)) is added"
    end
    set_objective_function(mod, obj)
    return mod, opf, pf, oplim, Pc, Pcc, Pccx
end

function add_contingencies!(opf::OPFsystem, oplim::Oplimits, mod::Model, ptdf::AbstractMatrix{<:Real}, c::Integer)
    p = JuMP.@variable(mod, [n in 1:length(opf.nodes)], base_name = @sprintf("p%s", c))
    @constraint(mod, mod[:inj_p0] .== p)
    ptdf_p = @expression(mod, ptdf * p)
    @constraint(mod, ptdf_p .>= -oplim.branch_rating * oplim.short_term_multi)
    @constraint(mod, ptdf_p .<= oplim.branch_rating * oplim.short_term_multi)
end

function add_contingencies!(Pc::Dict{<:Integer,ExprC}, opf::OPFsystem, oplim::Oplimits, mod::Model, obj::AbstractJuMPScalar, islands::Vector,
    island::Integer, ptdf::AbstractMatrix{<:Real}, c::Integer
)
    pc = JuMP.@variable(mod, [n in 1:length(opf.nodes)], base_name = @sprintf("pc%s", c))
    pgu, pgd, prc, lsc = init_P!(Pc, opf, oplim, mod, obj, islands, island, c)
    Pc[c].pc = pc

    inj_pc = @expression(mod, opf.mgx' * (mod[:pg0] + pgu - pgd))
    add_to_expression!.(inj_pc, opf.mrx' * (mod[:pr] - prc))
    add_to_expression!.(inj_pc, opf.mdcx' * mod[:pfdc0])
    add_to_expression!.(inj_pc, opf.mdx' * (lsc - mod[:pd]))
    @constraint(mod, inj_pc .== pc)

    ptdf_pc = @expression(mod, ptdf * pc)
    @constraint(mod, ptdf_pc .>= -oplim.branch_rating * oplim.short_term_multi)
    @constraint(mod, ptdf_pc .<= oplim.branch_rating * oplim.short_term_multi)
end

function add_contingencies!(Pcc::Dict{<:Integer,ExprCC}, opf::OPFsystem, oplim::Oplimits, mod::Model, obj::AbstractJuMPScalar, islands::Vector,
    island::Integer, ptdf::AbstractMatrix{<:Real}, c::Integer
)
    pcc = JuMP.@variable(mod, [n in 1:length(opf.nodes)], base_name = @sprintf("pcc%s", c))
    pgu, pgd, pfdccc, prcc, lscc = init_P!(Pcc, opf, oplim, mod, obj, islands, island, c)
    Pcc[c].pcc = pcc

    inj_pcc = @expression(mod, opf.mgx' * (mod[:pg0] + pgu - pgd))
    add_to_expression!.(inj_pcc, opf.mrx' * (mod[:pr] - prcc))
    add_to_expression!.(inj_pcc, opf.mdcx' * pfdccc)
    add_to_expression!.(inj_pcc, opf.mdx' * (lscc - mod[:pd]))
    @constraint(mod, inj_pcc .== pcc)

    ptdf_pcc = @expression(mod, ptdf * pcc)
    @constraint(mod, ptdf_pcc .>= -oplim.branch_rating * oplim.long_term_multi)
    @constraint(mod, ptdf_pcc .<= oplim.branch_rating * oplim.long_term_multi)
end

function add_contingencies!(Pccx::Dict{<:Integer,ExprCCX}, opf::OPFsystem, oplim::Oplimits, mod::Model, obj::AbstractJuMPScalar, islands::Vector,
    island::Integer, ptdf::AbstractMatrix{<:Real}, c::Integer
)
    pccx = JuMP.@variable(mod, [n in 1:length(opf.nodes)], base_name = @sprintf("pccx%s", c))
    pgd, prcc, lscc = init_P!(Pccx, opf, oplim, mod, obj, islands, island, c)
    Pccx[c].pccx = pccx

    inj_pccx = @expression(mod, opf.mgx' * (mod[:pg0] - pgd))
    add_to_expression!.(inj_pccx, opf.mrx' * (mod[:pr] - prcc))
    add_to_expression!.(inj_pccx, opf.mdcx' * mod[:pfdc0])
    add_to_expression!.(inj_pccx, opf.mdx' * (lscc - mod[:pd]))
    @constraint(mod, inj_pccx .== pccx)

    @constraint(mod, ptdf * pccx .>= -oplim.branch_rating)
    @constraint(mod, ptdf * pccx .<= oplim.branch_rating)
end


function init_P!(Pc::Dict{<:Integer,ExprC}, opf::OPFsystem, oplim::Oplimits, mod::Model, obj::JuMP.AbstractJuMPScalar, 
    islands::Vector, island::Integer, c::Integer
)
    pgu = JuMP.@variable(mod, [g in 1:length(opf.ctrl_generation)], base_name = @sprintf("pgu%s", c),
        lower_bound = 0.0, upper_bound = oplim.rampup[g] * 0.0)
    # active power variables for the generators in contingencies ramp up 
    pgd = JuMP.@variable(mod, [g in 1:length(opf.ctrl_generation)], base_name = @sprintf("pgd%s", c),
        lower_bound = 0.0) #, upper_bound = rampdown[g] * 0.0)
    # and ramp down
    prc = JuMP.@variable(mod, [d in 1:length(opf.renewables)], base_name = @sprintf("prc%s", c),
        lower_bound = 0.0, upper_bound = oplim.pr_lim[d] * (typeof(oplim.max_curtail) <: Real ? 1.0 : oplim.max_curtail[d]))
    lsc = JuMP.@variable(mod, [d in 1:length(opf.demands)], base_name = @sprintf("lsc%s", c),
        lower_bound = 0.0, upper_bound = oplim.pd_lim[d] * (typeof(oplim.max_shed) <: Real ? 1.0 : oplim.max_shed[d]))
    Pc[c] = ExprC(JuMP.VariableRef[], pgu, pgd, prc, lsc)

    p_survive = 1.0 - oplim.p_failure
    add_to_expression!(obj, opf.prob[c], sum(opf.voll' * lsc))
    add_to_expression!(obj, opf.prob[c], sum(oplim.ramp_mult * 1 * prc))
    for (cost, gu, gd) in zip(opf.cost_ctrl_gen, pgu, pgd)
        add_to_expression!.(obj, opf.prob[c], p_survive * (cost.ramp * gu))
        add_to_expression!.(obj, opf.prob[c], p_survive * (cost.ramp * gd))
        # add_to_expression!.(obj, opf.prob[c], p_survive * oplim.ramp_mult * (cost[1] * g^2 + cost[2] * g))
    end

    # Add new constraints that limit the corrective variables within operating limits
    balance_pc = JuMP.@expression(mod, sum(lsc))
    add_to_expression!.(balance_pc, -prc)
    add_to_expression!.(balance_pc, pgu)
    add_to_expression!.(balance_pc, -pgd)
    JuMP.@constraint(mod, balance_pc == 0.0)
    if isempty(islands)
        JuMP.@constraint(mod, mod[:pg0] .+ pgu .- pgd .>= oplim.pg_lim_min)
        if typeof(oplim.max_shed) <: Real
            @constraint(mod, sum(lsc) <= oplim.max_shed)
        end
    else
        itr = length(islands[island]) < 2 ? Int[] : islands[island]
        pgc = @expression(mod, opf.mgx[:, itr]' * (mod[:pg0] + pgu - pgd))
        for n in itr
            JuMP.@constraints(mod, begin
                [g = opf.mgx[:,n].nzind], mod[:pg0][g] + pgu[g] - pgd[g] >= oplim.pg_lim_min[g]
                [g = opf.mgx[:,n].nzind], mod[:pg0][g] + pgu[g] - pgd[g] <= oplim.pg_lim_max[g]
            end)
        end
        itr = isempty(itr) ? (1:length(opf.nodes)) : islands[1:end.!=island]
        for in_vec in itr
            for n in in_vec
                for g in opf.mgx[:,n].nzind
                    JuMP.set_upper_bound(pgu[g], oplim.pg_lim_max[g] - oplim.pg_lim_min[g])
                    JuMP.@constraint(mod, mod[:pg0][g] + pgu[g] - pgd[g] == 0.0)
                end
                fix!(prc, oplim.pr_lim, opf.mrx[:,n].nzind)
                fix!(lsc, oplim.pd_lim, opf.mdx[:,n].nzind)
            end
        end
    end
    return pgu, pgd, prc, lsc
end

function init_P!(Pcc::Dict{<:Integer,ExprCC}, opf::OPFsystem, oplim::Oplimits, mod::Model, obj::JuMP.AbstractJuMPScalar, 
    islands::Vector, island::Integer, c::Integer
)
    # Add corrective variables
    pgu = JuMP.@variable(mod, [g in 1:length(opf.ctrl_generation)], base_name = @sprintf("pgu%s", c),
        lower_bound = 0.0, upper_bound = oplim.rampup[g] * oplim.ramp_minutes)
    # active power variables for the generators in contingencies ramp up 
    pgd = JuMP.@variable(mod, [g in 1:length(opf.ctrl_generation)], base_name = @sprintf("pgd%s", c),
        lower_bound = 0.0) #, upper_bound = rampdown[g] * ramp_minutes)
    # and ramp down
    pfdccc = JuMP.@variable(mod, [d in 1:length(opf.dc_branches)], base_name = @sprintf("pfdccc%s", c),
        lower_bound = oplim.dc_lim_min[d], upper_bound = oplim.dc_lim_max[d])
    prcc = JuMP.@variable(mod, [r in 1:length(opf.renewables)], base_name = @sprintf("prcc%s", c),
        lower_bound = 0.0, upper_bound = oplim.pr_lim[r] * (typeof(oplim.max_curtail) <: Real ? 1.0 : oplim.max_curtail[d]))
    lscc = JuMP.@variable(mod, [d in 1:length(opf.demands)], base_name = @sprintf("lscc%s", c),
        lower_bound = 0.0, upper_bound = oplim.pd_lim[d] * (typeof(oplim.max_shed) <: Real ? 1.0 : oplim.max_shed[d]))
    # load curtailment variables in in contingencies
    Pcc[c] = ExprCC(JuMP.VariableRef[], pgu, pgd, pfdccc, prcc, lscc)

    p_survive = 1.0 - oplim.p_failure
    # Extend the objective with the corrective variables
    # obj = objective_function(mod)
    add_to_expression!.(obj, opf.prob[c],
        # (1.0 - p_failure) * (sum(opf.voll' * lscc) + sum(60 * (pgu + pgd)) # + # TODO: remove 60 and uncomment next lines for non-4-area analysis!!!!
        p_survive * (sum(opf.voll' * lscc) # + sum(opf.cost_ctrl_gen' * ramp_mult * (pgu + pgd))
        # (sum(opf.voll' * lscc) +
        # sum(opf.cost_ctrl_gen' * ramp_mult * pgu) # +
        # sum(opf.cost_ctrl_gen' * pgd)
        ))
    add_to_expression!(obj, opf.prob[c], sum(p_survive * oplim.ramp_mult * 1 * prcc))
    for (cost, gu, gd) in zip(opf.cost_ctrl_gen, pgu, pgd)
        add_to_expression!.(obj, opf.prob[c], p_survive * (cost.ramp * gu))
        add_to_expression!.(obj, opf.prob[c], p_survive * (cost.ramp * gd))
        # add_to_expression!.(obj, opf.prob[c], p_survive * oplim.ramp_mult * (cost[1] * gu^2 + cost[2] * gu))
        # add_to_expression!.(obj, opf.prob[c], p_survive * oplim.ramp_mult * (cost[1] * gd^2 + cost[2] * gd))
    end

    # Add new constraints that limit the corrective variables within operating limits
    balance_pcc = JuMP.@expression(mod, sum(pgu))
    add_to_expression!.(balance_pcc, lscc)
    add_to_expression!.(balance_pcc, -pgd)
    add_to_expression!.(balance_pcc, -prcc)
    JuMP.@constraint(mod, balance_pcc == 0.0)
    if isempty(islands)
        JuMP.@constraints(mod, begin
            mod[:pg0] .+ pgu .- pgd .>= oplim.pg_lim_min
            mod[:pg0] .+ pgu .- pgd .<= oplim.pg_lim_max
        end)
        if typeof(oplim.max_shed) <: Real
            @constraint(mod, sum(lscc) <= oplim.max_shed)
        end
    else
        itr = length(islands[island]) < 2 ? Int[] : islands[island]
        for n in itr
            JuMP.@constraints(mod, begin
                [g = opf.mgx[:,n].nzind], mod[:pg0][g] + pgu[g] - pgd[g] >= oplim.pg_lim_min[g]
                [g = opf.mgx[:,n].nzind], mod[:pg0][g] + pgu[g] - pgd[g] <= oplim.pg_lim_max[g]
            end)
        end
        itr = isempty(itr) ? (1:length(opf.nodes)) : islands[1:end.!=island]
        for in_vec in itr
            for n in in_vec
                for g in opf.mgx[:,n].nzind
                    JuMP.set_upper_bound(pgu[g], oplim.pg_lim_max[g] - oplim.pg_lim_min[g])
                    JuMP.@constraint(mod, mod[:pg0][g] + pgu[g] - pgd[g] == 0.0)
                end
                fix!(pfdccc, opf.mdcx[:,n].nzind)
                fix!(prcc, oplim.pr_lim, opf.mrx[:,n].nzind)
                fix!(lscc, oplim.pd_lim, opf.mdx[:,n].nzind)
            end
        end
    end
    return pgu, pgd, pfdccc, prcc, lscc
end

function init_P!(Pccx::Dict{<:Integer,ExprCCX}, opf::OPFsystem, oplim::Oplimits, mod::Model, obj::JuMP.AbstractJuMPScalar, 
    islands::Vector, island::Integer, c::Integer
)
    # Add corrective variables
    pgd = JuMP.@variable(mod, [g in 1:length(opf.ctrl_generation)], base_name = @sprintf("pgdx%s", c),
        lower_bound = 0.0)#, upper_bound = oplim.rampdown[g] * oplim.ramp_minutes)
    # and ramp down
    prcc = JuMP.@variable(mod, [r in 1:length(opf.renewables)], base_name = @sprintf("prccx%s", c),
        lower_bound = 0.0, upper_bound = oplim.pr_lim[d] * (typeof(oplim.max_curtail) <: Real ? 1.0 : oplim.max_curtail[d]))
    lscc = JuMP.@variable(mod, [d in 1:length(opf.demands)], base_name = @sprintf("lsccx%s", c),
        lower_bound = 0.0, upper_bound = oplim.pd_lim[d] * (typeof(oplim.max_shed) <: Real ? 1.0 : oplim.max_shed[d]))
    # load curtailment variables in in contingencies
    Pccx[c] = ExprCCX(JuMP.VariableRef[], pgd, prcc, lscc)

    # Extend the objective with the corrective variables
    add_to_expression!(obj, opf.prob[c] * oplim.p_failure, sum(opf.voll' * lscc))
    add_to_expression!(obj, opf.prob[c], sum(oplim.p_failure .* 1 .* prcc))
    add_to_expression!(obj, opf.prob[c] * oplim.p_failure, sum(60 * pgd))

    # Add new constraints that limit the corrective variables within operating limits
    balance_pccx = JuMP.@expression(mod, sum(lscc))
    add_to_expression!.(balance_pccx, -pgd)
    JuMP.@constraint(mod, balance_pccx == 0.0)
    if isempty(islands)
        JuMP.@constraint(mod, mod[:pg0] .- pgd .>= 0.0)
    else
        itr = length(islands[island]) < 2 ? Int[] : islands[island]
        for n in itr
            JuMP.@constraint(mod, [g = opf.mgx[:,n].nzind], mod[:pg0][g] - pgd[g] >= 0.0)
        end
        itr = isempty(itr) ? (1:length(opf.nodes)) : islands[1:end.!=island]
        for in_vec in itr
            for n in in_vec
                JuMP.@constraint(mod, [g = opf.mgx[:,n].nzind], mod[:pg0][g] - pgd[g] == 0.0)
                fix!(prcc, oplim.pr_lim, opf.mrx[:,n].nzind)
                fix!(lscc, oplim.pd_lim, opf.mdx[:,n].nzind)
            end
        end
    end
    return pgd, prcc, lscc
end

""" Force a variable to equal a value"""
fix!(var::AbstractVector{JuMP.VariableRef}, val::AbstractVector{<:Real}, vec::AbstractVector) =
    JuMP.fix.(var[vec], val[vec]; force=true)

""" Force a variable to equal zero"""
fix!(var::AbstractVector{JuMP.VariableRef}, vec::AbstractVector) =
    JuMP.fix.(var[vec], 0.0; force=true)

""" Force a variable to equal its current value in the model """
fix_values!(mod::Model, symb::Symbol) = JuMP.fix.(mod[symb], get_value(mod, symb), force=true)
fix_values!(mod::Model, var::AbstractVector{JuMP.VariableRef}) = JuMP.fix.(var, get_value(mod, var), force=true)

""" Fix all base case varibles to its current values in the model """
function fix_base_case!(mod::Model)
    fix_values!(mod, :pg0)
    fix_values!(mod, :pfdc0)
    fix_values!(mod, :ls0)
    fix_values!(mod, :pr0)
end

""" Fix all short term varibles to its current values in the model """
function fix_short_term!(mod::Model, Pc::Dict{<:Integer,ExprC})
    for (_, c) in Pc
        fix_values!(mod, c.pgu)
        fix_values!(mod, c.pgd)
        fix_values!(mod, c.prc)
        fix_values!(mod, c.lsc)
    end
end

""" Fix all long term varibles to its current values in the model """
function fix_long_term!(mod::Model, Pcc::Dict{<:Integer,ExprCC})
    for (_, c) in Pcc
        fix_values!(mod, c.pgu)
        fix_values!(mod, c.pgd)
        fix_values!(mod, c.pfdccc)
        fix_values!(mod, c.prcc)
        fix_values!(mod, c.lscc)
    end
end

""" Calculate the cost of VOLL for the provided variables """
calc_cens(mod::Model, opf::OPFsystem, var::AbstractVector{JuMP.VariableRef}) = sum(opf.voll .* get_value(mod, var))

""" Calculate the cost of generation for the provided variables """
calc_ctrl_cost(mod::Model, opf::OPFsystem, var::AbstractVector{JuMP.VariableRef}, symb=:var) =
    sum(getproperty(c, symb) * g for (c, g) in zip(opf.cost_ctrl_gen, get_value(mod, var)))
# sum(c[1] * g^2 + c[2] * g for (c, g) in zip(opf.cost_ctrl_gen, get_value(mod, var)))

""" Calculate the cost of the modelled slack for all contingencies. """
function calc_slack_cost(mod::Model, opf::OPFsystem, pf::DCPowerFlow; dist_slack=Float64[])
    costs = zeros(length(opf.contingencies))
    pg0 = get_value(mod, :pg0)
    if isempty(dist_slack)
        slack_cost = sum(last.(get_generator_cost.(opf.ctrl_generation[opf.mgx[:,pf.slack].nzind])))
    else
        dist_slack = dist_slack / sum(dist_slack)
        slack_cost = opf.mgx' * (getproperty.(opf.cost_ctrl_gen, :ramp) .* dist_slack .* pg0)
    end
    pg0 = opf.mgx' * pg0
    pr = opf.mrx' * (get_value(mod, :pr) - get_value(mod, :pr0))
    pr_cost = pr
    pd = opf.mdx' * (get_value(mod, :pd) - get_value(mod, :ls0))
    ls_cost = opf.mdx' * (opf.voll' * (get_value(mod, :pd) - get_value(mod, :ls0)))
    for (i, c_obj) in enumerate(opf.contingencies)
        cont = typesort_component(c_obj, opf)
        if is_islanded(pf, cont[2], cont[1])
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont[2], cont[1], pf.slack)
            nodes = reduce(vcat, islands[1:end.!=island], init=[])
            ΔP = sum(pg0[nodes]) + sum(pr[nodes]) - sum(pd[nodes])
            costs[i] = sum(slack_cost * ΔP) + sum(pr_cost[nodes]) + sum(ls_cost[nodes])
        elseif typeof(c_obj) <: StaticInjection
            costs[i] = sum(slack_cost * JuMP.value(mod[:pg0][cont[1]]))
        end
    end
    return costs
end

""" Calculate the base case cost """
calc_objective(mod::Model, opf::OPFsystem) =
    calc_ctrl_cost(mod, opf, mod[:pg0]) + calc_cens(mod, opf, mod[:ls0])

""" Calculate the short term state cost """
calc_objective(mod::Model, opf::OPFsystem, expr::ExprC) =
    calc_ctrl_cost(mod, opf, expr.pgu, :ramp) + calc_ctrl_cost(mod, opf, expr.pgd, :ramp) + calc_cens(mod, opf, expr.lsc)

""" Calculate the long term state cost """
calc_objective(mod::Model, opf::OPFsystem, expr::ExprCC) =
    calc_ctrl_cost(mod, opf, expr.pgu, :ramp) + calc_ctrl_cost(mod, opf, expr.pgd, :ramp) + calc_cens(mod, opf, expr.lscc)

""" Calculate the long term with corrective failure state cost """
calc_objective(mod::Model, opf::OPFsystem, expr::ExprCCX) =
    60 * calc_ctrl_cost(mod, opf, expr.pgdx, :ramp) + calc_cens(mod, opf, expr.lsccx)

""" Calculate a state cost """
calc_objective(mod::Model, opf::OPFsystem, P::Dict) =
    sum((opf.prob[i] * calc_objective(mod, opf, c) for (i, c) in P), init=0.0)

""" Calculate the total cost """
calc_objective(mod::Model, opf::OPFsystem, Pc::Dict{<:Integer,ExprC},
    Pcc::Dict{<:Integer,ExprCC}, Pccx::Dict{<:Integer,ExprCCX}
) = calc_objective(mod, opf) + calc_objective(mod, opf, Pc) +
    calc_objective(mod, opf, Pcc) + calc_objective(mod, opf, Pccx)

function print_costs(mod::Model, opf::OPFsystem, pf::DCPowerFlow; dist_slack=Float64[]
)
    base_cost = calc_objective(mod, opf)
    slack_cost = sum(opf.prob .* calc_slack_cost(mod, opf, pf, dist_slack=dist_slack))
    @printf "Base cost %.2f; " base_cost
    @printf "Slack cost %.2f; " slack_cost
    @printf "Total cost %.2f \n" base_cost + slack_cost
    return base_cost, slack_cost
end

function get_costs(mod::Model, opf::OPFsystem, pf::DCPowerFlow, Pc::Dict{<:Integer,ExprC},
    Pcc::Dict{<:Integer,ExprCC}, Pccx::Dict{<:Integer,ExprCCX}, dist_slack=Float64[]
)
    costs = SparseArrays.spzeros(length(opf.contingencies), 5)
    costs[:,1] .= calc_objective(mod, opf)
    for c in keys(Pc)
        costs[c,2] = calc_objective(mod, opf, Pc[c])
    end
    for c in keys(Pcc)
        costs[c,3] = calc_objective(mod, opf, Pcc[c])
    end
    for c in keys(Pccx)
        costs[c,4] = calc_objective(mod, opf, Pccx[c])
    end
    costs[:,5] = calc_slack_cost(mod, opf, pf, dist_slack=dist_slack)
    return costs, [:Base, :Pc, :Pcc, :Pccx, :distslack]
end