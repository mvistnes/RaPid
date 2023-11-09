""" Operational limits type """
struct Oplimits{TR<:Real}
    max_shed::TR
    max_curtail::TR
    ramp_mult::TR
    ramp_minutes::TR
    p_failure::TR
    short_term_multi::TR
    long_term_multi::TR

    branch_rating::Vector{TR}
    pg_lim_min::Vector{TR}
    pg_lim_max::Vector{TR}
    pr_lim::Vector{TR}
    dc_lim_min::Vector{TR}
    dc_lim_max::Vector{TR}
    pd_lim::Vector{TR}
    rampup::Vector{TR}
    rampdown::Vector{TR}
end

""" Constructor for Oplimits """
function oplimits(
    opf::OPFsystem,
    max_shed::TR,
    max_curtail::TR,
    ramp_mult::TR,
    ramp_minutes::TR,
    p_failure::TR,
    short_term_multi::TR,
    long_term_multi::TR
) where {TR<:Real}
    branch_rating = get_rate.(opf.branches)
    (pg_lim_min, pg_lim_max) = split_pair(get_active_power_limits.(opf.ctrl_generation))
    pr_lim::Vector{TR} = get_active_power.(opf.renewables)
    (dc_lim_min::Vector{TR}, dc_lim_max::Vector{TR}) = split_pair(get_active_power_limits_from.(opf.dc_branches))
    pd_lim = get_active_power.(opf.demands)
    (rampup, rampdown) = split_pair(get_ramp_limits.(opf.ctrl_generation))
    return Oplimits{TR}(max_shed, max_curtail, ramp_mult, ramp_minutes, p_failure, short_term_multi, long_term_multi,
        branch_rating, pg_lim_min, pg_lim_max, pr_lim, dc_lim_min, dc_lim_max, pd_lim, rampup, rampdown)
end

""" Run an OPF of a power system """
function opf_base(type::OPF, system::System, optimizer;
    voll=Float64[],
    contingencies=Component[],
    prob=Float64[],
    time_limit_sec::Int64=600,
    unit_commit::Bool=false,
    ramp_minutes::Real=10.0,
    ramp_mult::Real=10.0,
    max_shed::Real=1.0,
    max_curtail::Real=1.0,
    short_term_multi::Real=1.5,
    long_term_multi::Real=1.2,
    p_failure=0.0,
    silent=true,
    debug=false
)
    @assert type.Base
    assert(type)
    @assert short_term_multi >= long_term_multi
    @assert !type.C2F | !iszero(p_failure)
    mod = create_model(optimizer, time_limit_sec=time_limit_sec, silent=silent, debug=debug)
    opf = isempty(voll) ? opfsystem(system) : opfsystem(system, voll, contingencies, prob, ramp_mult)
    idx = get_nodes_idx(opf.nodes)
    pf = DCPowerFlow(opf.nodes, opf.branches, idx)
    opf.dc_branches = DCBranch[]

    list = make_list(opf, idx, opf.nodes)
    Pc = Dict{Int,ExprC}()
    Pcc = Dict{Int,ExprCC}() 
    Pccx = Dict{Int,ExprCCX}() 

    oplim = oplimits(opf, max_shed, max_curtail, ramp_mult, ramp_minutes, p_failure, short_term_multi, long_term_multi)

    @variables(mod, begin
        p0[n in 1:length(opf.nodes)]
        # active power injection on each node
        oplim.pg_lim_min[g] <= pg0[g in 1:length(opf.ctrl_generation)] <= oplim.pg_lim_max[g]
        # active power variables for the generators
        oplim.dc_lim_min[l] <= pfdc0[l in 1:length(opf.dc_branches)] <= oplim.dc_lim_max[l]
        # power flow on DC branches
        pd[d in 1:length(opf.demands)]
        0.0 <= ls0[d in 1:length(opf.demands)] <= oplim.pd_lim[d] * oplim.max_shed
        # demand curtailment variables
        pr[d in 1:length(opf.renewables)]
        0.0 <= pr0[d in 1:length(opf.renewables)] <= oplim.pr_lim[d] * oplim.max_curtail
        # renewable curtailment variables
    end)

    @objective(mod, Min, sum(c.var * g for (c, g) in zip(opf.cost_ctrl_gen, pg0)) + sum(opf.voll' * ls0) + sum(pr0 * 1))
    # @objective(mod, Min, sum(c[1] * g^2 + c[2] * g for (c, g) in zip(opf.cost_ctrl_gen, pg0)) + sum(opf.voll' * ls0) + sum(pr0 * 1))

    JuMP.fix.(pd, get_active_power.(opf.demands))
    JuMP.fix.(pr, get_active_power.(opf.renewables))

    @expression(mod, inj_p0, -p0)
    for n = 1:length(opf.nodes)
        add_to_expression!(inj_p0[n], sum((pg0[g] for g in list[n].ctrl_generation), init=0.0))
        add_to_expression!(inj_p0[n], sum((pr[d] for d in list[n].renewables), init=0.0))
        add_to_expression!(inj_p0[n], -1, sum((pd[d] for d in list[n].demands), init=0.0))
        add_to_expression!(inj_p0[n], sum((beta(opf.nodes[n], opf.dc_branches[l]) * pfdc0[l] for l in list[n].dc_branches), init=0.0))
        add_to_expression!(inj_p0[n], -1, sum((pr0[d] for d in list[n].renewables), init=0.0))
        add_to_expression!(inj_p0[n], sum((ls0[d] for d in list[n].demands), init=0.0))
    end
    @constraint(mod, inj_p0 .== 0.0)

    @expression(mod, ptdf0, pf.ϕ * p0)
    @constraint(mod, branch_lim_n, -oplim.branch_rating .<= ptdf0)
    @constraint(mod, branch_lim_p, ptdf0 .<= oplim.branch_rating)
    @expression(mod, balance, sum(pg0, init=0.0))
    add_to_expression!.(balance, ls0)
    add_to_expression!.(balance, pr)
    add_to_expression!.(balance, -pd)
    add_to_expression!.(balance, -pr0)
    @constraint(mod, power_balance, balance == 0.0)

    if unit_commit
        add_unit_commit!(opf)
    end

    if any([type.P, type.C1, type.C2, type.C2F])
        return add_all_contingencies!(type, opf, oplim, mod, list, pf, idx, Pc, Pcc, Pccx)
    end
    return mod, opf, pf, oplim, Pc, Pcc, Pccx
end

function add_all_contingencies!(type::OPF, opf::OPFsystem, oplim::Oplimits, mod::Model, list::Vector{<:CTypes{Int}},
    pf::DCPowerFlow, idx::Dict{<:Int,<:Int},
    Pc::Dict{<:Integer,ExprC}, Pcc::Dict{<:Integer,ExprCC}, Pccx::Dict{<:Integer,ExprCCX}
)
    obj = objective_function(mod)
    for (i, c_obj) in enumerate(opf.contingencies)
        (typelist, c, cont) = typesort_component(c_obj, opf, idx)
        if is_islanded(pf, cont, c)
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont, c, pf.slack)
            ptdf = get_isf(pf, cont, c, islands, island, island_b)
        else
            ptdf = get_isf(pf, cont, c)
            islands = Vector{Vector{Int}}[]
            island = 0
        end
        type.P && add_contingencies!(opf, oplim, mod, ptdf, list, i)
        type.C1 && add_contingencies!(Pc, opf, oplim, mod, obj, islands, island, ptdf, list, i)
        type.C2 && add_contingencies!(Pcc, opf, oplim, mod, obj, islands, island, ptdf, list, i)
        type.C2F && add_contingencies!(Pccx, opf, oplim, mod, obj, islands, island, ptdf, list, i)

        @debug "Contingency $(get_name(c_obj)) is added"
    end
    set_objective_function(mod, obj)
    return mod, opf, pf, oplim, Pc, Pcc, Pccx
end

function add_contingencies!(opf::OPFsystem, oplim::Oplimits, mod::Model, ptdf::AbstractMatrix{<:Real}, list::Vector{<:CTypes{Int}}, c::Integer)
    p = JuMP.@variable(mod, [n in 1:length(opf.nodes)], base_name = @sprintf("p%s", c))

    inj_p = @expression(mod, -p)
    for n = 1:length(opf.nodes)
        add_to_expression!(inj_p[n], sum((mod[:pg0][g] for g in list[n].ctrl_generation), init=0.0))
        add_to_expression!(inj_p[n], sum((beta(opf.nodes[n], opf.dc_branches[l]) * mod[:pfdc0][l] for l in list[n].dc_branches), init=0.0))
        add_to_expression!(inj_p[n], sum((mod[:pr][d] for d in list[n].renewables), init=0.0))
        add_to_expression!(inj_p[n], -1, sum((mod[:pd][d] for d in list[n].demands), init=0.0))
        add_to_expression!(inj_p[n], -1, sum((mod[:pr0][d] for d in list[n].renewables), init=0.0))
        add_to_expression!(inj_p[n], sum((mod[:ls0][d] for d in list[n].demands), init=0.0))
    end
    @constraint(mod, inj_p .== 0.0)
    ptdf_p = @expression(mod, ptdf * p)
    @constraint(mod, ptdf_p .>= -oplim.branch_rating * oplim.short_term_multi)
    @constraint(mod, ptdf_p .<= oplim.branch_rating * oplim.short_term_multi)
end

function add_contingencies!(Pc::Dict{<:Integer,ExprC}, opf::OPFsystem, oplim::Oplimits, mod::Model, obj::AbstractJuMPScalar, islands::Vector,
    island::Integer, ptdf::AbstractMatrix{<:Real}, list::Vector{<:CTypes{Int}}, c::Integer
)
    pc = JuMP.@variable(mod, [n in 1:length(opf.nodes)], base_name = @sprintf("pc%s", c))
    pgc, prc, lsc = init_P!(Pc, opf, oplim, mod, obj, list, islands, island, c)
    Pc[c].pc = pc

    inj_pc = @expression(mod, -pc)
    for n = 1:length(opf.nodes)
        add_to_expression!(inj_pc[n], sum((mod[:pg0][g] for g in list[n].ctrl_generation), init=0.0))
        add_to_expression!(inj_pc[n], sum((beta(opf.nodes[n], opf.dc_branches[l]) * mod[:pfdc0][l] for l in list[n].dc_branches), init=0.0))
        add_to_expression!(inj_pc[n], sum((mod[:pr][d] for d in list[n].renewables), init=0.0))
        add_to_expression!(inj_pc[n], -1, sum((mod[:pd][d] for d in list[n].demands), init=0.0))
        add_to_expression!(inj_pc[n], -1, sum((pgc[g] for g in list[n].ctrl_generation), init=0.0))
        add_to_expression!(inj_pc[n], -1, sum((prc[d] for d in list[n].renewables), init=0.0))
        add_to_expression!(inj_pc[n], sum((lsc[d] for d in list[n].demands), init=0.0))
    end
    @constraint(mod, inj_pc .== 0.0)
    ptdf_pc = @expression(mod, ptdf * pc)
    @constraint(mod, ptdf_pc .>= -oplim.branch_rating * oplim.short_term_multi)
    @constraint(mod, ptdf_pc .<= oplim.branch_rating * oplim.short_term_multi)
end

function add_contingencies!(Pcc::Dict{<:Integer,ExprCC}, opf::OPFsystem, oplim::Oplimits, mod::Model, obj::AbstractJuMPScalar, islands::Vector,
    island::Integer, ptdf::AbstractMatrix{<:Real}, list::Vector{<:CTypes{Int}}, c::Integer
)
    pcc = JuMP.@variable(mod, [n in 1:length(opf.nodes)], base_name = @sprintf("pcc%s", c))
    pgu, pgd, pfdccc, prcc, lscc = init_P!(Pcc, opf, oplim, mod, obj, list, islands, island, c)
    Pcc[c].pcc = pcc

    inj_pcc = @expression(mod, -pcc)
    for n = 1:length(opf.nodes)
        add_to_expression!(inj_pcc[n], sum((mod[:pg0][g] for g in list[n].ctrl_generation), init=0.0))
        add_to_expression!(inj_pcc[n], sum((mod[:pr][d] for d in list[n].renewables), init=0.0))
        add_to_expression!(inj_pcc[n], -1, sum((mod[:pd][d] for d in list[n].demands), init=0.0))
        add_to_expression!(inj_pcc[n], sum((pgu[g] - pgd[g] for g in list[n].ctrl_generation), init=0.0))
        add_to_expression!(inj_pcc[n], sum((beta(opf.nodes[n], opf.dc_branches[l]) * pfdccc[l] for l in list[n].dc_branches), init=0.0))
        add_to_expression!(inj_pcc[n], -1, sum((prcc[d] for d in list[n].renewables), init=0.0))
        add_to_expression!(inj_pcc[n], sum((lscc[d] for d in list[n].demands), init=0.0))
    end
    @constraint(mod, inj_pcc .== 0.0)
    ptdf_pcc = @expression(mod, ptdf * pcc)
    @constraint(mod, ptdf_pcc .>= -oplim.branch_rating * oplim.long_term_multi)
    @constraint(mod, ptdf_pcc .<= oplim.branch_rating * oplim.long_term_multi)
end

function add_contingencies!(Pccx::Dict{<:Integer,ExprCCX}, opf::OPFsystem, oplim::Oplimits, mod::Model, obj::AbstractJuMPScalar, islands::Vector,
    island::Integer, ptdf::AbstractMatrix{<:Real}, list::Vector{<:CTypes{Int}}, c::Integer
)
    pccx = JuMP.@variable(mod, [n in 1:length(opf.nodes)], base_name = @sprintf("pccx%s", c))
    pgd, prcc, lscc = init_P!(Pccx, opf, oplim, mod, obj, list, islands, island, c)
    Pccx[c].pccx = pccx

    inj_pccx = @expression(mod, -pccx)
    for n = 1:length(opf.nodes)
        add_to_expression!(inj_pccx[n], sum((mod[:pg0][g] for g in list[n].ctrl_generation), init=0.0))
        add_to_expression!(inj_pccx[n], sum((mod[:pr][d] for d in list[n].renewables), init=0.0))
        add_to_expression!(inj_pccx[n], -1, sum((mod[:pd][d] for d in list[n].demands), init=0.0))
        add_to_expression!(inj_pccx[n], sum((-pgd[g] for g in list[n].ctrl_generation), init=0.0))
        add_to_expression!(inj_pccx[n], -1, sum((prcc[d] for d in list[n].renewables), init=0.0))
        add_to_expression!(inj_pccx[n], sum((lscc[d] for d in list[n].demands), init=0.0))
    end
    @constraint(mod, inj_pccx .== 0.0)
    @constraint(mod, ptdf * pccx .>= -oplim.branch_rating)
    @constraint(mod, ptdf * pccx .<= oplim.branch_rating)
end


function init_P!(Pc::Dict{<:Integer,ExprC}, opf::OPFsystem, oplim::Oplimits, mod::Model, obj::JuMP.AbstractJuMPScalar, list::Vector{<:CTypes{Int}},
    islands::Vector, island::Integer, c::Integer
)
    pgc = JuMP.@variable(mod, [g in 1:length(opf.ctrl_generation)], base_name = @sprintf("pgc%s", c),
        lower_bound = 0.0)
    prc = JuMP.@variable(mod, [d in 1:length(opf.renewables)], base_name = @sprintf("prc%s", c),
        lower_bound = 0.0, upper_bound = oplim.pr_lim[d] * oplim.max_curtail)
    lsc = JuMP.@variable(mod, [d in 1:length(opf.demands)], base_name = @sprintf("lsc%s", c),
        lower_bound = 0.0, upper_bound = oplim.pd_lim[d] * oplim.max_shed)
    Pc[c] = ExprC(JuMP.VariableRef[], pgc, prc, lsc)

    p_survive = 1.0 - oplim.p_failure
    add_to_expression!(obj, opf.prob[c], sum(opf.voll' * lsc))
    add_to_expression!(obj, opf.prob[c], sum(oplim.ramp_mult * 1 * prc))
    for (cost, g) in zip(opf.cost_ctrl_gen, pgc)
        add_to_expression!.(obj, opf.prob[c], p_survive * oplim.ramp_mult * (cost.var * g))
        # add_to_expression!.(obj, opf.prob[c], p_survive * oplim.ramp_mult * (cost[1] * g^2 + cost[2] * g))
    end

    # Add new constraints that limit the corrective variables within operating limits
    balance_pc = JuMP.@expression(mod, sum(lsc))
    add_to_expression!.(balance_pc, -prc)
    add_to_expression!.(balance_pc, -pgc)
    JuMP.@constraint(mod, balance_pc == 0.0)
    if isempty(islands)
        JuMP.@constraint(mod, mod[:pg0] .- pgc .>= 0)
    else
        itr = length(islands[island]) < 2 ? Int[] : islands[island]
        for n in itr
            JuMP.@constraint(mod, [g = list[n].ctrl_generation], mod[:pg0][g] - pgc[g] >= 0.0)
        end
        itr = isempty(itr) ? (1:length(opf.nodes)) : islands[1:end.!=island]
        for in_vec in itr
            for n in in_vec
                JuMP.@constraint(mod, [g = list[n].ctrl_generation], mod[:pg0][g] - pgc[g] == 0.0)
                fix!(prc, oplim.pr_lim, list[n].renewables)
                fix!(lsc, oplim.pd_lim, list[n].demands)
            end
        end
    end
    return pgc, prc, lsc
end

function init_P!(Pcc::Dict{<:Integer,ExprCC}, opf::OPFsystem, oplim::Oplimits, mod::Model, obj::JuMP.AbstractJuMPScalar, list::Vector{<:CTypes{Int}},
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
        lower_bound = 0.0, upper_bound = oplim.pr_lim[r] * oplim.max_curtail)
    lscc = JuMP.@variable(mod, [d in 1:length(opf.demands)], base_name = @sprintf("lscc%s", c),
        lower_bound = 0.0, upper_bound = oplim.pd_lim[d] * oplim.max_shed)
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
            mod[:pg0] .+ pgu .- pgd .>= 0.0
            mod[:pg0] .+ pgu .- pgd .<= oplim.pg_lim_max
        end)
    else
        itr = length(islands[island]) < 2 ? Int[] : islands[island]
        for n in itr
            JuMP.@constraints(mod, begin
                [g = list[n].ctrl_generation], mod[:pg0][g] + pgu[g] - pgd[g] >= 0.0
                [g = list[n].ctrl_generation], mod[:pg0][g] + pgu[g] - pgd[g] <= oplim.pg_lim_max[g]
            end)
        end
        itr = isempty(itr) ? (1:length(opf.nodes)) : islands[1:end.!=island]
        for in_vec in itr
            for n in in_vec
                JuMP.@constraint(mod, [g = list[n].ctrl_generation], mod[:pg0][g] - pgd[g] == 0.0)
                fix!(pgu, list[n].ctrl_generation)
                fix!(pfdccc, list[n].dc_branches)
                fix!(prcc, oplim.pr_lim, list[n].renewables)
                fix!(lscc, oplim.pd_lim, list[n].demands)
            end
        end
    end
    return pgu, pgd, pfdccc, prcc, lscc
end

function init_P!(Pccx::Dict{<:Integer,ExprCCX}, opf::OPFsystem, oplim::Oplimits, mod::Model, obj::JuMP.AbstractJuMPScalar, list::Vector{<:CTypes{Int}},
    islands::Vector, island::Integer, c::Integer
)
    # Add corrective variables
    pgd = JuMP.@variable(mod, [g in 1:length(opf.ctrl_generation)], base_name = @sprintf("pgdx%s", c),
        lower_bound = 0.0)#, upper_bound = oplim.rampdown[g] * oplim.ramp_minutes)
    # and ramp down
    prcc = JuMP.@variable(mod, [r in 1:length(opf.renewables)], base_name = @sprintf("prccx%s", c),
        lower_bound = 0.0, upper_bound = oplim.pr_lim[d] * oplim.max_curtail)
    lscc = JuMP.@variable(mod, [d in 1:length(opf.demands)], base_name = @sprintf("lsccx%s", c),
        lower_bound = 0.0, upper_bound = oplim.pd_lim[d] * oplim.max_shed)
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
            JuMP.@constraint(mod, [g = list[n].ctrl_generation], mod[:pg0][g] - pgd[g] >= 0.0)
        end
        itr = isempty(itr) ? (1:length(opf.nodes)) : islands[1:end.!=island]
        for in_vec in itr
            for n in in_vec
                JuMP.@constraint(mod, [g = list[n].ctrl_generation], mod[:pg0][g] - pgd[g] == 0.0)
                JuMP.fix.(prcc, oplim.pr_lim; force=true)
                JuMP.fix.(lscc, oplim.pd_lim; force=true)
            end
        end
    end
    return pgd, prcc, lscc
end

fix!(var::AbstractVector{JuMP.VariableRef}, val::AbstractVector{<:Real}, vec::AbstractVector) =
    JuMP.fix.(var[vec], val[vec]; force=true)
fix!(var::AbstractVector{JuMP.VariableRef}, vec::AbstractVector) =
    JuMP.fix.(var[vec], 0.0; force=true)

fix_values!(mod::Model, symb::Symbol) = JuMP.fix.(mod[symb], get_value(mod, symb), force=true)
fix_values!(mod::Model, var::AbstractVector{JuMP.VariableRef}) = JuMP.fix.(var, get_value(mod, var), force=true)

function fix_base_case!(mod::Model)
    fix_values!(mod, :pg0)
    fix_values!(mod, :pfdc0)
    fix_values!(mod, :ls0)
    fix_values!(mod, :pr0)
end

function fix_short_term!(mod::Model, Pc::Dict{<:Integer,ExprC})
    for (_, c) in Pc
        fix_values!(mod, c.pgc)
        fix_values!(mod, c.prc)
        fix_values!(mod, c.lsc)
    end
end

function fix_long_term!(mod::Model, Pcc::Dict{<:Integer,ExprCC})
    for (_, c) in Pcc
        fix_values!(mod, c.pgu)
        fix_values!(mod, c.pgd)
        fix_values!(mod, c.pfdccc)
        fix_values!(mod, c.prcc)
        fix_values!(mod, c.lscc)
    end
end

calc_cens(mod::Model, opf::OPFsystem, var::AbstractVector{JuMP.VariableRef}) = sum(opf.voll .* get_value(mod, var))
calc_ctrl_cost(mod::Model, opf::OPFsystem, var::AbstractVector{JuMP.VariableRef}) =
    sum(c[2] * g for (c, g) in zip(opf.cost_ctrl_gen, get_value(mod, var)))
    # sum(c[1] * g^2 + c[2] * g for (c, g) in zip(opf.cost_ctrl_gen, get_value(mod, var)))

calc_objective(mod::Model, opf::OPFsystem) =
    calc_ctrl_cost(mod, opf, mod[:pg0]) + calc_cens(mod, opf, mod[:ls0])
calc_objective(mod::Model, opf::OPFsystem, expr::ExprC, ramp_mult=1.0) =
    ramp_mult * calc_ctrl_cost(mod, opf, expr.pgc) + calc_cens(mod, opf, expr.lsc)
calc_objective(mod::Model, opf::OPFsystem, expr::ExprCC, ramp_mult=1.0) =
    ramp_mult * (calc_ctrl_cost(mod, opf, expr.pgu) + calc_ctrl_cost(mod, opf, expr.pgd)) + calc_cens(mod, opf, expr.lscc)
calc_objective(mod::Model, opf::OPFsystem, expr::ExprCCX, ramp_mult=1.0) =
    ramp_mult * 60 * calc_ctrl_cost(mod, opf, expr.pgdx) + calc_cens(mod, opf, expr.lsccx)
calc_objective(mod::Model, opf::OPFsystem, P::Dict, ramp_mult=1.0) =
    sum((opf.prob[i] * calc_objective(mod, opf, c, ramp_mult) for (i, c) in P), init=0.0)
calc_objective(mod::Model, opf::OPFsystem, Pc::Dict{<:Integer,ExprC},
    Pcc::Dict{<:Integer,ExprCC}, Pccx::Dict{<:Integer,ExprCCX}, ramp_mult=1.0
) = calc_objective(mod, opf) + calc_objective(mod, opf, Pc, ramp_mult) +
    calc_objective(mod, opf, Pcc, ramp_mult) + calc_objective(mod, opf, Pccx, ramp_mult)
