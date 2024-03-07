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
    max_curtail::Union{TR,Vector{TR}}
    dc_lim_min::Vector{TR}
    dc_lim_max::Vector{TR}
    pd_lim::Vector{TR}
    max_shed::Union{TR,Vector{TR}}
end

""" Constructor for Oplimits """
function oplimits(
    system::System,
    max_shed::Union{TR,Vector{TR}},
    max_curtail::Union{TR,Vector{TR}},
    dist_slack::Vector{TR},
    ramp_mult::TR,
    ramp_minutes::TR,
    p_failure::TR,
    short_term_multi::Union{TR,Vector{TR}},
    long_term_multi::Union{TR,Vector{TR}}
) where {TR<:Real}
    generation = sort_components!(get_generation(system))
    branches = sort_components!(get_branches(system))
    dc_branches = sort_components!(get_dc_branches(system))
    demands = sort_components!(get_demands(system))

    branch_rating::Vector{TR} = get_rate.(branches)
    (pg_lim_min::Vector{TR}, pg_lim_max::Vector{TR}) = split_pair(get_active_power_limits.(generation))
    (rampup::Vector{TR}, rampdown::Vector{TR}) = split_pair(get_ramp_limits.(generation))
    (dc_lim_min::Vector{TR}, dc_lim_max::Vector{TR}) = split_pair(get_active_power_limits_from.(dc_branches))
    pd_lim::Vector{TR} = get_active_power.(demands)
    if !isempty(dist_slack)
        dist_slack = dist_slack / sum(dist_slack)
    end
    return Oplimits{TR}(ramp_mult, ramp_minutes, p_failure, branch_rating, short_term_multi, long_term_multi, 
        pg_lim_min, pg_lim_max, rampup, rampdown, dist_slack, max_curtail, dc_lim_min, dc_lim_max, pd_lim, max_shed)
end

mutable struct Case{TR<:Real,TI<:Integer}
    model::Model
    opf::OPFsystem{TR}
    pf::DCPowerFlow{TR,TI}
    oplim::Oplimits{TR}
    brc_up::Dict{TI,ConstraintRef}
    brc_down::Dict{TI,ConstraintRef}
    Pc::Dict{TI,ExprC}
    Pcc::Dict{TI,ExprCC}
    Pccx::Dict{TI,ExprCCX}
end

""" Initialize an OPF of a power system """
function opf_base(type::OPF, system::System, optimizer;
    voll::Vector{<:Number}=Float64[],
    prob::Vector{<:Number}=Float64[],
    contingencies::Vector{Pair{String, Int64}}=Pair{String,Int64}[],
    dist_slack::Vector{<:Real}=Float64[],
    time_limit_sec::Int64=600,
    ramp_minutes::Real=10.0,
    ramp_mult::Real=10.0,
    renew_cost::Number=0.0,
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
    opf = isempty(voll) ? opfsystem(system) : opfsystem(system, voll, contingencies, prob; renew_cost=renew_cost)
    pf = DCPowerFlow(system)
    oplim = oplimits(system, max_shed, max_curtail, dist_slack, ramp_mult, ramp_minutes, p_failure, short_term_multi, long_term_multi)

    return add_base_constraints!(type, opf, oplim, m, pf)
end

function add_base_constraints!(type::OPF, opf::OPFsystem, oplim::Oplimits, m::Model, pf::DCPowerFlow)
    brc_up = Dict{Int64,ConstraintRef}()
    brc_down = Dict{Int64,ConstraintRef}()
    Pc = Dict{Int64,ExprC}()
    Pcc = Dict{Int64,ExprCC}()
    Pccx = Dict{Int64,ExprCCX}()


    @variables(m, begin
        p0[n in 1:size(opf.mgx, 2)]
        # active power injection on each node
        oplim.pg_lim_min[i] <= pg0[i in 1:size(opf.mgx, 1)] <= oplim.pg_lim_max[i]
        # active power variables for the generators
        oplim.dc_lim_min[i] <= pfdc0[i in 1:size(opf.mdcx, 1)] <= oplim.dc_lim_max[i]
        # power flow on DC branches
        pd[d in 1:size(opf.mdx, 1)]
        0.0 <= ls0[i in 1:size(opf.mdx, 1)] <= oplim.pd_lim[i] * (typeof(oplim.max_shed) <: Real ? 1.0 : oplim.max_shed[i])
        # demand curtailment variables
    end)

    @objective(m, Min, sum(c.var * g for (c, g) in zip(opf.cost_ctrl_gen, pg0)) + sum(opf.voll' * ls0))
    # @objective(m, Min, sum(c[1] * g^2 + c[2] * g for (c, g) in zip(opf.cost_ctrl_gen, pg0)) + sum(opf.voll' * ls0) 

    JuMP.fix.(pd, oplim.pd_lim)
    if typeof(oplim.max_shed) <: Real
        @constraint(m, sum_max_shed, sum(ls0) <= oplim.max_shed)
    end

    @expression(m, inj_p0, opf.mgx' * pg0)
    add_to_expression!.(inj_p0, opf.mdcx' * pfdc0)
    add_to_expression!.(inj_p0, opf.mdx' * (ls0 - pd))
    @constraint(m, injected_power, inj_p0 .== p0)

    # add_branch_constraints!(m, pf.ϕ, p0, oplim.branch_rating)
    @expression(m, balance, sum(pg0, init=0.0))
    add_to_expression!.(balance, -pd)
    add_to_expression!.(balance, ls0)
    @constraint(m, power_balance, balance == 0.0)

    if any([type.P, type.C1, type.C2, type.C2F])
        return add_all_contingencies!(type, opf, oplim, m, pf, brc_up, brc_down, Pc, Pcc, Pccx)
    end
    return m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx
end

" Add all constraints to the model "
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
    return m
end

" Add a branch constraint for branch to the model "
function add_branch_constraint!(m::Model, pf::DCPowerFlow, p::AbstractVector{VariableRef}, 
    brc_up::Dict{<:Integer, ConstraintRef}, brc_down::Dict{<:Integer, ConstraintRef}, branch::Integer, rating::Real
)
    # ptdf = calc_isf_vec(pf, branch)
    ptdf = view(pf.ϕ, branch, :)
    # ptdf0 = GenericAffExpr(0.0, Pair.(p, ptdf[i,:])) 
    ptdf0 = @expression(m, AffExpr())
    for (i,j) in zip(ptdf, p)
        add_to_expression!(ptdf0, i, j)
    end
    brc_down[branch] = @constraint(m, ptdf0 + rating >= 0.0)
    brc_up[branch] = @constraint(m, ptdf0 - rating <= 0.0)
    return m
end

" Add branch limits to overloaded branches "
function constrain_branches!(m::Model, pf::DCPowerFlow, oplim::Oplimits, brc_up::Dict{<:Integer, ConstraintRef}, brc_down::Dict{<:Integer, ConstraintRef}, 
    total_solve_time::Real, atol::Real=1e-6
)
    # if !has_values(m)
        # Note: While using a direct_model, this check fails after the model is modified for some solvers
        total_solve_time = update_model!(m, pf, total_solve_time)
    # end
    while true
        ol_br = find_overloaded_branches(pf.F, oplim.branch_rating, atol)
        isempty(ol_br) && break 
        termination_status(m) != JuMP.OPTIMAL && break
        for br in ol_br
            add_branch_constraint!(m, pf, m[:p0], brc_up, brc_down, br, oplim.branch_rating[br])
            @info "Branch $br added"
        end
        total_solve_time = update_model!(m, pf, total_solve_time)
    end
    return total_solve_time
end
constrain_branches!(case::Case, total_solve_time::Real, atol::Real=1e-6) = 
    constrain_branches!(case.model, case.pf, case.oplim, case.brc_up, case.brc_down, total_solve_time, atol)

""" Solve model and update the power flow object """
function update_model!(m::Model, pf::DCPowerFlow, total_solve_time::Real)
    solve_model!(m)
    total_solve_time += solve_time(m)
    Pᵢ = get_value(m, :p0)
    calc_θ!(pf, Pᵢ)
    calc_Pline!(pf)
    return total_solve_time
end

function add_all_contingencies!(type::OPF, opf::OPFsystem, oplim::Oplimits, m::Model,
    pf::DCPowerFlow, brc_up::Dict{<:Integer, ConstraintRef}, brc_down::Dict{<:Integer, ConstraintRef}, 
    Pc::Dict{<:Integer,ExprC}, Pcc::Dict{<:Integer,ExprCC}, Pccx::Dict{<:Integer,ExprCCX}
)

    obj = objective_function(m)
    !isempty(oplim.dist_slack) && set_dist_slack!(pf.ϕ, opf.mgx, oplim.dist_slack)
    for (c, c_obj) in enumerate(opf.contingencies)
        c_i = c_obj.second
        c_n = c_obj.first == "branch" ? Tuple(opf.mbx[c_i,:].nzind) : opf.mgx[c_i,:].nzind[1]
        if is_islanded(pf, c_n, c_i)
            islands = island_detection_thread_safe(pf.B, c_n)
            for (n,nodes) in enumerate(islands)
                if length(nodes) < 2
                    deleteat!(islands, n)
                end
            end
            islands_b = [find_island_branches(island, pf.DA, c_i) for island in islands]
            ptdf = [calc_isf(pf, c_n, c_i, islands, n, ibranches) for (n, ibranches) in enumerate(islands_b)]
        else
            islands = [Vector{Int64}()]
            ptdf = [calc_isf(pf, c_n, c_i)]
            set_tol_zero!.(ptdf)
        end
        type.P && add_contingency!(opf, oplim, m, brc_up, brc_down, ptdf, c)
        type.C1 && add_contingency!(Pc, opf, oplim, m, brc_up, brc_down, obj, islands, ptdf, c)
        type.C2 && add_contingency!(Pcc, opf, oplim, m, brc_up, brc_down, obj, islands, ptdf, c)
        type.C2F && add_contingency!(Pccx, opf, oplim, m, brc_up, brc_down, obj, islands, ptdf, c)

        @debug "Contingency $(get_name(c_obj)) is added"
    end
    set_objective_function(m, obj)
    return m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx
end
add_all_contingencies!(type::OPF, case::Case) = 
    add_all_contingencies!(type, case.opf, case.opflim, case.model, case.pf, case.brc_up, case.brc_down, case.Pc, case.Pcc, case.Pccx)

function add_contingency!(opf::OPFsystem, oplim::Oplimits, m::Model, 
    brc_up::Dict{<:Integer, ConstraintRef}, brc_down::Dict{<:Integer, ConstraintRef}, 
    ptdf::Vector{<:AbstractMatrix{<:Real}}, c::Integer
)
    p = @variable(m, [n in 1:size(opf.mgx, 2)], base_name = "p"*string(c))
    @constraint(m, m[:inj_p0] .== p)
    for x in ptdf 
        add_branch_constraints!(m, x, p, brc_up, brc_down, oplim.branch_rating * oplim.short_term_multi)
    end
    # add_branch_constraints!(m, opf, pf, p, oplim.branch_rating * oplim.short_term_multi, c)
end

function add_contingency!(Pc::Dict{<:Integer,ExprC}, opf::OPFsystem, oplim::Oplimits, m::Model, 
    brc_up::Dict{<:Integer, ConstraintRef}, brc_down::Dict{<:Integer, ConstraintRef}, obj::AbstractJuMPScalar, islands::Vector,
    ptdf::Vector{<:AbstractMatrix{<:Real}}, c::Integer
)
    pc = @variable(m, [n in 1:size(opf.mgx, 2)], base_name = "pc"*string(c))
    pgu, pgd, lsc = init_P!(Pc, opf, oplim, m, obj, islands, c)

    Pc[c].pc = pc

    inj_pc = @expression(m, opf.mgx' * (m[:pg0] + pgu - pgd))
    add_to_expression!.(inj_pc, opf.mdcx' * m[:pfdc0])
    add_to_expression!.(inj_pc, opf.mdx' * (lsc - m[:pd]))
    @constraint(m, inj_pc .== pc)

    for x in ptdf 
        add_branch_constraints!(m, x, pc, brc_up, brc_down, oplim.branch_rating * oplim.short_term_multi)
    end
    # add_branch_constraints!(m, opf, pf, pc, oplim.branch_rating * oplim.short_term_multi, c)
end

function add_contingency!(Pcc::Dict{<:Integer,ExprCC}, opf::OPFsystem, oplim::Oplimits, m::Model, 
    brc_up::Dict{<:Integer, ConstraintRef}, brc_down::Dict{<:Integer, ConstraintRef}, obj::AbstractJuMPScalar, islands::Vector,
    ptdf::Vector{<:AbstractMatrix{<:Real}}, c::Integer
)
    pcc = @variable(m, [n in 1:size(opf.mgx, 2)], base_name = "pcc"*string(c))
    pgu, pgd, pfdccc, lscc = init_P!(Pcc, opf, oplim, m, obj, islands, c)
    Pcc[c].pcc = pcc

    inj_pcc = @expression(m, opf.mgx' * (m[:pg0] + pgu - pgd))
    add_to_expression!.(inj_pcc, opf.mdcx' * pfdccc)
    add_to_expression!.(inj_pcc, opf.mdx' * (lscc - m[:pd]))
    @constraint(m, inj_pcc .== pcc)

    for x in ptdf 
        add_branch_constraints!(m, x, pcc, brc_up, brc_down, oplim.branch_rating * oplim.long_term_multi)
    end
    # add_branch_constraints!(m, opf, pf, pcc, oplim.branch_rating * oplim.long_term_multi, c)
end

function add_contingency!(Pccx::Dict{<:Integer,ExprCCX}, opf::OPFsystem, oplim::Oplimits, m::Model, 
    brc_up::Dict{<:Integer, ConstraintRef}, brc_down::Dict{<:Integer, ConstraintRef}, obj::AbstractJuMPScalar, islands::Vector,
    ptdf::Vector{<:AbstractMatrix{<:Real}}, c::Integer
)
    pccx = @variable(m, [n in 1:size(opf.mgx, 2)], base_name = "pccx"*string(c))
    pgd, lscc = init_P!(Pccx, opf, oplim, m, obj, islands, c)
    Pccx[c].pccx = pccx

    inj_pccx = @expression(m, opf.mgx' * (m[:pg0] - pgd))
    add_to_expression!.(inj_pccx, opf.mdcx' * m[:pfdc0])
    add_to_expression!.(inj_pccx, opf.mdx' * (lscc - m[:pd]))
    @constraint(m, inj_pccx .== pccx)

    for x in ptdf 
        add_branch_constraints!(m, x, pccx, brc_up, brc_down, oplim.branch_rating * oplim.long_term_multi)
    end
    # add_branch_constraints!(m, opf, pf, pccx, oplim.branch_rating * oplim.long_term_multi, c)
end


function init_P!(Pc::Dict{<:Integer,ExprC}, opf::OPFsystem, oplim::Oplimits, m::Model, obj::AbstractJuMPScalar, 
    islands::Vector, c::Integer
)
    pgu = @variable(m, [g in 1:length(m[:pg0])], base_name = "pguc"*string(c),
        lower_bound = 0.0, upper_bound = oplim.rampup[g] * 0.0)
    # active power variables for the generators in contingencies ramp up 
    pgd = @variable(m, [g in 1:length(m[:pg0])], base_name = "pgdc"*string(c),
        lower_bound = 0.0) #, upper_bound = rampdown[g] * 0.0)
    # and ramp down
    lsc = @variable(m, [d in 1:length(m[:pd])], base_name = "lsc"*string(c),
        lower_bound = 0.0, upper_bound = oplim.pd_lim[d] * (typeof(oplim.max_shed) <: Real ? 1.0 : oplim.max_shed[d]))
    Pc[c] = ExprC(VariableRef[], pgu, pgd, lsc)

    p_survive = 1.0 - oplim.p_failure
    add_to_expression!(obj, opf.prob[c], sum(opf.voll' * lsc))
    for (cost, gu, gd) in zip(opf.cost_ctrl_gen, pgu, pgd)
        add_to_expression!(obj, opf.prob[c] * p_survive * cost.var * oplim.ramp_mult, gu)
        add_to_expression!(obj, opf.prob[c] * p_survive * cost.var * oplim.ramp_mult, gd)
        # add_to_expression!.(obj, opf.prob[c], p_survive * oplim.ramp_mult * (cost[1] * g^2 + cost[2] * g))
    end

    # Add new constraints that limit the corrective variables within operating limits
    balance_pc = @expression(m, sum(lsc))
    add_to_expression!.(balance_pc, pgu)
    add_to_expression!.(balance_pc, -pgd)
    @constraint(m, balance_pc == 0.0)
    if isempty(islands[1])
        @constraint(m, m[:pg0] .+ pgu .- pgd .>= oplim.pg_lim_min)
        # @constraint(m, m[:pg0] .+ pgu .- pgd .<= oplim.pg_lim_max)
        if typeof(oplim.max_shed) <: Real
            @constraint(m, sum(lsc) <= oplim.max_shed)
        end
    else
        for itr in islands, n in itr, g = opf.mgx[:,n].nzind
            @constraints(m, begin
                m[:pg0][g] + pgu[g] - pgd[g] >= oplim.pg_lim_min[g]
                # m[:pg0][g] + pgu[g] - pgd[g] <= oplim.pg_lim_max[g]
            end)
        end
        if typeof(oplim.max_shed) <: Real
            expr = AffExpr()
            for itr in islands, n in itr, d in opf.mdx[:,n].nzind
                add_to_expression!(expr, lsc[d])
            end
            @constraint(m, expr <= oplim.max_shed)
        end
        itr = 1:size(opf.mgx, 2)
        for island in islands
            itr = setdiff(itr, island)
        end
        for in_vec in itr, n in in_vec
            for g in opf.mgx[:,n].nzind
                set_upper_bound(pgu[g], oplim.pg_lim_max[g] - oplim.pg_lim_min[g])
                @constraint(m, m[:pg0][g] + pgu[g] - pgd[g] == 0.0)
            end
            fix!(lsc, oplim.pd_lim, opf.mdx[:,n].nzind)
        end
    end
    return pgu, pgd, lsc
end

function init_P!(Pcc::Dict{<:Integer,ExprCC}, opf::OPFsystem, oplim::Oplimits, m::Model, obj::AbstractJuMPScalar, 
    islands::Vector, c::Integer
)
    # Add corrective variables
    pgu = @variable(m, [g in 1:length(m[:pg0])], base_name = @sprintf("pgucc%s", c),
        lower_bound = 0.0, upper_bound = oplim.rampup[g] * oplim.ramp_minutes)
    # active power variables for the generators in contingencies ramp up 
    pgd = @variable(m, [g in 1:length(m[:pg0])], base_name = @sprintf("pgdcc%s", c),
        lower_bound = 0.0) #, upper_bound = rampdown[g] * ramp_minutes)
    # and ramp down
    pfdccc = @variable(m, [d in 1:length(m[:pfdc0])], base_name = @sprintf("pfdccc%s", c),
        lower_bound = oplim.dc_lim_min[d], upper_bound = oplim.dc_lim_max[d])
    lscc = @variable(m, [d in 1:length(m[:pd])], base_name = @sprintf("lscc%s", c),
        lower_bound = 0.0, upper_bound = oplim.pd_lim[d] * (typeof(oplim.max_shed) <: Real ? 1.0 : oplim.max_shed[d]))
    # load curtailment variables in in contingencies
    Pcc[c] = ExprCC(VariableRef[], pgu, pgd, pfdccc, lscc)

    p_survive = 1.0 - oplim.p_failure
    # Extend the objective with the corrective variables
    add_to_expression!.(obj, opf.prob[c],
        # (1.0 - p_failure) * (sum(opf.voll' * lscc) + sum(60 * (pgu + pgd)) # + # TODO: remove 60 and uncomment next lines for non-4-area analysis!!!!
        p_survive * (sum(opf.voll' * lscc) # + sum(opf.cost_ctrl_gen' * ramp_mult * (pgu + pgd))
        # (sum(opf.voll' * lscc) +
        # sum(opf.cost_ctrl_gen' * ramp_mult * pgu) # +
        # sum(opf.cost_ctrl_gen' * pgd)
        ))
    for (cost, gu, gd) in zip(opf.cost_ctrl_gen, pgu, pgd)
        add_to_expression!(obj, opf.prob[c] * p_survive * cost.var * oplim.ramp_mult, gu)
        add_to_expression!(obj, opf.prob[c] * p_survive * cost.var * oplim.ramp_mult, gd)
        # add_to_expression!.(obj, opf.prob[c], p_survive * oplim.ramp_mult * (cost[1] * gu^2 + cost[2] * gu))
        # add_to_expression!.(obj, opf.prob[c], p_survive * oplim.ramp_mult * (cost[1] * gd^2 + cost[2] * gd))
    end

    # Add new constraints that limit the corrective variables within operating limits
    balance_pcc = @expression(m, sum(pgu))
    add_to_expression!.(balance_pcc, lscc)
    add_to_expression!.(balance_pcc, -pgd)
    @constraint(m, balance_pcc == 0.0)
    if isempty(islands[1])
        @constraints(m, begin
            m[:pg0] .+ pgu .- pgd .>= oplim.pg_lim_min
            m[:pg0] .+ pgu .- pgd .<= oplim.pg_lim_max
        end)
        if typeof(oplim.max_shed) <: Real
            @constraint(m, sum(lscc) <= oplim.max_shed)
        end
    else
        for itr in islands, n in itr, g in opf.mgx[:,n].nzind
            expr = @expression(m, m[:pg0][g] + pgu[g] - pgd[g])
            @constraint(m, expr >= oplim.pg_lim_min[g])
            @constraint(m, expr <= oplim.pg_lim_max[g])
        end
        if typeof(oplim.max_shed) <: Real
            expr = AffExpr()
            for itr in islands, n in itr, d in opf.mdx[:,n].nzind
                add_to_expression!(expr, lscc[d])
            end
            @constraint(m, expr <= oplim.max_shed)
        end
        itr = 1:size(opf.mgx, 2)
        for island in islands
            itr = setdiff(itr, island)
        end
        for in_vec in itr, n in in_vec
            for g in opf.mgx[:,n].nzind
                set_upper_bound(pgu[g], oplim.pg_lim_max[g] - oplim.pg_lim_min[g])
                @constraint(m, m[:pg0][g] + pgu[g] - pgd[g] == 0.0)
            end
            fix!(pfdccc, opf.mdcx[:,n].nzind)
            fix!(lscc, oplim.pd_lim, opf.mdx[:,n].nzind)
        end
    end
    return pgu, pgd, pfdccc, lscc
end

function init_P!(Pccx::Dict{<:Integer,ExprCCX}, opf::OPFsystem, oplim::Oplimits, m::Model, obj::AbstractJuMPScalar, 
    islands::Vector, c::Integer
)
    # Add corrective variables
    pgd = @variable(m, [g in 1:length(m[:pg0])], base_name = @sprintf("pgdx%s", c),
        lower_bound = 0.0)#, upper_bound = oplim.rampdown[g] * oplim.ramp_minutes)
    # and ramp down
    lscc = @variable(m, [d in 1:length(m[:pd])], base_name = @sprintf("lsccx%s", c),
        lower_bound = 0.0, upper_bound = oplim.pd_lim[d] * (typeof(oplim.max_shed) <: Real ? 1.0 : oplim.max_shed[d]))
    # load curtailment variables in in contingencies
    Pccx[c] = ExprCCX(VariableRef[], pgd, lscc)

    # Extend the objective with the corrective variables
    add_to_expression!(obj, opf.prob[c] * oplim.p_failure, sum(opf.voll' * lscc))
    add_to_expression!(obj, opf.prob[c] * oplim.p_failure, sum(60 * pgd))

    # Add new constraints that limit the corrective variables within operating limits
    balance_pccx = @expression(m, sum(lscc))
    add_to_expression!.(balance_pccx, -pgd)
    @constraint(m, balance_pccx == 0.0)
    if isempty(islands[1])
        @constraint(m, m[:pg0] .- pgd .>= 0.0)
    else
        for itr in islands, n in itr, g = opf.mgx[:,n].nzind
            @constraint(m, m[:pg0][g] - pgd[g] >= 0.0)
        end
        if typeof(oplim.max_shed) <: Real
            expr = AffExpr()
            for itr in islands, n in itr, d in opf.mdx[:,n].nzind
                add_to_expression!(expr, lscc[d])
            end
            @constraint(m, expr <= oplim.max_shed)
        end
        itr = 1:size(opf.mgx, 2)
        for island in islands
            itr = setdiff(itr, island)
        end
        for in_vec in itr, n in in_vec
            for g = opf.mgx[:,n].nzind
                @constraint(m, m[:pg0][g] - pgd[g] == 0.0)
            end
            fix!(lscc, oplim.pd_lim, opf.mdx[:,n].nzind)
        end
    end
    return pgd, lscc
end

""" Force a variable to equal a value"""
fix!(var::AbstractVector{VariableRef}, val::AbstractVector{<:Real}, vec::AbstractVector) =
    JuMP.fix.(var[vec], val[vec]; force=true)

""" Force a variable to equal zero"""
fix!(var::AbstractVector{VariableRef}, vec::AbstractVector) =
    JuMP.fix.(var[vec], 0.0; force=true)

""" Force a variable to equal its current value in the model """
fix_values!(m::Model, symb::Symbol) = JuMP.fix.(m[symb], get_value(m, symb), force=true)
fix_values!(m::Model, var::AbstractVector{VariableRef}) = JuMP.fix.(var, get_value(m, var), force=true)

""" Fix all base case varibles to its current values in the model """
function fix_base_case!(m::Model)
    fix_values!(m, :pg0)
    fix_values!(m, :pfdc0)
    fix_values!(m, :ls0)
end

""" Fix all contingency varibles to its current values in the model """
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
    if isempty(case.oplim.dist_slack)
        slack_cost = sum(getproperty.(opf.cost_ctrl_gen[opf.mgx[:,pf.slack].nzind], :var) * case.oplim.ramp_mult)
    else
        slack_cost = opf.mgx' * (getproperty.(opf.cost_ctrl_gen, :var) * case.oplim.ramp_mult .* case.oplim.dist_slack .* pg0)
    end
    pg0 = opf.mgx' * pg0
    pd = opf.mdx' * (get_value(m, :pd) - get_value(m, :ls0))
    ls_cost = opf.mdx' * (opf.voll' * (get_value(m, :pd) - get_value(m, :ls0)))
    for (i, c_obj) in enumerate(opf.contingencies)
        cont = typesort_component(c_obj, opf)
        if is_islanded(pf, cont[2], cont[1])
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont[2], cont[1], pf.slack)
            nodes = reduce(vcat, islands[1:end.!=island], init=[])
            ΔP = sum(pg0[nodes]) - sum(pd[nodes])
            costs[i] = sum(slack_cost * abs(ΔP)) + sum(ls_cost[nodes])
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
calc_objective(m::Model, opf::OPFsystem, expr::ExprC, ramp_mult::Real) =
    (calc_ctrl_cost(m, opf, expr.pgu, :var) + calc_ctrl_cost(m, opf, expr.pgd, :var)) * ramp_mult + calc_cens(m, opf, expr.lsc)

""" Calculate the long term state cost """
calc_objective(m::Model, opf::OPFsystem, expr::ExprCC, ramp_mult::Real) =
    (calc_ctrl_cost(m, opf, expr.pgu, :var) + calc_ctrl_cost(m, opf, expr.pgd, :var)) * ramp_mult + calc_cens(m, opf, expr.lscc)

""" Calculate the long term with corrective failure state cost """
calc_objective(m::Model, opf::OPFsystem, expr::ExprCCX, ramp_mult::Real) =
    60 * calc_ctrl_cost(m, opf, expr.pgdx, :var) * ramp_mult + calc_cens(m, opf, expr.lsccx)

""" Calculate a state cost """
calc_objective(m::Model, opf::OPFsystem, P::Dict) =
    sum((opf.prob[i] * calc_objective(m, opf, c) for (i, c) in P), init=0.0)

""" Calculate the total cost """
calc_objective(m::Model, opf::OPFsystem, Pc::Dict{<:Integer,ExprC},
    Pcc::Dict{<:Integer,ExprCC}, Pccx::Dict{<:Integer,ExprCCX}, ramp_mult::Real
) = calc_objective(m, opf) + calc_objective(m, opf, Pc, ramp_mult) +
    calc_objective(m, opf, Pcc, ramp_mult) + calc_objective(m, opf, Pccx, ramp_mult)

function print_costs(case::Case)
    base_cost = calc_objective(case.model, case.opf)
    slack_cost = sum(case.opf.prob .* calc_slack_cost(case))
    @printf "Base cost %.2f; " base_cost
    @printf "Slack cost %.2f; " slack_cost
    @printf "Total cost %.2f \n" base_cost + slack_cost
    return base_cost, slack_cost
end

function get_costs(case::Case)
    costs = SparseArrays.spzeros(length(case.opf.contingencies), 5)
    costs[:,1] .= calc_objective(case.model, case.opf)
    for c in keys(case.Pc)
        costs[c,2] = calc_objective(case.model, case.opf, case.Pc[c], case.oplim.ramp_mult)
    end
    for c in keys(case.Pcc)
        costs[c,3] = calc_objective(case.model, case.opf, case.Pcc[c], case.oplim.ramp_mult)
    end
    for c in keys(case.Pccx)
        costs[c,4] = calc_objective(case.model, case.opf, case.Pccx[c], case.oplim.ramp_mult)
    end
    costs[:,5] = calc_slack_cost(case)
    return costs, [:Base, :Pc, :Pcc, :Pccx, :distslack]
end