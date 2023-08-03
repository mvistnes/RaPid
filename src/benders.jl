# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022 

""" Benders type """
mutable struct Benders
    idx::Dict{<:Int, <:Int}
    list::Vector{CTypes}

    linerating::Vector
    pg_lim_min::Vector
    pg_lim_max::Vector
    pr_lim::Vector
    dc_lim_min::Vector
    dc_lim_max::Vector
    pd_lim::Vector
    rampup::Vector
    rampdown::Vector

    obj::AffExpr
    Pc::Dict{Int, NTuple{3, Any}} # Holds the short term variables for contingencies
    Pcc::Dict{Int, NTuple{5, Any}} # Holds the long term variables for contingencies
    Pccx::Dict{Int, NTuple{2, Any}} # Holds the long term variables for contingencies, no ramp up allowed

    Pᵢ::Vector{<:Real}
    Pg::Vector{<:Real}
    Pd::Vector{<:Real}
end

""" Constructor for OPFmodel """
function benders(opfm::OPFmodel)
    idx = get_nodes_idx(opfm.nodes)
    list = make_list(opfm, idx, opfm.nodes)
    # cgen = connectivitymatrix(system, length(nodes), idx)
    linerating = get_rate.(opfm.branches)
    (pg_lim_min, pg_lim_max) = split_pair(get_active_power_limits.(opfm.ctrl_generation))
    pr_lim = get_active_power.(opfm.renewables)
    (dc_lim_min, dc_lim_max) = split_pair(get_active_power_limits_from.(opfm.dc_branches))
    pd_lim = get_active_power.(opfm.demands)
    (rampup, rampdown) = split_pair(get_ramp_limits.(opfm.ctrl_generation))

    obj = objective_function(opfm.mod)
    Pc = Dict{Int, NTuple{3, Any}}() 
    Pcc = Dict{Int, NTuple{5, Any}}() 
    Pccx = Dict{Int, NTuple{2, Any}}() 
    # Pᵢ = get_net_Pᵢ(opfm, idx)
    Pᵢ = get_value(opfm.mod, :p0)
    Pg = get_controllable(opfm, idx)
    Pd = get_Pd(opfm, idx) # Fixed injection
    @assert isapprox(Pg, (Pᵢ - Pd); atol=1e-5) string(Pg - (Pᵢ - Pd))
    return Benders(idx, list, linerating, pg_lim_min, pg_lim_max, pr_lim, dc_lim_min, dc_lim_max, pd_lim, 
        rampup, rampdown, obj, Pc, Pcc, Pccx, Pᵢ, Pg, Pd)
end

""" 
Solve the optimization model using Benders decomposition.

Function creates extra production constraints on generators in the system
based on the power transfer distribution factors of the generators in the
system. 
"""
function run_benders(type::OPF, system::System, optimizer, voll=nothing, prob=nothing, contingencies=nothing; 
        ramp_minutes = 10, ramp_mult = 10, max_shed = 1.0, lim = 1e-6, 
        branch_short_term_limit_multi::Real = 1.5, branch_long_term_limit_multi::Real = 1.2, 
        gen_short_term_limit_multi::Real = 2.0, gen_long_term_limit_multi::Real = 1.0, p_failure = 0.0, branch_c = nothing, rate_c = 0.0, debug = false)
    total_time = time()
    # opfm = scopf(SC, system, optimizer, voll=voll, prob=prob, contingencies=contingencies, ramp_minutes=ramp_minutes, max_shed=max_shed, 
    #     renewable_prod = 1.0, debug = debug)
    opfm, pf, _, _ = SCOPF.opf(SC, system, optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed)
    # if !isnothing(branch_c) 
    #     @info "Flow constraint on branch $branch_c is $rate_c."
    #     JuMP.@constraint(opfm.mod, opfm.mod[:pf0][findfirst(x -> x == branch_c, opfm.branches)] == rate_c)
    # end
    MOI.set(opfm.mod, MOI.Silent(), true) # supress output from the solver
    solve_model!(opfm.mod)
    termination_status(opfm.mod) != MOI.OPTIMAL && return opfm, [], [], [], []
    total_solve_time = solve_time(opfm.mod)
    @debug "lower_bound = $(objective_value(opfm.mod))"

    # Set variables
    bd = benders(opfm)
    calc_θ!(pf, bd.Pᵢ)
    calc_Pline!(pf)
    overloads = zeros(length(opfm.contingencies))
    # overloads = Dict{Int, Vector{Tuple{Int64, Float64}}}()

    # pf = SCOPF.DCPowerFlow(opfm.nodes, opfm.branches, bd.Pᵢ, bd.idx)
    ΔPc = zeros(length(bd.Pg))
    ΔPcc = zeros(length(bd.Pg))
    ΔPccx = zeros(length(bd.Pg))
    pre = 0
    corr = 0
    # print_contingency_overflow(opfm, pf, Pcc, 1.0)

    contids = get_branch_bus_idx(opfm.branches, opfm.contingencies, bd.idx)
    bufs = [similar(pf.ϕ) for _ in 1:Threads.nthreads()]
    bufs_F = [similar(pf.F) for _ in 1:Threads.nthreads()]
    bufs_X = [similar(pf.X) for _ in 1:Threads.nthreads()]
    # num_islands = is_islanded.([pf], getindex.(contids, 2), getindex.(contids, 1))
    # if sum(num_islands) > 0.1 * length(opfm.contingencies) 
    mod_lock = Threads.ReentrantLock();
    # for contid in contids 
    Threads.@threads for contid in contids 
        ptdf = bufs[Threads.threadid()]
        flow = bufs_F[Threads.threadid()]
        X = bufs_X[Threads.threadid()]
        (c, cont) = contid
        # cont  = get_bus_idx(opfm.contingencies[c], bd.idx)
        islands = Vector{Vector{Int64}}() 
        island = 0
        if !is_islanded(pf, cont, c)
            calculate_line_flows!(flow, pf, cont, c)
        else
            # ptdf = similar(pf.ϕ)
            islands, island = find_system_state(ptdf, X, pf, cont, c)
            LinearAlgebra.mul!(flow, ptdf, bd.Pᵢ)
        end

        # Calculate the power flow with the new outage and find if there are any overloads
        overload = filter_overload(flow, bd.linerating * branch_short_term_limit_multi)

        if !isempty(overload)
            overloads[c] = maximum(x->abs(x[2]), overload)
        end
        
        if !isempty(islands)
            Threads.lock(mod_lock) do
                init_Pc(opfm, bd, islands, island, max_shed, cont, c)
                if type == PCSC::OPF
                    init_Pcc(opfm, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c)
                    if p_failure > 0.0
                        init_Pccx(opfm, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c)
                    end
                end
            end
            @info "Island: Contingency line $(cont[1])-$(cont[2])-i_$(c)"
        end
    end
    set_objective_function(opfm.mod, bd.obj)
    total_solve_time = update_model!(opfm, pf, bd, total_solve_time)
    termination_status(opfm.mod) != MOI.OPTIMAL && return opfm, pf, bd.Pc, bd.Pcc, bd.Pccx
    GC.safepoint()
    # end

    ptdf = copy(pf.ϕ)
    X = copy(pf.X)
    flow = copy(pf.F)
    permutation = sortperm(overloads, rev=true)
    cut_added = 1
    for iterations in 1:length(opfm.contingencies)
        if cut_added == 0 # loops until no new cuts are added for the contingencies
            @printf "END: Total solve time %.4f. Total time %.4f.\n" total_solve_time (time() - total_time)
            return opfm, pf, bd.Pc, bd.Pcc, bd.Pccx
        end
        cut_added = 0
        @info "Iteration $iterations"

        for contid in contids[permutation]
            (c, cont) = contid
            islands, island = find_system_state(ptdf, X, pf, cont, c)

            # if !isempty(islands) && get(bd.Pc, c, 0) == 0
            #     _ = get_Pc(opfm, bd, islands, island, max_shed, cont, c)
            #     if type == PCSC::OPF
            #         _ = get_Pcc(opfm, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c)
            #         if p_failure > 0.0
            #             _ = get_Pccx(opfm, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c)
            #         end
            #     end
            #     @info "Island: Contingency line $(cont[1])-$(cont[2])-i_$(c)"
            #     total_solve_time = update_model!(opfm, pf, bd, total_solve_time)
            # end

            # Calculate the power flow with the new outage and find if there are any overloads
            ΔPc = get_ΔPc(opfm.mod, bd.list, ΔPc, bd.Pc, c)
            LinearAlgebra.mul!(flow, ptdf, (bd.Pᵢ .+ ΔPc))
            overloads_c = filter_overload(flow, bd.linerating * branch_short_term_limit_multi)
            if type == PCSC::OPF
                ΔPcc = get_ΔPcc(opfm.mod, opfm.nodes, opfm.dc_branches, bd.list, ΔPcc, bd.Pcc, c)
                LinearAlgebra.mul!(flow, ptdf, (bd.Pᵢ .+ ΔPcc))
                overloads_cc = filter_overload(flow, bd.linerating * branch_long_term_limit_multi)
                if p_failure > 0.0
                    ΔPccx = get_ΔPccx(opfm.mod, bd.list, ΔPccx, bd.Pccx, c)
                    LinearAlgebra.mul!(flow, ptdf, (bd.Pᵢ .+ ΔPccx))
                    overloads_ccx = filter_overload(flow, bd.linerating * branch_long_term_limit_multi)
                end
            end

            # Cannot change the model before all data is exctracted!
            if !isempty(overloads_c)
                # overloads[c] = overload
                cut_added, pre = set_Pc(opfm, bd, ΔPc, ptdf, overloads_c, islands, island, max_shed, cont, c, cut_added, lim, pre)
            end
            if type == PCSC::OPF 
                if !isempty(overloads_cc) 
                    cut_added, corr = set_Pcc(opfm, bd, ΔPcc, ptdf, overloads_cc, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c, cut_added, lim, corr)
                end
                if p_failure > 0.0 && !isempty(overloads_ccx)
                    cut_added = set_Pccx(opfm, bd, ΔPccx, ptdf, overloads_ccx, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c, cut_added, lim, corr)
                end
            end
            if cut_added > 1
                total_solve_time = update_model!(opfm, pf, bd, total_solve_time)
                cut_added = 1
            end
            termination_status(opfm.mod) != MOI.OPTIMAL && return opfm, pf, bd.Pc, bd.Pcc, bd.Pccx
        end
        
    end
    @warn "Reached $(iterations) iterations without a stable solution."
    return opfm, pf, bd.Pc, bd.Pcc, bd.Pccx
end         

""" Solve model and update the power flow object """
function update_model!(opfm, pf, bd, total_solve_time)
    solve_model!(opfm.mod)
    total_solve_time += solve_time(opfm.mod)

    # Pᵢ = calc_Pᵢ(pf)
    # Pᵢ = get_net_Pᵢ(opfm, idx)
    bd.Pᵢ = get_value(opfm.mod, :p0)
    @. bd.Pg = bd.Pᵢ - bd.Pd
    calc_θ!(pf, bd.Pᵢ)
    calc_Pline!(pf)
    return total_solve_time
end

function find_overload(pf::DCPowerFlow, cont, c, linerating, branch_short_term_limit_multi, branch_long_term_limit_multi)
    lmt = min(branch_short_term_limit_multi, branch_long_term_limit_multi)
    flow = calculate_line_flows(pf, cont, c)
    return flow == zero(eltype(pf.X)) ? Vector(undef,1) : filter_overload(flow, linerating * lmt)
end

function find_system_state(ptdf::Matrix{T}, X::Matrix{T}, pf::DCPowerFlow, cont, c) where T<:Real
    islands = Vector{Vector{Int64}}() 
    island = 0
    if !is_islanded(pf, cont, c)
        get_isf!(ptdf, X, pf.X, pf.B, pf.DA, cont, c) # ptdf calculation is more computational expensive than line flow
    else
        islands = island_detection_thread_safe(pf.B, cont[1], cont[2])
        island = find_ref_island(islands, pf.slack)
        island_b = find_island_branches(islands[island], pf.DA, c)
        fill!(ptdf, zero(T))
        ptdf[island_b, islands[island]] = get_isf(pf.DA, pf.B, cont, c, pf.slack, islands[island], island_b)
    end
    
    return islands, island
end

""" Return the short term power injection change at each node. """
function get_ΔPc(mod::Model, list, ΔP, Pc, c)
    fill!(ΔP, zero(eltype(ΔP)))
    x = get(Pc, c, 0)
    if x != 0
        pgc = get_value(mod, x[1])
        prc = get_value(mod, x[2])
        lsc = get_value(mod, x[3])
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

function init_Pc(opfm::OPFmodel, bd, islands, island, max_shed, cont, c)
    pgc = JuMP.@variable(opfm.mod, [g in 1:length(opfm.ctrl_generation)], base_name = @sprintf("pgc%s", c), lower_bound = 0)
    prc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.renewables)], base_name = @sprintf( "prc%s", c), lower_bound = 0)
    lsc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.demands)], base_name = @sprintf( "lsc%s", c), lower_bound = 0)
    bd.Pc[c] = (pgc, prc, lsc)
    # obj = objective_function(opfm.mod)
    add_to_expression!.(bd.obj, opfm.prob[c], sum(opfm.voll' * lsc))
    # set_objective_function(opfm.mod, bd.obj)

    # Add new constraints that limit the corrective variables within operating limits
    JuMP.@constraint(opfm.mod, sum(lsc) == sum(prc) + sum(pgc))
    if isempty(islands)
        JuMP.@constraints(opfm.mod, begin
            0 .<= opfm.mod[:pg0] .- pgc
            prc .<= bd.pr_lim .* max_shed
            lsc .<= bd.pd_lim .* max_shed
        end)
    else
        for n in islands[island]
            JuMP.@constraints(opfm.mod, begin 
                [g = bd.list[n].ctrl_generation],  0 <= opfm.mod[:pg0][g] - pgc[g]
                [g = bd.list[n].renewables], prc[g] <= bd.pr_lim[g] * max_shed
                [d = bd.list[n].demands], lsc[d] <= bd.pd_lim[d] * max_shed
            end)
        end
        for in_vec in islands[1:end .!= island]
            for n in in_vec
                JuMP.@constraints(opfm.mod, begin 
                    [g = bd.list[n].ctrl_generation], opfm.mod[:pg0][g] == pgc[g]
                    [g = bd.list[n].renewables], prc[g] == bd.pr_lim[g]
                    [d = bd.list[n].demands], lsc[d] == bd.pd_lim[d]
                end)
            end
        end
    end
    return pgc, prc, lsc
end

function get_Pc(opfm::OPFmodel, bd, islands, island, max_shed, cont, c)
    if get(bd.Pc, c, 0) == 0  # If the contingency is not run before, a set of corrective variables is added
        init_Pc(opfm, bd, islands, island, max_shed, cont, c) 
        set_objective_function(opfm.mod, bd.obj)
    end
    return get(bd.Pc, c, 0)
end

function set_Pc(opfm::OPFmodel, bd, ΔPc, ptdf, overloads, islands, island, max_shed, cont, c, cut_added, lim, id)
    if !isempty(overloads)
        (pgc, prc, lsc) = get_Pc(opfm, bd, islands, island, max_shed, cont, c)
        for (i,ol) in overloads
            expr = JuMP.@expression(opfm.mod, sum(
                    (ptdf[i, inode] * (bd.Pg[inode] + ΔPc[inode] - 
                    sum((opfm.mod[:pg0][ctrl] - pgc[ctrl] for ctrl in sublist.ctrl_generation), init=0.0) -
                    sum((beta(sublist.node, opfm.dc_branches[d]) * opfm.mod[:pfdc0][d] for d in sublist.dc_branches), init=0.0) +
                    sum((prc[r] for r in sublist.renewables), init=0.0) -
                    sum((lsc[d] for d in sublist.demands), init=0.0)
                    ) for (inode, sublist) in enumerate(bd.list)), init=0.0
                ) )
            
            id += 1
            @info @sprintf "Pre %d: Contingency line %d-%d-i_%d; overload on %s of %.4f" id cont[1] cont[2] c opfm.branches[i].name ol
            @debug "Cut added: $(sprint_expr(expr,lim))\n"
            # set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
            if ol < 0
                pre_cut = JuMP.@constraint(opfm.mod, expr <= ol)
            else
                pre_cut = JuMP.@constraint(opfm.mod, expr >= ol)
            end
            cut_added = 2
        end
    end
    return cut_added, id
end

""" Return the long term power injection change at each node. """
function get_ΔPcc(mod::Model, nodes::AbstractVector{Bus}, dc_branches::AbstractVector{DCBranch}, list, ΔP, Pcc, c)
    fill!(ΔP, zero(eltype(ΔP)))
    x = get(Pcc, c, 0)
    if x != 0
        pgu = get_value(mod, x[1])
        pgd = get_value(mod, x[2])
        pfdccc = get_value(mod, x[3])
        prcc = get_value(mod, x[4])
        lscc = get_value(mod, x[5])
        for (i, n) in enumerate(list)
            for g in n.ctrl_generation
                ΔP[i] += (pgu[g] - pgd[g])
            end
            for g in n.dc_branches
                ΔP[i] += beta(nodes[i], dc_branches[g]) * pfdccc[g]
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

function init_Pcc(opfm::OPFmodel, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c)
    # Add corrective variables
    pgu = JuMP.@variable(opfm.mod, [g in 1:length(opfm.ctrl_generation)], base_name = @sprintf( "pgu%s", c), lower_bound = 0)
        # active power variables for the generators in contingencies ramp up 
    pgd = JuMP.@variable(opfm.mod, [g in 1:length(opfm.ctrl_generation)], base_name = @sprintf( "pgd%s", c), lower_bound = 0)
            # and ramp down
    pfdccc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.dc_branches)], base_name = @sprintf( "pfdccc%s", c))
    prcc = JuMP.@variable(opfm.mod, [r in 1:length(opfm.renewables)], base_name = @sprintf( "prcc%s", c), lower_bound = 0)
    lscc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.demands)], base_name = @sprintf( "lscc%s", c), lower_bound = 0)
            # load curtailment variables in in contingencies
    bd.Pcc[c] = (pgu, pgd, pfdccc, prcc, lscc)

    # Extend the objective with the corrective variables
    # obj = objective_function(opfm.mod)
    add_to_expression!.(bd.obj, opfm.prob[c], 
        # (1.0 - p_failure) * (sum(opfm.voll' * lscc) + sum(60 * (pgu + pgd)) # + # TODO: remove 60 and uncomment next lines for non-4-area analysis!!!!
        (1.0 - p_failure) * (sum(opfm.voll' * lscc) # + sum(opfm.cost_ctrl_gen' * ramp_mult * (pgu + pgd))
        # (sum(opfm.voll' * lscc) +
        # sum(opfm.cost_ctrl_gen' * ramp_mult * pgu) # +
        # sum(opfm.cost_ctrl_gen' * pgd)
        ) )
    add_to_expression!.(bd.obj, opfm.prob[c], (1.0 - p_failure) * (sum(opfm.cost_ctrl_gen' * ramp_mult * pgu)) )
    add_to_expression!.(bd.obj, opfm.prob[c], (1.0 - p_failure) * (sum(opfm.cost_ctrl_gen' * ramp_mult * pgd)) )
    # set_objective_function(opfm.mod, bd.obj)

    # Add new constraints that limit the corrective variables within operating limits
    JuMP.@constraint(opfm.mod, sum(pgu) + sum(lscc) == sum(pgd) + sum(prcc))
    if isempty(islands)
        JuMP.@constraints(opfm.mod, begin
            pgu .<= bd.rampup * ramp_minutes
            pgd .<= bd.rampdown * ramp_minutes
            0 .<= opfm.mod[:pg0] .+ pgu .- pgd
            opfm.mod[:pg0] .+ pgu .- pgd .<= bd.pg_lim_max
            bd.dc_lim_min .<= pfdccc
            pfdccc .<= bd.dc_lim_max
            prcc .<= bd.pr_lim .* max_shed
            lscc .<= bd.pd_lim .* max_shed
        end)
    else
        for n in islands[island]
            JuMP.@constraints(opfm.mod, begin
                [g = bd.list[n].ctrl_generation], pgu[g] <= bd.rampup[g] * ramp_minutes
                [g = bd.list[n].ctrl_generation], pgd[g] <= bd.rampdown[g] * ramp_minutes
                [g = bd.list[n].ctrl_generation], 0 <= opfm.mod[:pg0][g] + pgu[g] - pgd[g]
                [g = bd.list[n].ctrl_generation], opfm.mod[:pg0][g] + pgu[g] - pgd[g] <= bd.pg_lim_max[g]
                [g = bd.list[n].dc_branches], bd.dc_lim_min[g] <= pfdccc[g]
                [g = bd.list[n].dc_branches], pfdccc[g] <= bd.dc_lim_max[g]
                [g = bd.list[n].renewables], prcc[g] <= bd.pr_lim[g] * max_shed
                [d = bd.list[n].demands], lscc[d] <= bd.pd_lim[d] * max_shed
            end)
        end
        for in_vec in islands[1:end .!= island]
            for n in in_vec
                JuMP.@constraints(opfm.mod, begin
                    [g = bd.list[n].ctrl_generation], opfm.mod[:pg0][g] == pgd[g] - pgu[g]
                    [g = bd.list[n].dc_branches], pfdccc[g] == 0
                    [g = bd.list[n].renewables], prcc[g] == bd.pr_lim[g]
                    [d = bd.list[n].demands], lscc[d] == bd.pd_lim[d]
                end)
            end
        end
    end
end

function get_Pcc(opfm::OPFmodel, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c)
    if get(bd.Pcc, c, 0) == 0 # If the contingency is not run before, a set of corrective variables is added
        init_Pcc(opfm, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c) 
        set_objective_function(opfm.mod, bd.obj)
    end
    return get(bd.Pcc, c, 0)
end

function set_Pcc(opfm::OPFmodel, bd, ΔPcc, ptdf, overloads, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c, cut_added, lim, id)
    if !isempty(overloads)
        (pgu, pgd, pfdccc, prcc, lscc) = get_Pcc(opfm, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c)

        # sort!(overloads, rev = true, by = x -> abs(x[2]))
        # (i,ol) = first(overloads)
        for (i,ol) in overloads
            # Finding and adding the Benders cut
            expr = JuMP.@expression(opfm.mod, sum((ptdf[i, inode] * (
                    bd.Pg[inode] + ΔPcc[inode] -
                    sum((opfm.mod[:pg0][ctrl] + pgu[ctrl] - pgd[ctrl] 
                        for ctrl in sublist.ctrl_generation), init=0.0) -
                    sum((beta(sublist.node, opfm.dc_branches[d]) * pfdccc[d] for d in sublist.dc_branches), init=0.0) +
                    sum((prcc[r] for r in sublist.renewables), init=0.0) -
                    sum((lscc[d] for d in sublist.demands), init=0.0)
                    ) for (inode, sublist) in enumerate(bd.list)), init=0.0
                ) )
            
            id += 1
            @info @sprintf "Corr %d: Contingency line %d-%d-i_%d; overload on %s of %.4f" id cont[1] cont[2] c opfm.branches[i].name ol
            @debug "Cut added: $(sprint_expr(expr,lim))\n"
            # set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
            if ol < 0
                corr_cut = JuMP.@constraint(opfm.mod, expr <= ol)
            else
                corr_cut = JuMP.@constraint(opfm.mod, expr >= ol)
            end
            cut_added = 3
        end
    end
    return cut_added, id
end

""" Return the long term power injection change at each node. """
function get_ΔPccx(mod::Model, list, ΔP, Pccx, c)
    fill!(ΔP, zero(eltype(ΔP)))
    x = get(Pccx, c, 0)
    if x != 0
        pgd = get_value(mod, x[1])
        lscc = get_value(mod, x[2])
        for (i, n) in enumerate(list)
            for g in n.ctrl_generation
                ΔP[i] -= pgd[g]
            end
            for g in n.demands
                ΔP[i] += lscc[g]
            end
        end
    end
    return ΔP
end

function init_Pccx(opfm::OPFmodel, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c)
    # Add corrective variables
    pgd = JuMP.@variable(opfm.mod, [g in 1:length(opfm.ctrl_generation)], base_name = @sprintf( "pgdx%s", c), lower_bound = 0)
            # and ramp down
    lscc = JuMP.@variable(opfm.mod, [d in 1:length(opfm.demands)], base_name = @sprintf( "lsccx%s", c), lower_bound = 0)
            # load curtailment variables in in contingencies
    bd.Pccx[c] = (pgd, lscc)

    # Extend the objective with the corrective variables
    # obj = objective_function(opfm.mod)
    add_to_expression!.(bd.obj, opfm.prob[c] * p_failure, sum(opfm.voll' * lscc))
    add_to_expression!.(bd.obj, opfm.prob[c] * p_failure, sum(60 * pgd))
    # set_objective_function(opfm.mod, bd.obj)

    # Add new constraints that limit the corrective variables within operating limits
    JuMP.@constraint(opfm.mod, sum(lscc) == sum(pgd))
    if isempty(islands)
        JuMP.@constraints(opfm.mod, begin
            pgd .<= bd.rampdown * ramp_minutes
            0 .<= opfm.mod[:pg0] .- pgd
            lscc .<= bd.pd_lim .* max_shed
        end)
    else
        for n in islands[island]
            JuMP.@constraints(opfm.mod, begin
                [g = bd.list[n].ctrl_generation], pgd[g] <= bd.rampdown[g] * ramp_minutes
                [g = bd.list[n].ctrl_generation], 0 <= opfm.mod[:pg0][g] - pgd[g]
                [d = bd.list[n].demands], lscc[d] <= bd.pd_lim[d] * max_shed
            end)
        end
        for in_vec in islands[1:end .!= island]
            for n in in_vec
                JuMP.@constraint(opfm.mod, [g = bd.list[n].ctrl_generation], 
                    opfm.mod[:pg0][g] == pgd[g])
                JuMP.@constraint(opfm.mod, [d = bd.list[n].demands], 
                    lscc[d] == bd.pd_lim[d])
            end
        end
    end
end

function get_Pccx(opfm::OPFmodel, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c)
    if get(bd.Pccx, c, 0) == 0 # If the contingency is not run before, a set of corrective variables is added
        init_Pccx(opfm, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c) 
        set_objective_function(opfm.mod, bd.obj)
    end
    return get(bd.Pccx, c, 0)
end

function set_Pccx(opfm::OPFmodel, bd, ΔPccx, ptdf, overloads, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c, cut_added, lim, id)
    if !isempty(overloads)
        (pgd, lscc) = get_Pccx(opfm, bd, islands, island, max_shed, ramp_mult, ramp_minutes, p_failure, cont, c)
        # sort!(overloads, rev = true, by = x -> abs(x[2]))
        # (i,ol) = first(overloads)
        for (i,ol) in overloads
            # Finding and adding the Benders cut
            expr = JuMP.@expression(opfm.mod, sum((ptdf[i, inode] * (
                    bd.Pg[inode] + ΔPccx[inode] -
                    sum((opfm.mod[:pg0][ctrl] - pgd[ctrl] 
                        for ctrl in sublist.ctrl_generation), init=0.0)  -
                    sum((lscc[d] for d in sublist.demands), init=0.0)
                    ) for (inode, sublist) in enumerate(bd.list)), init=0.0
                ) )
            
            @info @sprintf "Corr.x %d: Contingency line %d-%d-i_%d; overload on %s of %.4f" id cont[1] cont[2] c opfm.branches[i].name ol
            @debug "Cut added: $(sprint_expr(expr,lim))\n"
            # set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
            if ol < 0
                corr_cut = JuMP.@constraint(opfm.mod, expr <= ol)
            else
                corr_cut = JuMP.@constraint(opfm.mod, expr >= ol)
            end
            cut_added = 3
        end
    end
    return cut_added
end

" An AffExpr nicely formatted to a string "
sprint_expr(expr::AffExpr, lim = 1e-6) = 
    join(Printf.@sprintf("%s%5.2f %s ", (x[2] > 0 ? "+" : "-"), abs(x[2]), x[1]) 
            for x in expr.terms if abs(x[2]) > lim) * 
        Printf.@sprintf("<= %s%.2f", (expr.constant > 0 ? "-" : " "), abs(expr.constant)
    )

function print_benders_results(opfm::OPFmodel, Pc::Dict{Int, NTuple{3, Any}} = Dict{Int, NTuple{3, Any}}(), 
        Pcc::Dict{Int, NTuple{5, Any}} = Dict{Int, NTuple{5, Any}}(), 
        Pccx::Dict{Int, NTuple{2, Any}} = Dict{Int, NTuple{2, Any}}(), lim::Real = 1e-6)
    function print_c(itr, val, symb, i_g, lim)
        for i in 1:length(opfm.contingencies)
            c = get(itr, i, 0)
            if c != 0 && JuMP.value(c[val][i_g]) > lim
                @printf("          c %12s: %s: %.3f\n", opfm.contingencies[i].name, symb, JuMP.value(c[val][i_g]))
            end
        end
    end
    for (i_g,g) in enumerate(opfm.ctrl_generation)
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(opfm.mod[:pg0][i_g]), get_active_power_limits(g).max)
        print_c(Pc, 1, "pgc", i_g, lim)
        print_c(Pcc, 1, "pgu", i_g, lim)
        print_c(Pcc, 2, "pgd", i_g, lim)
        print_c(Pccx, 1, "pgdx", i_g, lim)
    end
    for (i_g,g) in enumerate(opfm.dc_branches)
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(opfm.mod[:pfdc0][i_g]), get_active_power_limits(g).max)
        print_c(Pcc, 3, "pfdccc", i_g, lim)
    end
    for (i_g,g) in enumerate(opfm.renewables)
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(opfm.mod[:pr0][i_g]), get_active_power_limits(g).max)
        print_c(Pc, 2, "prc", i_g, lim)
        print_c(Pcc, 4, "prcc", i_g, lim)
    end
    for (i_g,g) in enumerate(opfm.demands)
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(opfm.mod[:ls0][i_g]), get_active_power(g))
        print_c(Pc, 3, "lsc", i_g, lim)
        print_c(Pcc, 5, "lscc", i_g, lim)
        print_c(Pccx, 2, "lsccx", i_g, lim)
    end
end
