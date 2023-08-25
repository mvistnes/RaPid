# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022 

""" 
Solve the optimization model using contingency selection.
"""
function run_contingency_select(type::OPF, system::System, optimizer, voll=nothing, prob=nothing, contingencies=nothing;
    ramp_minutes=10, ramp_mult=10, max_shed=1.0, lim=1e-6,
    branch_short_term_limit_multi::Real=1.5, branch_long_term_limit_multi::Real=1.2,
    gen_short_term_limit_multi::Real=2.0, gen_long_term_limit_multi::Real=1.0, p_failure=0.0, branch_c=nothing, rate_c=0.0, debug=false
)
    @assert branch_short_term_limit_multi >= branch_long_term_limit_multi
    total_time = time()
    # LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
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
    pre = 0
    corr = 0

    contids = get_branch_bus_idx(opfm.branches, opfm.contingencies, bd.idx)
    flow = similar(pf.F)
    θ = similar(pf.θ)
    ptdf = copy(pf.ϕ)
    # bufs_F = [similar(pf.F) for _ in 1:Threads.nthreads()]
    # bufs_θ = [similar(pf.θ) for _ in 1:Threads.nthreads()]
    # num_islands = is_islanded.([pf], getindex.(contids, 2), getindex.(contids, 1))
    # if sum(num_islands) > 0.1 * length(opfm.contingencies) 
    # LinearAlgebra.BLAS.set_num_threads(1)
    # mod_lock = Threads.ReentrantLock()
    for contid in contids 
    # Threads.@threads for contid in contids
        # flow = bufs_F[Threads.threadid()]
        # θ = bufs_θ[Threads.threadid()]
        (c, cont) = contid
        # cont  = get_bus_idx(opfm.contingencies[c], bd.idx)
        if !is_islanded(pf, cont, c)
            calculate_line_flows!(flow, pf, cont, c)
        else
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont, c, pf.slack)
            get_isf!(ptdf, pf.DA, pf.B, cont, c, pf.slack, islands[island], island_b)
            add_short_term_contingencies(opfm, bd.obj, bd.Pc, islands, island, ptdf, bd.list, bd.pr_lim, bd.pd_lim, max_shed, bd.linerating * branch_short_term_limit_multi, cont, c)
            add_long_term_contingencies(opfm, bd.obj, bd.Pcc, islands, island, ptdf, bd.list, bd.pr_lim, bd.pd_lim, max_shed, bd.linerating * branch_long_term_limit_multi, ramp_mult, ramp_minutes, bd.rampup, bd.rampdown, bd.dc_lim_min, bd.dc_lim_max, bd.pg_lim_max, cont, c)
            pre += 1
            corr += 1
            @debug "Island: Contingency line $(cont[1])-$(cont[2])-i_$(c)"
            calculate_line_flows!(flow, θ, pf.DA, pf.B, bd.Pᵢ, cont, c, pf.slack, islands[island], island_b)
        end

        # Calculate the power flow with the new outage and find if there are any overloads
        overload = filter_overload(flow, bd.linerating * branch_short_term_limit_multi)

        if !isempty(overload)
            overloads[c] = maximum(x -> abs(x[2]), overload)
        end
    end
    set_objective_function(opfm.mod, bd.obj)
    total_solve_time = update_model!(opfm, pf, bd, total_solve_time)
    termination_status(opfm.mod) != MOI.OPTIMAL && return opfm, pf, bd.Pc, bd.Pcc, bd.Pccx
    GC.safepoint()
    # end
    # LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())

    ΔPc = zeros(length(bd.Pg))
    ΔPcc = zeros(length(bd.Pg))
    ΔPccx = zeros(length(bd.Pg))

    X = copy(pf.X)
    B = copy(pf.B)
    θ = copy(pf.θ)
    flow = copy(pf.F)
    overloads_c = Vector{Tuple{Int,Float64}}[]
    overloads_cc = Vector{Tuple{Int,Float64}}[]
    overloads_ccx = Vector{Tuple{Int,Float64}}[]
    islands = Vector{Vector{Int}}[]
    island = 0
    island_b = Int[]

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
            empty!(overloads_c)
            empty!(overloads_cc)
            empty!(overloads_ccx)
            if is_islanded(pf, cont, c)
                islands, island, island_b = handle_islands(pf.B, pf.DA, cont, c, pf.slack) 
            else
                empty!(islands)
                island = 0
                empty!(island_b)
            end
            
            ΔPc, ΔPcc, ΔPccx, overload = calculate_contingency_line_flows!(ΔPc, ΔPcc, ΔPccx, flow, θ, B, Val(type), p_failure, opfm, pf, bd, cont, c, islands, island, island_b, branch_short_term_limit_multi, branch_long_term_limit_multi)
            if overload
                if get(bd.Pc, c, 0) == 0
                    ΔPc, ΔPcc, ΔPccx = get_ΔP!(ΔPc, ΔPcc, ΔPccx, Val(type), p_failure, opfm, bd, c)
                end
                # Calculate the power flow with the new outage and find if there are any overloads
                if !is_islanded(pf, cont, c)
                    get_isf!(ptdf, X, pf.X, pf.B, pf.DA, cont, c) # ptdf calculation is more computational expensive than line flow
                else
                    get_isf!(ptdf, pf.DA, pf.B, cont, c, pf.slack, islands[island], island_b)
                end
                overloads_c, overloads_cc, overloads_ccx = find_overloads(Val(type), p_failure, flow, ptdf, bd, ΔPc, ΔPcc, ΔPccx, branch_short_term_limit_multi, branch_long_term_limit_multi)
            end
            # (overload || !isempty(overloads_c) || !isempty(overloads_cc)) && println(c, " ", Int(overload), Int(!isempty(overloads_c)), Int(!isempty(overloads_cc)))

            # Cannot change the model before all data is exctracted!
            if !isempty(overloads_c)
                add_short_term_contingencies(opfm, bd.obj, bd.Pc, islands, island, ptdf, bd.list, bd.pr_lim, bd.pd_lim, max_shed, bd.linerating * branch_short_term_limit_multi, cont, c)
                pre += 1
                @info @sprintf "Pre %d: Contingency line %d-%d-i_%d" pre cont[1] cont[2] c 
                cut_added = 2
            end
            if !isempty(overloads_cc)
                add_long_term_contingencies(opfm, bd.obj, bd.Pcc, islands, island, ptdf, bd.list, bd.pr_lim, bd.pd_lim, max_shed, bd.linerating * branch_long_term_limit_multi, ramp_mult, ramp_minutes, bd.rampup, bd.rampdown, bd.dc_lim_min, bd.dc_lim_max, bd.pg_lim_max, cont, c)
                corr += 1
                @info @sprintf "Corr %d: Contingency line %d-%d-i_%d" corr cont[1] cont[2] c
                cut_added = 2
            end
            if cut_added > 1
                set_objective_function(opfm.mod, bd.obj)
                total_solve_time = update_model!(opfm, pf, bd, total_solve_time)
                cut_added = 1
            end
            termination_status(opfm.mod) != MOI.OPTIMAL && return opfm, pf, bd.Pc, bd.Pcc, bd.Pccx
        end

    end
    @warn "Reached $(length(opfm.contingencies)) iterations without a stable solution."
    return opfm, pf, bd.Pc, bd.Pcc, bd.Pccx
end
