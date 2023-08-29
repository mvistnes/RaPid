# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022 

""" 
Solve the optimization model using contingency selection.
"""
function run_contingency_select(type::OPF, system::System, optimizer, voll=nothing, prob=nothing, contingencies=nothing;
    time_limit_sec::Int64=600, ramp_minutes=10, ramp_mult=10, max_shed=1.0, max_curtail=1.0, lim=1e-6,
    short_term_multi::Real=1.5, long_term_multi::Real=1.2,
    p_failure=0.0, branch_c=nothing, rate_c=0.0, debug=false
)
    @assert short_term_multi >= long_term_multi
    total_time = time()
    mod, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SC, system, optimizer, voll=voll, contingencies=contingencies, prob=prob,
        time_limit_sec=time_limit_sec, ramp_minutes=ramp_minutes, ramp_mult=ramp_mult, max_shed=max_shed, max_curtail=max_curtail,
        short_term_multi=short_term_multi, long_term_multi=long_term_multi, p_failure=p_failure, debug=debug)
    MOI.set(mod, MOI.Silent(), true) # supress output from the solver
    solve_model!(mod)
    termination_status(mod) != MOI.OPTIMAL && return mod, opf, pf, oplim, Pc, Pcc, Pccx
    total_solve_time = solve_time(mod)
    @debug "lower_bound = $(objective_value(mod))"

    # Set variables
    bd = benders(opf, mod)
    calc_θ!(pf, bd.Pᵢ)
    calc_Pline!(pf)
    overloads = zeros(length(opf.contingencies))
    pre = 0
    corr = 0

    contids = get_branch_bus_idx(opf.branches, opf.contingencies, bd.idx)
    flow = similar(pf.F)
    θ = similar(pf.θ)
    ptdf = copy(pf.ϕ)
    for contid in contids 
        (c, cont) = contid
        if !is_islanded(pf, cont, c)
            calculate_line_flows!(flow, pf, cont, c)
            # Calculate the power flow with the new outage and find if there are any overloads
            overload = filter_overload(flow, oplim.branch_rating * oplim.short_term_multi)
            if !isempty(overload)
                overloads[c] = maximum(x -> abs(x[2]), overload)
            end
        else
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont, c, pf.slack)
            get_isf!(ptdf, pf.DA, pf.B, cont, c, pf.slack, islands[island], island_b)
            add_short_term_contingencies!(Pc, opf, oplim, mod, bd.obj, islands, island, ptdf, bd.list, c)
            type == PCSC::OPF && add_long_term_contingencies!(Pcc, opf, oplim, mod, bd.obj, islands, island, ptdf, bd.list, c)
            pre += 1
            corr += 1
            @debug "Island: Contingency line $(cont[1])-$(cont[2])-i_$(c)"
        end
    end
    set_objective_function(mod, bd.obj)
    total_solve_time = update_model!(mod, pf, bd, total_solve_time)
    termination_status(mod) != MOI.OPTIMAL && return mod, opf, pf, oplim, Pc, Pcc, Pccx
    GC.safepoint()

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
        
        ΔPc, ΔPcc, ΔPccx, overload = calculate_contingency_line_flows!(ΔPc, ΔPcc, ΔPccx, flow, θ, B, Val(type), Pc, Pcc, Pccx, opf, oplim, mod, pf, bd, cont, c, islands, island, island_b)
        if overload
            if get(Pc, c, 0) == 0
                ΔPc, ΔPcc, ΔPccx = get_ΔP!(ΔPc, ΔPcc, ΔPccx, Pc, Pcc, Pccx, Val(type), opf, mod, bd.list, oplim.p_failure, c)
            end
            # Calculate the power flow with the new outage and find if there are any overloads
            if !is_islanded(pf, cont, c)
                get_isf!(ptdf, X, pf.X, pf.B, pf.DA, cont, c) # ptdf calculation is more computational expensive than line flow
            else
                get_isf!(ptdf, pf.DA, pf.B, cont, c, pf.slack, islands[island], island_b)
            end
            overloads_c, overloads_cc, overloads_ccx = find_overloads(Val(type), oplim.p_failure, flow, ptdf, bd.Pᵢ, oplim, ΔPc, ΔPcc, ΔPccx)
        end
        # (overload || !isempty(overloads_c) || !isempty(overloads_cc)) && println(c, " ", Int(overload), Int(!isempty(overloads_c)), Int(!isempty(overloads_cc)))

        # Cannot change the model before all data is exctracted!
        if !isempty(overloads_c)
            add_short_term_contingencies!(Pc, opf, oplim, mod, bd.obj, islands, island, ptdf, bd.list, c)
            pre += 1
            @info @sprintf "Pre %d: Contingency line %d-%d-i_%d" pre cont[1] cont[2] c 
            cut_added = 2
        end
        if type == PCSC::OPF && !isempty(overloads_cc)
            add_long_term_contingencies!(Pcc, opf, oplim, mod, bd.obj, islands, island, ptdf, bd.list, c)
            corr += 1
            @info @sprintf "Corr %d: Contingency line %d-%d-i_%d" corr cont[1] cont[2] c
            cut_added = 2
        end
        if cut_added > 1
            set_objective_function(mod, bd.obj)
            total_solve_time = update_model!(mod, pf, bd, total_solve_time)
            cut_added = 1
        end
        termination_status(mod) != MOI.OPTIMAL && return mod, opf, pf, oplim, Pc, Pcc, Pccx
    end
    # @warn "Reached $(length(opf.contingencies)) iterations without a stable solution."
    return mod, opf, pf, oplim, Pc, Pcc, Pccx
end
