# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022 

""" 
Solve the optimization model using contingency selection.
"""
function run_contingency_select(
    type::OPF, 
    system::System, 
    optimizer, 
    voll=nothing, 
    prob=nothing, 
    contingencies=nothing;
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
    branch_c=nothing, 
    rate_c=0.0, 
    silent=true,
    debug=false
)
    mod, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SC, system, optimizer, voll=voll, contingencies=contingencies, prob=prob,
        time_limit_sec=time_limit_sec, ramp_minutes=ramp_minutes, ramp_mult=ramp_mult, max_shed=max_shed, max_curtail=max_curtail,
        short_term_multi=short_term_multi, long_term_multi=long_term_multi, p_failure=p_failure, silent=silent, debug=debug)
    return run_contingency_select(SC, type, mod, opf, pf, oplim, Pc, Pcc, Pccx, lim, max_itr, branch_c, rate_c, debug)
end
        
function run_contingency_select(
    basetype::OPF, 
    type::OPF, 
    mod::Model, 
    opf::OPFsystem, 
    pf::DCPowerFlow, 
    oplim::Oplimits, 
    Pc::Dict{<:Integer, Main.SCOPF.ExprC}, 
    Pcc::Dict{<:Integer, Main.SCOPF.ExprCC}, 
    Pccx::Dict{<:Integer, Main.SCOPF.ExprCCX},
    lim=1e-6,
    max_itr=length(opf.contingencies),
    branch_c=nothing, 
    rate_c=0.0, 
    debug=false
)
    @assert isempty(Pc) || basetype == PSC::OPF
    @assert isempty(Pcc) && isempty(Pccx)
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
            basetype == SC::OPF && add_contingencies!(Pc, opf, oplim, mod, bd.obj, islands, island, ptdf, bd.list, c)
            type == PCSC::OPF && add_contingencies!(Pcc, opf, oplim, mod, bd.obj, islands, island, ptdf, bd.list, c)
            type == PCFSC::OPF && add_contingencies!(Pccx, opf, oplim, mod, bd.obj, islands, island, ptdf, bd.list, c)
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
    olc = Vector{Tuple{Int,Float64}}[]
    olcc = Vector{Tuple{Int,Float64}}[]
    olccx = Vector{Tuple{Int,Float64}}[]
    islands = Vector{Vector{Int}}[]
    island = 0
    island_b = Int[]

    permutation = sortperm(overloads, rev=true)
    contids = contids[permutation]
    ix = [i for (i, (c, cont)) in enumerate(contids) if !is_islanded(pf, cont, c)]
    contids = contids[ix]

    i = 1
    iterations = 1
    cut_added = 0
    while !isempty(contids)
        (c, cont) = contids[i]
        empty!(olc)
        empty!(olcc)
        empty!(olccx)
        if is_islanded(pf, cont, c)
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont, c, pf.slack) 
        else
            empty!(islands)
            island = 0
            empty!(island_b)
        end
        
        olc = basetype == SC::OPF ? calculate_contingency_overload!(ΔPc, flow, (oplim.branch_rating * oplim.short_term_multi), θ, B, Pc, opf, mod, pf, bd, cont, c, islands, island, island_b) : Tuple{Int64,Float64}[]
        olcc = type >= PCSC::OPF ? calculate_contingency_overload!(ΔPcc, flow, (oplim.branch_rating * oplim.long_term_multi), θ, B, Pccx, opf, mod, pf, bd, cont, c, islands, island, island_b) : Tuple{Int64,Float64}[]
        olccx = type >= PCFSC::OPF ? calculate_contingency_overload!(ΔPccx, flow, (oplim.branch_rating * oplim.long_term_multi), θ, B, Pccx, opf, mod, pf, bd, cont, c, islands, island, island_b) : Tuple{Int64,Float64}[]
        if !isempty(olc) || !isempty(olcc) || !isempty(olccx)
            # Calculate the power flow with the new outage and find if there are any overloads
            if !is_islanded(pf, cont, c)
                get_isf!(ptdf, X, pf.X, pf.B, pf.DA, cont, c) # ptdf calculation is more computational expensive than line flow
            else
                get_isf!(ptdf, pf.DA, pf.B, cont, c, pf.slack, islands[island], island_b)
            end
        end
        # (overload || !isempty(olc) || !isempty(olcc)) && println(c, " ", Int(overload), Int(!isempty(olc)), Int(!isempty(olcc)))

        # Cannot change the model before all data is exctracted!
        if !isempty(olc)
            add_contingencies!(Pc, opf, oplim, mod, bd.obj, islands, island, ptdf, bd.list, c)
            pre += 1
            @info @sprintf "Pre %d: Contingency line %d-%d-i_%d" pre cont[1] cont[2] c 
            cut_added = 2
        end
        if !isempty(olcc)
            add_contingencies!(Pcc, opf, oplim, mod, bd.obj, islands, island, ptdf, bd.list, c)
            corr += 1
            @info @sprintf "Corr %d: Contingency line %d-%d-i_%d" corr cont[1] cont[2] c
            cut_added = 2
        end
        if cut_added > 1
            set_objective_function(mod, bd.obj)
            total_solve_time = update_model!(mod, pf, bd, total_solve_time)
            deleteat!(contids, findfirst(x -> x[1] == c, contids))
            cut_added = 1
        else
            i += 1
        end
        termination_status(mod) != MOI.OPTIMAL && return mod, opf, pf, oplim, Pc, Pcc, Pccx
        if i > length(contids)
            iterations > max_itr && break
            if cut_added == 0 # loops until no new cuts are added for the contingencies
                @printf "END: Total solve time %.4f.\n" total_solve_time
                return mod, opf, pf, oplim, Pc, Pcc, Pccx
            else
                cut_added = 0
            end
            iterations += 1
            @info "Iteration $iterations"
            i = 1
        end
    end
    # @warn "Reached $(length(opf.contingencies)) iterations without a stable solution."
    return mod, opf, pf, oplim, Pc, Pcc, Pccx
end
