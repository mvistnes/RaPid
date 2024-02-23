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
    dist_slack=Float64[],
    time_limit_sec::Int64=600, 
    ramp_minutes=10, 
    ramp_mult=10, 
    max_shed=1.0, 
    max_curtail=1.0, 
    lim=1e-14,
    short_term_multi::Real=1.5, 
    long_term_multi::Real=1.2, 
    max_itr=max(length(contingencies), 5),
    p_failure=0.0, 
    silent=true,
    debug=false
)
    case = Case(opf_base(OPF(true, false, false, false, false), system, optimizer, voll=voll, contingencies=contingencies, prob=prob,
    dist_slack=dist_slack, time_limit_sec=time_limit_sec, ramp_minutes=ramp_minutes, ramp_mult=ramp_mult, max_shed=max_shed, max_curtail=max_curtail,
    short_term_multi=short_term_multi, long_term_multi=long_term_multi, p_failure=p_failure, silent=silent, debug=debug)...)

    total_solve_time = constrain_branches!(case, 0.0)
    if !type.P & !type.C1 & !type.C2 & !type.C2F
        return case, total_solve_time
    end
    return run_contingency_select!(type, case, lim, max_itr)
end

run_contingency_select!(type::OPF, case::Case, lim=1e-14, max_itr=max(length(case.opf.contingencies), 5)) = 
    run_contingency_select!(type, case.model, case.opf, case.pf, case.oplim, case.brc_up, case.brc_down, case.Pc, case.Pcc, case.Pccx, lim, max_itr)

function run_contingency_select!(
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
    lim=1e-14,
    max_itr=max(length(contingencies), 5)
)
    assert(type)
    @assert !type.C1 || isempty(Pc)
    @assert !type.C2 || isempty(Pcc)
    @assert !type.C2F || isempty(Pccx)
    
    total_solve_time = solve_time(m)
    termination_status(m) != MOI.OPTIMAL && return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time
    @debug "lower_bound = $(objective_value(m))"

    # Set variables
    !isempty(oplim.dist_slack) && set_dist_slack!(pf.ϕ, opf.mgx, oplim.dist_slack)
    bd = benders(opf, m)
    
    overloads = zeros(length(opf.contingencies))
    pre = 0
    corr1 = 0
    corr2 = 0
    corr2f = 0

    pg = get_value(m, :pg0)
    for (i, c_obj) in enumerate(opf.contingencies)
        cont = typesort_component(c_obj, opf)
        if is_islanded(pf, cont[2], cont[1])
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont[2], cont[1], pf.slack)
            ptdf = calc_isf(pf, cont[2], cont[1], islands, island, island_b)
            set_tol_zero!(ptdf)
            if type.P 
                add_contingency!(opf, pf, oplim, m, brc_up, brc_down, ptdf, i)
                pre += 1
            end
            if type.C1 
                add_contingency!(Pc, opf, pf, oplim, m, brc_up, brc_down, bd.obj, islands, island, ptdf, i)
                corr1 += 1
            end
            if type.C2  
                add_contingency!(Pcc, opf, pf, oplim, m, brc_up, brc_down, bd.obj, islands, island, ptdf, i)
                corr2 += 1
            end
            if type.C2F 
                add_contingency!(Pccx, opf, pf, oplim, m, brc_up, brc_down, bd.obj, islands, island, ptdf, i)
                corr2f += 1
            end
            @debug "Island: Contingency $(string(typeof(c_obj))) $(get_name(c_obj))"
            overloads[i] = -1.0
        else
            if typeof(c_obj) <: ACBranch
                flow = calculate_contingency_line_flows(m, pf, bd.Pᵢ, cont, i, c_obj, Int[], Int[])
            else
                flow = calculate_contingency_line_flows(m, pf, bd.Pᵢ, cont, i, c_obj, Int[], Int[], pg[cont[1]])
            end
            # Calculate the power flow with the new outage and find if there are any overloads
            overload = filter_overload(flow, oplim.branch_rating * oplim.short_term_multi)
            if !isempty(overload)
                overloads[i] = maximum(x -> abs(x[2]), overload)
            end
        end
    end
    
    total_solve_time = update_model!(m, pf, oplim, brc_up, brc_down, bd, total_solve_time)
    termination_status(m) != MOI.OPTIMAL && return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time

    ΔPc = zeros(length(bd.Pg))
    ΔPcc = zeros(length(bd.Pg))
    ΔPccx = zeros(length(bd.Pg))

    olc = Vector{Tuple{Int,Float64}}()
    olcc = Vector{Tuple{Int,Float64}}()
    olccx = Vector{Tuple{Int,Float64}}()
    islands = Vector{Vector{Int}}()
    island = 0
    island_b = Int[]
    
    brst = oplim.branch_rating * oplim.short_term_multi
    brlt = oplim.branch_rating * oplim.long_term_multi

    permutation = sortperm(overloads, rev=true)
    cont = typesort_component(opf.contingencies[permutation[end]], opf)
    i = length(permutation)
    while is_islanded(pf, cont[2], cont[1])
        pop!(permutation)
        cont = typesort_component(opf.contingencies[permutation[end]], opf)
    end

    i = 1
    iterations = 1
    cut_added = 0
    while !isempty(permutation)
        i_c = permutation[i]
        c_obj = opf.contingencies[i_c]
        cont = typesort_component(c_obj, opf)
        if is_islanded(pf, cont[2], cont[1])
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont[2], cont[1], pf.slack)
            inodes = islands[island]
        else
            empty!(islands)
            island = 0
            empty!(island_b)
            inodes = Int[]
        end

        if type.P
            olc = calculate_contingency_overload(brst, m, pf, bd, cont, i_c, c_obj, inodes, island_b)
        end
        if type.C1
            olc = calculate_contingency_overload!(ΔPc, brst, Pc, opf, m, pf, bd, cont, i_c, c_obj, inodes, island_b)
        end
        if type.C2
            olcc = calculate_contingency_overload!(ΔPcc, brlt, Pcc, opf, m, pf, bd, cont, i_c, c_obj, inodes, island_b)
        end
        if type.C2F
            olccx = calculate_contingency_overload!(ΔPccx, brlt, Pccx, opf, m, pf, bd, cont, i_c, c_obj, inodes, island_b)
        end
        if !isempty(olc) || !isempty(olcc) || !isempty(olccx) # ptdf calculation is more computational expensive than line flow
            if is_islanded(pf, cont[2], cont[1])
                ptdf = calc_isf(pf, cont[2], cont[1], islands, island, island_b)
            else
                ptdf = calc_isf(pf, cont[2], cont[1])
            end
            set_tol_zero!(ptdf)
        end

        # Cannot change the model before all data is exctracted!
        if !isempty(olc)
            if type.P
                add_contingency!(opf, pf, oplim, m, brc_up, brc_down, ptdf, i_c)
                pre += 1
                @info @sprintf "Pre %d: Contingency %s %s" pre string(typeof(c_obj)) PowerSystems.get_name(c_obj)
            else
                add_contingency!(Pc, opf, pf, oplim, m, brc_up, brc_down, bd.obj, islands, island, ptdf, i_c)
                corr1 += 1
                @info @sprintf "Corr1 %d: Contingency %s %s" corr1 string(typeof(c_obj)) PowerSystems.get_name(c_obj)
            end
            cut_added = 2
        end
        if !isempty(olcc)
            add_contingency!(Pcc, opf, pf, oplim, m, brc_up, brc_down, bd.obj, islands, island, ptdf, i_c)
            corr2 += 1
            @info @sprintf "Corr2 %d: Contingency %s %s" corr2 string(typeof(c_obj)) PowerSystems.get_name(c_obj)
            cut_added = 2
        end
        if !isempty(olccx)
            add_contingency!(Pccx, opf, pf, oplim, m, brc_up, brc_down, bd.obj, islands, island, ptdf, i_c)
            corr2f += 1
            @info @sprintf "Corr2 %d: Contingency %s %s" corr2f string(typeof(c_obj)) PowerSystems.get_name(c_obj)
            cut_added = 2
        end
        if cut_added > 1
            total_solve_time = update_model!(m, pf, oplim, brc_up, brc_down, bd, total_solve_time)
            deleteat!(permutation, i)
            cut_added = 1
        else
            i += 1
        end
        termination_status(m) != MOI.OPTIMAL && return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time
        if i > length(permutation)
            iterations > max_itr && break
            if cut_added == 0 # loops until no new cuts are added for the contingencies
                # @printf "END: Total solve time %.4f.\n" total_solve_time
                print_cuts(type, pre, corr1, corr2, corr2f)
                return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time
            else
                cut_added = 0
            end
            iterations += 1
            @info "\n------------------\nIteration $iterations"
            i = 1
        end
    end
    # @warn "Reached $(length(opf.contingencies)) iterations without a stable solution."
    return Case(m, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx), total_solve_time
end
