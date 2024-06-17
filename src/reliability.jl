function run_type!(result, i, c, type, goal, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    case = Case(opf_base(type, sys, optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, time_limit_sec=time_limit_sec)...);
    add_branch_constraints!(case.model, case.pf.ϕ, case.model[:p0], case.brc_up, case.brc_down, case.oplim.branch_rating)
    solve_model!(case.model)
    MOI.get(case.model, MOI.ResultCount()) < 1 && return
    if type != goal
        fix_base_case!(case.model)
        # type.C1 && fix_contingencies!(case.model, case.Pc)
        # type.C2 && fix_contingencies!(case.model, case.Pcc)
        add_all_contingencies!(goal - type, case)
        solve_model!(case.model);
    end
    push!(result, extract_results(case))
    return result
end

function run_types!(result, i, types, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    for (c, type) in enumerate(types)
        run_type!(result, i, c, type, types[end], optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=p_failure, time_limit_sec=time_limit_sec)
    end
    return result
end

function run_corrective_with_base(case::Case, vals)
    JuMP.fix.(case.model[:pg0], vals[:pg0], force=true)
    JuMP.fix.(case.model[:pfdc0], vals[:pfdc0], force=true)
    JuMP.fix.(case.model[:ls0], vals[:ls0], force=true)
    JuMP.fix.(case.model[:pr0], vals[:pr0], force=true)
    solve_model!(case.model)
    case, tot_t = run_benders!(PC2_SCOPF - Base_SCOPF, case, 0.0)
    
    if JuMP.is_solved_and_feasible(case.model)
        return extract_results(case)
    else
        for (k,x) in case.Pc
            JuMP.set_upper_bound.(x.pgu, case.oplim.rampup)
        end
        # JuMP.delete_upper_bound.(case.model[:pgrd])
        @time solve_model!(case.model)
        if JuMP.is_solved_and_feasible(case.model)
            return extract_results(case)
        else
            vals[:obj_val] = NaN
            return vals
        end
    end
end

function run_benders_type!(result, type, goal, optimizer, optimizer2, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    case, tot_t = run_benders(type, sys, optimizer, voll, prob, cont, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, time_limit_sec=time_limit_sec)
    !is_solved_and_feasible(case.model) && return
    if type != goal
        vals = extract_results(case)
        case_i, _ = run_benders(Base_SCOPF, sys, optimizer2, voll, prob, cont, max_shed=max_shed, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, time_limit_sec=time_limit_sec);
        result[type] = run_corrective_with_base(case_i, vals)
    else
        result[type] = extract_results(case)
    end
    return result
end

function run_benders_types!(result, types, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    for type in types
        run_benders_type!(result, type, types[end], optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=p_failure, time_limit_sec=time_limit_sec)
    end
    return result
end

function run_reliability_calculation_benders(types, optimizer, system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    # demands = sort_components!(get_demands(system))
    # pd = demands .|> get_active_power
    # hours = read_x_data("data\\ieee_std_load_profile.txt")
    # hours = 0.5:0.1:1.3
    hours = 0.8:0.1:1.1

    result = Dict()
    Logging.disable_logging(Logging.Info)
    # buffer_sys = [deepcopy(system) for _ in 1:Threads.nthreads()]
    # Threads.@threads for i in eachindex(hours)
    for i in eachindex(hours)
        h = hours[i]
        # sys = buffer_sys[Threads.threadid()]
        set_active_power!.(get_components(StaticLoad, system), (get_components(StaticLoad, system) .|> get_max_active_power) * h)
        # cont = sort_components!(get_branches(sys))

        run_benders_types!(result, types, optimizer, system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=p_failure, time_limit_sec=time_limit_sec)

        print(i, " ")
    end
    Logging.disable_logging(Logging.Debug)
    # set_active_power!.(get_components(StaticLoad, system), (get_components(StaticLoad, system) .|> get_max_active_power))
    return result
end

read_x_data(fname) = vec(DelimitedFiles.readdlm(fname,Float64))

get_objective(mod::Model) = termination_status(mod) != MOI.OPTIMAL ? -1.0 : objective_value(mod)
