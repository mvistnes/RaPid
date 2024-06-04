
function set_time_series_value!(vec::Vector{<:Component}, t::Dates.DateTime)
    for x in vec
        v = get_time_series_values(SingleTimeSeries, x, "max_active_power", 
            start_time = t, len=1)
        set_active_power!(x, v[1])
    end
end

function run_corrective_with_base(case::SCOPF.Case, vals)
    JuMP.fix.(case.model[:pg0], vals[:pg0], force=true)
    JuMP.fix.(case.model[:pfdc0], vals[:pfdc0], force=true)
    JuMP.fix.(case.model[:ls0], vals[:ls0], force=true)
    SCOPF.solve_model!(case.model)
    @time SCOPF.run_benders!(SCOPF.PC2_SCOPF - SCOPF.Base_SCOPF, case)
    
    if JuMP.is_solved_and_feasible(case.model)
        new_vals = SCOPF.extract_results(case)
        return new_vals
    else
        JuMP.delete_upper_bound.(case.model[:pgru])
        # JuMP.delete_upper_bound.(case.model[:pgrd])
        @time SCOPF.solve_model!(case.model)
        if JuMP.is_solved_and_feasible(case.model)
            new_vals = SCOPF.extract_results(case)
            return new_vals
        else
            vals[:obj_val] = NaN
            return vals
        end
    end
end

function run_corrective_with_base(vals, system, scenarioes, voll, max_shed, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long)
    res = []
    for (i,s) in enumerate(scenarioes)
        println("\n", i)
        case_i, _ = SCOPF.run_benders(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(), voll, s[2], s[1], max_shed=max_shed, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10);
        push!(res, run_corrective_with_base(case_i, vals))
    end
    return res
end

function run_cases(system, scenarioes, voll, prob, contingencies, max_shed, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim)
    Logging.disable_logging(Logging.Info)
    res = []
    for (i,s) in enumerate(scenarioes)
        println("\nScenario prob ", i)
        @time case, tot_t = SCOPF.run_benders(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(), voll, s[2], s[1], max_shed=max_shed, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10) 
        push!(res, SCOPF.extract_results(case))
    end

    println("\nN-1 prob")
    @time case, tot_t = SCOPF.run_benders(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(), voll, prob/8760, contingencies, max_shed=max_shed, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10)
    vals = SCOPF.extract_results(case)
    res_n1 = run_corrective_with_base(vals, system, scenarioes, voll, max_shed, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long)

    res_prev = []
    for (i,s) in enumerate(scenarioes)
        println("\nScenario preventive ", i)
        @time case, tot_t = SCOPF.run_benders(SCOPF.P_SCOPF, system, Gurobi.Optimizer(), voll, s[2], s[1], max_shed=max_shed, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=prev_lim, long_term_multi=long, p_failure=0.00, max_itr=10) 
        vals = SCOPF.extract_results(case)
        push!(res_prev, run_corrective_with_base(case, vals))
    end

    println("\nN-1 preventive")
    @time case, tot_t = SCOPF.run_benders(SCOPF.P_SCOPF, system, Gurobi.Optimizer(), voll, prob/8760, contingencies, max_shed=max_shed, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=prev_lim, long_term_multi=long, p_failure=0.00, max_itr=10)
    vals = SCOPF.extract_results(case)
    res_prev_n1 = run_corrective_with_base(vals, system, scenarioes, voll, max_shed, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long)
    
    println("\nN-0 ")
    case = SCOPF.Case(SCOPF.opf_base(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(), voll=voll, renew_cost=renew_cost, renew_ramp=renew_ramp, max_shed=max_shed)...)
    SCOPF.constrain_branches!(case, 0.0)
    vals = SCOPF.extract_results(case)
    res_base = run_corrective_with_base(vals, system, scenarioes, voll, max_shed, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long)
    
    return res, res_prev, res_prev_n1, res_n1, res_base
end

function plot_contingency_data(results::Vector{<:Vector}, state::Symbol, property::Symbol, labels::Matrix{String}, forced_ls::Vector{<:Dict})
    data = [sort(reduce(vcat, [[sum(i[property]) - (get(forced_ls[t], c, 0) == 0 ? 0 : forced_ls[t][c]) for (c,i) in x[state]] for (t,x) in enumerate(r)]), rev=true) for r in results]
    xmax = maximum(filter(x->!isnothing(x), findfirst.(x->x<1e-3, data)); init=0)
    pl = Plots.plot(data, labels=labels, leg=:topright, xaxis=[1,xmax], xlabel="Contingency in scenario", ylabel="Load shed", thickness_scaling=1.2)
    return pl, data
end