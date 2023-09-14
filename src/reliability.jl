
function run_reliability_calculation(system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=0.00)
    # demands = sort_components!(get_demands(system))
    # pd = demands .|> get_active_power
    hours = read_x_data("data\\ieee_std_load_profile.txt")

    result = [[],[],[],[],[],[],[]]
    Logging.disable_logging(Logging.Info)
    for h in 0.8:0.05:1.5 # hours[1:50]
        set_active_power!.(get_components(StaticLoad, system), (get_components(StaticLoad, system) .|> get_max_active_power) * h)

        model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.SC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed,
            ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
        SCOPF.solve_model!(model);
        push!(result[2], get_objective(model))
        last(result[2]) < 0.0 && pop!(result[2]) < 0.0 && continue
        # oplim.max_shed = 1.0
        idx = get_nodes_idx(opf.nodes)
        list = make_list(opf, idx, opf.nodes)
        SCOPF.fix_base_case(model)
        model, opf, pf, oplim, Pc, Pcc, Pccx = add_all_contingencies!(SCOPF.SC, SCOPF.PCSC, Pc, Pcc, Pccx, opf, oplim, model, list, pf, idx)
        SCOPF.solve_model!(model);
        push!(result[3], get_objective(model))

        model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.PSC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed,
            ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
        SCOPF.solve_model!(model);
        push!(result[4], get_objective(model))
        # oplim.max_shed = 1.0
        idx = get_nodes_idx(opf.nodes)
        list = make_list(opf, idx, opf.nodes)
        SCOPF.fix_base_case(model)
        SCOPF.fix_short_term(model, Pc)
        model, opf, pf, oplim, Pc, Pcc, Pccx = add_all_contingencies!(SCOPF.PSC, SCOPF.PCSC, Pc, Pcc, Pccx, opf, oplim, model, list, pf, idx)
        SCOPF.solve_model!(model);
        push!(result[5], get_objective(model))

        model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.PCSC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
        SCOPF.solve_model!(model);
        push!(result[6], get_objective(model))

        model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.run_benders(SCOPF.PCSC, system, Gurobi.Optimizer, voll, prob, contingencies, max_shed=max_shed,
            ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure);
        push!(result[7], get_objective(model))
        # push!(result, (h, obj, obj_tot, objP, objP_tot, objPC, objB))
        push!(result[1], h)
    end
    Logging.disable_logging(Logging.Debug)
    set_active_power!.(get_components(StaticLoad, system), (get_components(StaticLoad, system) .|> get_max_active_power))
    return result
end

read_x_data(fname) = vec(readdlm(fname,Float64))

get_objective(mod::Model) = termination_status(mod) != MOI.OPTIMAL ? -1.0 : objective_value(mod)