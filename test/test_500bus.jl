# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

# julia.NumThreads": "4" # On laptop with 4 physical cores

function test_g500()
    system = SCOPF.System("data\\matpower\\ACTIVSg500.m")
    # system = System("data\\matpower\\ACTIVSg2000.m")
    # system = System("data\\matpower\\case_ACTIVSg10k.m")
    voll, prob, contingencies = SCOPF.setup(system, 100., 400.);
    # SCOPF.fix_generation_cost!(system);

    # nodes = SCOPF.sort_components!(SCOPF.get_nodes(system))
    # active_capacities = zeros(length(nodes))
    # reactive_capacities = zeros(length(nodes))
    # cost = zeros(length(nodes))
    # for g in SCOPF.get_generation(system)
    #     if typeof(g) == HydroDispatch
    #         SCOPF.set_operation_cost!(g, 30)
    #     end
    # end
    # for g in SCOPF.get_ctrl_generation(system)
    #     n = g.bus.number
    #     active_capacities[n] += g.active_power_limits.max
    #     reactive_capacities[n] += g.reactive_power_limits.max
    #     cost[n] = max(SCOPF.get_generator_cost(g)[2], cost[n])
    # end
    # PowerSystems.remove_component!.([system], SCOPF.get_ctrl_generation(system))
    # for (i, (ac, rc, c)) in enumerate(zip(active_capacities, reactive_capacities, cost))
    #     if ac > 0.0
    #         PowerSystems.add_component!(system, PowerSystems.ThermalStandard(string(i), true, true, 
    #             nodes[i], ac, rc, sqrt(ac^2+rc^2), (min=0.0, max=ac), (min=-rc, max=rc), (up=0.01*ac, down=0.01*ac), 
    #             PowerSystems.ThreePartCost(PowerSystems.VariableCost((0.0, c)), 0.0, 0.0, 0.0), 100.0))
    #     end
    # end
    SCOPF.set_ramp_limits!(system, 0.01);
    SCOPF.set_renewable_prod!(system, 0.5)
    # set_rate!.(SCOPF.get_branches(system), get_rate.(SCOPF.get_branches(system))*0.8);

    # voll = fill(6000, length(SCOPF.get_demands(system)))
    # branches = SCOPF.sort_components!(SCOPF.get_branches(system));
    # contingencies = branches[1:10:end]
    prob = fill(0.01, length(contingencies))
    short = 1.5
    long = 1.25
    ramp_minutes = 10.
    max_shed = 0.5
    ramp_mult = 10.

    # @time opfm_norm = SCOPF.scopf(SCOPF.SC, system, SCOPF.Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=1.0,
    #     ramp_minutes=ramp_minutes, short_term_limit_multi=short, debug=true);
    # @time SCOPF.solve_model!(opfm_norm.mod);
    @time mod_ptdf, opf_ptdf, pf_ptdf, oplim_ptdf, Pc_ptdf, Pcc_ptdf, Pccx_ptdf = SCOPF.opf_base(SCOPF.SC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    @time SCOPF.solve_model!(mod_ptdf);

    println("Start Contingency select")
    @time mod_cont, opf_cont, pf_cont, oplim_cont, Pc_cont, Pcc_cont, Pccx_cont, tot_t = SCOPF.run_contingency_select(SCOPF.PCSC, system, SCOPF.Gurobi.Optimizer, voll, prob, contingencies, 
        max_shed=max_shed, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00);
    println("Solver time: ", tot_t, " Objective value Contingency select: ", JuMP.objective_value(mod_cont))

    println("Start Benders")
    @time model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = SCOPF.run_benders(SCOPF.PCSC, system, Gurobi.Optimizer, voll, prob, contingencies, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00);
    println("Solver time: ", tot_t, " Objective value Benders: ", JuMP.objective_value(model))

    println("Start PTDF")
    @time mod_ptdf, opf_ptdf, pf_ptdf, oplim_ptdf, Pc_ptdf, Pcc_ptdf, Pccx_ptdf = SCOPF.opf(SCOPF.PCSC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    @time SCOPF.solve_model!(mod_ptdf);
    println("Objective value PTDF: ", JuMP.objective_value(mod_ptdf))
    
    # opfm_norm = SCOPF.scopf(SCOPF.PCSC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=1.0, 
    #     ramp_minutes=10, short_term_limit_multi=short, debug=true);
    # SCOPF.solve_model!(opfm_norm.mod);

    return (opfm, pf, Pc, Pcc, Pccx), (opfm_cont, pf_cont, Pc_cont, Pcc_cont, Pccx_cont), (opfm_ptdf, pf_ptdf, Pc_ptdf, Pcc_ptdf)
end