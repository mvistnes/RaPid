# julia.NumThreads": "4" # On laptop with 4 physical cores

function test_g500()
    # system, voll, contingencies, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec = setup_system(joinpath("data","matpower","RTS_GMLC.m"))
    system, voll, contingencies, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec = setup_system(joinpath("data","matpower","ACTIVSg500.m"))
    # system, voll, contingencies, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec = setup_system(joinpath("data","matpower","ACTIVSg2000.m"))
    # system, voll, contingencies, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec = setup_system(joinpath("data","matpower","case_ACTIVSg10k.m"))

    # @time opfm_norm = SCOPF.scopf(SCOPF.SC, system, SCOPF.Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=1.0,
    #     ramp_minutes=ramp_minutes, short_term_limit_multi=short, debug=true);
    # @time SCOPF.solve_model!(opfm_norm.mod);
    @time mod_ptdf, opf_ptdf, pf_ptdf, oplim_ptdf, Pc_ptdf, Pcc_ptdf, Pccx_ptdf = SCOPF.opf_base(SCOPF.P_SCOPF, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    @time SCOPF.solve_model!(mod_ptdf);

    println("Start Contingency select")
    @time mod_cont, opf_cont, pf_cont, oplim_cont, Pc_cont, Pcc_cont, Pccx_cont, tot_t = SCOPF.run_contingency_select(SCOPF.PC_SCOPF, system, Gurobi.Optimizer, voll, prob, contingencies, 
        max_shed=max_shed, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00);
    println("Solver time: ", tot_t, " Objective value Contingency select: ", JuMP.objective_value(mod_cont))

    println("Start Benders")
    @time model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = SCOPF.run_benders(SCOPF.PC_SCOPF, system, Gurobi.Optimizer, voll, prob, contingencies, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00);
    println("Solver time: ", tot_t, " Objective value Benders: ", JuMP.objective_value(model))

    println("Start PTDF")
    @time mod_ptdf, opf_ptdf, pf_ptdf, oplim_ptdf, Pc_ptdf, Pcc_ptdf, Pccx_ptdf = SCOPF.opf(SCOPF.PC_SCOPF, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    @time SCOPF.solve_model!(mod_ptdf);
    println("Objective value PTDF: ", JuMP.objective_value(mod_ptdf))
    
    # opfm_norm = SCOPF.scopf(SCOPF.PCSC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=1.0, 
    #     ramp_minutes=10, short_term_limit_multi=short, debug=true);
    # SCOPF.solve_model!(opfm_norm.mod);

    return (opfm, pf, Pc, Pcc, Pccx), (opfm_cont, pf_cont, Pc_cont, Pcc_cont, Pccx_cont), (opfm_ptdf, pf_ptdf, Pc_ptdf, Pcc_ptdf)
end