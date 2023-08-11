

function run_benders(fname::String)
    system = System(fname)
    set_rate!.(SCOPF.get_branches(system), get_rate.(SCOPF.get_branches(system))*0.8);
    SCOPF.fix_generation_cost!(system);
    voll, prob, contingencies = SCOPF.setup(system, 1, 4);

    @time opfm_sc = SCOPF.scopf(SCOPF.SC, system, Gurobi.Optimizer, voll=voll);
    @time SCOPF.solve_model!(opfm_sc.mod);
    SCOPF.print_results(opfm_sc)
    SCOPF.print_power_flow(opfm_sc)

    @time opfm, pf, Pc, Pcc = SCOPF.run_benders(SCOPF.PCSC, system, voll, prob, short_term_limit_multi=1.2);
    SCOPF.print_benders_results(opfm, Pc, Pcc)
    SCOPF.print_power_flow(opfm)

    @time opfm_norm = SCOPF.scopf(SCOPF.PCSC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, short_term_limit_multi=1.2);
    @time SCOPF.solve_model!(opfm_norm.mod);
    SCOPF.print_corrective_results(opfm_norm)
    SCOPF.print_power_flow(opfm_norm)
end

run_benders("data\\matpower\\case5.m")
run_benders("data\\matpower\\IEEE_RTS.m")
run_benders("data\\matpower\\RTS_GMLC.m")
# run_benders("data\\matpower\\ACTIVSg2000.m")
# run_benders("data\\matpower\\caseACTIVSg10k.m")