# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
include("C:\\Users\\matiaskv\\OneDrive - NTNU\\Documents\\PhD\\SCOPF\\src\\SCOPF.jl")
using PowerSystems
import JuMP
import Gurobi
import Logging
import Random
Random.seed!(42)
Logging.disable_logging(Logging.Info)

function comparison(fname::String)
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

function run_timed_cases(sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long)
    run_time = []
    t = @timed begin model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.SC, sys, Gurobi.Optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    SCOPF.solve_model!(model);
    end
    push!(run_time, (t.time, SCOPF.solve_time(model)))
    idx = SCOPF.get_nodes_idx(opf.nodes);
    list = SCOPF.make_list(opf, idx, opf.nodes);

    t = @timed begin model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.SC, sys, Gurobi.Optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.add_all_contingencies!(SCOPF.PSC, SCOPF.PCSC, opf, oplim, model, list, pf, idx, Pc, Pcc, Pccx)
    SCOPF.solve_model!(model);
    end
    push!(run_time, (t.time, SCOPF.solve_time(model)))

    t = @timed begin model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.PSC, sys, Gurobi.Optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    SCOPF.solve_model!(model);
    end
    push!(run_time, (t.time, SCOPF.solve_time(model)))

    t = @timed begin model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.PCSC, sys, Gurobi.Optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    SCOPF.solve_model!(model);
    end
    push!(run_time, (t.time, SCOPF.solve_time(model)))
    return run_time
end

function run_timed_contingency_select_cases(sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=0.00)
    run_time = []
    t = @timed begin model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.SC, sys, Gurobi.Optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    SCOPF.solve_model!(model);
    end
    push!(run_time, (t.time, SCOPF.solve_time(model)))
    
    t = @timed begin model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.SC, sys, Gurobi.Optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = SCOPF.run_contingency_select!(SCOPF.PSC, SCOPF.PCSC, model, opf, pf, oplim, Pc, Pcc, Pccx)
    end
    push!(run_time, (t.time, tot_t))

    t = @timed begin model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = SCOPF.run_contingency_select(SCOPF.PSC, sys, Gurobi.Optimizer, voll, prob, cont, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure);
    end
    push!(run_time, (t.time, tot_t))

    t = @timed begin model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = SCOPF.run_contingency_select(SCOPF.PCSC, sys, Gurobi.Optimizer, voll, prob, cont, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure);
    end
    push!(run_time, (t.time, tot_t))
        
end

function run_timed_benders_cases(sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=0.00)
    run_time = []
    t = @timed begin model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.SC, sys, Gurobi.Optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    SCOPF.solve_model!(model);
    end
    push!(run_time, (t.time, SCOPF.solve_time(model)))

    t = @timed begin model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.SC, sys, Gurobi.Optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = SCOPF.run_benders!(SCOPF.PSC, SCOPF.PCSC, model, opf, pf, oplim, Pc, Pcc, Pccx)
    end
    push!(run_time, (t.time, tot_t))
    
    t = @timed begin model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = SCOPF.run_benders(SCOPF.PSC, sys, Gurobi.Optimizer, voll, prob, cont, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure);
    end
    push!(run_time, (t.time, tot_t))

    t = @timed begin model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = SCOPF.run_benders(SCOPF.PCSC, sys, Gurobi.Optimizer, voll, prob, cont, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure);
    end
    push!(run_time, (t.time, tot_t))
    return run_time
end
function run_all(systems)
    result = SCOPF.SystemRunData(4, length(systems))
    Logging.disable_logging(Logging.Info)
    # Threads.@threads for i in eachindex(systems)
    for i in eachindex(systems)
        sys, voll, cont, prob, short, long, ramp_minutes, ramp_mult, max_shed = systems[i]
        SCOPF.run_cases!(result, i, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long)
        println(i)
    end
    Logging.disable_logging(Logging.Debug)
    return result
end

function run_all_contingency_select(systems)
    result = SCOPF.SystemRunData(4, length(systems))
    Logging.disable_logging(Logging.Info)
    # Threads.@threads for i in eachindex(systems)
    for i in eachindex(systems)
        sys, voll, cont, prob, short, long, ramp_minutes, ramp_mult, max_shed = systems[i]
        SCOPF.run_contingency_select_cases!(result, i, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long)
        println(i)
    end
    Logging.disable_logging(Logging.Debug)
    return result
end

function run_all_benders(systems)
    result = SCOPF.SystemRunData(4, length(systems))
    Logging.disable_logging(Logging.Info)
    # Threads.@threads for i in eachindex(systems)
    for i in eachindex(systems)
        sys, voll, cont, prob, short, long, ramp_minutes, ramp_mult, max_shed = systems[i]
        SCOPF.run_benders_cases!(result, i, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long)
        println(i)
    end
    Logging.disable_logging(Logging.Debug)
    return result
end

function setup_ieee_rts(fname::String)
    system = SCOPF.System(fname)
    SCOPF.set_rate!.(SCOPF.get_branches(system), SCOPF.get_rate.(SCOPF.get_branches(system)) * 0.8);
    SCOPF.set_operation_cost!.(SCOPF.get_gens_h(system), [15.0, 16.0, 17.0, 18.0, 19.0, 20.0])

    voll = [4304., 5098., 5245., 5419., 4834., 5585., 5785., 5192., 4575., 5244., 4478., 5698., 4465., 4859., 5032., 5256., 4598.]
    branches = SCOPF.sort_components!(SCOPF.get_branches(system));
    contingencies = branches
    prob =
        [ # spesified for the RTS-96
            0.24, 0.51, 0.33, 0.39, 0.48, 0.38, 0.02, 0.36, 0.34, 0.33, 0.30, 0.44, 0.44,
            0.02, 0.02, 0.02, 0.02, 0.40, 0.39, 0.40, 0.52, 0.49, 0.38, 0.33, 0.41, 0.41,
            0.41, 0.35, 0.34, 0.32, 0.54, 0.35, 0.35, 0.38, 0.38, 0.34, 0.34, 0.45 # branches
        ]
    # prob /= 8760
    prob /= 100
    short = 1.2
    long = 1.0
    ramp_minutes = 10.
    max_shed = 0.1
    ramp_mult = 10.
    return system, voll, contingencies, prob, short, long, ramp_minutes, ramp_mult, max_shed
end

function setup_system(fname::String)
    system = SCOPF.System(fname)
    voll, prob, contingencies = SCOPF.setup(system, 100., 400.);
    SCOPF.set_ramp_limits!(system, 0.01);
    SCOPF.set_renewable_prod!(system, 0.5)
    prob = fill(0.01, length(contingencies))
    short = 1.5
    long = 1.25
    ramp_minutes = 10.
    max_shed = 0.5
    ramp_mult = 10.
    return system, voll, contingencies, prob, short, long, ramp_minutes, ramp_mult, max_shed
end

function run_test_benders()
    systems = []
    push!(systems, setup_system(joinpath("data","matpower","case5.m")))
    push!(systems, setup_ieee_rts(joinpath("data","matpower","IEEE_RTS.m")))
    push!(systems, setup_system(joinpath("data","matpower","RTS_GMLC.m")))
    # push!(systems, setup_system(joinpath("data","matpower","ACTIVSg500.m")))
    # push!(systems, setup_system(joinpath("data","matpower","ACTIVSg2000.m")))
    # push!(systems, setup_system(joinpath("data","matpower","caseACTIVSg10k.m")))
    results1 = run_all_benders(systems);
    results2 = run_all_contingency_select(systems);
    results3 = run_all(systems);
    println("Benders")
    SCOPF.print_data(results1)
    println("Cont")
    SCOPF.print_data(results2)
    println("PTDF")
    SCOPF.print_data(results3)

    rt = []
    for x in systems
        system, voll, contingencies, prob, short, long, ramp_minutes, ramp_mult, max_shed = x
        t = run_timed_cases(system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long)
        push!(rt, t)
        t = run_timed_contingency_select_cases(system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long)
        push!(rt, t)
        t = run_timed_benders_cases(system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long)
        push!(rt, t)
    end
    return rt
end
