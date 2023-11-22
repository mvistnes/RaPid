using PowerSystems
import JuMP
import Gurobi
import Logging
import Random
Random.seed!(42)

function run_timed_cases(cases, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    run_time = []
    for case in cases
        t = @timed begin model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf_base(case, sys, optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
            ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, time_limit_sec=time_limit_sec);
        SCOPF.solve_model!(model);
        end
        push!(run_time, (t.time, SCOPF.solve_time(model)))
    end
    return run_time
end

function run_timed_contingency_select_cases(cases, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    run_time = []
    for case in cases
        t = @timed model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = SCOPF.run_contingency_select(case, sys, optimizer, voll, prob, cont, max_shed=max_shed,
            ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure, time_limit_sec=time_limit_sec);
        push!(run_time, (t.time, tot_t))
    end
    return run_time        
end

function run_timed_benders_cases(cases, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    run_time = []
    for case in cases
        t = @timed model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = SCOPF.run_benders(case, sys, optimizer, voll, prob, cont, max_shed=max_shed,
            ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure, time_limit_sec=time_limit_sec);
        push!(run_time, (t.time, tot_t))
    end
    return run_time
end

function run_all(systems, optimizer, cases)
    result = SCOPF.SystemRunData(4, length(systems))
    Logging.disable_logging(Logging.Info)
    # Threads.@threads for i in eachindex(systems)
    for i in eachindex(systems)
        sys, voll, cont, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec = systems[i]
        SCOPF.run_cases!(result, i, cases, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, time_limit_sec=time_limit_sec)
        println(i)
    end
    Logging.disable_logging(Logging.Debug)
    return result
end

function run_all_contingency_select(systems, optimizer, cases)
    result = SCOPF.SystemRunData(4, length(systems))
    Logging.disable_logging(Logging.Info)
    # Threads.@threads for i in eachindex(systems)
    for i in eachindex(systems)
        sys, voll, cont, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec = systems[i]
        SCOPF.run_contingency_select_cases!(result, i, cases, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, time_limit_sec=time_limit_sec)
        println(i)
    end
    Logging.disable_logging(Logging.Debug)
    return result
end

function run_all_benders(systems, optimizer, cases)
    result = SCOPF.SystemRunData(4, length(systems))
    Logging.disable_logging(Logging.Info)
    # Threads.@threads for i in eachindex(systems)
    for i in eachindex(systems)
        sys, voll, cont, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec = systems[i]
        SCOPF.run_benders_cases!(result, i, cases, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, time_limit_sec=time_limit_sec)
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
    ramp_mult = 2.
    time_limit_sec = 100
    return system, voll, contingencies, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec
end

function setup_system(fname::String)
    system = SCOPF.System(fname)
    voll, prob, contingencies = SCOPF.setup(system, 100., 400.);
    SCOPF.set_ramp_limits!(system);
    SCOPF.set_renewable_prod!(system, 0.5)
    prob = fill(0.01, length(contingencies))
    short = 1.5
    long = 1.25
    ramp_minutes = 10.
    max_shed = 0.5
    ramp_mult = 2.
    time_limit_sec = length(contingencies)^2 + 10
    return system, voll, contingencies, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec
end

function run_test_benders()
    optimizer = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "Threads" => Threads.nthreads())
    # optimizer = JuMP.optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.GLP_MSG_ON)
    systems = []
    push!(systems, setup_system(joinpath("data","matpower","case5.m")))
    push!(systems, setup_ieee_rts(joinpath("data","matpower","IEEE_RTS.m")))
    # push!(systems, setup_system(joinpath("data","matpower","RTS_GMLC.m")))
    # push!(systems, setup_system(joinpath("data","matpower","ACTIVSg500.m")))
    # push!(systems, setup_system(joinpath("data","matpower","ACTIVSg2000.m")))
    # push!(systems, setup_system(joinpath("data","matpower","case_ACTIVSg10k.m")))
    
    cases = [SCOPF.Base_SCOPF, SCOPF.P_SCOPF, SCOPF.OPF(true, false, true, false, false), SCOPF.PC2_SCOPF]
    results1 = run_all_benders(systems, optimizer, cases);
    # results2 = run_all_contingency_select(systems, optimizer, cases);
    results3 = run_all(systems, optimizer, cases);
    println("Benders")
    SCOPF.print_data(results1)
    # println("Cont")
    # SCOPF.print_data(results2)
    println("PTDF")
    SCOPF.print_data(results3)

    rt = []
    for x in systems
        system, voll, contingencies, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec = x
        t = run_timed_cases(cases, optimizer, system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, time_limit_sec=time_limit_sec)
        push!(rt, "ptdf"=>t)
        # t = run_timed_contingency_select_cases(cases, optimizer, system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, time_limit_sec=time_limit_sec)
        # push!(rt, "cont"=>t)
        t = run_timed_benders_cases(cases, optimizer, system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, time_limit_sec=time_limit_sec)
        push!(rt, "bend"=>t)
    end
    return rt
end
