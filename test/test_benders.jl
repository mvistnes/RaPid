using PowerSystems
using DataFrames
import JuMP
import Gurobi
import Logging
import Random
Random.seed!(42)

function run_timed_types(types, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    run_time = []
    for type in types
        t = @timed begin
            model, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx = SCOPF.opf_base(type, sys, optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
                ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, time_limit_sec=time_limit_sec)
            SCOPF.solve_model!(model)
        end
        push!(run_time, (t.time, SCOPF.solve_time(model)))
    end
    return run_time
end

function run_timed_contingency_select_types(types, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    run_time = []
    for type in types
        t = @timed case, tot_t = SCOPF.run_contingency_select(type, sys, optimizer, voll, prob, cont, max_shed=max_shed,
            ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure, time_limit_sec=time_limit_sec)
        push!(run_time, (t.time, tot_t))
    end
    return run_time
end

function run_timed_benders_types(types, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    run_time = []
    for type in types
        t = @timed case, tot_t = SCOPF.run_benders(type, sys, optimizer, voll, prob, cont, max_shed=max_shed,
            ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure, time_limit_sec=time_limit_sec)
        push!(run_time, (t.time, tot_t))
    end
    return run_time
end

" func= run_types! or run_contingency_select_types! or run_benders_types!"
function run_all(systems, optimizer, types, func)
    result = SCOPF.SystemRunData(4, length(systems))
    Logging.disable_logging(Logging.Info)
    # Threads.@threads for i in eachindex(systems)
    for i in eachindex(systems)
        sys, voll, cont, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec = systems[i]
        func(result, i, types, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, time_limit_sec=time_limit_sec)
        println(i)
    end
    Logging.disable_logging(Logging.Debug)
    return result
end

function setup_ieee_rts(fname::String)
    system = SCOPF.System(fname)
    costh = [15.0, 16.0, 17.0, 18.0, 19.0, 20.0]
    for (i,g) in enumerate(SCOPF.get_gens_h(system))
        SCOPF.set_operation_cost!(g, costh[mod(i,length(costh))+1])
    end

    # voll = [4304.0, 5098.0, 5245.0, 5419.0, 4834.0, 5585.0, 5785.0, 5192.0, 4575.0, 5244.0, 4478.0, 5698.0, 4465.0, 4859.0, 5032.0, 5256.0, 4598.0]
    base_voll = [6.20 4.89 5.30 5.62 6.11 5.50 5.41 5.40 2.30 4.14 5.39 3.41 3.01 3.54 3.75 2.29 3.64] * 100
    if length(SCOPF.get_demands(system)) == length(base_voll)
        voll = vec(base_voll)
        SCOPF.set_rate!.(SCOPF.get_branches(system), SCOPF.get_rate.(SCOPF.get_branches(system)) * 0.8)
    else
        voll = vec([base_voll... base_voll... base_voll...])
        SCOPF.set_renewable_prod!(system, 0.2)
    end
    branches = SCOPF.sort_components!(SCOPF.get_branches(system))
    contingencies = branches
    prob =
        [ # spesified for the RTS-96
            0.24, 0.51, 0.33, 0.39, 0.48, 0.38, 0.02, 0.36, 0.34, 0.33, 0.30, 0.44, 0.44,
            0.02, 0.02, 0.02, 0.02, 0.40, 0.39, 0.40, 0.52, 0.49, 0.38, 0.33, 0.41, 0.41,
            0.41, 0.35, 0.34, 0.32, 0.54, 0.35, 0.35, 0.38, 0.38, 0.34, 0.34, 0.45 # branches
        ]
    # prob /= 8760
    prob /= 100
    if length(contingencies) > length(prob)
        prob = [prob[mod(i,length(prob))+1] for i in 1:length(contingencies)]
    end
    short = 1.2
    long = 1.0
    ramp_minutes = 10.0
    max_shed = sum(get_active_power.(SCOPF.get_demands(system)))/100*10
    ramp_mult = 2.0
    time_limit_sec = 100
    return system, voll, contingencies, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec
end

function setup_system(fname::String)
    system = SCOPF.System(fname)
    voll, prob, contingencies = SCOPF.setup(system, 100.0, 400.0)
    SCOPF.set_ramp_limits!(system)
    SCOPF.set_renewable_prod!(system, 0.5)
    prob = fill(1/length(contingencies), length(contingencies))
    short = 1.5
    long = 1.25
    ramp_minutes = 10.0
    max_shed = sum(get_active_power.(SCOPF.get_demands(system)))/100*10
    ramp_mult = 2.0
    time_limit_sec = length(contingencies)^2 + 10
    return system, voll, contingencies, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec
end

function run_test_benders()
    optimizer = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "Threads" => Threads.nthreads())
    # optimizer = JuMP.optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.GLP_MSG_ON)
    systems = []
    push!(systems, setup_system(joinpath("data", "matpower", "case5.m")))
    push!(systems, setup_ieee_rts(joinpath("data", "matpower", "IEEE_RTS.m")))
    # push!(systems, setup_system(joinpath("data", "matpower", "RTS_GMLC.m")))
    push!(systems, setup_system(joinpath("data","matpower","ACTIVSg500.m")))
    push!(systems, setup_system(joinpath("data","matpower","ACTIVSg2000.m")))
    push!(systems, setup_system(joinpath("data", "matpower7", "case2383wp.m")))
    push!(systems, setup_system(joinpath("data", "matpower7", "case2746wp.m")))
    # push!(systems, setup_system(joinpath("data","matpower","case_ACTIVSg10k.m")))

    # types = [SCOPF.Base_SCOPF, SCOPF.P_SCOPF, SCOPF.OPF(true, false, true, false, false), SCOPF.PC2_SCOPF]
    # results1 = run_all(systems, optimizer, types, SCOPF.run_benders_types!)
    # # results2 = run_all_contingency_select(systems, optimizer, types);
    # results3 = run_all(systems, optimizer, types, SCOPF.run_types!)
    # println("Benders")
    # SCOPF.print_data(results1)
    # println("Cont")
    # SCOPF.print_data(results2)
    # println("PTDF")
    # SCOPF.print_data(results3)

    rt = Dict()
    for x in systems
        system, voll, contingencies, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec = x
        # t = run_timed_types(types, optimizer, system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, time_limit_sec=time_limit_sec)
        # push!(rt, "ptdf" => t)
        # t = run_timed_contingency_select_types(types, optimizer, system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, time_limit_sec=time_limit_sec)
        # push!(rt, "cont"=>t)
        # t = run_timed_benders_types(types, optimizer, system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, time_limit_sec=time_limit_sec)
        rt[length(get_components(ACBus, system))] = run_test_scopf(system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, time_limit_sec)
    end
    return rt
end

function run_test_scopf(system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, dist_slack, time_limit_sec)
    result = Dict()
    for type in [SCOPF.Base_SCOPF, SCOPF.P_SCOPF, SCOPF.PC2_SCOPF]
        println(type)
        SCOPF.run_benders_type!(result, type, SCOPF.PC2_SCOPF, Gurobi.Optimizer(GRB_ENV), Gurobi.Optimizer(GRB_ENV), 
            system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, dist_slack, time_limit_sec)
    end
    return result
end

function run_test_scopf_simple(system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, dist_slack, time_limit_sec; all_post_c=true)
    result = Dict()
    for type in [SCOPF.Base_SCOPF, SCOPF.P_SCOPF, SCOPF.OPF(true, false, true, false, false), SCOPF.OPF(true, false, false, true, false), SCOPF.PC_SCOPF, SCOPF.PC2_SCOPF]
        println(type)
        case, _ = SCOPF.run_benders(type, system, Gurobi.Optimizer(GRB_ENV), voll, prob, contingencies, max_shed=max_shed, 
            ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, dist_slack=dist_slack, 
            time_limit_sec=time_limit_sec, all_post_c=all_post_c);
        result[type] = SCOPF.extract_results(case)
    end
    return result
end

function run_tests(system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, dist_slack, time_limit_sec; all_post_c=true)
    result = []
    logrange(start,stepmul,length) = start .* stepmul .^ (0:(length-1))
    for x in logrange(0.1, 2, 10)
        push!(result, x => run_test_scopf_simple(system, voll*x, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, 
            long, dist_slack, time_limit_sec, all_post_c=all_post_c))
    end
    return result
end

function run_range(system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, dist_slack, time_limit_sec; all_post_c=true)
    result = []
    logrange(start,stepmul,length) = start .* stepmul .^ (0:(length-1))
    for x in logrange(0.1, 2, 10)
        case, tot_t = SCOPF.run_benders(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll*x, prob, contingencies, 
            max_shed=max_shed, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, 
            dist_slack=dist_slack, time_limit_sec=time_limit_sec, all_post_c=all_post_c);
        push!(result, x => SCOPF.extract_results(case))
    end
    return result
end

function run_load_profile(system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, dist_slack, time_limit_sec; all_post_c=true)
    demands = SCOPF.sort_components!((SCOPF.get_demands(system)))
    pd0 = get_active_power.(demands)
    result = []
    profile = parse.(Float64, readlines("data/ieee_std_load_profile.txt"))
    for x in profile
        set_active_power!.(demands, pd0*x)
        push!(result, x => run_test_scopf_simple(system, voll*x, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, dist_slack, 
            time_limit_sec, all_post_c=all_post_c))
    end
    return result
end

const GRB_ENV = Gurobi.Env()

function run_random(system, voll, prob, contingencies, all_post_c, iterations)
    result = []
    for i in 1:iterations
        x = rand(0.1:0.01:10)
        y = rand(0.1:0.01:10)
        max_shed = rand(0.1:0.01:50)
        ramp_mult = rand(1:0.01:10)
        ramp_minutes = 10.
        short = rand(1.25:0.01:1.5)
        long = rand(1.0:0.01:1.25)
        reset_timer!(Main.to)
        case, tot_t = SCOPF.run_benders(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll*x, prob*y, contingencies, max_shed=max_shed, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, all_post_c=all_post_c);
        vals = SCOPF.extract_results(case)
        vals[:t] = tot_t
        push!(result, (i, x, y, max_shed, ramp_mult, ramp_minutes, short, long) => vals)
    end
    return result
end

function run_random_p(system, voll, prob, contingencies, all_post_c, iterations)
    result = []
    for i in 1:iterations
        x = rand(0.1:0.01:10)
        y = 1.0
        max_shed = 1.0
        ramp_mult = 1.0
        ramp_minutes = 10.
        short = rand(1.0:0.01:1.25)
        long = 1.0
        reset_timer!(Main.to)
        case, tot_t = SCOPF.run_benders(SCOPF.P_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll*x, prob*y, contingencies, max_shed=max_shed, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, all_post_c=all_post_c);
        vals = SCOPF.extract_results(case)
        vals[:t] = tot_t
        push!(result, (i, x, y, max_shed, ramp_mult, ramp_minutes, short, long) => vals)
    end
    return result
end

make_dataframe(result) = DataFrames.DataFrame(
    :VOLL => getindex.(getindex.(result, 1), 2), 
    :p => getindex.(getindex.(result, 1), 3), 
    :max_shed => getindex.(getindex.(result, 1), 4), 
    :ramp_mult => getindex.(getindex.(result, 1), 5), 
    :short => getindex.(getindex.(result, 1), 7), 
    :long => getindex.(getindex.(result, 1), 8), 
    :pg0 => [sum(getproperty.(case.opf.cost_ctrl_gen, :var) .* y[2][:pg0]) for y in result], 
    :lsc => [sum((sum(x[:lsc], init=0.0) for (k,x) in y[2][:Pc]), init=0.0) for y in result], 
    :lscc => [sum((sum(x[:lscc], init=0.0) for (k,x) in y[2][:Pcc]), init=0.0) for y in result], 
    :cens_c1 => [sum((sum((case.opf.voll .* x[:lsc]), init=0.0) for (k,x) in y[2][:Pc]), init=0.0) for y in result], 
    :cens_c2 => [sum((sum((case.opf.voll .* x[:lscc]), init=0.0) for (k,x) in y[2][:Pcc]), init=0.0) for y in result], 
    :obj_val => getindex.(getindex.(result, 2), :obj_val))

function make_dataframe_big(result) 
    df = DataFrames.DataFrame(
        type = SCOPF.OPF[],
        pd0 = Float64[],
        obj_val = Float64[],
        pg0 = Float64[],
        lsc = Float64[],
        lscc = Float64[])
    for i in result
        for (type,x) in i[2]
            push!(df, [type, i[1], x[:obj_val], sum(getproperty.(case.opf.cost_ctrl_gen, :var) .* x[:pg0]),
            sum((sum((case.opf.voll .* y[:lsc]), init=0.0) for (k,y) in x[:Pc]), init=0.0),
            sum((sum((case.opf.voll .* y[:lscc]), init=0.0) for (k,y) in x[:Pcc]), init=0.0)])
        end
    end
    return df
end

af_categorical_colours = ["#12436D" "#28A197" "#801650" "#F46A25" "#3D3D3D" "#A285D1"]

function scatter_plot(df::DataFrames.DataFrame, x::Symbol, y::Symbol, mark::Symbol=:none)
    if mark == :none
        Plots.scatter(df[!,x], df[!,y],  
            xlabel=string(x), ylabel=string(y), 
            leg=:none)
    else
        Plots.scatter(df[!,x], df[!,y], 
            markersize=(df[!,mark] .- minimum(df[!,mark]) )/(maximum(df[!,mark])-minimum(df[!,mark]))*8. .+ 2., 
            xlabel=string(x), ylabel=string(y), title=string(mark), 
            leg=:none)
    end
end