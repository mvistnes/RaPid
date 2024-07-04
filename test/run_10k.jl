# const to = TimerOutput()
const GRB_ENV = Gurobi.Env()

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

function run_random(system, voll, prob, contingencies, all_post_c, iterations, fname)
    df = DataFrames.DataFrame(
	    VOLL = Float64[],
        p = Float64[],
	    max_shed = Float64[],
	    ramp_mult = Float64[],
	    short = Float64[],
	    long = Float64[],
        pg0 = Float64[],
        lsc = Float64[],
        lscc = Float64[],
        cens_c1 = Float64[],
        cens_c2 = Float64[],
        obj_val = Float64[])
    result = []
    for i in 1:iterations
        x = rand(0.1:0.01:10)
        y = rand(0.1:0.01:10)
        max_shed = rand(0.1:0.01:50)
        ramp_mult = rand(1:0.01:10)
        ramp_minutes = 10.
        short = rand(1.25:0.01:1.5)
        long = rand(1.0:0.01:1.25)
        # reset_timer!(Main.to)
        case, tot_t = SCOPF.run_benders(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll*x, prob*y, contingencies, max_shed=max_shed, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, all_post_c=all_post_c);
        vals = SCOPF.extract_results(case)
        vals[:t] = tot_t
        push!(result, (i, x, y, max_shed, ramp_mult, ramp_minutes, short, long) => vals)
	push!(df, (x, y, max_shed, ramp_mult, short, long, 
	sum(getproperty.(case.opf.cost_ctrl_gen, :var) .* vals[:pg0]),
	sum((sum((y[:lsc]), init=0.0) for (k,y) in vals[:Pc]), init=0.0), 
	sum((sum((y[:lscc]), init=0.0) for (k,y) in vals[:Pcc]), init=0.0),
	sum((sum((case.opf.voll .* y[:lsc]), init=0.0) for (k,y) in vals[:Pc]), init=0.0), 
	sum((sum((case.opf.voll .* y[:lscc]), init=0.0) for (k,y) in vals[:Pcc]), init=0.0), vals[:obj_val]))
	CSV.write(fname, df)
    end
    return result, df
end


make_dataframe(result, case) = DataFrames.DataFrame(
    :VOLL => getindex.(getindex.(result, 1), 2), 
    :p => getindex.(getindex.(result, 1), 3), 
    :max_shed => getindex.(getindex.(result, 1), 4), 
    :ramp_mult => getindex.(getindex.(result, 1), 5), 
    :short => getindex.(getindex.(result, 1), 7), 
    :long => getindex.(getindex.(result, 1), 8), 
    :obj_val => getindex.(getindex.(result, 2), :obj_val), 
    :pg0 => [sum(getproperty.(case.opf.cost_ctrl_gen, :var) .* y[2][:pg0]) for y  in result], 
    :lsc => [sum((sum((case.opf.voll .* x[:lsc]), init=0.0) for (k,x) in y[2][:Pc]), init=0.0) for y in result], 
    :lscc => [sum((sum((case.opf.voll .* x[:lscc]), init=0.0) for (k,x) in y[2][:Pcc]), init=0.0) for y in result])


system, voll, contingencies, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec = setup_system(joinpath("data","matpower","IEEE_RTS.m"))
result, df = run_random(system, voll, prob, contingencies, true, 10, "IEEE_RTS_bus_random.csv")
println("Mean time ", mean([y[2][:t] for y  in result]))
# CSV.write("IEEE_RTS_bus_random.csv", df)