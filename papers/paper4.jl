

system = include("data/RTS_GMLC/config.jl")
base_costs, corrective_costs = calc_costs(results, scenarioes_n1)
write_to_file("s2_"*name, results)
write_to_file("s2_"*name, base_costs, corrective_costs)
plb, plc, plt = plotsave_weather(name, labels, base_costs, corrective_costs)
plb, plc, plt = plotsave_weather(name, base_costs, corrective_costs)  

fname = joinpath("data", "matpower", "ACTIVSg2000.m")
n = "g2000"
case, tot_t = SCOPF.run_benders(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, prob, contingencies, max_shed=max_shed, reserve=reserve, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10); 

costs = distributed_costs(results, case.opf, scenarioes_n1, [5])
Plots.plot!(yaxis=[900,1000])
Plots.savefig("results/"*name*"_distributed_costs_hour5_zoom.pdf")
Plots.plot!(yaxis=[0,1000])
Plots.savefig("results/"*name*"_distributed_costs_hour5.pdf")

max_shed = sum(get_active_power.(SCOPF.get_demands(system))) * 0.01  
system = include("data/RTS_GMLC/config.jl")
case, tot_t = SCOPF.run_benders(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, prob, contingencies, max_shed=max_shed, reserve=reserve, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10); 
costs2 = distributed_costs(results2, case.opf, scenarioes2_n1, [26])  
costs = distributed_costs([results_w2[1:end-1]..., results2[end-1:end]...], case.opf, scenarioes2_n1, [26])
Plots.savefig("results/"*name*"_distributed_costs_w2_hour26.pdf")

set_parameters(demands, demand_mult, renewables, renew_mult, t)
results_w2 = FileIO.load("results/s2"*name*".jld2", "results");
results2 = FileIO.load("results/w2_"*name*".jld2", "results");
costs2 = distributed_costs(results2, case.opf, scenarioes2_n1, [20]) 
costs = distributed_costs([results_w2[1:end-1]..., results2[end-1:end]...], case.opf, scenarioes2_n1, [20])
Plots.savefig("results/"*name*"_distributed_costs_w2_hour20.pdf")
results = FileIO.load("results/"*name*".jld2", "results");
