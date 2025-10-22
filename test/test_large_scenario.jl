import Dates, DataFrames, CSV, Plots, StatsPlots, Measures, Gurobi, Printf, FileIO, JLD2, Random, Logging, JuMP
using PowerSystems, LaTeXStrings
import RaPidSCOPF as SCOPF
const GRB_ENV = Gurobi.Env()
Random.seed!(53715)

for (n, fname) in [
        "g500"=>joinpath("data", "matpower", "ACTIVSg500.m"), 
        "g2000"=>joinpath("data", "matpower", "ACTIVSg2000.m"), 
        "g10k"=>joinpath("data", "matpower", "case_ACTIVSg10k.m")
        ]

    system = SCOPF.System(fname)
    voll, prob, contingencies = SCOPF.setup(system, 100.0, 400.0)
    SCOPF.set_ramp_limits!(system)
    SCOPF.set_renewable_prod!(system, 0.5)
    prob = fill(0.01, length(contingencies))
    short = 1.5
    long = 1.25
    ramp_minutes = 10.0
    max_shed = sum(get_active_power.(SCOPF.get_demands(system))) * 0.05
    ramp_mult = 2.0
    time_limit_sec = length(contingencies)^2 + 10
    branches = SCOPF.sort_components!(SCOPF.get_branches(system));
    ia_br = SCOPF.get_interarea(branches)
    if isempty(ia_br)
        ia_br = sort(branches, by=x->x.rate, rev=true)[1:50]
    end

    w = [branches[i] ∈ ia_br ? Storm(x/100/8760, 0.01*(0.5+rand()), rand(3:4), 1, rand(1:2), rand(1:2)) : storm() for (i,x) in enumerate(prob)]
    p = generate_p_from_weather(prob/8760, 7, w)
    scenarioes = generate_contingency_probabilities(p, 5, 10000)
    scenarioes_n1 = scenario_plus_n1(scenarioes, length(branches))

    SCOPF.set_ramp_limits!(system)
    prev_lim = long
    reserve = 0.0
    ramp_mult = 2.0
    renew_cost = 0.00
    renew_ramp = 0
    labels = ["NkC" "NkCz" "NkCz2" "N1C" "N1P"]
    # labels = ["NkC" "NkCz" "NkCz2" "NkP" "N1C" "N1P" "N0"]
    # xguidefont=fsize, yguidefont=fsize, xtickfont=fsize, ytickfont=fsize, legendfont=fsize
    name = Printf.@sprintf "%s_prev%d_ls%d_ramp%d_r%d" n prev_lim*100 max_shed*100 renew_ramp reserve*100
    println("\n\n"*name*"\n")
    try
        case, tot_t = SCOPF.run_decomposition(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, scenarioes_n1[2][2], scenarioes_n1[2][1], max_shed=max_shed, reserve=reserve, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10);
        write_to_file(name*"_h2", SCOPF.extract_results(case))
        case, tot_t = SCOPF.run_decomposition(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, scenarioes_n1[5][2], scenarioes_n1[5][1], max_shed=max_shed, reserve=reserve, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10);
        write_to_file(name*"_h5", SCOPF.extract_results(case))
        results = run_plotsave_weather(name, system, labels, scenarioes_n1, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim);
    catch e
        println(e)
    end

end