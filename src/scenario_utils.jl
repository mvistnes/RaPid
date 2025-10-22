
function set_time_series_value!(vec::Vector{<:Component}, t::Dates.DateTime)
    for x in vec
        v = get_time_series_values(SingleTimeSeries, x, "max_active_power", 
            start_time = t, len=1)
        set_active_power!(x, v[1])
    end
end

function run_corrective_with_base(case::SCOPF.Case, vals)
    JuMP.fix.(case.model[:pg0], vals[:pg0], force=true)
    # JuMP.fix.(case.model[:pgru], vals[:pgru], force=true)
    # JuMP.fix.(case.model[:pgrd], vals[:pgrd], force=true)
    JuMP.fix.(case.model[:pfdc0], vals[:pfdc0], force=true)
    JuMP.fix.(case.model[:ls0], vals[:ls0], force=true)
    SCOPF.solve_model!(case.model)
    @time SCOPF.run_decomposition!(SCOPF.PC2_SCOPF - SCOPF.Base_SCOPF, case)
    
    if JuMP.is_solved_and_feasible(case.model)
        new_vals = SCOPF.extract_results(case)
        return new_vals
    else
        # JuMP.delete_upper_bound.(case.model[:pgru])
        # # JuMP.delete_upper_bound.(case.model[:pgrd])
        # @time SCOPF.solve_model!(case.model)
        # if JuMP.is_solved_and_feasible(case.model)
        #     new_vals = SCOPF.extract_results(case)
        #     return new_vals
        # else
            vals[:obj_val] = NaN
            return vals
        # end
    end
end

function run_corrective_with_base(vals, system, scenarioes, voll, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long)
    res = []
    for (i,s) in enumerate(scenarioes)
        println("\n", i)
        case_i, _ = SCOPF.run_decomposition(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10);
        push!(res, run_corrective_with_base(case_i, vals))
    end
    return res
end

function run_cases(system, scenarioes, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim)
    Logging.disable_logging(Logging.Info)
    res = []
    for (i,s) in enumerate(scenarioes)
        println("\nScenario prob ", i)
        @time case, tot_t = SCOPF.run_decomposition(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10) 
        push!(res, SCOPF.extract_results(case))
    end

    res_zc = []
    for (i,s) in enumerate(scenarioes)
        println("\nScenario prob zc ", i)
        @time case, tot_t = SCOPF.run_decomposition(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, 
            ramp_mult=0.0, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10) 
        vals = SCOPF.extract_results(case)
        @time case, tot_t = SCOPF.run_decomposition(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, 
            ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10) 
        push!(res_zc, run_corrective_with_base(case, vals))
    end
    
    res_zvoll = []
    for (i,s) in enumerate(scenarioes)
        println("\nScenario prob zvoll ", i)
        @time case, tot_t = SCOPF.run_decomposition(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, 
            ramp_mult=0.0, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, zero_c_cost=true, max_itr=10) 
        vals = SCOPF.extract_results(case)
        @time case, tot_t = SCOPF.run_decomposition(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, 
            ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10) 
        push!(res_zvoll, run_corrective_with_base(case, vals))
    end

    println("\nN-1 prob")
    @time case, tot_t = SCOPF.run_decomposition(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, prob/8760, contingencies, max_shed=max_shed, reserve=reserve, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10)
    vals = SCOPF.extract_results(case)
    res_n1 = run_corrective_with_base(vals, system, scenarioes, voll, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long)

    res_prev = []
    for (i,s) in enumerate(scenarioes)
        println("\nScenario preventive ", i)
        @time case, tot_t = SCOPF.run_decomposition(SCOPF.P_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=prev_lim, long_term_multi=long, p_failure=0.00, max_itr=10) 
        vals = SCOPF.extract_results(case)
        push!(res_prev, run_corrective_with_base(case, vals))
    end

    println("\nN-1 preventive")
    @time case, tot_t = SCOPF.run_decomposition(SCOPF.P_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, prob/8760, contingencies, max_shed=max_shed, reserve=reserve, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=prev_lim, long_term_multi=long, p_failure=0.00, max_itr=10)
    vals = SCOPF.extract_results(case)
    res_prev_n1 = run_corrective_with_base(vals, system, scenarioes, voll, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long)
    
    println("\nN-0 ")
    case = SCOPF.Case(SCOPF.opf_base(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll=voll, renew_cost=renew_cost, renew_ramp=renew_ramp, max_shed=max_shed)...)
    JuMP.@constraint(case.model, sum(case.model[:pgru]) >= sum(case.oplim.pd_lim)*reserve)
    JuMP.@constraint(case.model, sum(case.model[:pgrd]) >= sum(case.oplim.pd_lim)*reserve)
    SCOPF.constrain_branches!(case, 0.0)
    vals = SCOPF.extract_results(case)
    res_base = run_corrective_with_base(vals, system, scenarioes, voll, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long)
    Logging.disable_logging(Logging.Debug)

    return res, res_zc, res_zvoll, res_prev, res_n1, res_prev_n1, res_base
end

function run_zc_case(system, scenarioes, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long)
    Logging.disable_logging(Logging.Info)
    res = []
    for (i,s) in enumerate(scenarioes)
        println("\nScenario prob ", i)
        @time case, tot_t = SCOPF.run_decomposition(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, 
            ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10) 
        push!(res, SCOPF.extract_results(case))
    end

    res_zc = []
    for (i,s) in enumerate(scenarioes)
        println("\nScenario prob zc ", i)
        @time case, tot_t = SCOPF.run_decomposition(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, 
            ramp_mult=0.0, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10) 
        vals = SCOPF.extract_results(case)
        @time case, tot_t = SCOPF.run_decomposition(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, 
            ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10) 
        push!(res_zc, run_corrective_with_base(case, vals))
    end
    
    res_zvoll = []
    for (i,s) in enumerate(scenarioes)
        println("\nScenario prob zvoll ", i)
        @time case, tot_t = SCOPF.run_decomposition(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, 
            ramp_mult=0.0, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, zero_c_cost=true, max_itr=10) 
        vals = SCOPF.extract_results(case)
        @time case, tot_t = SCOPF.run_decomposition(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, 
            ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10) 
        push!(res_zvoll, run_corrective_with_base(case, vals))
    end
    Logging.disable_logging(Logging.Debug)

    return res, res_zc, res_zvoll
end

function set_parameters(demands, demand_mult, renewables, renew_mult, t)
    # t = Dates.DateTime("2020-01-01T"*string(i)*":00:00")
    set_time_series_value!(demands, t)
    set_active_power!.(demands, get_active_power.(demands)*demand_mult)
    set_time_series_value!(renewables, t)
    set_active_power!.(renewables, get_active_power.(renewables)*renew_mult)
    return nothing
end
    # res, res_prev, res_n1, res_prev_n1, res_base = FileIO.load("results/"*name*".jld2", "Scen_Prob","Scen_Prev", "N-1_Prob", "N-1_Prev", "N-0");

function run_cases_weather2(results, scenarioes2, system, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim)
    Logging.disable_logging(Logging.Info)
    res2 = []
    res = results[findfirst(x->x[1]=="NkC",results)][2]
    for (i,(r, s)) in enumerate(zip(res, scenarioes2))
        println("\nScen.2 ", i, " prob")
        case, _ = SCOPF.run_decomposition(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10)
        push!(res2, run_corrective_with_base(case, r))
    end
    
    res2_zc = []
    res = results[findfirst(x->x[1]=="NkCz",results)][2]
    for (i,(r, s)) in enumerate(zip(res, scenarioes2))
        println("\nScenario prob zc ", i)
        @time case, tot_t = SCOPF.run_decomposition(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, 
            ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10) 
        push!(res2_zc, run_corrective_with_base(case, r))
    end
    
    res2_zvoll = []
    res = results[findfirst(x->x[1]=="NkCz2",results)][2]
    for (i,(r, s)) in enumerate(zip(res, scenarioes2))
        println("\nScenario prob zvoll ", i)
        @time case, tot_t = SCOPF.run_decomposition(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, 
            ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10) 
        push!(res2_zvoll, run_corrective_with_base(case, r))
    end

    println("\nN-1 prob")
    @time case, tot_t = SCOPF.run_decomposition(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, prob/8760, contingencies, max_shed=max_shed, reserve=reserve, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10)
    vals = SCOPF.extract_results(case)
    res2_n1 = run_corrective_with_base(vals, system, scenarioes2, voll, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long)

    res2_prev = []
    res = results[findfirst(x->x[1]=="NkP",results)][2]
    for (i,(r, s)) in enumerate(zip(res, scenarioes2))
        println("\nScen.2 ", i, " preventive")
        case, _ = SCOPF.run_decomposition(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10)
        push!(res2_prev, run_corrective_with_base(case, r))
    end
    
    res2_prev_n1 = []
    for (i,s) in enumerate(scenarioes2)
        println("\nScenario preventive ", i)
        @time case, tot_t = SCOPF.run_decomposition(SCOPF.P_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, reserve=reserve, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=prev_lim, long_term_multi=long, p_failure=0.00, max_itr=10) 
        vals = SCOPF.extract_results(case)
        push!(res2_prev_n1, run_corrective_with_base(case, vals))
    end
    Logging.disable_logging(Logging.Debug)

    return res2, res2_zc, res2_zvoll, res2_prev, res2_n1, res2_prev_n1
end

function calc_costs(results, scenarioes)
    base_costs = [n=>[isnan(x[:obj_val]) ? NaN : x[:costs][1,:Base] for x in r] for (n,r) in results]
    corrective_costs = [n=>[isnan(x[:obj_val]) ? NaN : (x[:costs][1,:r] + sum(p[2] .* (x[:costs][:,:Pc] + x[:costs][:,:Pcc]))) for (x,p) in zip(r,scenarioes)] for (n,r) in results]
    return base_costs, corrective_costs
end
   
function write_to_file(name, base_costs, corrective_costs)
    CSV.write("results/"*name*"_base_cost.csv", DataFrames.DataFrame(base_costs...))
    CSV.write("results/"*name*"_corrective_cost.csv", DataFrames.DataFrame(corrective_costs...)) 
end
function write_to_file(name, results)
    FileIO.save("results/"*name*".jld2", "results", results)
end

function plotsave_weather(name, base_costs, corrective_costs)
    plb = Plots.plot(last.(base_costs)./10, labels=reshape(getindex.(base_costs, 1), 1, :), 
        line=:steppre, xlabel="Time [h]", ylabel="Cost [k\$/h]", palette = Plots.palette(:seaborn_colorblind), thickness_scaling=1.2)
    Plots.savefig(plb, "results/"*name*"_base_cost.pdf")

    plc = Plots.plot(last.(corrective_costs)./10, labels=reshape(getindex.(base_costs, 1), 1, :), 
        line=:steppre, xlabel="Time [h]", ylabel="Cost [k\$/h]", palette = Plots.palette(:seaborn_colorblind), thickness_scaling=1.2)
    Plots.savefig(plc, "results/"*name*"_corrective_cost.pdf")

    plt = Plots.plot((last.(base_costs) + last.(corrective_costs))./10, labels=reshape(getindex.(base_costs, 1), 1, :), 
        line=:steppre, xlabel="Time [h]", ylabel="Cost [k\$/h]", palette = Plots.palette(:seaborn_colorblind), thickness_scaling=1.2)
    Plots.savefig(plt, "results/"*name*"_objective.pdf")
    return plb, plc, plt
end

function plotsave_weather2(name, base_costs, corrective_costs, base_costs2, corrective_costs2)
    plc = Plots.plot(last.(corrective_costs).*100, labels=reshape(getindex.(base_costs, 1), 1, :), line=(:steppre, :dash), alpha=0.5, xlabel="Time [h]", ylabel="Cost [\$/h]", palette = Plots.palette(:seaborn_colorblind)[1:5])
    Plots.plot!(plc, last.(corrective_costs2).*100, labels=reshape(getindex.(base_costs, 1), 1, :).*" S2", line=:steppre, xlabel="Time [h]", ylabel="Cost [\$/h]", palette = Plots.palette(:seaborn_colorblind)[1:5], leg=:topright, thickness_scaling=1.2)
    Plots.savefig(plc, "results/"*name*"_corrective_cost_s2.pdf")

    plt = Plots.plot((last.(base_costs) + last.(corrective_costs)).*100, labels=reshape(getindex.(base_costs, 1), 1, :), line=(:steppre, :dash), alpha=0.5, xlabel="Time [h]", ylabel="Cost [\$/h]", palette = Plots.palette(:seaborn_colorblind)[1:5])
    Plots.plot!(plt, (last.(base_costs2) + last.(corrective_costs2)).*100, labels=reshape(getindex.(base_costs, 1), 1, :).*" S2", line=:steppre, xlabel="Time [h]", ylabel="Cost [\$/h]", palette = Plots.palette(:seaborn_colorblind)[1:5], leg=:topright, thickness_scaling=1.2)
    Plots.savefig(plt, "results/"*name*"_objective_s2.pdf")
    return plc, plt
end

function plotsave_renewable(name, results, system)
    max_renewable = sum((get_active_power(g) for g in SCOPF.sort_components!(SCOPF.get_generation(system)) if typeof(g) <: RenewableGen), init=0.0)
    Plots.savefig(Plots.plot([[SCOPF.get_all_type_prod(system, r)[3,:] for (_,r) in results]... fill(max_renewable, length(results[1][2]))].*100, 
        label=[reshape(getindex.(results, 1), 1, :)... "Max renewable"], 
        line=:steppre, xlabel="Time [h]", ylabel="Generation [MWh]", palette = Plots.palette(:seaborn_colorblind), leg=:topright, thickness_scaling=1.2),
        "results/"*name*"_renewable.pdf")
end
   
function run_plotsave_weather(name, system, labels, scenarioes, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim)
    results = vec(labels) .=> run_cases(system, scenarioes, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim)
    base_costs, corrective_costs = calc_costs(results, scenarioes)
    write_to_file(name, results)
    write_to_file(name, base_costs, corrective_costs)
    plotsave_weather(name, base_costs, corrective_costs)
    plotsave_renewable(name, results, system)
    return results
end
function run_plotsave_weather(name, system, labels, scenarioes, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim, demands, demand_mult, renewables, renew_mult, t)
    set_parameters(demands, demand_mult, renewables, renew_mult, t)
    return run_plotsave_weather(name, system, labels, scenarioes, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim)
end
function run_plotsave_weather_zc(name, system, labels, scenarioes, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long)
    results = vec(labels) .=> run_zc_case(system, scenarioes, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long)
    base_costs, corrective_costs = calc_costs(results, scenarioes)
    write_to_file(name, results)
    write_to_file(name, base_costs, corrective_costs)
    plotsave_weather(name, base_costs, corrective_costs)
    plotsave_renewable(name, results, system)
    return results
end
function run_plotsave_weather_zc(name, system, labels, scenarioes, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, demands, demand_mult, renewables, renew_mult, t)
    set_parameters(demands, demand_mult, renewables, renew_mult, t)
    return run_plotsave_weather_zc(name, system, labels, scenarioes, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long)
end
function run_plotsave_weather2(results, name, system, labels, scenarioes2, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim)
    results = vec(labels) .=> run_cases_weather2(results, scenarioes2, system, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim)
    base_costs, corrective_costs = calc_costs(results, scenarioes2)
    write_to_file(name, results)
    write_to_file(name, base_costs, corrective_costs)
    plotsave_weather(name, base_costs, corrective_costs)
    return results
end

function plot_contingency_data(results::Vector, state::Symbol, property::Symbol, labels::Matrix{String}, forced_ls::Vector{<:Dict})
    data = [sort(reduce(vcat, [[sum(i[property]) - (get(forced_ls[t], c, 0) == 0 ? 0 : forced_ls[t][c]) for (c,i) in x[state]] for (t,x) in enumerate(r)]), rev=true) for r in results]
    xmax = maximum(filter(x->!isnothing(x), findfirst.(x->x<1e-3, data)); init=0)
    pl = Plots.plot(data.*100, labels=reshape(getindex.(results, 1), 1, :), leg=:topright, xaxis=[1,xmax], xlabel="Contingency", ylabel="Load shed [MWh/h]", thickness_scaling=1.2)
    return pl, data
end

function plot_probability_contingency_data(results::Vector, state::Symbol, property::Symbol, labels::Matrix{String}, forced_ls::Vector{<:Dict}, scenarioes)
    xax = reduce(vcat, last.(scenarioes))
    data = [[zeros(length(s[2])) for s in scenarioes] for _ in 1:length(results)]
    for (i,(_,r)) in enumerate(results), (j,t) in enumerate(scenarioes), (k,c) in enumerate(t[2])
        if get(r[j][state], k, 0) != 0
            if get(forced_ls[j], k, 0) != 0
                data[i][j][k] = sum(r[j][state][k][property]) - forced_ls[j][k]
            else
                data[i][j][k] = sum(r[j][state][k][property])
            end
        end
    end
    data = reduce.(vcat, data)
    # data = [reduce(vcat, [[sum(i[property]) - (get(forced_ls[t], c, 0) == 0 ? 0 : forced_ls[t][c]) for (c,i) in x[state]] for (t,x) in enumerate(r)]) for r in results]
    perm = sortperm.(data, rev=true)
    xax = cumsum.(getindex.([xax], perm))./sum(xax).*100
    xmax = maximum(filter(x->!isnothing(x), findfirst.(x->x<1e-3, data)); init=0)
    pl = Plots.plot(xax, getindex.(data, perm).*100, labels=reshape(getindex.(results, 1), 1, :), leg=:topright, xlabel="Total cumulative probability [%]", ylabel="Load shed [MWh/h]")
    return pl, xax, getindex.(data, perm).*100
end

function plot_contingency_data(results::Vector, state::Symbol, property::Symbol, labels::Matrix{String}, forced_ls::Vector{<:Dict}, scenarioes)
    data = [[sum(i[property]) - (get(forced_ls[t], c, 0) == 0 ? 0 : forced_ls[t][c]) for (c,i) in x[state]] for (t,x) in enumerate(results)]
    perm = sortperm.(data, rev=true)
    xax = reduce(vcat, last.(scenarioes))
    xax = cumsum.(getindex.([xax], perm))
    xmax = maximum(filter(x->!isnothing(x), findfirst.(x->x<1e-3, data)); init=0)
    pl = Plots.plot(xax, getindex.(data, perm).*100, leg=:none, xlabel="Contingency", ylabel="Load shed [MWh/h]", thickness_scaling=1.2)
    return pl, data
end

function table_results(results, system, scenarioes)
    base_costs, corrective_costs = calc_costs(results, scenarioes)
    forced_ls = [SCOPF.calc_forced_ls(case, x[2], x[1]) for x in scenarioes];
    _, xax_st, data_st = plot_probability_contingency_data(results, :Pc, :lsc, reshape(getindex.(results, 1), 1, :), forced_ls, scenarioes)
    _, xax_lt, data_lt = plot_probability_contingency_data(results, :Pcc, :lscc, reshape(getindex.(results, 1), 1, :), forced_ls, scenarioes)
    max_renewable = sum(get_active_power(g) for g in SCOPF.sort_components!(SCOPF.get_generation(system)) if typeof(g) <: RenewableGen)
    return DataFrames.DataFrame([
        sum.(last.(base_costs))/10, getindex.(last.(base_costs), 26)/10, 
        sum.(last.(corrective_costs))/10, getindex.(last.(corrective_costs), 26)/10,
        sum.(last.(base_costs) + last.(corrective_costs))/10, getindex.(last.(base_costs) + last.(corrective_costs), 26)/10,
        sum.([SCOPF.get_all_type_prod(system, r)[3,:] for (_,r) in results])./10, getindex.([SCOPF.get_all_type_prod(system, r)[3,:] for (_,r) in results], 26)./10,
        (sum(fill(max_renewable, size(scenarioes,1))) .- sum.([SCOPF.get_all_type_prod(system, r)[3,:] for (_,r) in results]))./10,
        (max_renewable .- getindex.([SCOPF.get_all_type_prod(system, r)[3,:] for (_,r) in results], 26))./10,
        findnext.(x->x<1e-6, data_st, 1) ./ sum(length.(getproperty.(scenarioes, 1))),
        findnext.(x->x<1e-6, data_lt, 1) ./ sum(length.(getproperty.(scenarioes, 1))),
        getindex.(xax_st, findnext.(x->x<1e-6, data_st, 1)),
        getindex.(xax_lt, findnext.(x->x<1e-6, data_lt, 1)),
        [sum(sum(scenarioes[j][2][i]*sum(x[:lsc]) for (i,x) in y[:Pc]) for (j,y) in enumerate(r)) for (_,r) in results].*100,
        [sum(sum(scenarioes[j][2][i]*sum(x[:lscc]) for (i,x) in y[:Pcc]) for (j,y) in enumerate(r)) for (_,r) in results].*100],
        ["Base cost [k\$]", "Base cost max [k\$]", "Post-c cost [k\$]", "Post-c cost max [k\$]", "Tot. cost [k\$]", "Tot. cost max [k\$]", 
        "Tot. VRE dispatched [GWh]", "Tot. VRE dispatched max [GWh]", "Tot. VRE curtail [GWh]", "Tot. VRE curtail max [GWh]",
        "Load shed st [%]", "Load shed lt [%]", "Prob. of ls st", "Prob. of ls lt", "Ls st [MW]", "Ls lt [MW]"])
end

function distributed_costs(opf::SCOPF.OPFsystem, scenario::Vector, x)
    return (
        sum(c.var * g for (c, g) in zip(opf.cost_gen, x[:pg0])),
        sum(opf.voll' * x[:ls0]),
        sum(c.ramp * g for (c, g) in zip(opf.cost_gen, x[:pgru])),
        sum(c.ramp * g for (c, g) in zip(opf.cost_gen, x[:pgrd])),
        # sum(scenario[i] * sum(c.ramp * g for (c, g) in zip(opf.cost_gen, c[:pgu])) for (i,c) in x[:Pc]),
        # sum(scenario[i] * sum(c.ramp * g for (c, g) in zip(opf.cost_gen, c[:pgd])) for (i,c) in x[:Pc]),
        maximum(opf.voll) * 10 * x[:pen][1],
        sum(scenario[i] * sum(opf.voll' * c[:lsc]) for (i,c) in x[:Pc]),
        sum(scenario[i] * sum(c.ramp * g for (c, g) in zip(opf.cost_gen, c[:pgu])) for (i,c) in x[:Pcc]),
        sum(scenario[i] * sum(c.ramp * g for (c, g) in zip(opf.cost_gen, c[:pgd])) for (i,c) in x[:Pcc]),
        sum(scenario[i] * sum(opf.voll' * c[:lscc]) for (i,c) in x[:Pcc])
    )
end

function distributed_costs(results::Vector, opf::SCOPF.OPFsystem, scenarioes::Vector, hours::Vector{Int64}=[x for x in 1:length(results[1][2])])
    costs = DataFrames.DataFrame(pg0=Float64[], ls0=Float64[], pgru=Float64[], pgrd=Float64[], pen=Float64[], lsc=Float64[], pgu=Float64[], pgd=Float64[], lscc=Float64[])
    for (r,(_,rval)) in enumerate(results)
        c = zeros(9)
        for (t,x) in enumerate(rval)
            if t ∈ hours
                c .+= distributed_costs(opf, scenarioes[t][2],x)
            end
        end
        push!(costs, c)
    end
    return costs
end