using PowerSystems
const PSY = PowerSystems
using JuMP
using Ipopt
using Printf
using PowerSystemCaseBuilder

function add_system_data_to_json(;
        file_name="system_data.json",
        data_name=joinpath("matpower","case5_re_uc_pwl.m"),
        time_series=joinpath("forecasts","5bus_ts","timeseries_pointers_da.json"),
        DATA_DIR="data"
    )
    # if !isdir(DATA_DIR)
    #     download(PowerSystems.UtilsData.TestData, folder = DATA_DIR)
    # end
    system_data = System(joinpath(DATA_DIR, data_name))
    add_time_series!(system_data, joinpath(DATA_DIR,time_series))
    to_json(system_data, file_name, force=true)
end

function p_scopf(system::System, optimizer; pen = 10000, time_limit_sec = 60)
    opf_m = Model(optimizer)
    set_time_limit_sec(opf_m, time_limit_sec)

    gen_names = get_name.(get_components(ThermalStandard, system))
    line_names = get_name.(get_components(Line, system))
    bus_names = get_name.(get_components(Bus, system))
    demand_names = get_name.(get_components(StaticLoad, system))
    contingencies = line_names

    # active power variables for the generators
    @variable(opf_m, pg0[g in gen_names] >= 0) 
    # power flow on lines in base case and contingencies
    @variable(opf_m, f0[l in line_names])
    @variable(opf_m, fc[l in line_names, c in contingencies])
    # voltage angle at a node in base case and contingencies
    @variable(opf_m, va0[b in bus_names])
    @variable(opf_m, vac[b in bus_names, c in contingencies])
    # load curtailment variables
    @variable(opf_m, ls0[b in demand_names] >= 0)

    # find value of beta, beta = 1 if from bus is the bus, beta = -1 if to bus is the bus, 0 else
    beta(bus, line) = bus == get_arc(line).from ? 1 : bus == get_arc(line).to ? -1 : 0

    # minimize cost of generation
    @objective(opf_m, Min, sum(
        # pg0[get_name(g)] * get_cost(get_variable(get_operation_cost(g)))[1][2]
        pg0[get_name(g)] * get_cost(get_variable(get_operation_cost(g)))[2]
        for g in get_components(ThermalGen, system)
        ) + pen * sum(ls0[d] for d in demand_names)
    )

    # incerted power at each bus for the base case and contingencies
    for b in get_components(Bus, system)
        @constraint(opf_m, sum(get_bus(g) == b ? pg0[get_name(g)] : 0 for g in get_components(ThermalStandard, system)) - 
                           sum(beta(b,l) * f0[get_name(l)] for l in get_components(Line, system)) ==
                           sum(get_active_power(l) - ls0[get_name(l)] for l in get_components(StaticLoad, system))
        )
        for c in contingencies
            @constraint(opf_m, sum(get_bus(g) == b ? pg0[get_name(g)] : 0 for g in get_components(ThermalStandard, system)) -
                                sum(beta(b,l) * fc[get_name(l),c] for l in get_components(Line, system)) ==
                                sum(get_active_power(l) - ls0[get_name(l)] for l in get_components(StaticLoad, system))
            )
        end
    end

    # power flow on line and line limits for the base case and contingencies
    for l in get_components(Line, system)
        l_name = get_name(l)
        @constraint(opf_m, f0[l_name] - sum(beta(b,l) * va0[get_name(b)] for b in get_components(Bus, system)) / get_x(l) == 0)
        @constraint(opf_m, -get_rate(l) <= f0[l_name] <= get_rate(l))
        for c in contingencies
            l_name = get_name(l)
            a = l != c # zero value if l is unavailable under contingency c
            @constraint(opf_m, fc[l_name,c] - a * sum(beta(b,l) * vac[get_name(b),c] for b in get_components(Bus, system)) / get_x(l) == 0)
            @constraint(opf_m, -get_rate(l) <= a * fc[l_name,c] <= get_rate(l))
        end
    end

    # restrict active power generation to min and max values
    for g in get_components(ThermalStandard, system)
        # set_lower_bound(pg[get_name(g)], get_active_power_limits(g).min) 
        set_upper_bound(pg0[get_name(g)], get_active_power_limits(g).max)
    end

    # restrict load shedding to max load per node
    for l in get_components(StaticLoad, system)
        set_upper_bound(ls0[get_name(l)], get_active_power(l))
    end

    optimize!(opf_m)
    # for g in get_components(ThermalStandard, system)
    #     @printf("%s = %.4f <= %.2f \n", get_name(g), value(pg0[get_name(g)]), get_active_power_limits(g).max)
    # end
    # for l in demand_names
    #     value(ls0[l]) > 0.00001 ? @printf("%s = %.4f \n", l,value(ls0[l])) : print()
    # end
    # for l in get_components(Line, system)
    #     @printf("%s = %.4f <= %.2f \n", get_name(l), value(f0[get_name(l)]), get_rate(l))
    # end

    return opf_m
end



function pc_scopf(system::System, optimizer; pen = 10000, time_limit_sec = 60)
    opf_m = Model(optimizer)
    set_time_limit_sec(opf_m, time_limit_sec)
    ramp_minutes = 10
    short_term_limit_multi = 1.5

    gen_names = get_name.(get_components(ThermalStandard, system))
    line_names = get_name.(get_components(Line, system))
    bus_names = get_name.(get_components(Bus, system))
    demand_names = get_name.(get_components(StaticLoad, system))
    contingencies = line_names # [bus_names; line_names]

    # active power variables for the generators in base case and contingencies ramp up and ramp down
    @variable(opf_m, pg0[g in gen_names] >= 0) 
    @variable(opf_m, pgu[g in gen_names, c in contingencies] >= 0) 
    @variable(opf_m, pgd[g in gen_names, c in contingencies] >= 0) 
    # power flow on lines in base case and in contingencies before and after corrective actions
    @variable(opf_m, f0[l in line_names])
    @variable(opf_m, fc[l in line_names, c in contingencies])
    @variable(opf_m, fcc[l in line_names, c in contingencies])
    # voltage angle at a node in base case and in contingencies before and after corrective actions
    @variable(opf_m, va0[b in bus_names])
    @variable(opf_m, vac[b in bus_names, c in contingencies])
    @variable(opf_m, vacc[b in bus_names, c in contingencies])
    # load curtailment variables in base case and in contingencies
    @variable(opf_m, ls0[d in demand_names] >= 0)
    @variable(opf_m, lsc[d in demand_names, c in contingencies] >= 0)

    prob = [ # spesified for the RTS-96
            # 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
            # 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, # generators
            0.24, 0.51, 0.33, 0.39, 0.48, 0.38, 0.02, 0.36, 0.34, 0.33, 0.30, 0.44, 0.44, 0.44, 
            0.02, 0.02, 0.02, 0.02, 0.40, 0.39, 0.40, 0.52, 0.49, 0.47, 0.38, 0.33, 0.41, 0.41, 
            0.41, 0.35, 0.34, 0.32, 0.54, 0.35, 0.35, 0.38, 0.38, 0.34, 0.34, 0.45, 0.46 # lines
    ]
    prob /= 8760

    # find value of beta, beta = 1 if from bus is the bus, beta = -1 if to bus is the bus, 0 else
    beta(bus, line) = bus == get_arc(line).from ? 1 : bus == get_arc(line).to ? -1 : 0

    # minimize cost of generation
    @objective(opf_m, Min, sum(
        # pg0[get_name(g)] * get_cost(get_variable(get_operation_cost(g)))[1][2]
        get_cost(get_variable(get_operation_cost(g)))[2] * (pg0[get_name(g)] + 
        sum(prob[i] * pgu[get_name(g),c] for (i,c) in enumerate(contingencies)))
        for g in get_components(ThermalGen, system)
        ) + 
        pen * sum(ls0[d] + sum(prob[i] * lsc[d,c] for (i,c) in enumerate(contingencies)) for d in demand_names) # removed "|CN1|", added "prob[c] * "
    )

    # incerted power at each bus for the base case and contingencies
    for b in get_components(Bus, system)
        @constraint(opf_m, sum(get_bus(g) == b ? pg0[get_name(g)] : 0 for g in get_components(ThermalStandard, system)) - 
                           sum(beta(b,l) * f0[get_name(l)] for l in get_components(Line, system)) ==
                           sum(get_active_power(l) - ls0[get_name(l)] for l in get_components(StaticLoad, system))
        )
        for c in contingencies
            @constraint(opf_m, sum(get_bus(g) == b ? pg0[get_name(g)] : 0 for g in get_components(ThermalStandard, system)) -
                                sum(beta(b,l) * fc[get_name(l),c] for l in get_components(Line, system)) ==
                                sum(get_active_power(l) - ls0[get_name(l)] for l in get_components(StaticLoad, system))
            )
            @constraint(opf_m, sum(get_bus(g) == b ? pg0[get_name(g)] + pgu[get_name(g),c] - pgd[get_name(g),c] : 0 
                                    for g in get_components(ThermalStandard, system)) -
                                sum(beta(b,l) * fcc[get_name(l),c] for l in get_components(Line, system)) ==
                                sum(get_active_power(l) - ls0[get_name(l)] - lsc[get_name(l),c] 
                                    for l in get_components(StaticLoad, system))
            )
        end
    end

    # power flow on line and line limits for the base case and contingencies
    for l in get_components(Line, system)
        l_name = get_name(l)
        @constraint(opf_m, f0[l_name] - sum(beta(b,l) * va0[get_name(b)] for b in get_components(Bus, system)) / get_x(l) == 0)
        @constraint(opf_m, -get_rate(l) <= f0[l_name] <= get_rate(l))
        for c in contingencies
            l_name = get_name(l)
            a = l != c # zero value if l is unavailable under contingency c
            @constraint(opf_m, fc[l_name,c] - a * sum(beta(b,l) * vac[get_name(b),c] for b in get_components(Bus, system)) / get_x(l) == 0)
            @constraint(opf_m, -get_rate(l)*short_term_limit_multi <= a * fc[l_name,c] <= get_rate(l)*short_term_limit_multi)
            @constraint(opf_m, fcc[l_name,c] - a * sum(beta(b,l) * vacc[get_name(b),c] for b in get_components(Bus, system)) / get_x(l) == 0)
            @constraint(opf_m, -get_rate(l) <= a * fcc[l_name,c] <= get_rate(l))
        end
    end

    # restrict active power generation to min and max values
    for g in get_components(ThermalStandard, system)
        g_name = get_name(g)
        # set_lower_bound(pg[get_name(g)], get_active_power_limits(g).min) 
        set_upper_bound(pg0[g_name], get_active_power_limits(g).max)
        for c in contingencies
            @constraint(opf_m, 0 <= pg0[g_name] + (pgu[g_name,c] - pgd[g_name,c]) <= get_active_power_limits(g).max)
            set_upper_bound(pgu[g_name,c], get_ramp_limits(g).up * ramp_minutes)
            set_upper_bound(pgd[g_name,c], get_ramp_limits(g).down * ramp_minutes)
        end
    end

    # restrict load shedding to load per node
    for l in get_components(StaticLoad, system)
        for c in contingencies
            @constraint(opf_m, ls0[get_name(l)] + lsc[get_name(l),c] <= get_active_power(l))
        end
    end

    optimize!(opf_m)
    return opf_m
end

function e_opf(system::System, optimizer, outage, pg0, ls0; time_limit_sec = 60)
    opf_m = Model(optimizer)
    set_time_limit_sec(opf_m, time_limit_sec)

    gen_names = get_name.(get_components(ThermalStandard, system))
    line_names = get_name.(get_components(Line, system))
    bus_names = get_name.(get_components(Bus, system))
    demand_names = get_name.(get_components(StaticLoad, system))

    # power flow on lines in base case and emergency
    @variable(opf_m, f0[l in line_names])
    @variable(opf_m, fe[l in line_names])
    # voltage angle at a node in emergency
    @variable(opf_m, vae[b in bus_names])
    # load curtailment variables in emergency
    @variable(opf_m, lse[d in demand_names] >= 0)

    # find value of beta, beta = 1 if from bus is the bus, beta = -1 if to bus is the bus, 0 else
    beta(bus, line) = bus == get_arc(line).from ? 1 : bus == get_arc(line).to ? -1 : 0

    # generation partition factors
    @expression(opf_m, e[gen], pgo[gen] / sum(pg0[g] for g in gen_names))

    # minimize cost of generation
    @objective(opf_m, Min, sum(lse))

    # incerted power at each bus for the base case
    for b in get_components(Bus, system)
        @constraint(opf_m, sum(get_bus(g) == b ? pg0[get_name(g)] - e[g] * sum(lse) : 0 for g in get_components(ThermalStandard, system)) - 
                           sum(beta(b,l) * f0[get_name(l)] for l in get_components(Line, system)) ==
                           sum(get_active_power(l) - ls0[get_name(l)] - lse[get_name(l)] for l in get_components(StaticLoad, system))
        )
    end

    # power flow on line and line limits for the base case
    for l in get_components(Line, system)
        l_name = get_name(l)
        @constraint(opf_m, fe[l_name] - outage[l_name] * sum(beta(b,l) * vae[get_name(b)] for b in get_components(Bus, system)) / get_x(l) == 0)
        @constraint(opf_m, -get_rate(l) <= fe[l_name] <= get_rate(l))
    end

    optimize!(opf_m)
    for d in demand_names
        value(lse[d]) > 0.00001 ? @printf("%s = %.4f \n", d,value(lse[d])) : print()
    end
    for l in get_components(Line, system)
        @printf("%s = %.4f <= %.2f \n", get_name(l), value(f0[get_name(l)]), get_rate(l))
    end

    return opf_m
end

# # corrective control failure probability
# phi(p, n) = sum((-1)^k * p^k * binomial(n,k) for k in 1:n)

# # socio economic cost function
# @expression(opf_m, sec, 
#     sum(cg[g] * (pg0[g] + sum(prob[c] * (pgu[g,c] - pgd[g,c]) for c in contingencies)) for g in gen_names) +
#     sum(voll[d] * (ls[d] + sum(prob[c] * lsc[d,c] for c in contingencies) for d in demand_names))
# )

# # severity function
# @expression(opf_m, severity, sum(voll[d] * lse[d] for d in demand_names))


# for l in get_name.(get_components(Line, ieee_rts))
#     for g in get_name.(get_components(ThermalStandard, ieee_rts))
#         if value.(results[:pgu])[g,l] > 0.00001 
#             @printf("%s: \t%s \t= %.5f \n", l, g, value.(results[:pgu])[g,l])
#         end
#     end
# end
        
# add_system_data_to_json()
# system_data = System("system_data.json")
# results = opf_model(system_data, Ipopt.Optimizer)
# value.(results[:pl])

function hot_start(model)
    x = all_variables(model)
    x_solution = value.(x)
    set_start_value.(x, x_solution)
    optimize!(model)
end