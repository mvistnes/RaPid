using PowerSystems
const PSY = PowerSystems
using JuMP
using Ipopt

function sl_scopf(system::System, optimizer)
    ramp_minutes = 10
    short_term_limit_multi = 1.5
    gen_names       = get_name.(get_components(Generator, system))
    branch_names    = get_name.(get_components(ACBranch, system))
    bus_names       = get_name.(get_components(Bus, system))
    demand_names    = get_name.(get_components(StaticLoad, system))
    contingencies = branch_names # [bus_names;branch_names]
    voll = JuMP.Containers.DenseAxisArray(
        [9144, 6496, 13884, 7185, 8574, 11703, 9062, 5013, 13498, 14698, 7812, 10102, 
        10402, 10691, 8190, 9384, 4875, 3992, 6368, 8203, 4854, 3062, 11128, 10774], 
        demand_names)
    
    p_opf_m = p_scopf(system, optimizer, gen_names, branch_names, bus_names, demand_names, voll, contingencies, short_term_limit_multi)
    optimize!(p_opf_m)
    pc_opf_m = pc_scopf(system, optimizer, p_opf_m, gen_names, branch_names, bus_names, demand_names, voll, contingencies, ramp_minutes)
    optimize!(pc_opf_m)
    return p_opf_m, pc_opf_m
end

function p_scopf(system, optimizer, gen_names, branch_names, bus_names, demand_names, voll, contingencies, short_term_limit_multi)
    opf_m = Model(optimizer)
    if GLPK.Optimizer == optimizer 
        set_optimizer_attribute(opf_m, "msg_lev", GLPK.GLP_MSG_ON)
    end
    
    @variables(opf_m, begin
        pg0[g in gen_names] >= 0     # active power variables for the generators
        pf0[l in branch_names]       # power flow on branches in base case
        pfc[l in branch_names, c in contingencies] # and contingencies
        va0[b in bus_names]          # voltage angle at a node in base case
        vac[b in bus_names, c in contingencies] # and contingencies
        ls0[b in demand_names] >= 0  # load curtailment variables
    end)

    # find value of beta, beta = 1 if from-bus is the bus, beta = -1 if to-bus is the bus, 0 else
    beta(bus, branch) = bus == get_arc(branch).from ? 1 : bus == get_arc(branch).to ? -1 : 0

    # incerted power at each bus for the base case
    @expression(opf_m, generators[b = get_components(Bus, system)], 
        sum(get_bus(g) == b ? pg0[get_name(g)] : 0 for g in get_components(Generator, system)))
    for b in get_components(Bus, system)
        demand = sum(get_bus(d) == b ? get_active_power(d) - ls0[get_name(d)] : 0 for d in get_components(StaticLoad, system))
        @constraint(opf_m, generators[b] - sum(beta(b,l) * pf0[get_name(l)] for l in get_components(ACBranch, system)) == demand)
        for c in contingencies
            @constraint(opf_m, generators[b] - sum(beta(b,l) * pfc[get_name(l),c] for l in get_components(ACBranch, system)) == demand)
        end
    end

    # power flow on branch and branch limits for the base case
    for l in get_components(ACBranch, system)
        l_name = get_name(l)
        @constraint(opf_m, pf0[l_name] - sum(beta(b,l) * va0[get_name(b)] for b in get_components(Bus, system)) / get_x(l) == 0)
        @constraint(opf_m, -get_rate(l) <= pf0[l_name] <= get_rate(l))
        for c in contingencies
            a = l_name != c # zero value if l is unavailable under contingency c
            @constraint(opf_m, pfc[l_name,c] - a * sum(beta(b,l) * vac[get_name(b),c] for b in get_components(Bus, system)) / get_x(l) == 0)
            @constraint(opf_m, -get_rate(l)*short_term_limit_multi <= a * pfc[l_name,c] <= get_rate(l)*short_term_limit_multi)
        end
    end

    # restrict active power generation to min and max values
    for g in get_components(Generator, system)
        # set_lower_bound(pg0[get_name(g)], get_active_power_limits(g).min) 
        set_upper_bound(pg0[get_name(g)], get_active_power_limits(g).max)
    end

    # restrict load shedding to max load per node
    for l in get_components(StaticLoad, system)
        set_upper_bound(ls0[get_name(l)], get_active_power(l))
    end

    # minimize socio-economic cost
    @objective(opf_m, Min, 
        sum(pg0[get_name(g)] * get_cost(get_variable(get_operation_cost(g)))[2]
        for g in get_components(ThermalStandard, system)) + 
        sum(pg0[get_name(g)] * 5 for g in get_components(HydroDispatch, system)) + 
        sum(voll[d] * ls0[d] for d in demand_names)
    )

    return opf_m
end


function pc_scopf(system::System, optimizer, p_opf_m, gen_names, branch_names, bus_names, demand_names, voll, contingencies, ramp_minutes)
    opf_m = Model(optimizer)
    if GLPK.Optimizer == optimizer 
        set_optimizer_attribute(opf_m, "msg_lev", GLPK.GLP_MSG_ON)
    end
    
    @variables(opf_m, begin
        pgu[g in gen_names, c in contingencies] >= 0    # active power variables for the generators in contingencies ramp up 
        pgd[g in gen_names, c in contingencies] >= 0        # and ramp down
        pfcc[l in branch_names, c in contingencies]     # power flow on branches in in contingencies after corrective actions
        vacc[b in bus_names, c in contingencies]        # voltage angle at a node in in contingencies after corrective actions
        lsc[d in demand_names, c in contingencies] >= 0 # load curtailment variables in in contingencies
    end)

    prob = [ # spesified for the RTS-96
            # 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
            # 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, # generators
            0.24, 0.51, 0.33, 0.39, 0.48, 0.38, 0.02, 0.36, 0.34, 0.33, 0.30, 0.44, 0.44, 0.44, 
            0.02, 0.02, 0.02, 0.02, 0.40, 0.39, 0.40, 0.52, 0.49, 0.47, 0.38, 0.33, 0.41, 0.41, 
            0.41, 0.35, 0.34, 0.32, 0.54, 0.35, 0.35, 0.38, 0.38, 0.34, 0.34, 0.45, 0.46 # branches
    ]
    prob /= 8760

    # find value of beta, beta = 1 if from bus is the bus, beta = -1 if to bus is the bus, 0 else
    beta(bus, branch) = bus == get_arc(branch).from ? 1 : bus == get_arc(branch).to ? -1 : 0

    # minimize socio-economic cost
    @objective(opf_m, Min, 
        sum(get_cost(get_variable(get_operation_cost(g)))[2] * ( 
            sum(prob[i] * (pgu[get_name(g),c] #=+ pgd[get_name(g),c]*0.1=#) for (i,c) in enumerate(contingencies)))
            for g in get_components(ThermalStandard, system)
        ) + 
        sum(5 * (sum(prob[i] * (pgu[get_name(g),c] #=+ pgd[get_name(g),c]*0.1=#)
            for (i,c) in enumerate(contingencies))) for g in get_components(HydroDispatch, system)
        ) +
        sum(voll[d] * (sum(prob[i] * lsc[d,c] for (i,c) in enumerate(contingencies))) for d in demand_names)
    )

    # incerted power at each bus for the base case and contingencies
    for b in get_components(Bus, system)
        for c in contingencies
            @constraint(opf_m, 
                sum(get_bus(g) == b ? value.(p_opf_m[:pg0][get_name(g)]) + pgu[get_name(g),c] - pgd[get_name(g),c] : 0 for g in get_components(Generator, system)) -
                sum(beta(b,l) * pfcc[get_name(l),c] for l in get_components(ACBranch, system)) == 
                sum(get_bus(d) == b ? get_active_power(d) - value.(p_opf_m[:ls0][get_name(d)]) - lsc[get_name(d),c] : 0 for d in get_components(StaticLoad, system))
            )
        end
    end
    
    # power flow on branch and branch limits for the base case and contingencies
    for l in get_components(ACBranch, system)
        l_name = get_name(l)
        for c in contingencies
            l_name = get_name(l)
            a = l_name != c # zero value if l is unavailable under contingency c
            @constraints(opf_m, begin
                pfcc[l_name,c] - a * sum(beta(b,l) * vacc[get_name(b),c] for b in get_components(Bus, system)) / get_x(l) == 0
                -get_rate(l) <= a * pfcc[l_name,c] <= get_rate(l)
            end)
        end
    end

    # restrict active power generation to min and max values
    for g in get_components(Generator, system)
        g_name = get_name(g)
        for c in contingencies
            @constraint(opf_m, 0 <= value.(p_opf_m[:pg0][g_name]) + (pgu[g_name,c] - pgd[g_name,c]) <= get_active_power_limits(g).max)
            set_upper_bound(pgu[g_name,c], get_ramp_limits(g).up * ramp_minutes)
            set_upper_bound(pgd[g_name,c], get_ramp_limits(g).down * ramp_minutes)
        end
    end

    # restrict load shedding to load per node
    for l in get_components(StaticLoad, system), c in contingencies
        @constraint(opf_m, value.(p_opf_m[:ls0][get_name(l)]) + lsc[get_name(l),c] <= get_active_power(l))
    end

    return opf_m
end

        
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