using PowerSystems
const PSY = PowerSystems
using JuMP
using Ipopt

function sl_scopf(system::System, optimizer)
    ramp_minutes = 10
    short_term_limit_multi = 1.5
    gens_t     = get_components(ThermalGen, system)
    gens_h     = get_components(HydroGen, system)
    branches   = get_components(ACBranch, system)
    buses      = get_components(Bus, system)
    demands    = get_components(StaticLoad, system)
    renewables = get_components(RenewableGen, system) # Renewable modelled as negative demand
    voll = JuMP.Containers.DenseAxisArray(
        [rand(1000:3000, length(demands)); rand(1:30, length(renewables))], 
        [get_name.(demands); get_name.(renewables)]
    )
    contingencies = get_name.(branches)[1:10] # [get_name.(buses);get_name.(branches)]
    
    p_opf_m = p_scopf(system, optimizer, gens_t, gens_h, branches, buses, demands, renewables, voll, contingencies, short_term_limit_multi)
    optimize!(p_opf_m)
    pc_opf_m = pc_scopf(system, optimizer, p_opf_m, gens_t, gens_h, branches, buses, demands, renewables, voll, contingencies, ramp_minutes)
    optimize!(pc_opf_m)
    return p_opf_m, pc_opf_m
end

function p_scopf(system, optimizer, gens_t, gens_h, branches, buses, demands, renewables, voll, contingencies, short_term_limit_multi)
    opf_m = Model(optimizer)
    if GLPK.Optimizer == optimizer 
        set_optimizer_attribute(opf_m, "msg_lev", GLPK.GLP_MSG_ON)
    end
    
    @variables(opf_m, begin
        pg0[g in [get_name.(gens_t); get_name.(gens_h)]] >= 0     # active power variables for the generators
        pf0[l in get_name.(branches)]       # power flow on branches in base case
        pfc[l in get_name.(branches), c in contingencies] # and contingencies
        va0[b in get_name.(buses)]          # voltage angle at a node in base case
        vac[b in get_name.(buses), c in contingencies] # and contingencies
        ls0[b in [get_name.(demands); get_name.(renewables)]] >= 0  # load curtailment variables
    end)

    # find value of beta, beta = 1 if from-bus is the bus, beta = -1 if to-bus is the bus, 0 else
    beta(bus, branch) = bus == get_arc(branch).from ? 1 : bus == get_arc(branch).to ? -1 : 0

    # incerted power at each bus for the base case and contingencies
    @expression(opf_m, generators[b = buses], 
        sum(get_bus(g) == b ? pg0[get_name(g)] : 0 for g in Iterators.flatten((gens_t, gens_h))))
    for b in buses
        demand = sum(get_bus(d) == b ? get_active_power(d) - ls0[get_name(d)] : 0 for d in demands)
        renewable = sum((get_bus(d) == b ? -get_active_power(d) + ls0[get_name(d)] : 0 for d in renewables), init = 0.0)
        c = @constraint(opf_m, generators[b] - sum(beta(b,l) * pf0[get_name(l)] for l in branches) == demand + renewable)
        set_name(c, @sprintf("inj_p[%s]", get_name(b)))
        for c in contingencies
            cc = @constraint(opf_m, opf_m[:generators][b] - sum(beta(b,l) * pfc[get_name(l),c] for l in branches) == demand + renewable)
            set_name(cc, @sprintf("inj_pc[%s,%s]", get_name(b), c))
        end
    end

    # power flow on branch and branch limits for the base case and contingencies
    for l in branches
        l_name = get_name(l)
        @constraint(opf_m, pf0[l_name] - sum(beta(b,l) * va0[get_name(b)] for b in buses) / get_x(l) == 0)
        @constraint(opf_m, -get_rate(l) <= pf0[l_name] <= get_rate(l))
        for c in contingencies
            a = l_name != c # zero value if l is unavailable under contingency c
            @constraint(opf_m, pfc[l_name,c] - a * sum(beta(b,l) * vac[get_name(b),c] for b in buses) / get_x(l) == 0)
            @constraint(opf_m, -get_rate(l) <= a * pfc[l_name,c] <= get_rate(l))
        end
    end

    # restrict active power generation to min and max values
    for g in Iterators.flatten((gens_t, gens_h))
        # set_lower_bound(pg0[get_name(g)], get_active_power_limits(g).min) 
        set_upper_bound(pg0[get_name(g)], get_active_power_limits(g).max)
    end

    # restrict load shedding to max load per node
    for d in Iterators.flatten((demands, renewables))
        set_upper_bound(ls0[get_name(d)], get_active_power(d))
    end

    # minimize socio-economic cost
    @objective(opf_m, Min, 
        sum(pg0[get_name(g)] * get_cost(get_variable(get_operation_cost(g)))[2]
        for g in gens_t) + 
        sum(pg0[get_name(g)] * 5 for g in gens_h) + 
        sum(voll[d] * ls0[d] for d in [get_name.(demands); get_name.(renewables)])
    )

    return opf_m
end


function pc_scopf(system::System, optimizer, p_opf_m, gens_t, gens_h, branches, buses, demands, renewables, voll, contingencies, ramp_minutes)
    opf_m = Model(optimizer)
    if GLPK.Optimizer == optimizer 
        set_optimizer_attribute(opf_m, "msg_lev", GLPK.GLP_MSG_ON)
    end
    
    @variables(opf_m, begin
        pgu[g in [get_name.(gens_t); get_name.(gens_h)], c in contingencies] >= 0    # active power variables for the generators in contingencies ramp up 
        pgd[g in [get_name.(gens_t); get_name.(gens_h)], c in contingencies] >= 0       # and ramp down
        pfcc[l in get_name.(branches), c in contingencies]         # and after corrective actions
        vacc[b in get_name.(buses), c in contingencies]            # and after corrective actions
        lsc[d in [get_name.(demands); get_name.(renewables)], c in contingencies] >= 0 # load curtailment variables in in contingencies
    end)

    prob = JuMP.Containers.DenseAxisArray((rand(length(contingencies)).*(0.5-0.02).+0.02)./8760, contingencies)
    # prob = [ # spesified for the RTS-96
    #         # 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
    #         # 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, # generators
    #         0.24, 0.51, 0.33, 0.39, 0.48, 0.38, 0.02, 0.36, 0.34, 0.33, 0.30, 0.44, 0.44, 0.44, 
    #         0.02, 0.02, 0.02, 0.02, 0.40, 0.39, 0.40, 0.52, 0.49, 0.47, 0.38, 0.33, 0.41, 0.41, 
    #         0.41, 0.35, 0.34, 0.32, 0.54, 0.35, 0.35, 0.38, 0.38, 0.34, 0.34, 0.45, 0.46 # branches
    # ]
    # prob /= 8760

    # find value of beta, beta = 1 if from bus is the bus, beta = -1 if to bus is the bus, 0 else
    beta(bus, branch) = bus == get_arc(branch).from ? 1 : bus == get_arc(branch).to ? -1 : 0

    # minimize socio-economic cost
    @objective(opf_m, Min, 
        sum(get_cost(get_variable(get_operation_cost(g)))[2] * (value.(p_opf_m[:pg0][get_name(g)]) + 
            sum(prob[c] * (pgu[get_name(g),c] #=+ pgd[get_name(g),c]*0.1=#) for c in contingencies))
            for g in gens_t
        ) + 
        sum(5 * (value.(p_opf_m[:pg0][get_name(g)]) + sum(prob[c] * (pgu[get_name(g),c] #=+ pgd[get_name(g),c]*0.1=#)
            for c in contingencies)) for g in gens_h
        ) +
        sum(voll[d] * (value.(p_opf_m[:ls0][d]) + sum(prob[c] * lsc[d,c] for c in contingencies)) 
            for d in [get_name.(demands); get_name.(renewables)])
    )

    # incerted power at each bus for the base case and contingencies
    for b in buses
        for c in contingencies
            ccc = @constraint(opf_m, 
                sum(get_bus(g) == b ? value.(p_opf_m[:pg0][get_name(g)]) + pgu[get_name(g),c] - pgd[get_name(g),c] : 0 for g in Iterators.flatten((gens_t, gens_h))) -
                sum(beta(b,l) * pfcc[get_name(l),c] for l in branches) == 
                sum(get_bus(d) == b ? get_active_power(d) - lsc[get_name(d),c] : 0 for d in demands) +
                sum(get_bus(d) == b ? -get_active_power(d) + lsc[get_name(d),c] : 0 for d in renewables)
            )
            set_name(ccc, @sprintf("inj_pcc[%s,%s]", get_name(b), c))
        end
    end
    
    # power flow on branch and branch limits for the base case and contingencies
    for l in branches
        l_name = get_name(l)
        for c in contingencies
            a = l_name != c # zero value if l is unavailable under contingency c
            @constraints(opf_m, begin
                pfcc[l_name,c] - a * sum(beta(b,l) * vacc[get_name(b),c] for b in buses) / get_x(l) == 0
                -get_rate(l) <= a * pfcc[l_name,c] <= get_rate(l)
            end)
        end
    end

    # restrict active power generation to min and max values
    for g in Iterators.flatten((gens_t, gens_h))
        g_name = get_name(g)
        for c in contingencies
            @constraint(opf_m, 0 <= value.(p_opf_m[:pg0][g_name]) + (pgu[g_name,c] - pgd[g_name,c]) <= get_active_power_limits(g).max)
            set_upper_bound(pgu[g_name,c], get_ramp_limits(g).up * ramp_minutes)
            set_upper_bound(pgd[g_name,c], get_ramp_limits(g).down * ramp_minutes)
        end
    end

    # restrict load shedding to load per node
    for l in Iterators.flatten((demands, renewables)), c in contingencies
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