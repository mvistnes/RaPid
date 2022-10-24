using PowerSystems
const PSY = PowerSystems
using JuMP
using Ipopt
using Printf
include("N-1_SCOPF.jl")

function sl_scopf(system::System, optimizer; voll=nothing, contingencies = nothing, prob=nothing)
    ramp_minutes = 10
    short_term_limit_multi = 1.5
    voll = isnothing(voll) ? make_voll(system) : voll
    contingencies = isnothing(contingencies) ? get_name.(branches(system)) : contingencies
    prob = isnothing(prob) ? make_prob(contingencies) : prob
    
    p_opf_m = sl_p_scopf(system, optimizer, contingencies, short_term_limit_multi, voll)
    solve_model!(p_opf_m)
    @assert termination_status(p_opf_m) == MOI.OPTIMAL
    pc_opf_m = sl_pc_scopf(system, optimizer, p_opf_m, contingencies, ramp_minutes, voll, prob)
    solve_model!(pc_opf_m)
    @assert termination_status(pc_opf_m) == MOI.OPTIMAL
    return p_opf_m, pc_opf_m
end

function sl_p_scopf(system::System, optimizer, contingencies, short_term_limit_multi, voll)
    opf_m = Model(optimizer)
    if GLPK.Optimizer == optimizer 
        set_optimizer_attribute(opf_m, "msg_lev", GLPK.GLP_MSG_ON)
    end
    
    @variables(opf_m, begin
        pg0[g in get_name.(get_ctrl_generation(system))] >= 0     # active power variables for the generators
        pf0[l in get_name.(branches(system))]       # power flow on branches in base case
        pfc[l in get_name.(branches(system)), c in contingencies] # and contingencies
        va0[b in get_name.(nodes(system))]          # voltage angle at a node in base case
        vac[b in get_name.(nodes(system)), c in contingencies] # and contingencies
        ls0[b in get_name.(get_nonctrl_generation(system))] >= 0  # load curtailment variables
    end)

    # find value of beta, beta = 1 if from-bus is the bus, beta = -1 if to-bus is the bus, 0 else
    beta(bus, branch) = bus == get_arc(branch).from ? 1 : bus == get_arc(branch).to ? -1 : 0

    # incerted power at each bus for the base case and contingencies
    @expression(opf_m, generators[b = nodes(system)], 
        sum(get_bus(g) == b ? pg0[get_name(g)] : 0 for g in get_ctrl_generation(system)))
    for b in nodes(system)
        demand = sum(get_bus(d) == b ? get_active_power(d) - ls0[get_name(d)] : 0 for d in demands(system))
        renewable = sum((get_bus(d) == b ? -get_active_power(d) + ls0[get_name(d)] : 0 for d in renewables(system)), init = 0.0)
        @constraint(opf_m, generators[b] - sum(beta(b,l) * pf0[get_name(l)] for l in branches(system)) == demand + renewable)
        # @constraint(opf_m, -π/2 <= va0[get_name(b)] <= π/2) # Not really needed, could be implemented with spesific line angle limits
        for c in contingencies
            @constraint(opf_m, opf_m[:generators][b] - sum(beta(b,l) * pfc[get_name(l),c] for l in branches(system)) == demand + renewable)
            # @constraint(opf_m, -π/2 <= vac[get_name(b),c] <= π/2) # Not really needed, could be implemented with spesific line angle limits
        end
    end
    
    slack = nothing
    for x in nodes(system)
        if x.bustype == BusTypes.REF 
            slack = x
            break
        end
    end
    @constraint(opf_m, va0[get_name(slack)] == 0)
    @constraint(opf_m, [c = contingencies], vac[get_name(slack),c] == 0)

    # power flow on branch and branch limits for the base case and contingencies
    for l in branches(system)
        l_name = get_name(l)
        @constraint(opf_m, pf0[l_name] - sum(beta(b,l) * va0[get_name(b)] for b in nodes(system)) / get_x(l) == 0)
        @constraint(opf_m, -get_rate(l) <= pf0[l_name] <= get_rate(l))
        for c in contingencies
            a = l_name != c # zero value if l is unavailable under contingency c
            @constraint(opf_m, pfc[l_name,c] - a * sum(beta(b,l) * vac[get_name(b),c] for b in nodes(system)) / get_x(l) == 0)
            @constraint(opf_m, a * -get_rate(l) * short_term_limit_multi <= pfc[l_name,c] <= a * get_rate(l) * short_term_limit_multi)
        end
    end

    # restrict active power generation to min and max values
    for g in get_ctrl_generation(system)
        # set_lower_bound(pg0[get_name(g)], get_active_power_limits(g).min) 
        set_upper_bound(pg0[get_name(g)], get_active_power_limits(g).max)
    end

    # restrict load shedding to max load per node
    for d in get_nonctrl_generation(system)
        set_upper_bound(ls0[get_name(d)], get_active_power(d))
    end

    # minimize socio-economic cost
    @objective(opf_m, Min, 
        sum(pg0[get_name(g)] * get_cost(get_variable(get_operation_cost(g)))[2]
        for g in gens_t(system)) + 
        sum(pg0[g] * 5 for g in get_name.(gens_h(system))) + 
        sum(voll[d] * ls0[d] for d in get_name.(get_nonctrl_generation(system)))
    )

    return opf_m
end


function sl_pc_scopf(system::System, optimizer, p_opf_m, contingencies, ramp_minutes, voll, prob)
    opf_m = Model(optimizer)
    if GLPK.Optimizer == optimizer 
        set_optimizer_attribute(opf_m, "msg_lev", GLPK.GLP_MSG_ON)
    end
    
    @variables(opf_m, begin
        pgu[g in get_name.(get_ctrl_generation(system)), c in contingencies] >= 0    # active power variables for the generators in contingencies ramp up 
        pgd[g in get_name.(get_ctrl_generation(system)), c in contingencies] >= 0       # and ramp down
        pfcc[l in get_name.(branches(system)), c in contingencies]         # and after corrective actions
        vacc[b in get_name.(nodes(system)), c in contingencies]            # and after corrective actions
        lsc[d in get_name.(get_nonctrl_generation(system)), c in contingencies] >= 0 # load curtailment variables in in contingencies
    end)

    # find value of beta, beta = 1 if from bus is the bus, beta = -1 if to bus is the bus, 0 else
    beta(bus, branch) = bus == get_arc(branch).from ? 1 : bus == get_arc(branch).to ? -1 : 0

    # minimize socio-economic cost
    @objective(opf_m, Min, 
        sum(get_cost(get_variable(get_operation_cost(g)))[2] * (value.(p_opf_m[:pg0][get_name(g)]) + 
            sum(prob[c] * (pgu[get_name(g),c] #=+ pgd[get_name(g),c]*0.1=#) for c in contingencies))
            for g in gens_t(system)
        ) + 
        sum(5 * (value.(p_opf_m[:pg0][get_name(g)]) + sum(prob[c] * (pgu[get_name(g),c] #=+ pgd[get_name(g),c]*0.1=#)
            for c in contingencies)) for g in gens_h(system)
        ) +
        sum(voll[d] * (value.(p_opf_m[:ls0][d]) + sum(prob[c] * lsc[d,c] for c in contingencies)) 
            for d in get_name.(get_nonctrl_generation(system))
        )
    )

    # incerted power at each bus for the base case and contingencies
    for b in nodes(system)
        for c in contingencies
            @constraint(opf_m, 
                sum(get_bus(g) == b ? value(p_opf_m[:pg0][get_name(g)]) + pgu[get_name(g),c] - pgd[get_name(g),c] : 0 for g in get_ctrl_generation(system)) -
                sum(beta(b,l) * pfcc[get_name(l),c] for l in branches(system)) == 
                sum(get_bus(d) == b ? get_active_power(d) - lsc[get_name(d),c] : 0 for d in demands(system)) +
                sum(get_bus(d) == b ? -get_active_power(d) + lsc[get_name(d),c] : 0 for d in renewables(system))
            )
            # @constraint(opf_m, -π/2 <= vacc[get_name(b),c] <= π/2) # Not really needed, could be implemented with spesific line angle limits
        end
    end
    
    slack = nothing
    for x in nodes(system)
        if x.bustype == BusTypes.REF 
            slack = x
            break
        end
    end
    @constraint(opf_m, [c = contingencies], vacc[get_name(slack),c] == 0)
    
    # power flow on branch and branch limits for the base case and contingencies
    for l in branches(system)
        l_name = get_name(l)
        for c in contingencies
            a = l_name != c # zero value if l is unavailable under contingency c
            @constraints(opf_m, begin
                pfcc[l_name,c] - a * sum(beta(b,l) * vacc[get_name(b),c] for b in nodes(system)) / get_x(l) == 0
                a * -get_rate(l) <= pfcc[l_name,c] <= a * get_rate(l)
            end)
        end
    end

    # restrict active power generation to min and max values
    for g in get_ctrl_generation(system)
        g_name = get_name(g)
        for c in contingencies
            @constraint(opf_m, 0 <= value.(p_opf_m[:pg0][g_name]) + (pgu[g_name,c] - pgd[g_name,c]) <= get_active_power_limits(g).max)
            set_upper_bound(pgu[g_name,c], get_ramp_limits(g).up * ramp_minutes)
            set_upper_bound(pgd[g_name,c], get_ramp_limits(g).down * ramp_minutes)
        end
    end

    # restrict load shedding to load per node
    for l in get_nonctrl_generation(system), c in contingencies
        @constraint(opf_m, value.(p_opf_m[:ls0][get_name(l)]) + lsc[get_name(l),c] <= get_active_power(l))
    end

    return opf_m
end
