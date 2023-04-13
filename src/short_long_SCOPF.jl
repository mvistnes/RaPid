# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
using PowerSystems
using JuMP
# include("N-1_SCOPF.jl")

function sl_scopf(system::System, optimizer; 
        voll = nothing, 
        contingencies = nothing, 
        prob = nothing,
        time_limit_sec = 600,
        unit_commit::Bool=false,
        max_shed = 0.1,
        max_curtail::Float64 = 1.0,
        ratio = 0.5, 
        short_term_limit_multi = 1.5,
        ramp_minutes::Int64 = 10,
        repair_time::Float64 = 1.0)
    system = set_renewable_prod!(system, ratio)
    voll = isnothing(voll) ? make_voll(system) : voll
    contingencies = isnothing(contingencies) ? get_name.(get_branches(system)) : contingencies
    prob = isnothing(prob) ? make_prob(contingencies) : prob
    
    p_opfm = scopf(PSC, system, optimizer, voll=voll, prob=prob, contingencies=contingencies, time_limit_sec=time_limit_sec,
        unit_commit=unit_commit, max_shed=max_shed, max_curtail=max_curtail, ratio=ratio, short_term_limit_multi=short_term_limit_multi, 
        ramp_minutes=ramp_minutes)
    p_opfm.mod = solve_model!(p_opfm.mod)
    c_mod = c_scopf(system, optimizer, voll, contingencies, prob, value.(p_opfm.mod[:pg0]), value.(p_opfm.mod[:ls0]), value.(p_opfm.mod[:lsc]), 
        unit_commit ? value.(p_opfm.mod[:u]) : nothing, time_limit_sec=time_limit_sec, unit_commit=unit_commit, max_shed=max_shed, max_curtail=max_curtail,
        ramp_minutes=ramp_minutes, repair_time=repair_time)
    c_mod = solve_model!(c_mod)

    return p_opfm, c_mod
end

function c_scopf(system::System, optimizer, voll, contingencies, prob, pg0, ls0, lsc, u = nothing; 
        time_limit_sec::Int64 = 600,
        unit_commit::Bool = false,
        max_shed::Float64 = 0.1,
        max_curtail::Float64 = 1.0,
        ramp_minutes::Int64 = 10,
        repair_time::Float64 = 1.0)

    model = Model(optimizer)
    cost = JuMP.Containers.DenseAxisArray(
        [[get_generator_cost(g)[2] for g in get_gens_t(system)]; [5 for _ in get_gens_h(system)]],
        get_name.(get_ctrl_generation(system))
    )

    
    @variables(model, begin
        0 <= pgu[g in get_name.(get_ctrl_generation(system)), c in contingencies]    # active power variables for the generators in contingencies ramp up 
        0 <= pgd[g in get_name.(get_ctrl_generation(system)), c in contingencies]       # and ramp down
        pfcc[l in get_name.(get_branches(system)), c in contingencies]         # and after corrective actions
        vacc[b in get_name.(get_nodes(system)), c in contingencies]            # and after corrective actions
        0 <= lscc[d in get_name.(get_nonctrl_generation(system)), c in contingencies] # load curtailment variables in in contingencies
    end)

    # minimize socio-economic cost
    @objective(model, Min, 
        sum(cost[g] * (pg0[g] + sum(prob[c] * repair_time * (pgu[g,c] #=+ pgd[get_name(g),c]*0.1=#) 
            for c in contingencies)) for g in get_name.(get_ctrl_generation(system))
        ) +
        sum(voll[d] * (ls0[d] + sum(prob[c] * (lsc[d,c] * ramp_minutes / 60 + lscc[d,c] * repair_time) for c in contingencies)) 
            for d in get_name.(get_nonctrl_generation(system))
        )
    )
    if unit_commit
        set_objective_function(model, objective_function(model) + 
        sum(u[get_name(g)] * get_generator_cost(g)[end] for g in gens_t(system)))
    end

    # incerted power at each bus for the base case and contingencies
    @constraint(model, inj_pcc[n = get_name.(get_nodes(system)), c = contingencies],
        sum(g.bus.name == n ? pg0[get_name(g)] + pgu[get_name(g),c] - pgd[get_name(g),c] : 0 for g in get_ctrl_generation(system)) -
        sum(beta(n,l) * pfcc[get_name(l),c] for l in get_branches(system)) == 
        sum(d.bus.name == n ? get_active_power(d) - lscc[get_name(d),c] : 0 for d in get_demands(system)) +
        sum(d.bus.name == n ? -get_active_power(d) + lscc[get_name(d),c] : 0 for d in get_renewables(system))
    )
    # @constraint(model, -π/2 .<= vacc .<= π/2) # Not really needed, could be implemented with spesific line angle limits
    
    i, slack = find_slack(system)
    @constraint(model, [c = contingencies], vacc[get_name(slack),c] == 0)
    
    # power flow on branch and branch limits for the base case and contingencies
    branch_rating = make_named_array(get_rate, get_branches(system))
    @constraint(model, pfcc_lim[l = get_name.(get_branches(system)), c = contingencies], 
        -branch_rating[l] .* a(l,c) .<= pfcc[l,c] .<= branch_rating[l] .* a(l,c)
    )
    x = make_named_array(get_x, get_branches(system))
    @constraint(model, pbcc[l = get_name.(get_branches(system)), c = contingencies],
        pfcc[l,c] .- a(l,c) .* sum(beta(system,l) .* vacc[:,c]) ./ x[l] .== 0
    )

    # restrict active power generation to min and max values
    for g in get_gens_t(system)
        g_name = get_name(g)
        for c in contingencies
            if unit_commit
                @constraint(model, 0 <= pg0[g_name] + (pgu[g_name,c] - pgd[g_name,c]) <= 
                    get_active_power_limits(g).max * u[g_name])
            else
                @constraint(model, 0 <= pg0[g_name] + (pgu[g_name,c] - pgd[g_name,c]) <= 
                    get_active_power_limits(g).max)
            end
            set_upper_bound(pgu[g_name,c], get_ramp_limits(g).up * ramp_minutes)
            set_upper_bound(pgd[g_name,c], get_ramp_limits(g).down * ramp_minutes)
        end
    end
    for g in get_gens_h(system)
        g_name = get_name(g)
        for c in contingencies
            @constraint(model, 0 <= pg0[g_name] + (pgu[g_name,c] - pgd[g_name,c]) <= 
                    get_active_power_limits(g).max)
            set_upper_bound(pgu[g_name,c], get_ramp_limits(g).up * ramp_minutes)
            set_upper_bound(pgd[g_name,c], get_ramp_limits(g).down * ramp_minutes)
        end
    end

    # restrict load shedding to load per node
    for l in get_demands(system), c in contingencies
        @constraint(model, lscc[get_name(l),c] <= get_active_power(l) * max_shed)
    end
    # restrict renewable shedding to renewable production
    for l in get_renewables(system), c in contingencies
        @constraint(model, lscc[get_name(l),c] <= get_active_power(l) * max_curtail)
    end

    return model
end

function d_scopf(system::System, optimizer, pg0; 
    time_limit_sec::Int64 = 600,
    max_shed::Float64 = 0.1,
    max_curtail::Float64 = 1.0,
    ramp_minutes::Int64 = 10)

    model = Model(optimizer)
    cost = JuMP.Containers.DenseAxisArray(
        [[get_generator_cost(g)[2] for g in get_gens_t(system)]; [5 for _ in get_gens_h(system)]],
        get_name.(get_ctrl_generation(system))
    )


    @variables(model, begin
        0 <= pgu[g in get_name.(get_ctrl_generation(system))]    # active power variables for the generators in contingencies ramp up 
        0 <= pgd[g in get_name.(get_ctrl_generation(system))]       # and ramp down
        pfcc[l in get_name.(get_branches(system))]         # and after corrective actions
        vacc[b in get_name.(get_nodes(system))]            # and after corrective actions
        0 <= lscc[d in get_name.(get_nonctrl_generation(system))] # load curtailment variables in in contingencies
        0 <= s[g in get_name.(get_ctrl_generation(system))]
    end)

    @objective(model, Min, sum(s))

    @constraint(model, inj_pcc[n = get_name.(get_nodes(system))],
        sum(g.bus.name == n ? pg0[get_name(g)] + pgu[get_name(g)] - pgd[get_name(g)] : 0 for g in get_ctrl_generation(system)) -
        sum(beta(n,l) * pfcc[get_name(l)] for l in get_branches(system)) == 
        sum(d.bus.name == n ? get_active_power(d) - lscc[get_name(d)] : 0 for d in get_demands(system)) +
        sum(d.bus.name == n ? -get_active_power(d) + lscc[get_name(d)] : 0 for d in get_renewables(system))
    )
    i, slack = find_slack(system)
    @constraint(model, vacc[get_name(slack)] == 0)

    branch_rating = make_named_array(get_rate, get_branches(system))
    @constraint(model, pfcc_lim[l = get_name.(get_branches(system))], 
        -branch_rating[l] .<= pfcc[l] .<= branch_rating[l]
    )
    x = make_named_array(get_x, get_branches(system))
    @constraint(model, pbcc[l = get_name.(get_branches(system))],
        pfcc[l] .- sum(beta(system,l) .* vacc[:]) ./ x[l] .== 0
    )

    for g in get_gens_t(system)
        g_name = get_name(g)
        @constraint(model, 0 <= pg0[g_name] + (pgu[g_name] - pgd[g_name]) <= 
            get_active_power_limits(g).max)
        # @constraint(model, pgu[g_name] <= get_ramp_limits(g).up * ramp_minutes + s[g_name])
        # @constraint(model, pgd[g_name] <= get_ramp_limits(g).down * ramp_minutes + s[g_name])
        @constraint(model, pgu[g_name] <= s[g_name])
        @constraint(model, pgd[g_name] <= s[g_name])
    end
    for g in get_gens_h(system)
        g_name = get_name(g)
        @constraint(model, 0 <= pg0[g_name] + (pgu[g_name] - pgd[g_name]) <= 
                get_active_power_limits(g).max)
        # @constraint(model, pgu[g_name] <= get_ramp_limits(g).up * ramp_minutes + s[g_name])
        # @constraint(model, pgd[g_name] <= get_ramp_limits(g).down * ramp_minutes + s[g_name])
        @constraint(model, pgu[g_name] <= s[g_name])
        @constraint(model, pgd[g_name] <= s[g_name])
    end

    for l in get_demands(system)
        @constraint(model, lscc[get_name(l)] <= get_active_power(l) * max_shed)
    end
    for l in get_renewables(system)
        @constraint(model, lscc[get_name(l)] <= get_active_power(l) * max_curtail)
    end

    return model
end
