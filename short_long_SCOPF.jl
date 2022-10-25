using PowerSystems
const PSY = PowerSystems
using JuMP
using Ipopt
using Printf
include("N-1_SCOPF.jl")

function sl_scopf(system::System, optimizer; 
        time_limit_sec = 60,
        voll = nothing, 
        contingencies = nothing, 
        unit_commit::Bool=false,
        max_shed = 0.1,
        ratio = 0.5, 
        short_term_limit_multi = 1.5,
        ramp_minutes = 10,
        prob = nothing)
    voll = isnothing(voll) ? make_voll(system) : voll
    contingencies = isnothing(contingencies) ? get_name.(branches(system)) : contingencies
    prob = isnothing(prob) ? make_prob(contingencies) : prob
    
    p_opfm = p_scopf(system, optimizer, contingencies=contingencies, 
        short_term_limit_multi=short_term_limit_multi, voll=voll)
    solve_model!(p_opfm.mod)
    
    set_renewable_prod!(system, ratio)
    model = Model(optimizer)
    
    @variables(model, begin
        pgu[g in get_name.(get_ctrl_generation(system)), c in contingencies] >= 0    # active power variables for the generators in contingencies ramp up 
        pgd[g in get_name.(get_ctrl_generation(system)), c in contingencies] >= 0       # and ramp down
        pfcc[l in get_name.(branches(system)), c in contingencies]         # and after corrective actions
        vacc[b in get_name.(nodes(system)), c in contingencies]            # and after corrective actions
        lsc[d in get_name.(get_nonctrl_generation(system)), c in contingencies] >= 0 # load curtailment variables in in contingencies
    end)

    # minimize socio-economic cost
    @objective(model, Min, 
        sum(get_generator_cost(g) * (value.(p_opfm.mod[:pg0][get_name(g)]) + 
            sum(prob[c] * (pgu[get_name(g),c] #=+ pgd[get_name(g),c]*0.1=#) for c in contingencies))
            for g in gens_t(system)
        ) + 
        sum(5 * (value.(p_opfm.mod[:pg0][get_name(g)]) + sum(prob[c] * (pgu[get_name(g),c] #=+ pgd[get_name(g),c]*0.1=#)
            for c in contingencies)) for g in gens_h(system)
        ) +
        sum(voll[d] * (value.(p_opfm.mod[:ls0][d]) + sum(prob[c] * lsc[d,c] for c in contingencies)) 
            for d in get_name.(get_nonctrl_generation(system))
        )
    )

    # incerted power at each bus for the base case and contingencies
    @constraint(model, inj_pcc[n = get_name.(nodes(system)), c = contingencies],
        sum(get_name(get_bus(g)) == n ? value(p_opfm.mod[:pg0][get_name(g)]) + pgu[get_name(g),c] - pgd[get_name(g),c] : 0 for g in get_ctrl_generation(system)) -
        sum(beta(n,l) * pfcc[get_name(l),c] for l in branches(system)) == 
        sum(get_name(get_bus(d)) == n ? get_active_power(d) - lsc[get_name(d),c] : 0 for d in demands(system)) +
        sum(get_name(get_bus(d)) == n ? -get_active_power(d) + lsc[get_name(d),c] : 0 for d in renewables(system))
    )
    # @constraint(model, -π/2 .<= vacc .<= π/2) # Not really needed, could be implemented with spesific line angle limits
    
    slack = find_slack(system)
    @constraint(model, [c = contingencies], vacc[get_name(slack),c] == 0)
    
    # power flow on branch and branch limits for the base case and contingencies
    branch_rating = make_named_array(get_rate, branches(system))
    @constraint(model, pfcc_lim[l = get_name.(branches(system)), c = contingencies], 
        -branch_rating[l] .* a(l,c) .<= pfcc[l,c] .<= branch_rating[l] .* a(l,c)
    )
    x = make_named_array(get_x, branches(system))
    @constraint(model, pbcc[l = get_name.(branches(system)), c = contingencies],
        pfcc[l,c] .- a(l,c) .* sum(beta(system,l) .* vacc[:,c]) ./ x[l] .== 0
    )

    # restrict active power generation to min and max values
    for g in get_ctrl_generation(system)
        g_name = get_name(g)
        for c in contingencies
            @constraint(model, 0 <= value.(p_opfm.mod[:pg0][g_name]) + (pgu[g_name,c] - pgd[g_name,c]) <= get_active_power_limits(g).max)
            set_upper_bound(pgu[g_name,c], get_ramp_limits(g).up * ramp_minutes)
            set_upper_bound(pgd[g_name,c], get_ramp_limits(g).down * ramp_minutes)
        end
    end

    # restrict load shedding to load per node
    for l in get_nonctrl_generation(system), c in contingencies
        @constraint(model, value.(p_opfm.mod[:ls0][get_name(l)]) + lsc[get_name(l),c] <= get_active_power(l))
    end

    solve_model!(model)
    return p_opfm, model
end
