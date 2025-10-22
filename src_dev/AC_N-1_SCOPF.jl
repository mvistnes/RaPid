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

function ac_scopf(type::OPF, system::System, optimizer; 
        voll = nothing, 
        contingencies = nothing, 
        prob = nothing,
        time_limit_sec::Integer = 600,
        unit_commit::Bool = false,
        max_shed::Real = 1.0,
        max_curtail::Real = 1.0,
        circuit_breakers::Bool=false,
        short_term_limit_multi::Real = 1.5,
        ramp_minutes::Real = 10,
        ramp_mult::Real = 10,
        debug=false
    )
    contingencies = isnothing(contingencies) ? sort_components!(get_branches(system)) : contingencies
    prob = isnothing(prob) ? make_prob(contingencies) : prob

    opfm = isnothing(voll) ? opfmodel(system, optimizer, time_limit_sec, debug=debug) : opfmodel(system, optimizer, time_limit_sec, voll, contingencies, prob, debug=debug)
    idx = get_nodes_idx(opfm.nodes)
    list = make_list(opfm, idx, opfm.nodes)
    slack = find_slack(opfm.nodes)

    vm0_limit::NamedTuple{(:min, :max), Tuple{Float64, Float64}} = (0.90, 1.05)
    p_lim = get_active_power_limits.(opfm.ctrl_generation)
    q_lim = get_reactive_power_limits.(opfm.ctrl_generation)
    pr = get_active_power.(opfm.renewables)
    qr = get_reactive_power.(opfm.renewables)
    pd = get_active_power.(opfm.demands)
    qd = get_reactive_power.(opfm.demands)
    @variables(opfm.mod, begin
        0 <= pg0[g in 1:length(opfm.ctrl_generation)] <= p_lim[g].max
            # active power variables for the generators
        q_lim[g].min <= qg0[g in 1:length(opfm.ctrl_generation)] <= q_lim[g].max
                # and reactive
        pf0_fr[l in 1:length(opfm.branches)]
            # active power flow from bus on lines in base case 
        qf0_fr[l in 1:length(opfm.branches)]
                # and reactive 
        pf0_to[l in 1:length(opfm.branches)]
            # active power flow to bus on lines in base case 
        qf0_to[l in 1:length(opfm.branches)]
                # and reactive 
        pfdc0[l in 1:length(opfm.dc_branches)]
            # power flow on DC branches
        qfdc0[l in 1:length(opfm.dc_branches)]
                # and reactive 
        vm0_limit.min <= vm0[b in 1:length(opfm.nodes)] <= vm0_limit.max
            # voltage magnitude at a node in base case 
        va0[b in 1:length(opfm.nodes)]
            # voltage angle at a node in base case 
        0 <= ls0[b in 1:length(opfm.demands)] <= pd[b] * max_shed
            # load curtailment variables
        0 <= pr0[b in 1:length(opfm.renewables)] <= pr[b] * max_curtail
            # renewable curtailment variables
    end)

    register(opfm.mod, :sum, 1, sum; autodiff = true)

    # minimize cost of generation
    @objective(opfm.mod, Min, # opfm.cost_ctrl_gen' * opfm.mod[:pg0] + opfm.voll' * opfm.mod[:ls0]
        sum(x * opfm.mod[:pg0][i] #= x[2] * opfm.mod[:pg0][i]^2 =# for (i,x) in enumerate(opfm.cost_ctrl_gen))+ 
        # opfm.cost_renewables' * opfm.mod[:pr0] + 
        sum(x * opfm.mod[:ls0][i] for (i,x) in enumerate(opfm.voll))
    )

    k = [idx[l.arc.from.number] for l in opfm.branches]
    m = [idx[l.arc.to.number] for l in opfm.branches]
    # incerted power at each bus for the base case and contingencies
    @constraint(opfm.mod, inj_p[n = 1:length(opfm.nodes)], 
        sum(isequal(k[l], n) * opfm.mod[:pf0_fr][l] - isequal(m[l], n) * opfm.mod[:pf0_to][l] for l in list[n].branches) == 
        sum(opfm.mod[:pfdc0][l] for l in list[n].dc_branches) +
        sum(opfm.mod[:pg0][g] for g in list[n].ctrl_generation) + 
        sum((pr[d] - opfm.mod[:pr0][d] for d in list[n].renewables), init = 0.0) - 
        sum((pd[d] - opfm.mod[:ls0][d] for d in list[n].demands), init = 0.0)
    )
    @constraint(opfm.mod, inj_q[n = 1:length(opfm.nodes)], 
        sum(isequal(k[l], n) * opfm.mod[:qf0_fr][l] - isequal(m[l], n) * opfm.mod[:qf0_to][l] for l in list[n].branches) == 
        sum(opfm.mod[:qfdc0][l] for l in list[n].dc_branches) +
        sum(opfm.mod[:qg0][g] for g in list[n].ctrl_generation) + 
        sum((qr[d] * (1 - opfm.mod[:pr0][d] / pr[d]) for d in list[n].renewables), init = 0.0) - 
        sum((qd[d] * (1 - opfm.mod[:ls0][d] / pd[d]) for d in list[n].demands), init = 0.0)
    )
    # @constraint(opfm.mod, va_lim, -π/2 <= opfm.mod[:va0] <= π/2) # Not really needed, could be implemented with spesific line angle limits
    @constraint(opfm.mod, opfm.mod[:va0][slack[1]] == 0) # Set voltage angle at reference bus
    
    # power flow on line and line limits
    for (l, branch) in enumerate(opfm.branches)
        branch_rating = get_rate(branch)
        g, b, B, tap, tr, ti = get_specs(branch)

        pf = opfm.mod[:pf0_fr][l]
        pt = opfm.mod[:pf0_to][l]
        qf = opfm.mod[:qf0_fr][l]
        qt = opfm.mod[:qf0_to][l]

        vm_fr = opfm.mod[:vm0][k[l]]
        vm_to = opfm.mod[:vm0][m[l]]
        va_fr = opfm.mod[:va0][k[l]]
        va_to = opfm.mod[:va0][m[l]]

        @NLconstraint(opfm.mod, pf^2 + qf^2 <= branch_rating^2)
        @NLconstraint(opfm.mod, pt^2 + qt^2 <= branch_rating^2)
        
        @NLconstraint(opfm.mod, pf == g * vm_fr^2 / tap^2 + 
            (-g * tr + b * ti) / tap^2 * (vm_fr * vm_to * cos(va_fr - va_to)) +
            (-b * tr - g * ti) / tap^2 * (vm_fr * vm_to * sin(va_fr - va_to))
        )
        @NLconstraint(opfm.mod, qf == -(b + B.from) * vm_fr^2 / tap^2 - 
            (-b * tr - g * ti) / tap^2 * (vm_fr * vm_to * cos(va_fr - va_to)) +
            (-g * tr + b * ti) / tap^2 * (vm_fr * vm_to * sin(va_fr - va_to))
        )
        @NLconstraint(opfm.mod, pt == g * vm_to^2 + 
            (-g * tr - b * ti) / tap^2 * (vm_to * vm_fr * cos(va_to - va_fr)) +
            (-b * tr + g * ti) / tap^2 * (vm_to * vm_fr * sin(va_to - va_fr))
        )
        @NLconstraint(opfm.mod, qt == -(b + B.to) * vm_to^2 - 
            (-b * tr + g * ti) / tap^2 * (vm_to * vm_fr * cos(va_to - va_fr)) +
            (-g * tr - b * ti) / tap^2 * (vm_to * vm_fr * sin(va_to - va_fr))
        )
    end
    # @NLconstraint(opfm.mod, pb0[l = 1:length(opfm.branches)], opfm.mod[:pf0][l] - opfm.mod[:vm0][k[l]] * 
    #     (opfm.mod[:vm0][m[l]] * (g[l] * cos(opfm.mod[:va0][k[l]] - opfm.mod[:va0][m[l]]) + 
    #     b[l] * sin(opfm.mod[:va0][k[l]] - opfm.mod[:va0][m[l]]))) == 0
    # )
    # @NLconstraint(opfm.mod, qb0[l = 1:length(opfm.branches)], opfm.mod[:qf0][l] - opfm.mod[:vm0][k[l]] * 
    #     (opfm.mod[:vm0][m[l]] * (g[l] * sin(opfm.mod[:va0][k[l]] - opfm.mod[:va0][m[l]]) - 
    #     b[l] * cos(opfm.mod[:va0][k[l]] - opfm.mod[:va0][m[l]]))) == 0
    # )

    if length(opfm.dc_branches) > 0 
        @error "Check constraints before run!"
        p_dc_lim_from = get_active_power_limits_from.(opfm.dc_branches)
        p_dc_lim_to = get_active_power_limits_to.(opfm.dc_branches)
        q_dc_lim_from  = get_reactive_power_limits_from.(opfm.dc_branches)
        q_dc_lim_to = get_reactive_power_limits_to.(opfm.dc_branches)
        @constraints(opfm.mod, begin
            pfdc0_lim_n[l = 1:length(opfm.dc_branches)], p_dc_lim_from[l].min <= opfm.mod[:pfdc0][l]
            pfdc0_lim_p[l = 1:length(opfm.dc_branches)], opfm.mod[:pfdc0][l] <= p_dc_lim_from[l].max
            qfdc0_lim_n[l = 1:length(opfm.dc_branches)], q_dc_lim_from[l].min <= opfm.mod[:qfdc0][l]
            qfdc0_lim_p[l = 1:length(opfm.dc_branches)], opfm.mod[:qfdc0][l] <= q_dc_lim_from[l].max
        end)
    end

    return opfm
end

function p_scopf(system::System, optimizer; pen = 10000, time_limit_sec = 120)
    opf_m = Model(optimizer)
    set_time_limit_sec(opf_m, time_limit_sec)
    vm0_limit::NamedTuple{(:min, :max), Tuple{Float64, Float64}} = (0.90, 1.05)
    vmc_limit::NamedTuple{(:min, :max), Tuple{Float64, Float64}} = (0.90, 1.10) # 0.85 for continental Europe, 0.9 for Nordic and UK

    gen_names       = get_name.(get_components(Generator, system))
    line_names      = get_name.(get_components(Branch, system))
    bus_names       = get_name.(get_components(Bus, system))
    demand_names    = get_name.(get_components(StaticLoad, system))
    contingencies   = line_names

    
    @variables(opf_m, begin
        pg0[g in gen_names] >= 0                    # active power variables for the generators
        qg0[g in gen_names]                             # and reactive
        pf0[l in line_names]                        # active power flow on lines in base case 
        pfc[l in line_names, c in contingencies]        # and contingencies
        qf0[l in line_names]                            # and reactive 
        qfc[l in line_names, c in contingencies]        # and contingencies
        vm0[b in bus_names] >= 0                    # voltage magnitude at a node in base case 
        vmc[b in bus_names, c in contingencies] >= 0    # and contingencies
        va0[b in bus_names]                         # voltage angle at a node in base case 
        vac[b in bus_names, c in contingencies]         # and contingencies
        ls0[b in demand_names] >= 0                 # load curtailment variables
    end)

    # find value of beta, beta = 1 if from bus is the bus, beta = -1 if to bus is the bus, 0 else
    beta = JuMP.Containers.DenseAxisArray{Int8}(undef, bus_names, line_names)
    for bus in bus_names
        for line in get_components(Branch, system)
            beta[bus, get_name(line)] = bus == get_name(get_arc(line).from) ? 1 : bus == get_name(get_arc(line).to) ? -1 : 0
        end
    end

    register(opf_m, :sum, 1, sum; autodiff = true)

    # minimize cost of generation
    @NLobjective(opf_m, Min, 
        sum(pg0[get_name(g)] * get_cost(get_variable(get_operation_cost(g)))[1] +
        pg0[get_name(g)]^2 * get_cost(get_variable(get_operation_cost(g)))[2]
        for g in get_components(ThermalGen, system)
        ) +  
        sum(pg0[get_name(g)] * 0.5 + pg0[get_name(g)]^2 * 5 for g in get_components(HydroDispatch, system)) + 
        pen * sum(ls0[d] for d in demand_names)
    )

    # incerted power at each bus for the base case and contingencies
    @expressions(opf_m,begin
        generators_p[b = bus_names], sum(get_name(get_bus(g)) == b ? pg0[get_name(g)] : 0 for g in get_components(Generator, system))
        generators_q[b = bus_names], sum(get_name(get_bus(g)) == b ? qg0[get_name(g)] : 0 for g in get_components(Generator, system))
    end)
    for b in bus_names
        demand_p = sum(get_name(get_bus(d)) == b ? get_active_power(d) - ls0[get_name(d)] : 0 for d in get_components(StaticLoad, system))
        demand_q = sum((1 - ls0[get_name(d)]/get_active_power(d))*get_reactive_power(d) for d in get_components(StaticLoad, system))
        @constraint(opf_m, generators_p[b] - sum(beta[b,l] * pf0[l] for l in line_names) == demand_p)
        @constraint(opf_m, generators_q[b] - sum(beta[b,l] * qf0[l] for l in line_names) == demand_q)
        for c in contingencies
            @constraint(opf_m, generators_p[b] - sum(beta[b,l] * pfc[l,c] for l in line_names) == demand_p)
            @constraint(opf_m, generators_q[b] - sum(beta[b,l] * qfc[l,c] for l in line_names) == demand_q)
        end
    end

    # register(opf_m, :real, 1, real; autodiff = true)
    # power flow on line and line limits for the base case and contingencies
    for l in get_components(Branch, system)
        y = 1/(get_r(l) + get_x(l)*im)
        (g,b) = (real(y), imag(y))
        l_rate = get_rate(l)
        l_name = get_name(l)
        @NLconstraint(opf_m, pf0[l_name] - sum(beta[k,l_name] * vm0[k] * 
            sum(vm0[m] * (g * cos(va0[k] - va0[m]) + b * sin(va0[k] - va0[m])) for m in bus_names) 
            for k in bus_names) == 0
        )
        @NLconstraint(opf_m, qf0[l_name] - sum(beta[k,l_name] * vm0[k] * 
            sum(vm0[m] * (g * sin(va0[k] - va0[m]) - b * cos(va0[k] - va0[m])) for m in bus_names) 
            for k in bus_names) == 0
        )
        @NLconstraint(opf_m, pf0[l_name]^2 + qf0[l_name]^2 <= l_rate^2)
        for c in contingencies
            a = l != c # zero value if l is unavailable under contingency c
            @NLconstraint(opf_m, pfc[l_name,c] - a * sum(beta[k,l_name] * vmc[k,c] * 
                sum(vmc[m,c] * (g * cos(vac[k,c] - vac[m,c]) + b * sin(vac[k,c] - vac[m,c])) for m in bus_names) 
                for k in bus_names) == 0
            )
            @NLconstraint(opf_m, qfc[l_name,c] - a * sum(beta[k,l_name] * vmc[k,c] * 
                sum(vmc[m,c] * (g * sin(vac[k,c] - vac[m,c]) - b * cos(vac[k,c] - vac[m,c])) for m in bus_names) 
                for k in bus_names) == 0
            )
            @NLconstraint(opf_m, pfc[l_name,c]^2 + qfc[l_name,c]^2 <= l_rate^2)
        end
    end
    @constraint(opf_m, va0[find_slack(opfm.nodes)[2].name] == 0)

    # restrict active and reactive power generation to min and max values
    for g in get_components(Generator, system)
        set_lower_bound(qg0[get_name(g)], get_reactive_power_limits(g).min) 
        set_upper_bound(pg0[get_name(g)], get_active_power_limits(g).max)
        set_upper_bound(qg0[get_name(g)], get_reactive_power_limits(g).max)
    end

    # restrict voltage magnitude to min and max values
    for b in bus_names
        set_lower_bound(vm0[b], vm0_limit.min) 
        set_upper_bound(vm0[b], vm0_limit.max)
        for c in contingencies
            set_lower_bound(vmc[b,c], vmc_limit.min)
            set_upper_bound(vmc[b,c], vmc_limit.max)
        end
    end

    # restrict load shedding to max load per node
    for d in get_components(StaticLoad, system)
        set_upper_bound(ls0[get_name(d)], get_active_power(d))
    end

    optimize!(opf_m)
    return opf_m
end



function pc_scopf(system::System, optimizer; pen = 10000, time_limit_sec = 100)
    opf_m = Model(optimizer)
    set_time_limit_sec(opf_m, time_limit_sec)
    vm0_limit::NamedTuple{(:min, :max), Tuple{Float64, Float64}} = (0.90, 1.05)
    vmc_limit::NamedTuple{(:min, :max), Tuple{Float64, Float64}} = (0.90, 1.10) # 0.85 for continental Europe, 0.9 for Nordic and UK
    ramp_minutes = 10
    short_term_limit_multi = 1.5

    gen_names = get_name.(get_components(Generator, system))
    line_names = get_name.(get_components(Branch, system))
    bus_names = get_name.(get_components(Bus, system))
    demand_names = get_name.(get_components(StaticLoad, system))
    contingencies = line_names # [bus_names; line_names]

    prob = [ # spesified for the RTS-96
            # 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
            # 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, # generators
            0.24, 0.51, 0.33, 0.39, 0.48, 0.38, 0.02, 0.36, 0.34, 0.33, 0.30, 0.44, 0.44, 0.44, 
            0.02, 0.02, 0.02, 0.02, 0.40, 0.39, 0.40, 0.52, 0.49, 0.47, 0.38, 0.33, 0.41, 0.41, 
            0.41, 0.35, 0.34, 0.32, 0.54, 0.35, 0.35, 0.38, 0.38, 0.34, 0.34, 0.45, 0.46 # lines
    ]
    prob /= 8760

    # active and reactive power variables for the generators in base case and contingencies ramp up and ramp down
    @variable(opf_m, pg0[g in gen_names] >= 0) 
    @variable(opf_m, pgu[g in gen_names, c in contingencies] >= 0) 
    @variable(opf_m, pgd[g in gen_names, c in contingencies] >= 0) 
    @variable(opf_m, qg0[g in gen_names]) 
    @variable(opf_m, qgc[g in gen_names, c in contingencies])
    # power flow on lines in base case and in contingencies before and after corrective actions
    @variable(opf_m, pf0[l in line_names])
    @variable(opf_m, pfc[l in line_names, c in contingencies]) # (c)ontingency
    @variable(opf_m, pfcc[l in line_names, c in contingencies]) # (c)ontingency (c)orrected
    @variable(opf_m, qf0[l in line_names])
    @variable(opf_m, qfc[l in line_names, c in contingencies])
    @variable(opf_m, qfcc[l in line_names, c in contingencies])
    # voltage magnitude at a node in base case and contingencies before and after corrective actions
    @variable(opf_m, vm0[b in bus_names] >= 0)
    @variable(opf_m, vmc[b in bus_names, c in contingencies] >= 0)
    @variable(opf_m, vmcc[b in bus_names, c in contingencies] >= 0)
    # voltage angle at a node in base case and in contingencies before and after corrective actions
    @variable(opf_m, va0[b in bus_names])
    @variable(opf_m, vac[b in bus_names, c in contingencies])
    @variable(opf_m, vacc[b in bus_names, c in contingencies])
    # load curtailment variables in base case and in contingencies
    @variable(opf_m, ls0[d in demand_names] >= 0)
    @variable(opf_m, lsc[d in demand_names, c in contingencies] >= 0)

    # find value of beta, beta = 1 if from bus is the bus, beta = -1 if to bus is the bus, 0 else
    beta = JuMP.Containers.DenseAxisArray{Int8}(undef, bus_names, line_names)
    for bus in bus_names
        for line in get_components(Branch, system)
            beta[bus, get_name(line)] = bus == get_name(get_arc(line).from) ? 1 : bus == get_name(get_arc(line).to) ? -1 : 0
        end
    end

    # minimize cost of generation
    @NLobjective(opf_m, Min, sum(
        get_cost(get_variable(get_operation_cost(g)))[1] * (pg0[get_name(g)] + 
        sum(prob[i] * pgu[get_name(g),c] for (i,c) in enumerate(contingencies)))
        for g in get_components(ThermalGen, system)
        ) + sum(
        get_cost(get_variable(get_operation_cost(g)))[2] * (pg0[get_name(g)]^2 + 
        sum(prob[i] * pgu[get_name(g),c]^2 for (i,c) in enumerate(contingencies)))
        for g in get_components(ThermalGen, system)
        ) + 
        pen * sum(ls0[d] + sum(prob[i] * lsc[d,c] for (i,c) in enumerate(contingencies)) for d in demand_names) # removed "|CN1|", added "prob[c] * "
    )

    # incerted power at each bus for the base case and contingencies
    for b in get_components(Bus, system)
        b_name = get_name(b)
        @constraint(opf_m, sum(get_bus(g) == b ? pg0[get_name(g)] : 0 for g in get_components(Generator, system)) - 
                           sum(beta[b_name,l] * pf0[l] for l in line_names) ==
                           sum(get_active_power(d) - ls0[get_name(d)] for d in get_components(StaticLoad, system))
        )
        @constraint(opf_m, sum(get_bus(g) == b ? qg0[get_name(g)] : 0 for g in get_components(Generator, system)) - 
                           sum(beta[b_name,l] * qf0[l] for l in line_names) == 
                           sum((1 - ls0[get_name(d)]/get_active_power(d))*get_reactive_power(d) 
                           for d in get_components(StaticLoad, system))
        )
        for c in contingencies
            @constraint(opf_m, sum(get_bus(g) == b ? pg0[get_name(g)] : 0 for g in get_components(Generator, system)) -
                                sum(beta[b_name,l] * pfc[l,c] for l in line_names) ==
                                sum(get_active_power(d) - ls0[get_name(d)] for d in get_components(StaticLoad, system))
            )
            @constraint(opf_m, sum(get_bus(g) == b ? qg0[get_name(g)] : 0 for g in get_components(Generator, system)) -
                                sum(beta[b_name,l] * qfc[l,c] for l in line_names) == 
                                sum((1 - ls0[get_name(d)]/get_active_power(d))*get_reactive_power(d) 
                                for d in get_components(StaticLoad, system))
            )
            @constraint(opf_m, sum(get_bus(g) == b ? pg0[get_name(g)] + pgu[get_name(g),c] - pgd[get_name(g),c] : 0 
                                    for g in get_components(Generator, system)) -
                                sum(beta[b_name,l] * pfcc[l,c] for l in line_names) ==
                                sum(get_active_power(d) - ls0[get_name(d)] - lsc[get_name(d),c] 
                                    for d in get_components(StaticLoad, system))
            )
            @constraint(opf_m, sum(get_bus(g) == b ? qg0[get_name(g)] + qgc[get_name(g),c] : 0 
                                    for g in get_components(Generator, system)) -
                                sum(beta[b_name,l] * qfcc[l,c] for l in line_names) == 
                                sum((1 - ls0[get_name(d)]/get_active_power(d) - lsc[get_name(d),c]/get_active_power(d)) * 
                                get_reactive_power(d) for d in get_components(StaticLoad, system))
            )
        end
    end

    # power flow on line and line limits for the base case and contingencies
    for l in get_components(Branch, system)
        y = 1/(get_r(l) + get_x(l)*im)
        (g,b) = (real(y), imag(y))
        l_rate = get_rate(l)
        l_name = get_name(l)
        @NLconstraint(opf_m, pf0[l_name] - sum(beta[k,l_name] * vm0[k] * 
            sum(vm0[m] * (g * cos(va0[k] - va0[m]) + b * sin(va0[k] - va0[m])) for m in bus_names) 
            for k in bus_names) == 0
        )
        @NLconstraint(opf_m, qf0[l_name] - sum(beta[k,l_name] * vm0[k] * 
            sum(vm0[m] * (g * sin(va0[k] - va0[m]) - b * cos(va0[k] - va0[m])) for m in bus_names) 
            for k in bus_names) == 0
        )
        @NLconstraint(opf_m, pf0[l_name]^2 + qf0[l_name]^2 <= l_rate^2)
        for c in contingencies
            l_name = get_name(l)
            a = l != c # zero value if l is unavailable under contingency c
            @NLconstraint(opf_m, pfc[l_name,c] - a * sum(beta[k,l_name] * vmc[k,c] * 
                sum(vmc[m,c] * (g * cos(vac[k,c] - vac[m,c]) + b * sin(vac[k,c] - vac[m,c])) for m in bus_names) 
                for k in bus_names) == 0
            )
            @NLconstraint(opf_m, qfc[l_name,c] - a * sum(beta[k,l_name] * vmc[k,c] * 
                sum(vmc[m,c] * (g * sin(vac[k,c] - vac[m,c]) - b * cos(vac[k,c] - vac[m,c])) for m in bus_names) 
                for k in bus_names) == 0
            )
            @NLconstraint(opf_m, pfc[l_name,c]^2 + qfc[l_name,c]^2 <= short_term_limit_multi^2 * l_rate^2)
            @NLconstraint(opf_m, pfcc[l_name,c] - a * sum(beta[k,l_name] * vmcc[k,c] * 
                sum(vmcc[m,c] * (g * cos(vacc[k,c] - vacc[m,c]) + b * sin(vacc[k,c] - vacc[m,c])) for m in bus_names) 
                for k in bus_names) == 0
            )
            @NLconstraint(opf_m, qfcc[l_name,c] - a * sum(beta[k,l_name] * vmcc[k,c] * 
                sum(vmcc[m,c] * (g * sin(vacc[k,c] - vacc[m,c]) - b * cos(vacc[k,c] - vacc[m,c])) for m in bus_names) 
                for k in bus_names) == 0
            )
            @NLconstraint(opf_m, pfcc[l_name,c]^2 + qfcc[l_name,c]^2 <= l_rate^2)
        end
    end

    # restrict active power generation to min and max values
    for g in get_components(Generator, system)
        g_name = get_name(g)
        set_lower_bound(qg0[g_name], get_reactive_power_limits(g).min) 
        set_upper_bound(pg0[g_name], get_active_power_limits(g).max)
        set_upper_bound(qg0[g_name], get_reactive_power_limits(g).max)
        for c in contingencies
            @constraint(opf_m, 0 <= pg0[g_name] + (pgu[g_name,c] - pgd[g_name,c]) <= get_active_power_limits(g).max)
            @constraint(opf_m, get_reactive_power_limits(g).min <= qg0[g_name] + qgc[g_name,c] <= get_reactive_power_limits(g).max)
            set_upper_bound(pgu[g_name,c], get_ramp_limits(g).up * ramp_minutes)
            set_upper_bound(pgd[g_name,c], get_ramp_limits(g).down * ramp_minutes)
        end
    end

    # restrict voltage magnitude to min and max values
    for b in bus_names
        set_lower_bound(vm0[b], vm0_limit.min) 
        set_upper_bound(vm0[b], vm0_limit.max)
        for c in contingencies
            set_lower_bound(vmc[b,c], vmc_limit.min)
            set_upper_bound(vmc[b,c], vmc_limit.max)
            set_lower_bound(vmcc[b,c], vm0_limit.min)
            set_upper_bound(vmcc[b,c], vm0_limit.max)
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

    gen_names = get_name.(get_components(Generator, system))
    line_names = get_name.(get_components(Branch, system))
    bus_names = get_name.(get_components(Bus, system))
    demand_names = get_name.(get_components(StaticLoad, system))

    # power flow on lines in base case and emergency
    @variable(opf_m, pf0[l in line_names])
    @variable(opf_m, fe[l in line_names])
    # voltage angle at a node in emergency
    @variable(opf_m, anglee[b in bus_names])
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
                           sum(beta[b,l] * pf0[get_name(l)] for l in get_components(Branch, system)) ==
                           sum(get_active_power(l) - ls0[get_name(l)] - lse[get_name(l)] for l in get_components(StaticLoad, system))
        )
    end

    # power flow on line and line limits for the base case
    for l in get_components(Branch, system)
        l_name = get_name(l)
        @constraint(opf_m, fe[l_name] - outage[l_name] * sum(beta[b,l] * anglee[get_name(b)] for b in get_components(Bus, system)) / get_x(l) == 0)
        @constraint(opf_m, -get_rate(l) <= fe[l_name] <= get_rate(l))
    end

    optimize!(opf_m)
    for d in demand_names
        value(lse[d]) > 0.00001 ? @printf("%s = %.4f \n", d,value(lse[d])) : print()
    end
    for l in get_components(Branch, system)
        @printf("%s = %.4f <= %.2f \n", get_name(l), value(pf0[get_name(l)]), get_rate(l))
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


# for l in get_name.(get_components(Branch, ieee_rts))
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