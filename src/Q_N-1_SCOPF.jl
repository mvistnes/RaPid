# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

function q_scopf(type::OPF, opfm::OPFmodel, optimizer; 
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
    qmod = create_model(optimizer, time_limit_sec, debug=debug)
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
    @variables(qmod, begin
        q_lim[g].min <= qg0[g in 1:length(opfm.ctrl_generation)] <= q_lim[g].max
        qf0_fr[l in 1:length(opfm.branches)]
        qf0_to[l in 1:length(opfm.branches)]
        qfdc0[l in 1:length(opfm.dc_branches)]
        vm0_limit.min <= vm0[b in 1:length(opfm.nodes)] <= vm0_limit.max
        va0[b in 1:length(opfm.nodes)]
    end)

    register(qmod, :sum, 1, sum; autodiff = true)

    # minimize cost of generation
    @objective(qmod, Min, qg0)

    k = [idx[l.arc.from.number] for l in opfm.branches]
    m = [idx[l.arc.to.number] for l in opfm.branches]
    # incerted power at each bus for the base case and contingencies
    @constraint(qmod, inj_q[n = 1:length(opfm.nodes)], 
        sum(isequal(k[l], n) * qf0_fr[l] - isequal(m[l], n) * qf0_to[l] for l in list[n].branches) == 
        sum(qfdc0[l] for l in list[n].dc_branches) +
        sum(qg0[g] for g in list[n].ctrl_generation) + 
        sum((qr[d] * (1 - value(opfm.mod[:pr0][d]) / pr[d]) for d in list[n].renewables), init = 0.0) - 
        sum((qd[d] * (1 - value(opfm.mod[:ls0][d]) / pd[d]) for d in list[n].demands), init = 0.0)
    )
    # @constraint(opfm.mod, va_lim, -π/2 <= opfm.mod[:va0] <= π/2) # Not really needed, could be implemented with spesific line angle limits
    @constraint(qmod, va0[slack[1]] == 0) # Set voltage angle at reference bus
    
    # power flow on line and line limits
    for (l, branch) in enumerate(opfm.branches)
        branch_rating = get_rate(branch)
        g, b, B, tap, tr, ti = get_specs(branch)

        pf = value(opfm.mod[:pf0_fr][l])
        pt = value(opfm.mod[:pf0_to][l])
        qf_lim = sqrt(branch_rating^2 - pf^2)
        qt_lim = sqrt(branch_rating^2 - pt^2)
        qf = qf0_fr[l]
        qt = qf0_to[l]

        vm_fr = vm0[k[l]]
        vm_to = vm0[m[l]]
        va_fr = va0[k[l]]
        va_to = va0[m[l]]

        @constraint(opfm.mod, qf <= qf_lim)
        @constraint(opfm.mod, qt <= qt_lim)
        
        @constraint(opfm.mod, pf == g * vm_fr^2 / tap^2 + 
            (-g * tr + b * ti) / tap^2 * (vm_fr * vm_to * cos(va_fr - va_to)) +
            (-b * tr - g * ti) / tap^2 * (vm_fr * vm_to * sin(va_fr - va_to))
        )
        @constraint(opfm.mod, qf == -(b + B.from) * vm_fr^2 / tap^2 - 
            (-b * tr - g * ti) / tap^2 * (vm_fr * vm_to * cos(va_fr - va_to)) +
            (-g * tr + b * ti) / tap^2 * (vm_fr * vm_to * sin(va_fr - va_to))
        )
        @constraint(opfm.mod, pt == g * vm_to^2 + 
            (-g * tr - b * ti) / tap^2 * (vm_to * vm_fr * cos(va_to - va_fr)) +
            (-b * tr + g * ti) / tap^2 * (vm_to * vm_fr * sin(va_to - va_fr))
        )
        @constraint(opfm.mod, qt == -(b + B.to) * vm_to^2 - 
            (-b * tr + g * ti) / tap^2 * (vm_to * vm_fr * cos(va_to - va_fr)) +
            (-g * tr - b * ti) / tap^2 * (vm_to * vm_fr * sin(va_to - va_fr))
        )
    end

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
