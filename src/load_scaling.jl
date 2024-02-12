function load_scale(type::OPF, system::System, optimizer;
    voll=Float64[],
    contingencies=Component[],
    prob=Float64[],
    dist_slack=Float64[],
    time_limit_sec::Int64=600,
    unit_commit::Bool=false,
    ramp_minutes::Real=10.0,
    ramp_mult::Real=10.0,
    max_shed=1.0,
    max_curtail=1.0,
    short_term_multi=1.5,
    long_term_multi=1.2,
    p_failure=0.0,
    silent=true,
    debug=false,
    atol=1e-10
)
    mod, opf, pf, oplim, Pc, Pcc, Pccx = opf_base(type, system, optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed,
    ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure)
    tot_t = constrain_branches!(mod, pf, oplim, 0.0)
    obj = objective_function(mod)
    areas = Dict(get_name.(PowerSystems.get_components(Area, system)) .=> Vector{StaticLoad}())
    for d in opf.demands
        areas[]
    pd = [get_active_power(d) for d in opf.demands if 
    while true
        ls0 = get_value(mod, mod[:ls0])
        ls = sum(ls0)
        ls < atol && break
        bit_ls = ls0 .> atol
        pd_ls = oplim.pd_lim[bit_ls]
        @constraint(mod, pd_ls, )
    end


    return Case(mod, opf, pf, oplim, Pc, Pcc, Pccx), tot_t
end