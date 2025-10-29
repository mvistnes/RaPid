"""
    Extentions to the SCOPF with risk constraints, circuit breakers 
    and unit commitment.
"""


""" Add unit commitment to thermal generation (not hydro) """
function add_unit_commit!(opfm::OPFmodel)
    @variable(opfm.mod, u[g in get_name.(get_gens_t(opfm.sys))], Bin)
    delete_lower_bound.(opfm.mod[:pg0])
    delete_upper_bound.(opfm.mod[:pg0])
    p_lim = make_named_array(get_active_power_limits, get_gens_t(opfm.sys))
    @constraint(opfm.mod, ucu[g in get_name.(get_gens_t(opfm.sys))], opfm.mod[:pg0][g] .<= getindex.(p_lim,2)[g] .* u[g]) 
    @constraint(opfm.mod, ucd[g in get_name.(get_gens_t(opfm.sys))], opfm.mod[:pg0][g] .>= getindex.(p_lim,1)[g] .* u[g]) 
    add_to_expression!(objective_function(opfm.mod),
        sum(u[get_name(g)] * get_generator_cost(g)[end] for g in get_gens_t(opfm.sys)))
    @info "Unit commitment added to the base case"
    return opfm
end
""" Add unit commitment """
function add_unit_commit_ccont!(opfm::OPFmodel)
    opfm = add_unit_commit!(opfm)
    delete.(opfm.mod, opfm.mod[:pg_lim])
    unregister(opfm.mod, :pg_lim)
    p_lim = make_named_array(get_active_power_limits, get_ctrl_generation(opfm.sys))
    @constraint(opfm.mod, pg_lim_d[g = get_name.(get_gens_t(opfm.sys)), c = opfm.contingencies], 
            0 <= opfm.mod[:pg0][g] + (opfm.mod[:pgu][g,c] - opfm.mod[:pgd][g,c])
        )
    @constraint(opfm.mod, pg_lim_u[g = get_name.(get_gens_t(opfm.sys)), c = opfm.contingencies], 
            opfm.mod[:pg0][g] + (opfm.mod[:pgu][g,c] - opfm.mod[:pgd][g,c]) <= p_lim[g].max * opfm.mod[:u][g]
        )
    @constraint(opfm.mod, pg_lim[g = get_name.(get_gens_h(opfm.sys)), c = opfm.contingencies], 
            0 <= opfm.mod[:pg0][g] + (opfm.mod[:pgu][g,c] - opfm.mod[:pgd][g,c]) <= p_lim[g].max
        )
    @info "- Constricted to contingencies"
    return opfm
end

""" Add circuit breakers, then power flow on branch and branch limits for the base case and short-term contingency"""
function add_circuit_breakers_cont!(
            opfm::OPFmodel, 
            x::JuMP.Containers.DenseAxisArray, 
            branch_rating::JuMP.Containers.DenseAxisArray, 
            short_term_limit_multi::Float64 = 1.0
        )
    @variable(opfm.mod, cbc[l in get_name.(get_branches(opfm.sys)), c in opfm.contingencies], Bin) # circuit breakers on get_branches in in contingencies before 
    @constraint(opfm.mod, pfc_lim_l[l = get_name.(get_branches(opfm.sys)), c = opfm.contingencies], 
            -branch_rating[l] .* a(l,c) .* cbc[l,c] .* short_term_limit_multi .<= opfm.mod[:pfc][l,c]
        )
    @constraint(opfm.mod, pfc_lim_u[l = get_name.(get_branches(opfm.sys)), c = opfm.contingencies], 
            opfm.mod[:pfc][l,c] .<= branch_rating[l] .* a(l,c) .* cbc[l,c] .* short_term_limit_multi
        )
    opfm = add_c_branch!(opfm, x, branch_rating)
    @constraint(opfm.mod, pbc[l = get_name.(get_branches(opfm.sys)), c = opfm.contingencies],
            opfm.mod[:pfc][l,c] .- a(l,c) .* cbc[l,c] .* sum(beta(opfm.sys,l) .* opfm.mod[:vac][:,c]) ./ x[l] .== 0
        )
    @info "- After contingency, before corrective actions, with circuit breakers"
    return opfm
end
""" Add circuit breakers, then power flow on branch and branch limits for the base case and short-term and long-term contingency """
function add_circuit_breakers_ccont!(
            opfm::OPFmodel, 
            x::JuMP.Containers.DenseAxisArray, 
            branch_rating::JuMP.Containers.DenseAxisArray, 
            short_term_limit_multi::Float64 = 1.0
        )
    @variable(opfm.mod, cbcc[l in get_name.(get_branches(opfm.sys)), c in opfm.contingencies], Bin) # and after corrective actions
    @constraint(opfm.mod, pfcc_lim_l[l = get_name.(get_branches(opfm.sys)), c = opfm.contingencies], 
            -branch_rating[l] .* a(l,c) .* cbcc[l,c] .<= opfm.mod[:pfcc][l,c]
        )
    @constraint(opfm.mod, pfcc_lim_u[l = get_name.(get_branches(opfm.sys)), c = opfm.contingencies], 
            opfm.mod[:pfcc][l,c] .<= branch_rating[l] .* a(l,c) .* cbcc[l,c]
        )
    opfm = add_circuit_breakers_cont!(opfm, x, branch_rating, short_term_limit_multi)
    @constraint(opfm.mod, pbcc[l = get_name.(get_branches(opfm.sys)), c = opfm.contingencies],
            opfm.mod[:pfcc][l,c] .- a(l,c) .* cbcc[l,c] .* sum(beta(opfm.sys,l) .* opfm.mod[:vacc][:,c]) ./ x[l] .== 0
        )
    @info "- After contingency and corrective actions, with circuit breakers"
    return opfm
end
add_circuit_breakers_cont!(opfm::OPFmodel, short_term_limit_multi::Float64 = 1.0) = 
    add_circuit_breakers_cont!(opfm, make_named_array(get_x, get_branches(opfm.sys)), make_named_array(get_rating, get_branches(opfm.sys)), short_term_limit_multi)
add_circuit_breakers_ccont!(opfm::OPFmodel, short_term_limit_multi::Float64 = 1.0) = 
    add_circuit_breakers_ccont!(opfm, make_named_array(get_x, get_branches(opfm.sys)), make_named_array(get_rating, get_branches(opfm.sys)), short_term_limit_multi)

function add_line_risk_constraints(opfm::OPFmodel; 
            short_term_limit_multi::Float64 = 1.5, 
            Kᵣ = 1, 
            risk = 1,
            lim = 0.9
        )
    mult = 1/(1-lim)
    rate = make_named_array(get_rating, get_branches(opfm.sys))
    @variable(opfm.mod, 0 <= Sev[l in get_name.(get_branches(opfm.sys)), c in opfm.contingencies])
    @constraint(opfm.mod, sev_plim[l = get_name.(get_branches(opfm.sys)), c = opfm.contingencies], 
            Sev[l,c] >= mult * (opfm.mod[:pfc][l,c]/(rate[l] * short_term_limit_multi) - lim)
        )
    @constraint(opfm.mod, sev_nlim[l = get_name.(get_branches(opfm.sys)), c = opfm.contingencies], 
            Sev[l,c] >= mult * (-opfm.mod[:pfc][l,c]/(rate[l] * short_term_limit_multi) - lim)
        )
    @constraint(opfm.mod, line_risk, 
            sum(opfm.prob[c] * sum(Sev[l,c] for l in get_name.(get_branches(opfm.sys)))
                for c in contingencies)
            <= Kᵣ * risk
        )
    @info "Line risk constraints added"
    return opfm
end

function add_load_shedding_risk_constraints(opfm::OPFmodel; 
        Kᵣ = 1, 
        risk = 100,
        shed_lim = 10
    )
    @expression(opfm.mod, shed[c = opfm.contingencies], 
            sum(opfm.mod[:ls0][d] - opfm.mod[:lsc][d,c] for d in get_name.(get_demands(opfm.sys)))
        )
    @constraint(opfm.mod, shed_lim, shed .<= shed_lim)
    @constraint(opfm.mod, shed_risk, sum(opfm.prob[c] * shed[c] for c in contingencies) <= Kᵣ * risk)
    @info "Load shedding risk constraints added"
    return opfm
end
