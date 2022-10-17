using PowerSystems
using JuMP
using Ipopt
using GLPK
using Printf

""" OPF model type """
mutable struct OPFmodel
    sys::System
    mod::Model # JuMP model object
    # beta::JuMP.Containers.SparseAxisArray
    contingencies::Union{Nothing, Vector{String}}
    voll::Union{Nothing, JuMP.Containers.DenseAxisArray}
    prob::Union{Nothing, JuMP.Containers.DenseAxisArray}
end

""" Constructor for OPFmodel """
function opfmodel(sys::System, optimizer, time_limit_sec)
    mod = Model(optimizer)
    if GLPK.Optimizer == optimizer 
        set_optimizer_attribute(mod, "msg_lev", GLPK.GLP_MSG_ON)
    end
    set_time_limit_sec(mod, time_limit_sec)
    
    # β = Dict{Tuple{String,String},Int8}()
    # for b in nodes(sys)
    #     for l in branches(sys)
    #         β[get_name(b), get_name(l)] = beta(b,l)
    #     end
    # end
    # β = JuMP.Containers.SparseAxisArray(β)

    return OPFmodel(sys, mod, nothing, nothing, nothing)
end

"""Find value of β, β = 1 if from-bus is the bus, β = -1 if to-bus is the bus, 0 else"""
beta(bus::Bus, branch::Branch) = bus == branch.arc.from ? 1 : bus == branch.arc.to ? -1 : 0
function beta(bus::String, sys::System) 
    b = get_component(Bus, sys, bus)
    return [beta(b, l) for l in branches(sys)]
end
function beta(sys::System, branch::String) 
    l = get_component(ACBranch, sys, branch)
    return [beta(b, l) for b in nodes(sys)]
end

""" An iterator to a type of power system component """
gens_t(sys::System) = get_components(ThermalGen, sys)
gens_h(sys::System) = get_components(HydroGen, sys)
branches(sys::System) = get_components(ACBranch, sys)
nodes(sys::System) = get_components(Bus, sys)
demands(sys::System) = get_components(StaticLoad, sys)
renewables(sys::System) = get_components(RenewableGen, sys) # Renewable modelled as negative demand
get_ctrl_generation(sys::System) = Iterators.flatten((gens_t(sys), gens_h(sys))) # An iterator of all controllable generators
get_nonctrl_generation(sys::System) = Iterators.flatten((demands(sys), renewables(sys))) # An iterator of all non-controllable load and generation

# it_name(::Type{T}, mos::Model) where {T <: Component} = get_name.(get_components(T,mod))

""" A Vector of the Value Of Lost Load for the demand and renewables """
make_voll(OPFm::OPFmodel) = JuMP.Containers.DenseAxisArray(
        [rand(1000:3000, length(demands(OPFm.sys))); rand(1:30, length(renewables(OPFm.sys)))], 
        [get_name.(demands(OPFm.sys)); get_name.(renewables(OPFm.sys))]
    )

make_prob(OPFm::OPFmodel) = JuMP.Containers.DenseAxisArray(
        (rand(length(OPFm.contingencies)).*(0.5-0.02).+0.02)./8760, 
        OPFm.contingencies
    )

""" Initialize variables for dc SCOPF """
function init_var_dc_SCOPF!(OPFm::OPFmodel)
    p_lim = JuMP.Containers.DenseAxisArray(
        [get_active_power_limits(g).max for g in get_ctrl_generation(OPFm.sys)],
        get_name.(get_ctrl_generation(OPFm.sys))
    )
    demand = JuMP.Containers.DenseAxisArray(
        [get_active_power(d) for d in demands(OPFm.sys)],
        get_name.(get_nonctrl_generation(OPFm.sys))
    )
    @variables(OPFm.mod, begin
        0 <= pg0[g in get_name.(get_ctrl_generation(OPFm.sys))] <= p_lim[g]     # active power variables for the generators
            # restricted to positive numbers (can be changed to minimum active power) and less than maximum active power
        pf0[l in get_name.(branches(OPFm.sys))]      # power flow on branches in base case
        va0[b in get_name.(nodes(OPFm.sys))]         # voltage angle at a node in base case
        0 <= ls0[d in get_name.(get_nonctrl_generation(OPFm.sys))] <= demand[d]    # demand curtailment variables
            # restrict load shedding to max load per node
    end)
    return OPFm
end
""" Initialize variables for dc P-SCOPF """
function init_var_dc_P_SCOPF!(OPFm::OPFmodel)
    OPFm = init_var_dc_SCOPF!(OPFm)
    # power flow on branches in contingencies
    @variable(OPFm.mod, pfc[l in get_name.(branches(OPFm.sys)), c in OPFm.contingencies])
    # voltage angle at a node in contingencies
    @variable(OPFm.mod, vac[b in get_name.(nodes(OPFm.sys)), c in OPFm.contingencies])
    return OPFm
end
""" Initialize variables for dc PC-SCOPF """
function init_var_dc_PC_SCOPF!(OPFm::OPFmodel)
    OPFm = init_var_dc_SCOPF!(OPFm)
    @variables(OPFm.mod, begin
        # pgc[g in get_name.(get_ctrl_generation(OPFm.sys)), c in OPFm.contingencies] >= 0    # active power variables for the generators in contingencies 
        pgu[g in get_name.(get_ctrl_generation(OPFm.sys)), c in OPFm.contingencies] >= 0    # active power variables for the generators in contingencies ramp up 
        pgd[g in get_name.(get_ctrl_generation(OPFm.sys)), c in OPFm.contingencies] >= 0       # and ramp down
        pfc[l in get_name.(branches(OPFm.sys)), c in OPFm.contingencies]      # power flow on branches in in contingencies before 
        pfcc[l in get_name.(branches(OPFm.sys)), c in OPFm.contingencies]         # and after corrective actions
        vac[b in get_name.(nodes(OPFm.sys)), c in OPFm.contingencies]         # voltage angle at a node in in contingencies before 
        vacc[b in get_name.(nodes(OPFm.sys)), c in OPFm.contingencies]            # and after corrective actions
        lsc[d in get_name.(get_nonctrl_generation(OPFm.sys)), c in OPFm.contingencies] >= 0 # load curtailment variables in in contingencies
        # cbc[l in get_name.(branches(OPFm.sys)), c in OPFm.contingencies], Bin      # circuit breakers on branches in in contingencies before 
        # cbcc[l in get_name.(branches(OPFm.sys)), c in OPFm.contingencies], Bin         # and after corrective actions
    end)
    return OPFm
end

function add_obj!(OPFm::OPFmodel)
    @objective(OPFm.mod, Min, 
        sum(OPFm.mod[:pg0][get_name(g)] * (get_cost ∘ get_variable ∘ get_operation_cost)(g)[2]
        for g in gens_t(OPFm.sys)) + 
        sum(OPFm.mod[:pg0][get_name(g)] * 5 for g in gens_h(OPFm.sys)) + 
        sum(OPFm.voll[d] * OPFm.mod[:ls0][d] for d in get_name.(get_nonctrl_generation(OPFm.sys)))
    )
    return OPFm
end
function add_obj_corrective!(OPFm::OPFmodel)
    @objective(OPFm.mod, Min, 
        sum((OPFm.mod[:pg0][get_name(g)] * (get_cost ∘ get_variable ∘ get_operation_cost)(g)[2] + 
            # sum(OPFm.prob[c] * (OPFm.mod[:pgc][get_name(g),c] - OPFm.mod[:pg0][get_name(g)]) for c in contingencies))
            sum(OPFm.prob[c] * (OPFm.mod[:pgu][get_name(g),c] #=+ OPFm.mod[:pgd][get_name(g),c]*0.1=#) for c in OPFm.contingencies))
            for g in gens_t(OPFm.sys)
        ) + 
        # sum(5 * (OPFm.mod[:pg0][get_name(g)] + sum(OPFm.prob[c] * (OPFm.mod[:pgc][get_name(g),c] - OPFm.mod[:pg0][get_name(g)])
        sum(5 * (OPFm.mod[:pg0][get_name(g)] + sum(OPFm.prob[c] * (OPFm.mod[:pgu][get_name(g),c] #=+ OPFm.mod[:pgd][get_name(g),c]*0.1=#)
            for c in OPFm.contingencies) ) for g in gens_h(OPFm.sys)
        ) +
        sum(OPFm.voll[d] * (OPFm.mod[:ls0][d] + sum(OPFm.prob[c] * OPFm.mod[:lsc][d,c] for c in OPFm.contingencies)) 
            for d in get_name.(get_nonctrl_generation(OPFm.sys)))
    )
    return OPFm
end

""" Set voltage angle at reference bus """
function set_ref_angle!(OPFm::OPFmodel)
    slack = find_slack(OPFm.sys)
    @constraint(OPFm.mod, OPFm.mod[:va0][get_name(slack)] == 0)
    return OPFm
end

""" Incerted power at each bus for the base case """
function add_c_bus!(OPFm::OPFmodel)
    @expression(OPFm.mod, ctrl[b = get_name.(nodes(OPFm.sys))], 
        sum(g.bus.name == b ? OPFm.mod[:pg0][get_name(g)] : 0 for g in get_ctrl_generation(OPFm.sys)))
    @expression(OPFm.mod, nonctrl[b = get_name.(nodes(OPFm.sys))], 
        sum((d.bus.name == b ? get_active_power(d) .- OPFm.mod[:ls0][get_name(d)] : 0 for d in demands(OPFm.sys)), init = 0.0) + 
        sum((d.bus.name == b ? -get_active_power(d) .+ OPFm.mod[:ls0][get_name(d)] : 0 for d in renewables(OPFm.sys)), init = 0.0))
    
    for n in nodes(OPFm.sys)
        inj = @constraint(OPFm.mod, ctrl[get_name(n)] - sum(beta(n, l) * OPFm.mod[:pf0][get_name(l)] for l in branches(OPFm.sys)) == nonctrl[get_name(n)])
        set_name(inj, @sprintf("inj[%s]", get_name(n)))
    end
    @constraint(OPFm.mod, va_lim, -π/2 .<= OPFm.mod[:va0] .<= π/2)
    return OPFm
end
function add_c_bus_cont!(OPFm::OPFmodel)
    OPFm = add_c_bus!(OPFm)
    for n in nodes(OPFm.sys)
        for c in OPFm.contingencies
            inj = @constraint(OPFm.mod, OPFm.mod[:ctrl][get_name(n)] - sum(beta(n,l) * OPFm.mod[:pfc][get_name(l),c] for l in branches(OPFm.sys)) == OPFm.mod[:nonctrl][get_name(n)])
            set_name(inj, @sprintf("inj_p[%s,%s]", get_name(n), c))
        end
    end
    @constraint(OPFm.mod, vac_lim, -π/2 .<= OPFm.mod[:vac] .<= π/2)
    return OPFm
end
function add_c_bus_ccont!(OPFm::OPFmodel)
    OPFm = add_c_bus_cont!(OPFm)
    for n in nodes(OPFm.sys)
        for c in OPFm.contingencies
            inj = @constraint(OPFm.mod, 
            # sum(get_bus(g) == n ? pgc[get_name(g),c] : 0 for g in get_ctrl_generation(OPFm)) -
            sum(get_bus(g) == n ? OPFm.mod[:pg0][get_name(g)] + OPFm.mod[:pgu][get_name(g),c] - OPFm.mod[:pgd][get_name(g),c] : 0 for g in get_ctrl_generation(OPFm.sys)) -
            sum(beta(n,l) * OPFm.mod[:pfcc][get_name(l),c] for l in branches(OPFm.sys)) == 
            sum(get_bus(d) == n ? get_active_power(d) - OPFm.mod[:lsc][get_name(d),c] : 0 for d in demands(OPFm.sys)) +
            sum(get_bus(d) == n ? -get_active_power(d) + OPFm.mod[:lsc][get_name(d),c] : 0 for d in renewables(OPFm.sys))
            )
            set_name(inj, @sprintf("inj_pcc[%s,%s]", get_name(n), c))
        end
    end
    @constraint(OPFm.mod, vacc_lim, -π/2 .<= OPFm.mod[:vacc] .<= π/2)
    return OPFm
end

""" Power flow on branch and branch limits for the base case """
function add_c_branch!(OPFm::OPFmodel)
    for l in branches(OPFm.sys)
        @constraint(OPFm.mod, OPFm.mod[:pf0][get_name(l)] - sum(beta(n, l) * OPFm.mod[:va0][get_name(n)] for n in nodes(OPFm.sys)) / get_x(l) == 0)
        @constraint(OPFm.mod, -get_rate(l) <= OPFm.mod[:pf0][get_name(l)] <= get_rate(l))
    end
    return OPFm
end
function add_c_branch_cont!(OPFm::OPFmodel, short_term_limit_multi::Float64 = 1.0)
    OPFm = add_c_branch!(OPFm)
    for l in branches(OPFm.sys)
        l_name = get_name(l)
        for c in OPFm.contingencies
            a = l_name != c # zero value if l is unavailable under contingency c
            @constraint(OPFm.mod, OPFm.mod[:pfc][l_name,c] - a * sum(beta(b,l) * OPFm.mod[:vac][get_name(b),c] for b in nodes(OPFm.sys)) / get_x(l) == 0)
            @constraint(OPFm.mod, -get_rate(l) * short_term_limit_multi <= a * OPFm.mod[:pfc][l_name,c] <= get_rate(l) * short_term_limit_multi)
        end
    end
    return OPFm
end
function add_c_branch_ccont!(OPFm::OPFmodel, short_term_limit_multi::Float64)
    OPFm = add_c_branch_cont!(OPFm, short_term_limit_multi)
    for l in branches(OPFm.sys)
        l_name = get_name(l)
        for c in OPFm.contingencies
            a = l_name != c # zero value if l is unavailable under contingency c
            @constraint(OPFm.mod, OPFm.mod[:pfcc][l_name,c] - a * sum(beta(b,l) * OPFm.mod[:vacc][get_name(b),c] for b in nodes(OPFm.sys)) / get_x(l) == 0)
            @constraint(OPFm.mod, -get_rate(l) <= a * OPFm.mod[:pfcc][l_name,c] <= get_rate(l))
        end
    end
    return OPFm
end

""" Restrict active power generation to a minimum level """
function add_min_P_gen!(OPFm::OPFmodel)
    for g in get_ctrl_generation(OPFm.sys)
        set_lower_bound(OPFm.mod.pg0[get_name(g)], get_active_power_limits(g).min) 
    end
    return OPFm
end

""" Restrict active power generation ramp to min and max values """
function add_lim_ramp_P_gen!(OPFm::OPFmodel, ramp_minutes)
    for g in get_ctrl_generation(OPFm.sys)
        g_name = get_name(g)
        for c in OPFm.contingencies
            # @constraint(OPFm.mod, 0 <= OPFm.mod[:pg0][g_name] <= get_active_power_limits(g).max)
            # @constraint(OPFm.mod, 0 <= OPFm.mod[:pgc][g_name,c] <= get_active_power_limits(g).max)
            # @constraint(OPFm.mod, -get_ramp_limits(g).down * ramp_minutes <= 
            #     OPFm.mod[:pgc][g_name,c] - OPFm.mod[:pg0][g_name] <= get_ramp_limits(g).up * ramp_minutes
            # )
            @constraint(OPFm.mod, 0 <= OPFm.mod[:pg0][g_name] + (OPFm.mod[:pgu][g_name,c] - OPFm.mod[:pgd][g_name,c]) <= get_active_power_limits(g).max)
            set_upper_bound(OPFm.mod[:pgu][g_name,c], get_ramp_limits(g).up * ramp_minutes)
            set_upper_bound(OPFm.mod[:pgd][g_name,c], get_ramp_limits(g).down * ramp_minutes)
        end
    end
    return OPFm
end

""" Restrict load shedding to load (or renewable production) per node """
function add_lim_nonctrl_shed!(OPFm::OPFmodel)
    for d in get_nonctrl_generation(OPFm.sys), c in OPFm.contingencies
        @constraint(OPFm.mod, OPFm.mod[:ls0][get_name(d)] + OPFm.mod[:lsc][get_name(d),c] <= get_active_power(d))
    end
    return OPFm
end

""" Add unit commitment """
function add_unit_commit!(OPFm::OPFmodel)
    @variable(OPFm.mod, u[g in get_name.(get_ctrl_generation(OPFm.sys))], Bin)
    delete_lower_bound.(OPFm.mod.pg0)
    delete_upper_bound.(OPFm.mod.pg0)
    p_lim = JuMP.Containers.DenseAxisArray(
        [get_active_power_limits(g) for g in get_ctrl_generation(OPFm.sys)],
        get_name.(get_ctrl_generation(OPFm.sys))
    )
    @constraint(uc, getindex.(p_lim,1) * u <= OPFm.mod[:pg0] <= getindex.(p_lim,2) * u) 
    JuMP.add_to_expression!(OPFm.mod.obj, sum(u[get_name(g)] * (get_cost ∘ get_variable ∘ get_operation_cost)(g)[3] for g in get_ctrl_generation(OPFm.sys)))
    return OPFm
end

""" Return the (first) slack bus in the system. """
function find_slack(sys::System)
    for x in nodes(sys)
        if x.bustype == BusTypes.REF
            return x
        end
    end
end


""" Run a SCOPF of a power system """
function scopf(system::System, optimizer; time_limit_sec = 60, voll = [0])
    opfm = opfmodel(system, optimizer, time_limit_sec)
    opfm.voll = length(voll) > 1 ? voll : make_voll(opfm)
    opfm = init_var_dc_SCOPF!(opfm) |> set_ref_angle! |> add_c_bus! |> add_c_branch! |> add_obj!
    optimize!(opfm.mod)
    @assert termination_status(opfm.mod) == MOI.OPTIMAL
    return opfm.mod
end

""" Run a P-SCOPF of a power system """
function p_scopf(system::System, optimizer; time_limit_sec = 60, voll = [0])
    opfm = opfmodel(system, optimizer, time_limit_sec)
    opfm.voll = length(voll) > 1 ? voll : make_voll(opfm)
    opfm.contingencies = get_name.(branches(opfm.sys)) # [1:10] # [get_name.(nodes(opfm.sys));get_name.(branches(opfm.sys))]
    opfm = init_var_dc_P_SCOPF!(opfm) |> set_ref_angle! |> add_c_bus_cont! |> add_c_branch_cont! |> add_obj!

    optimize!(opfm.mod)
    @assert termination_status(opfm.mod) == MOI.OPTIMAL
    return opfm.mod
end

""" Run a PC-SCOPF of a power system """
function pc_scopf(system::System, optimizer; time_limit_sec = 60, voll = [0], prob = [0])
    ramp_minutes = 10
    short_term_limit_multi = 1.5
    opfm = opfmodel(system, optimizer, time_limit_sec)
    opfm.contingencies = get_name.(branches(opfm.sys)) # [1:10] # [get_name.(nodes(opfm.sys));get_name.(branches(opfm.sys))]
    opfm.voll = length(voll) > 1 ? voll : make_voll(opfm)
    opfm.prob = length(prob) > 1 ? prob : make_prob(opfm)
    opfm = init_var_dc_PC_SCOPF!(opfm) |> set_ref_angle! |> add_c_bus_ccont! |> add_lim_nonctrl_shed! |> add_obj_corrective!
    opfm = add_c_branch_ccont!(opfm, short_term_limit_multi) 
    opfm = add_lim_ramp_P_gen!(opfm, ramp_minutes)

    optimize!(opfm.mod)
    @assert termination_status(opfm.mod) == MOI.OPTIMAL
    return opfm.mod
end

