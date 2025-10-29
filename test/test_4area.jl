import Plots

system = System(joinpath("data","4-area","4_area.raw"))
set_active_power!.(SCOPF.sort_components!(get_components(StaticLoad, system)), 
    [258.81*0.425, 0*0.425, 323.52*0.425, 95.819*0.17, 336.15*0.17, 208.64*0.17, 0*1.19, 257.34*1.19, 331.03*1.19, 68.269*1.7] ./100)
# SCOPF.fix_generation_cost!(system);
# SCOPF.set_ramp_limits!(system, 1.0);

# # Remove area 4
# PowerSystems.remove_component!(StandardLoad, system, "load400241")
# PowerSystems.remove_component!(Bus, system, "24")
# PowerSystems.remove_component!(Line, system, "19-24-i_1")

voll = fill(6000, length(SCOPF.get_demands(system)))
branches = SCOPF.sort_components!(SCOPF.get_branches(system));
c = [9,10,18,19]
contingencies = branches[c]
set_x!.(contingencies[2:end], get_x.(contingencies[2:end]) * 1.5)
for g in SCOPF.get_generation(system)
    if g.bus.area.name == "1"
        SCOPF.set_operation_cost!(g, 30)
    elseif g.bus.area.name == "2"
        SCOPF.set_operation_cost!(g, 40)
    elseif g.bus.area.name == "3"
        SCOPF.set_operation_cost!(g, 50)
    end
end
prob = fill(0.01, length(contingencies))
short = 1.25
long = 1.0

opfm_norm = SCOPF.scopf(SCOPF.SC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, short_term_limit_multi=short);
SCOPF.solve_model!(opfm_norm.mod);
opfm, pf, Pc, Pcc, Pccx = SCOPF.run_decomposed_optimization(SCOPF.PCSC, system, voll, prob, contingencies, max_shed=1.0, 
       ramp_minutes=0.5, branch_short_term_limit_multi=short, branch_long_term_limit_multi=long, p_failure=0.01);

i0 = get_rating(contingencies[2])
x = []
y = []
cost = []
for i in 0.5:0.05:1.0
    println("Rating ", i*i0)
    # set_rating!(contingencies[2], i*i0)
    global opfm, pf, Pc, Pcc, Pccx = SCOPF.run_decomposed_optimization(SCOPF.PCSC, system, voll, prob, contingencies, max_shed=1.0, 
        ramp_minutes=0.5, branch_short_term_limit_multi=short, branch_long_term_limit_multi=long, p_failure=0.00, branch_c=contingencies[2], rate_c=i*i0);
    # SCOPF.print_power_flow(opfm)
    # global opfm_norm = SCOPF.scopf(SCOPF.SC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, short_term_limit_multi=short);
    # set_rating!(contingencies[2], i0)
    # SCOPF.solve_model!(opfm_norm.mod);
    if JuMP.termination_status(opfm.mod) == JuMP.MOI.OPTIMAL
        mod = corrective_shedding(opfm, Gurobi.Optimizer, long_term_limit_multi=long)
        SCOPF.solve_model!(mod);
        push!(y, JuMP.objective_value(opfm.mod))
        # push!(y, JuMP.objective_value(opfm_norm.mod))
        push!(cost, JuMP.objective_value(mod))
        push!(x, JuMP.value(opfm.mod[:pf0][10]))
    end
end
println(y)
Plots.plot([(0.5:0.05:1.0) * i0],y, title="Rate vs objective value")
Plots.plot([(0.5:0.05:1.0) * i0], 1 .* cost .+ y, title="Rate vs objective value + corrective_shedding")
Plots.plot(x,y, title="Flow vs objective value")
set_rating!(contingencies[2], i0)

# x = []
# for i in [1, 2, 5, 10]    
#     println("Prob ", 0.01*i)
#     prob = fill(0.01*i, length(branches))
#     opfm, pf, Pc, Pcc = SCOPF.run_decomposed_optimization(SCOPF.PCSC, system, voll, prob, contingencies, max_shed=1.0, 
#        ramp_minutes=0.5, branch_short_term_limit_multi=1.2, branch_long_term_limit_multi=1.0);
#     push!(x, JuMP.value.(opfm.mod[:pf0]))
# end
# prob = fill(0.01, length(branches))
# for i in axes(x[1],1)
#     for j in axes(x,1)
#         SCOPF.@printf "%5.2f " x[j][i]
#     end
#     println("")
# end

function corrective_shedding(opfm::SCOPF.OPFmodel, optimizer;
        long_term_limit_multi::Real = 1.5,
        ramp_minutes::Real = 10
    )
    mod = JuMP.Model(optimizer)
    idx = SCOPF.get_nodes_idx(opfm.nodes)
    list = SCOPF.make_list(opfm, idx, opfm.nodes)

    slack = SCOPF.find_slack(opfm.nodes)
    islands = SCOPF.get_all_islands(opfm, slack[1]) 
    
    JuMP.@variables(mod, begin
        0 <= pgd[g in 1:length(opfm.ctrl_generation), c in 1:length(opfm.contingencies)]
            # and ramp down
        pfcc[l in 1:length(opfm.branches), c in 1:length(opfm.contingencies)]
            # power flow on branches in in contingencies after corrective actions
        vacc[b in 1:length(opfm.nodes), c in 1:length(opfm.contingencies)]
            # voltage angle at a node in in contingencies after corrective actions
        lscc[d in 1:length(opfm.demands), c in 1:length(opfm.contingencies)]
            # load curtailment variables in in contingencies
    end)
    JuMP.@objective(mod, Min,
        sum(opfm.prob[c] * (sum(opfm.voll' * lscc[:,c]) + sum(60 * pgd[:,c])) for c in 1:length(opfm.contingencies))
    )
    for (c,island) in enumerate(islands)
        isempty(island) && continue
        for n in island
            JuMP.@constraint(mod,
                    sum(SCOPF.beta(opfm.nodes[n], opfm.branches[l]) * pfcc[l,c] for l in list[n].branches) == 
                    sum(JuMP.value(opfm.mod[:pg0][g]) - pgd[g,c] for g in list[n].ctrl_generation) -
                    sum(PowerSystems.get_active_power(opfm.demands[d]) * 1.0 - lscc[d,c] for d in list[n].demands)
                )
        end
    end
    JuMP.@constraint(mod, ref_vacc, vacc[slack[1],:] .== 0) # Set voltage angle at reference bus
    branch_rating = PowerSystems.get_rating.(opfm.branches) * long_term_limit_multi
    JuMP.@constraint(mod, pfcc_lim[l = 1:length(opfm.branches), c = 1:length(opfm.contingencies)], 
        -branch_rating[l] * !isequal(opfm.branches[l],opfm.contingencies[c]) <= pfcc[l,c] <= branch_rating[l] * !isequal(opfm.branches[l],opfm.contingencies[c])
    )
    
    JuMP.@constraint(mod, pbcc[l = 1:length(opfm.branches), c = 1:length(opfm.contingencies)],
        pfcc[l,c] - !isequal(opfm.branches[l],opfm.contingencies[c]) * (vacc[idx[opfm.branches[l].arc.from.number],c] - 
        vacc[idx[opfm.branches[l].arc.to.number],c]) / get_x(opfm.branches[l]) == 0
    )

    listPd = PowerSystems.get_active_power.(opfm.demands)
    for (c,island) in enumerate(islands)
        for n in island
            JuMP.@constraint(mod, [d in list[n].demands], 0 <= lscc[d,c] <= listPd[d])
        end
        if length(island) < length(opfm.nodes)
            for n in setdiff(1:length(opfm.nodes), island)
                JuMP.@constraint(mod, [d in list[n].demands], lscc[d,c] == listPd[d])
            end
        end
    end
    
    p_lim = getindex.(PowerSystems.get_active_power_limits.(opfm.ctrl_generation), :max)
    (rampup, rampdown) = SCOPF.split_pair(SCOPF.get_ramp_limits.(opfm.ctrl_generation))
    for (c,island) in enumerate(islands)
        if length(island) < length(opfm.nodes)
            for n in island
                JuMP.@constraint(mod, [g in list[n].ctrl_generation], pgd[g,c] .<= rampdown * ramp_minutes)
                JuMP.@constraint(mod, [g in list[n].ctrl_generation], 
                    0 <= JuMP.value(opfm.mod[:pg0][g]) - pgd[g,c] <= p_lim[g])
            end
            for n in setdiff(1:length(opfm.nodes), island)
                JuMP.@constraint(mod, [g in list[n].ctrl_generation], 
                    JuMP.value(opfm.mod[:pg0][g]) == pgd[g,c])
            end
        else
            JuMP.@constraint(mod, pgd[:,c] .<= rampdown * ramp_minutes)
            JuMP.@constraint(mod, 
                0 .<= JuMP.value.(opfm.mod[:pg0]) .- pgd[:,c]
            )
        end
    end

    return mod
end
