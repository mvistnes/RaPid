# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

# julia.NumThreads": "4" # On laptop with 4 physical cores

system = SCOPF.System("data\\matpower\\ACTIVSg500.m")
# system = System("data\\matpower\\ACTIVSg2000.m")
voll, prob, contingencies = SCOPF.setup(system, 100, 400);
SCOPF.fix_generation_cost!(system);
# nodes = SCOPF.sort_components!(SCOPF.get_nodes(system))
# active_capacities = zeros(length(nodes))
# reactive_capacities = zeros(length(nodes))
# cost = zeros(length(nodes))
# for g in SCOPF.get_generation(system)
#     if typeof(g) == HydroDispatch
#         SCOPF.set_operation_cost!(g, 30)
#     end
# end
# for g in SCOPF.get_ctrl_generation(system)
#     n = g.bus.number
#     active_capacities[n] += g.active_power_limits.max
#     reactive_capacities[n] += g.reactive_power_limits.max
#     cost[n] = max(SCOPF.get_generator_cost(g)[2], cost[n])
# end
# PowerSystems.remove_component!.([system], SCOPF.get_ctrl_generation(system))
# for (i, (ac, rc, c)) in enumerate(zip(active_capacities, reactive_capacities, cost))
#     if ac > 0.0
#         PowerSystems.add_component!(system, PowerSystems.ThermalStandard(string(i), true, true, 
#             nodes[i], ac, rc, sqrt(ac^2+rc^2), (min=0.0, max=ac), (min=-rc, max=rc), (up=0.01*ac, down=0.01*ac), 
#             PowerSystems.ThreePartCost(PowerSystems.VariableCost((0.0, c)), 0.0, 0.0, 0.0), 100.0))
#     end
# end
SCOPF.set_ramp_limits!(system, 0.01);
# SCOPF.set_renewable_prod!(system, 0.5)
# set_rate!.(SCOPF.get_branches(system), get_rate.(SCOPF.get_branches(system))*0.8);

# voll = fill(6000, length(SCOPF.get_demands(system)))
# branches = SCOPF.sort_components!(SCOPF.get_branches(system));
# c = [7,12,13,21,22,23,27]
# contingencies = branches[1:10:end]
# prob = 
#     [ # spesified for the RTS-96
#     # 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
#     # 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, # generators
#     0.24, 0.51, 0.33, 0.39, 0.48, 0.38, 0.02, 0.36, 0.34, 0.33, 0.30, 0.44, 0.44,  
#     0.02, 0.02, 0.02, 0.02, 0.40, 0.39, 0.40, 0.52, 0.49, 0.38, 0.33, 0.41, 0.41, 
#     0.41, 0.35, 0.34, 0.32, 0.54, 0.35, 0.35, 0.38, 0.38, 0.34, 0.34, 0.45 # branches
#     ]
# prob /= 8760
# prob = prob[c]
prob = fill(0.01, length(contingencies))
short = 1.2
long = 1.0
ramp_minutes = 10
max_shed = 0.1
ramp_mult = 10

@time opfm_norm = SCOPF.scopf(SCOPF.SC, system, SCOPF.Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=1.0, 
    ramp_minutes=ramp_minutes, short_term_limit_multi=short, debug=true);
@time SCOPF.solve_model!(opfm_norm.mod);
@time opfm, pf, Pc, Pcc, Pccx = SCOPF.run_benders(SCOPF.PCSC, system, SCOPF.Gurobi.Optimizer, voll, prob, contingencies, max_shed=max_shed, 
    ramp_minutes=ramp_minutes, branch_short_term_limit_multi=short, branch_long_term_limit_multi=long, p_failure=0.00);
@time opfm_ptdf, Pc_ptdf, Pcc_ptdf = SCOPF.opf(SCOPF.PCSC, system, SCOPF.Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed, 
    ramp_minutes=ramp_minutes, short_term_limit_multi=short, long_term_limit_multi=long);
@time SCOPF.solve_model!(opfm_ptdf.mod);
# opfm_norm = SCOPF.scopf(SCOPF.PCSC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=1.0, 
#     ramp_minutes=10, short_term_limit_multi=short, debug=true);
# SCOPF.solve_model!(opfm_norm.mod);