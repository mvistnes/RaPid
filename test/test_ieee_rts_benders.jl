# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

using PowerSystems
import JuMP
import Gurobi

# open("output.txt", "w") do out
#     redirect_stdout(out) do

system = SCOPF.System("data\\matpower\\IEEE_RTS.m")
# SCOPF.fix_generation_cost!(system);
nodes = SCOPF.sort_components!(SCOPF.get_nodes(system))
# idx = SCOPF.get_nodes_idx(nodes)
# active_capacities = zeros(length(nodes))
# reactive_capacities = zeros(length(nodes))
# cost = zeros(length(nodes))
# for g in SCOPF.get_generation(system)
#     if typeof(g) == HydroDispatch
#         SCOPF.set_operation_cost!(g, 30)
#     end
# end
# for g in SCOPF.get_generation(system)
#     n = idx[g.bus.number]
#     active_capacities[n] += g.active_power_limits.max
#     reactive_capacities[n] += g.reactive_power_limits.max
#     cost[n] = max(SCOPF.get_generator_cost(g)[2], cost[n])
# end
# PowerSystems.remove_component!.([system], SCOPF.get_generation(system))
# for (i, (ac, rc, c)) in enumerate(zip(active_capacities, reactive_capacities, cost))
#     if ac > 0.0
#         PowerSystems.add_component!(system, PowerSystems.ThermalStandard(string(i), true, true, 
#             nodes[i], ac, rc, sqrt(ac^2+rc^2), (min=0.0, max=ac), (min=-rc, max=rc), (up=0.01*ac, down=0.01*ac), 
#             PowerSystems.ThreePartCost(PowerSystems.VariableCost((0.0, c)), 0.0, 0.0, 0.0), 100.0))
#     end
# end
# SCOPF.set_ramp_limits!(system, 0.01)
SCOPF.set_rate!.(SCOPF.get_branches(system), SCOPF.get_rate.(SCOPF.get_branches(system)) * 0.8);
SCOPF.set_operation_cost!.(SCOPF.get_gens_h(system), [15.0, 16.0, 17.0, 18.0, 19.0, 20.0])

# voll, prob, contingencies = SCOPF.setup(system, 10, 40);
voll = [4304, 5098, 5245, 5419, 4834, 5585, 5785, 5192, 4575, 5244, 4478, 5698, 4465, 4859, 5032, 5256, 4598]
branches = SCOPF.sort_components!(SCOPF.get_branches(system));
# c = [7,12,13,21,22,23,27]
contingencies = branches
prob =
    [ # spesified for the RTS-96
        # 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
        # 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, # generators
        0.24, 0.51, 0.33, 0.39, 0.48, 0.38, 0.02, 0.36, 0.34, 0.33, 0.30, 0.44, 0.44,
        0.02, 0.02, 0.02, 0.02, 0.40, 0.39, 0.40, 0.52, 0.49, 0.38, 0.33, 0.41, 0.41,
        0.41, 0.35, 0.34, 0.32, 0.54, 0.35, 0.35, 0.38, 0.38, 0.34, 0.34, 0.45 # branches
    ]
# prob /= 8760
prob /= 100
# prob = prob[c]
short = 1.2
long = 1.0
ramp_minutes = 10
max_shed = 0.1
ramp_mult = 10

println("Start PTDF")
@time opfm_ptdf, pf_ptdf, Pc_ptdf, Pcc_ptdf = SCOPF.opf(SCOPF.PCSC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed,
    ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_limit_multi=short, long_term_limit_multi=long);
@time SCOPF.solve_model!(opfm_ptdf.mod);
println("Objective value PTDF: ", JuMP.objective_value(opfm_ptdf.mod))
SCOPF.print_contingency_P(opfm_ptdf, Pc_ptdf, Pcc_ptdf)
# @time opfm_norm = SCOPF.scopf(SCOPF.PCSC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed, 
#     ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_limit_multi=short);
# @time SCOPF.solve_model!(opfm_norm.mod);

println("Start Benders")
@time opfm, pf, Pc, Pcc, Pccx = SCOPF.run_benders(SCOPF.PCSC, system, Gurobi.Optimizer, voll, prob, contingencies, max_shed=max_shed,
    ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, branch_short_term_limit_multi=short, branch_long_term_limit_multi=long, p_failure=0.00);
println("Objective value Benders: ", JuMP.objective_value(opfm.mod))
SCOPF.print_contingency_P(opfm, Pc, Pcc, Pccx)

# SCOPF.print_corrective_results(opfm_norm)
# SCOPF.print_benders_results(opfm_ptdf, Pc_ptdf, Pcc_ptdf)
# SCOPF.print_benders_results(opfm, Pc, Pcc)

function change_slack(system, nodes, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long)
    slack = SCOPF.find_slack(nodes)[1]
    slacktype = BusTypes.PV
    res = []
    for n in eachindex(nodes)
        nodes[slack].bustype = slacktype
        slacktype = nodes[n].bustype
        nodes[n].bustype = BusTypes.REF
        slack = n
        println(SCOPF.find_slack(nodes)[1])
        opfm, pf, Pc, Pcc, Pccx = SCOPF.run_benders(SCOPF.PCSC, system, Gurobi.Optimizer, voll, prob, contingencies, max_shed=max_shed,
            ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, branch_short_term_limit_multi=short, branch_long_term_limit_multi=long, p_failure=0.00)
        push!(res, JuMP.objective_value(opfm.mod))
    end
    for r in res
        SCOPF.@printf "%12.9f\t" r
    end
end

# end
# end