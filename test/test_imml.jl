using .SCOPF
using PowerSystems

# system = get_system("ELK14.json")
system = System("data\\ELK14\\A5.m")
opfm = scopf(SCOPF.SC, system, Gurobi.Optimizer)
solve_model!(opfm.mod)
nodes = get_sorted_nodes(opfm.sys)
branches = get_sorted_branches(opfm.sys)
idx = get_nodes_idx(nodes)
slack = find_slack(nodes)[1]

A = calc_A(branches, length(nodes), idx)
D = calc_D(branches)
B = calc_B(A, D)
X = calc_X(B, slack)
ϕ = calc_isf(A, D, X)
Pᵢ = get_net_Pᵢ(opfm, nodes, idx)
θ = run_pf(B, Pᵢ)

pf = SCOPF.DCPowerFlow(nodes, branches, idx, Pᵢ)
b1 = get_bus_idx(branches[1], idx)
@test calc_Pline(A, D, get_changed_angles(pf.X, pf.DA, pf.θ, b1[1], b1[2], slack)) ≈ 
    calculate_line_flows(calc_isf(D*A, get_changed_X(X, D*A, b1[1], b1[2])), Pᵢ) ≈
    calculate_line_flows(F, ptdf, X, B, angles, b1[1], b1[2], 1)

b3 = get_bus_idx(branches[3], idx)
@test get_changed_angles(pf.X, pf.DA, pf.θ, b2[1], b2[2], slack)
@test get_changed_X(X, D*A, b2[1], b2[2])
@test calculate_line_flows(F, ptdf, X, B, angles, b2[1], b2[2], 2)

@test θ ≈ pf.θ
@test D*A*θ ≈ calculate_line_flows(pf)

b5 = get_bus_idx(branches[5], idx)
@test_throws DivideError get_changed_anglespf.(X, pf.DA, pf.θ, b5[1], b5[2], slack)
@test_throws DivideError get_changed_X(X, D*A, b5[1], b5[2])
@test_throws DivideError calculate_line_flows(F, ptdf, X, B, angles, b5[1], b5[2], 2)
