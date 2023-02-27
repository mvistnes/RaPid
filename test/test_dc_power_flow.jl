using .SCOPF
using PowerSystems

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
@test θ ≈ pf.θ
@test D*A*θ ≈ calculate_line_flows(pf)