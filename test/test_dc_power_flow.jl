using Test
using PowerSystems
import JuMP
import Gurobi
import LinearAlgebra

# system = System("data\\ELK14\\A5.m"); c1 = 1; c2 = 5
# system = SCOPF.System("data\\matpower\\IEEE_RTS.m"); c1 = 1; c2 = 11
system = SCOPF.System("data\\matpower\\ACTIVSg500.m"); c1 = 2; c2 = 1
opfm = SCOPF.scopf(SCOPF.SC, system, Gurobi.Optimizer)
SCOPF.solve_model!(opfm.mod)
nodes = SCOPF.get_nodes(opfm.sys) |> SCOPF.sort_components!
branches = SCOPF.get_branches(opfm.sys) |> SCOPF.sort_components!
idx = SCOPF.get_nodes_idx(nodes)
bx = SCOPF.get_bus_idx.(branches, [idx])
slack = SCOPF.find_slack(nodes)[1]

A = SCOPF.calc_A(branches, idx)
D = SCOPF.calc_D(branches)
DA = D*A
B = SCOPF.calc_B(A, DA)
X = SCOPF.calc_X(B, slack)
ϕ = SCOPF.calc_isf(D, A, X)
Pᵢ = SCOPF.get_net_Pᵢ(opfm)
θ = SCOPF.run_pf(B, Pᵢ)
θ .-= θ[slack]

pf = SCOPF.DCPowerFlow(nodes, branches, Pᵢ, idx)
@test θ ≈ pf.θ
@test D*A*θ ≈ SCOPF.calculate_line_flows!(pf, Pᵢ)
@test D*A*θ ≈ SCOPF.calc_Pline!(pf)

flow1 = copy(pf.F)
flow2 = copy(pf.F)
flow3 = copy(pf.F)

cont = bx[c1]
@btime SCOPF.calculate_line_flows!(flow1, pf, cont, c1)
@btime SCOPF.calculate_line_flows!(flow2, θ, B, pf.DA, pf.B, Pᵢ, cont, c1, pf.slack)
@btime SCOPF.get_isf!(ϕ, X, pf.X, pf.B, pf.DA, cont, c1)
@btime LinearAlgebra.mul!(flow3, ϕ, Pᵢ)
@test flow1 ≈ flow2
@test flow2 ≈ flow3

ΔPc = rand(length(Pᵢ))
@btime SCOPF.calculate_line_flows!(flow2, θ, B, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont, c1, pf.slack)
@btime SCOPF.get_isf!(ϕ, X, pf.X, pf.B, pf.DA, cont, c1)
@btime LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc))
@test flow2 ≈ flow3

cont = bx[c2]
islands, island, island_b = SCOPF.handle_islands(pf.B, pf.DA, cont, c2, pf.slack)
@btime SCOPF.calculate_line_flows!(flow2, θ, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont, c2, pf.slack, islands[island], island_b)
@btime SCOPF.get_isf!(ϕ, pf.DA, pf.B, cont, c2, pf.slack, islands[island], island_b)
@btime LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc))
@test flow2 ≈ flow3
