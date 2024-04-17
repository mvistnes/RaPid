@testset "Test dc power flow" begin

# SETUP
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
# system = System(joinpath("data","ELK14","A5.m")); c1 = [1]; c2 = [5]
# system = SCOPF.System(joinpath("cases","IEEE_RTS.m")); c1 = [1]; c2 = [11]; c3 = [1,2]; c4 = [3,9]
system = SCOPF.System(joinpath("cases","RTS_GMLC.m")); c1 = [1]; c2 = [52]; c3 = [82, 106]; c4 = [90,111]
# system = SCOPF.System(joinpath("cases","ACTIVSg500.m")); c1 = [2]; c2 = [1]
# system = SCOPF.System(joinpath("cases","ACTIVSg2000.m")); c1 = [1]; c2 = [9]
# system = SCOPF.System(joinpath("cases","case_ACTIVSg10k.m")); c1 = [1]; c2 = [4]
# SCOPF.fix_generation_cost!(system);
voll = SCOPF.make_voll(system)
model, opf, pf, oplim, brc_up, brc_down, _, _, _ = SCOPF.opf_base(SCOPF.OPF(true, false, false, false, false), system, HiGHS.Optimizer(), voll=voll);
SCOPF.constrain_branches!(model, pf, oplim, brc_up, brc_down, 0.0)
nodes = SCOPF.sort_components!(SCOPF.get_nodes(system))
branches = SCOPF.sort_components!(SCOPF.get_branches(system))
bx = SCOPF.get_bus_idx.(branches, [SCOPF.get_nodes_idx(nodes)])
slack = SCOPF.find_slack(nodes)[1]

# CALC MATRICES DIRECTLY
A = SCOPF.calc_A(bx, size(opf.mbx,2))
D = SCOPF.calc_D(branches)
DA = D*A
B = SCOPF.calc_B(A, DA)
X = SCOPF.calc_X(B, slack)
ϕ = SCOPF.calc_isf(DA, X)
Pᵢ = SCOPF.get_value(model, :p0)
θ = B \ Pᵢ
θ .-= θ[slack]
F = DA*θ
@test θ ≈ SCOPF.calc_θ(X, Pᵢ)
@test F ≈ SCOPF.calculate_line_flows(ϕ, Pᵢ)
@test Pᵢ ≈ SCOPF.calc_Pᵢ(B, θ)

# CHECK DCPF-STRUCT CORRECTNESS
pf = SCOPF.DCPowerFlow(system, Pᵢ)
@test B ≈ pf.B
@test DA ≈ pf.DA
@test ϕ ≈ pf.ϕ
@test X ≈ pf.X
@test θ ≈ pf.θ
@test F ≈ SCOPF.calculate_line_flows!(pf, Pᵢ)
@test F ≈ SCOPF.calc_Pline!(pf)
@test θ ≈ SCOPF.calc_θ!(pf, Pᵢ)
@test Pᵢ ≈ SCOPF.calc_Pᵢ(pf)

K = SCOPF.calc_klu(B, pf.slack)   
# CONTAINERS FOR POWER FLOW
flow1 = copy(pf.F)
flow2 = copy(pf.F)
flow3 = copy(pf.F)
flow4 = copy(pf.F)

cont1 = bx[c1]
ΔPc = zeros(length(Pᵢ))
ΔPc[1] += 0.2
ΔPc[end] -= 0.3
θ₂ = pf.X * (Pᵢ .+ ΔPc)
F₂ = pf.DA * θ₂

# CONTINGENCY WITHOUT POWER INJECTION CHANGE OR ISLANDING
pfc = SCOPF.DCPowerFlow(nodes, branches[1:end .!= c1[1]], opf.idx)
SCOPF.calculate_line_flows!(flow1, θ, B, pf.DA, pf.B, Pᵢ, cont1, c1, pf.slack) # inverse with theta
@test (pfc.ϕ * Pᵢ) ≈ flow1[1:end .!= c1[1]]
@test flow1[c1[1]] ≈ 0.0
SCOPF.calc_isf!(ϕ, pf.DA, pf.B, cont1, c1, pf.slack); LinearAlgebra.mul!(flow2, ϕ, Pᵢ); # inverse with ptdf
@test flow1 ≈ flow2

# CONTINGENCY WITH POWER INJECTION CHANGE AND WITHOUT ISLANDING
SCOPF.calculate_line_flows!(flow1, θ, B, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont1, c1, pf.slack) # inverse with theta
@test (pfc.ϕ * (Pᵢ .+ ΔPc)) ≈ flow1[1:end .!= c1[1]]
@test flow1[c1[1]] ≈ 0.0
SCOPF.calc_isf!(ϕ, pf.DA, pf.B, cont1, c1, pf.slack); LinearAlgebra.mul!(flow2, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow1 ≈ flow2

# CONTINGENCY WITH POWER INJECTION CHANGE AND ISLANDING
cont2 = bx[c2]
islands, island_b = SCOPF.handle_islands(pf.B, pf.DA, cont2[1], c2[1])
pfc = SCOPF.DCPowerFlow(nodes[islands[1]], branches[island_b[1]], SCOPF.get_nodes_idx(nodes[islands[1]]))
SCOPF.calculate_line_flows!(flow1, θ, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont2, c2, pf.slack, islands[1], island_b[1]) # inverse with theta
@test (pfc.ϕ * (Pᵢ .+ ΔPc)[islands[1]]) ≈ flow1[island_b[1]]
@test flow1[c2[1]] ≈ 0.0
SCOPF.calc_isf!(ϕ, pf.DA, pf.B, cont2, c2, pf.slack, islands[1], island_b[1]); LinearAlgebra.mul!(flow2, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow1 ≈ flow2
SCOPF.calc_isf!(ϕ, pf.ϕ, islands[1], island_b[1]); LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc)); # hack
@test flow2 ≈ flow3 
SCOPF.calculate_line_flows!(flow4, ϕ, pf.ϕ, (Pᵢ .+ ΔPc), islands[1], island_b[1]); # hack2
@test flow3 ≈ flow4 

# CONTINGENCY WITH POWER INJECTION CHANGE, ISLANDING, AND DISTRIBUTED SLACK
SCOPF.set_dist_slack!(pfc.ϕ, opf.mgx[:,islands[1]], oplim.pg_lim_max)
SCOPF.set_dist_slack!(pf.ϕ, opf.mgx, oplim.pg_lim_max)
SCOPF.calculate_line_flows!(flow1, θ, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont2, c2, pf.slack, islands[1], island_b[1]) # inverse with theta
@test flow1[c2[1]] ≈ 0.0
SCOPF.calc_isf!(ϕ, pf.DA, pf.B, cont2, c2, pf.slack, islands[1], island_b[1]); 
    SCOPF.set_dist_slack!(ϕ[island_b[1],islands[1]], opf.mgx[:,islands[1]], oplim.pg_lim_max); 
    LinearAlgebra.mul!(flow2, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow1[island_b[1]] ≈ flow2[island_b[1]]
SCOPF.calc_isf!(ϕ, pf.ϕ, islands[1], island_b[1]); 
    SCOPF.set_dist_slack!(ϕ[island_b[1],islands[1]], opf.mgx[:,islands[1]], oplim.pg_lim_max); 
    LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc)); # hack
@test flow2[island_b[1]] ≈ flow3[island_b[1]] broken=true
SCOPF.calculate_line_flows!(flow4, ϕ, pf.ϕ, (Pᵢ .+ ΔPc), islands[1], island_b[1]); # hack2
@test (pfc.ϕ * (Pᵢ .+ ΔPc)[islands[1]]) ≈ flow1[island_b[1]] broken=true
@test flow3 ≈ flow4 

# CONTINGENCY WITHOUT POWER INJECTION CHANGE OR ISLANDING
cont3 = bx[c3]
pfc = SCOPF.DCPowerFlow(nodes, branches[setdiff(1:end, c3)], opf.idx)
SCOPF.calculate_line_flows!(flow1, θ, B, pf.DA, pf.B, Pᵢ, cont3, c3, pf.slack) # inverse with theta
@test pfc.ϕ * Pᵢ ≈ flow1[setdiff(1:end, c3)]
SCOPF.calc_isf!(ϕ, pf.DA, pf.B, cont3, c3, pf.slack); LinearAlgebra.mul!(flow2, ϕ, Pᵢ); # inverse with ptdf
@test flow1 ≈ flow2

# CONTINGENCY WITH POWER INJECTION CHANGE AND WITHOUT ISLANDING
SCOPF.calculate_line_flows!(flow1, θ, B, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont3, c3, pf.slack) # inverse with theta
SCOPF.calc_isf!(ϕ, pf.DA, pf.B, cont3, c3, pf.slack); LinearAlgebra.mul!(flow2, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow1 ≈ flow2

# CONTINGENCY WITH POWER INJECTION CHANGE AND ISLANDING
cont4 = bx[c4]
islands, islands_b = SCOPF.handle_islands(pf.B, pf.DA, cont4, c4)
pfc = SCOPF.DCPowerFlow(nodes[islands[1]], branches[islands_b[1]], SCOPF.get_nodes_idx(nodes[islands[1]]))
SCOPF.calculate_line_flows!(flow1, θ, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont4, c4, pf.slack, islands[1], islands_b[1]) # inverse with theta
@test pfc.ϕ * (Pᵢ .+ ΔPc)[islands[1]] ≈ flow1[islands_b[1]]
SCOPF.calc_isf!(ϕ, pf.DA, pf.B, cont4, c4, pf.slack, islands[1], islands_b[1]); LinearAlgebra.mul!(flow2, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow1[islands_b[1]] ≈ flow2[islands_b[1]]

end