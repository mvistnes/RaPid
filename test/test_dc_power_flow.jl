@testset "Test dc power flow" begin

# SETUP
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
# system = System(joinpath("data","ELK14","A5.m")); c1 = 1; c2 = 5
system = SCOPF.System(joinpath("cases","IEEE_RTS.m")); c1 = 1; c2 = 11
# system = SCOPF.System(joinpath("cases","RTS_GMLC.m")); c1 = 1; c2 = 52
# system = SCOPF.System(joinpath("cases","ACTIVSg500.m")); c1 = 2; c2 = 1
# system = SCOPF.System(joinpath("cases","ACTIVSg2000.m")); c1 = 2921; c2 = 9
# system = SCOPF.System(joinpath("cases","case_ACTIVSg10k.m")); c1 = 1; c2 = 4
SCOPF.fix_generation_cost!(system);
voll = SCOPF.make_voll(system)
model, opf, pf, oplim, brc_up, brc_down, _, _, _ = SCOPF.opf_base(SCOPF.OPF(true, false, false, false, false), system, HiGHS.Optimizer(), voll=voll);
SCOPF.constrain_branches!(model, pf, oplim, brc_up, brc_down, 0.0)
bx = SCOPF.get_bus_idx.(opf.branches, [opf.idx])
slack = SCOPF.find_slack(opf.nodes)[1]

# CALC MATRICES DIRECTLY
A = SCOPF.calc_A(opf.branches, length(opf.nodes), opf.idx)
D = SCOPF.calc_D(opf.branches)
DA = D*A
B = SCOPF.calc_B(A, DA)
X = SCOPF.calc_X(B, slack)
ϕ = SCOPF.calc_ptdf(DA, X)
Pᵢ = SCOPF.get_value(model, :p0)
θ = B \ Pᵢ
θ .-= θ[slack]
F = DA*θ
@test θ ≈ SCOPF.run_pf(pf.K, Pᵢ, pf.slack)

# CHECK DCPF-STRUCT CORRECTNESS
pf = SCOPF.DCPowerFlow(opf.nodes, opf.branches, opf.idx, Pᵢ)
@test B ≈ pf.B
@test DA ≈ pf.DA
@test ϕ ≈ pf.ϕ
@test X ≈ pf.X
@test θ ≈ pf.θ
@test F ≈ ϕ * Pᵢ
@test F ≈ SCOPF.calc_Pline!(pf)

K = SCOPF.calc_klu(B, pf.slack)   
# CONTAINERS FOR POWER FLOW
flow1 = copy(pf.F)
flow2 = copy(pf.F)
flow3 = copy(pf.F)
flow4 = copy(pf.F)

cont1 = bx[c1]
ΔPc = rand(length(Pᵢ))
θ₂ = SCOPF.run_pf!(similar(pf.θ), pf.K, (Pᵢ .+ ΔPc), pf.slack)
F₂ = pf.DA * θ₂

# CONTINGENCY WITHOUT POWER INJECTION CHANGE OR ISLANDING
SCOPF.calculate_line_flows!(flow1, θ, B, pf.DA, pf.B, Pᵢ, cont1, c1, pf.slack) # inverse with theta
SCOPF.calc_ptdf!(ϕ, K, pf.DA, pf.B, cont1, c1, pf.slack); LinearAlgebra.mul!(flow2, ϕ, Pᵢ); # inverse with ptdf
@test flow1 ≈ flow2

# CONTINGENCY WITH POWER INJECTION CHANGE AND WITHOUT ISLANDING
SCOPF.calculate_line_flows!(flow1, θ, B, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont1, c1, pf.slack) # inverse with theta
SCOPF.calc_ptdf!(ϕ, K, pf.DA, pf.B, cont1, c1, pf.slack); LinearAlgebra.mul!(flow2, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow1 ≈ flow2

# CONTINGENCY WITH POWER INJECTION CHANGE AND ISLANDING
cont2 = bx[c2]
islands, island, island_b = SCOPF.handle_islands(pf.B, pf.DA, cont2, c2, pf.slack)
SCOPF.calculate_line_flows!(flow1, θ, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont2, c2, pf.slack, islands[island], island_b) # inverse with theta
SCOPF.calc_ptdf!(ϕ, pf.DA, pf.B, cont2, c2, pf.slack, islands[island], island_b); LinearAlgebra.mul!(flow2, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow1 ≈ flow2
SCOPF.calc_ptdf!(ϕ, pf.ϕ, islands[island], island_b); LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc)); # hack
@test flow2 ≈ flow3 
SCOPF.calculate_line_flows!(flow4, ϕ, pf.ϕ, (Pᵢ .+ ΔPc), islands[island], island_b); # hack2
@test flow3 ≈ flow4 

SCOPF.set_dist_slack!(pf, opf)
SCOPF.calculate_line_flows!(flow1, θ, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont2, c2, pf.slack, islands[island], island_b) # inverse with theta
SCOPF.calc_ptdf!(ϕ, pf.DA, pf.B, cont2, c2, pf.slack, islands[island], island_b); LinearAlgebra.mul!(flow2, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow1 ≈ flow2
SCOPF.calc_ptdf!(ϕ, pf.ϕ, islands[island], island_b); LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc)); # hack
@test flow2 ≈ flow3 broken=true
SCOPF.calculate_line_flows!(flow4, ϕ, pf.ϕ, (Pᵢ .+ ΔPc), islands[island], island_b); # hack2
@test flow3 ≈ flow4 


# CONTINGENCY WITH POWER INJECTION CHANGE, ISLANDING, AND DISTRIBUTED SLACK
pfc = SCOPF.DCPowerFlow(opf.nodes[islands[island]], opf.branches[island_b], SCOPF.get_nodes_idx(opf.nodes[islands[island]]))
SCOPF.set_dist_slack!(pfc.ϕ, opf.mgx[:,islands[island]], oplim.pg_lim_max)
SCOPF.set_dist_slack!(pf.ϕ, opf.mgx, oplim.pg_lim_max)
SCOPF.calculate_line_flows!(flow1, θ, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont2, c2, pf.slack, islands[island], island_b) # inverse with theta
@test flow1[c2[1]] ≈ 0.0
SCOPF.calc_ptdf!(ϕ, pf.DA, pf.B, cont2, c2, pf.slack, islands[island], island_b); 
    SCOPF.set_dist_slack!(ϕ[island_b,islands[island]], opf.mgx[:,islands[island]], oplim.pg_lim_max); 
    LinearAlgebra.mul!(flow2, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow1[island_b] ≈ flow2[island_b]
SCOPF.calc_ptdf!(ϕ, pf.ϕ, islands[island], island_b); 
    SCOPF.set_dist_slack!(ϕ[island_b,islands[island]], opf.mgx[:,islands[island]], oplim.pg_lim_max); 
    LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc)); # hack
@test flow2[island_b] ≈ flow3[island_b] broken=true
SCOPF.calculate_line_flows!(flow4, ϕ, pf.ϕ, (Pᵢ .+ ΔPc), islands[island], island_b); # hack2
@test (pfc.ϕ * (Pᵢ .+ ΔPc)[islands[island]]) ≈ flow1[island_b] broken=true
@test flow3 ≈ flow4 

end