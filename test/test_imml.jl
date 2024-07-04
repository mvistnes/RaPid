@testset "Test IMML" begin

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
model, opf, pf, oplim, brc_up, brc_down, _, _, _ = SCOPF.opf_base(SCOPF.Base_SCOPF, system, HiGHS.Optimizer(), voll=voll);
SCOPF.constrain_branches!(model, pf, oplim, brc_up, brc_down, 0.0)
bx = SCOPF.get_bus_idx.(opf.branches, [opf.idx])
slack = SCOPF.find_slack(opf.nodes)[1]

Pᵢ = SCOPF.get_value(model, :p0)
pf = SCOPF.DCPowerFlow(opf.nodes, opf.branches, opf.idx, Pᵢ)
B = copy(pf.B)
X = copy(pf.X)
ϕ = copy(pf.ϕ)
SCOPF.update_model!(model, pf, 0.0)
θ = copy(pf.θ)
F = copy(pf.F)
K = SCOPF.calc_klu(B, pf.slack)
  
# CONTAINERS FOR POWER FLOW
flow1 = copy(pf.F)
flow2 = copy(pf.F)
flow3 = copy(pf.F)
flow4 = copy(pf.F)
flow5 = copy(pf.F)
flow6 = copy(pf.F)
flow7 = copy(pf.F)

cont1 = bx[c1]
ΔPc = rand(length(Pᵢ))
θ₂ = SCOPF.run_pf!(similar(pf.θ), pf.K, (Pᵢ .+ ΔPc), pf.slack)
F₂ = pf.DA * θ₂

# CONTINGENCY WITHOUT POWER INJECTION CHANGE OR ISLANDING
SCOPF.calculate_line_flows!(flow1, pf, cont1, c1) # IMML flow
SCOPF.calculate_line_flows!(flow2, θ, B, pf.DA, pf.B, Pᵢ, cont1, c1, pf.slack) # inverse with theta
@test flow1 ≈ flow2
SCOPF.calc_isf!(ϕ, X, pf.X, pf.B, pf.DA, cont1, c1); LinearAlgebra.mul!(flow3, ϕ, Pᵢ); # IMML ptdf
@test flow2 ≈ flow3
SCOPF.calc_isf!(ϕ, K, pf.DA, pf.B, cont1, c1, pf.slack); LinearAlgebra.mul!(flow4, ϕ, Pᵢ); # inverse with ptdf
@test flow3 ≈ flow4
SCOPF.calc_Pline!(flow5, θ, pf.X, pf.B, pf.DA, pf.θ, cont1, c1) # IMML theta
@test flow4 ≈ flow5 

# CONTINGENCY WITH POWER INJECTION CHANGE AND WITHOUT ISLANDING
SCOPF.calculate_line_flows!(flow1, pf, cont1, c1, Pᵢ=(Pᵢ .+ ΔPc)) # IMML flow
SCOPF.calculate_line_flows!(flow2, θ, B, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont1, c1, pf.slack) # inverse with theta
@test flow1 ≈ flow2
SCOPF.calc_Pline!(flow3, θ, pf.X, pf.B, pf.DA, θ₂, cont1, c1) # IMML theta
@test flow2 ≈ flow3
SCOPF.calc_isf!(ϕ, X, pf.X, pf.B, pf.DA, cont1, c1); LinearAlgebra.mul!(flow4, ϕ, (Pᵢ .+ ΔPc)); # IMML ptdf
@test flow3 ≈ flow4
SCOPF.calc_isf!(ϕ, K, pf.DA, pf.B, cont1, c1, pf.slack); LinearAlgebra.mul!(flow5, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow4 ≈ flow5

# CONTINGENCY WITH POWER INJECTION CHANGE AND ISLANDING
cont2 = bx[c2]
islands, island, island_b = SCOPF.handle_islands(pf.B, pf.DA, cont2, c2, pf.slack)
SCOPF.calculate_line_flows!(flow1, pf, cont2, c2, Pᵢ=(Pᵢ .+ ΔPc), nodes=islands[island], branches=island_b) # IMML flow
SCOPF.calculate_line_flows!(flow2, θ, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont2, c2, pf.slack, islands[island], island_b) # inverse with theta
@test flow1 ≈ flow2
SCOPF.calc_isf!(ϕ, pf.DA, pf.B, cont2, c2, pf.slack, islands[island], island_b); LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow2 ≈ flow3 
SCOPF.calc_Pline!(flow4, θ, pf, cont2, c2, Pᵢ=(Pᵢ .+ ΔPc), nodes=islands[island], branches=island_b) # IMML theta
@test flow3 ≈ flow4
SCOPF.calc_isf!(ϕ, X, pf.X, pf.B, pf.DA, cont2, c2); LinearAlgebra.mul!(flow5, ϕ, (Pᵢ .+ ΔPc)); # IMML ptdf
@test flow4 ≈ flow5 

# CONTINGENCY WITH POWER INJECTION CHANGE, ISLANDING, AND DISTRIBUTED SLACK
SCOPF.set_dist_slack!(pf.ϕ, opf.mgx, oplim.pg_lim_max)
SCOPF.calculate_line_flows!(flow1, pf, cont1, c1, Pᵢ=(Pᵢ .+ ΔPc)) # IMML flow
SCOPF.calculate_line_flows!(flow2, θ, B, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont1, c1, pf.slack) # inverse with theta
@test flow1 ≈ flow2
SCOPF.calc_isf!(ϕ, X, pf.X, pf.B, pf.DA, cont1, c1); LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc)); # IMML ptdf
@test flow2 ≈ flow3
SCOPF.calc_isf!(ϕ, pf.K, pf.DA, pf.B, cont1, c1, pf.slack); LinearAlgebra.mul!(flow4, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow3 ≈ flow4
dist_slack = getproperty.(get_active_power_limits.(opf.ctrl_generation), [:max])
slack_array = dist_slack / sum(dist_slack)
pf_ds = SCOPF.DCPowerFlow(opf.nodes, opf.branches, opf.idx, Pᵢ, slack_array=slack_array, mgx=opf.mgx)
SCOPF.calculate_line_flows!(flow5, pf_ds, cont1, c1, Pᵢ=(Pᵢ .+ ΔPc)) # IMML flow
@test flow4 ≈ flow5

end