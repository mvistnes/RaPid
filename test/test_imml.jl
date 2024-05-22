@testset "Test IMML" begin

# SETUP
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
# system = System(joinpath("data","ELK14","A5.m")); c1 = 1; c2 = 5
# system = SCOPF.System(joinpath("cases","IEEE_RTS.m")); c1 = [1]; c2 = [11]; c3 = [1,2]; c4 = [3,9]
system = SCOPF.System(joinpath("cases","RTS_GMLC.m")); c1 = [1]; c2 = [52]; c3 = [1,2]; c4 = [90,111]
# system = SCOPF.System(joinpath("cases","ACTIVSg500.m")); c1 = 2; c2 = 1
# system = SCOPF.System(joinpath("cases","ACTIVSg2000.m")); c1 = 1; c2 = 9
# system = SCOPF.System(joinpath("cases","case_ACTIVSg10k.m")); c1 = 1; c2 = 4
# SCOPF.fix_generation_cost!(system);
voll = SCOPF.make_voll(system)
model, opf, pf, oplim, brc_up, brc_down, _, _, _ = SCOPF.opf_base(SCOPF.Base_SCOPF, system, HiGHS.Optimizer(), voll=voll);
SCOPF.constrain_branches!(model, pf, oplim, brc_up, brc_down, 0.0)
bx = [Tuple(opf.mbx[i,:].nzind) for i in axes(opf.mbx,1)]
Pᵢ = SCOPF.get_value(model, :p0)

B = copy(pf.B)
X = copy(pf.X)
ϕ = copy(pf.ϕ)
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
θ₂ = pf.X * (Pᵢ .+ ΔPc)
F₂ = pf.DA * θ₂

# CONTINGENCY WITHOUT POWER INJECTION CHANGE OR ISLANDING
SCOPF.calculate_line_flows!(flow1, pf, cont1[1], c1[1]) # IMML flow
SCOPF.calculate_line_flows!(flow2, θ, B, pf.DA, pf.B, Pᵢ, cont1, c1, pf.slack) # inverse with theta
@test flow1 ≈ flow2
SCOPF.calc_isf!(ϕ, X, pf.X, pf.B, pf.DA, cont1[1], c1[1]); LinearAlgebra.mul!(flow3, ϕ, Pᵢ); # IMML ptdf
@test flow2 ≈ flow3
SCOPF.calc_isf!(ϕ, pf.DA, pf.B, cont1, c1, pf.slack); LinearAlgebra.mul!(flow4, ϕ, Pᵢ); # inverse with ptdf
@test flow3 ≈ flow4

# CONTINGENCY WITH POWER INJECTION CHANGE AND WITHOUT ISLANDING
SCOPF.calculate_line_flows!(flow1, pf, cont1[1], c1[1], (Pᵢ .+ ΔPc)) # IMML flow
SCOPF.calculate_line_flows!(flow2, θ, B, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont1, c1, pf.slack) # inverse with theta
@test flow1 ≈ flow2
SCOPF.calc_isf!(ϕ, X, pf.X, pf.B, pf.DA, cont1[1], c1[1]); LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc)); # IMML ptdf
@test flow2 ≈ flow3
SCOPF.calc_isf!(ϕ, pf.DA, pf.B, cont1, c1, pf.slack); LinearAlgebra.mul!(flow4, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow3 ≈ flow4

# CONTINGENCY WITH POWER INJECTION CHANGE AND ISLANDING
cont2 = bx[c2]
islands, island_b = SCOPF.handle_islands(pf.B, pf.DA, cont2[1], c2[1])
slack = SCOPF.find_slack(islands[1], pf.slack, oplim.pg_lim_max, opf.mgx)
SCOPF.calculate_line_flows!(flow1, θ, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont2, c2, slack, islands[1], island_b[1]) # inverse with theta
SCOPF.calc_isf!(ϕ, pf.DA, pf.B, cont2, c2, slack, islands[1], island_b[1]); LinearAlgebra.mul!(flow2, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow1 ≈ flow2
SCOPF.calc_isf!(ϕ, pf.ϕ, islands[1], island_b[1]); LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc)); # hack
@test flow2 ≈ flow3 
SCOPF.calculate_line_flows!(flow4, ϕ, pf.ϕ, (Pᵢ .+ ΔPc), islands[1], island_b[1]); # hack2
@test flow3 ≈ flow4

cont3 = bx[c3]

# MULTI-CONTINGENCY WITHOUT POWER INJECTION CHANGE OR ISLANDING
SCOPF.calculate_line_flows!(flow1, pf, cont3, c3, Pᵢ) # IMML flow
SCOPF.calculate_line_flows!(flow2, θ, B, pf.DA, pf.B, Pᵢ, cont3, c3, pf.slack) # inverse with theta
@test flow1 ≈ flow2
X = SCOPF.calc_X(pf.X, pf.A, pf.D, c3); LinearAlgebra.mul!(ϕ, pf.DA, X); LinearAlgebra.mul!(flow3, ϕ, Pᵢ); flow3[c3] .= 0.0; # IMML ptdf
@test flow2 ≈ flow3
SCOPF.calc_isf!(ϕ, pf.DA, pf.B, cont3, c3, pf.slack); LinearAlgebra.mul!(flow4, ϕ, Pᵢ); # inverse with ptdf
@test flow3 ≈ flow4

# MULTI-CONTINGENCY WITH POWER INJECTION CHANGE AND WITHOUT ISLANDING
SCOPF.calculate_line_flows!(flow1, pf, cont3, c3, (Pᵢ .+ ΔPc)) # IMML flow
SCOPF.calculate_line_flows!(flow2, θ, B, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont3, c3, pf.slack) # inverse with theta
@test flow1 ≈ flow2
ϕ .= SCOPF.calc_isf(pf, cont3, c3); LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc)); # IMML ptdf
@test flow2 ≈ flow3
SCOPF.calc_isf!(ϕ, pf.DA, pf.B, cont3, c3, pf.slack); LinearAlgebra.mul!(flow4, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow3 ≈ flow4

# MULTI-CONTINGENCY WITH POWER INJECTION CHANGE AND ISLANDING
cont4 = bx[c4]
islands, island_b = SCOPF.handle_islands(pf.B, pf.DA, cont4, c4)
slack = SCOPF.find_slack(islands[1], pf.slack, oplim.pg_lim_max, opf.mgx)
SCOPF.calculate_line_flows!(flow1, θ, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont4, c4, slack, islands[1], island_b[1]) # inverse with theta
SCOPF.calc_isf!(ϕ, pf.DA, pf.B, cont4, c4, slack, islands[1], island_b[1]); LinearAlgebra.mul!(flow2, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow1 ≈ flow2
SCOPF.calc_isf!(ϕ, pf.ϕ, islands[1], island_b[1]); LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc)); # hack
@test flow2 ≈ flow3 broken=true
SCOPF.calculate_line_flows!(flow4, ϕ, pf.ϕ, (Pᵢ .+ ΔPc), islands[1], island_b[1]); # hack2
@test flow3 ≈ flow4 broken=true

# CONTINGENCY WITH POWER INJECTION CHANGE, ISLANDING, AND DISTRIBUTED SLACK
SCOPF.set_dist_slack!(pf.ϕ, opf.mgx, oplim.pg_lim_max)
SCOPF.calculate_line_flows!(flow1, pf, cont1[1], c1[1], (Pᵢ .+ ΔPc)) # IMML flow
SCOPF.calculate_line_flows!(flow2, θ, B, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont1, c1, pf.slack) # inverse with theta
@test flow1 ≈ flow2 broken=true
SCOPF.calc_isf!(ϕ, X, pf.X, pf.B, pf.DA, cont1[1], c1[1]); LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc)); # IMML ptdf
@test flow2 ≈ flow3
SCOPF.calc_isf!(ϕ, pf.DA, pf.B, cont1, c1, pf.slack); LinearAlgebra.mul!(flow4, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow3 ≈ flow4

end