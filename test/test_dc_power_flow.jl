using Test
using Printf
using PowerSystems
import JuMP
import HiGHS
import LinearAlgebra
import Random
Random.seed!(42)

# SETUP
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
# system = System("data\\ELK14\\A5.m"); c1 = 1; c2 = 5
system = SCOPF.System("../cases/IEEE_RTS.m"); c1 = 1; c2 = 11
# system = SCOPF.System("data\\matpower\\RTS_GMLC.m"); c1 = 1; c2 = 52
# system = SCOPF.System("data\\matpower\\ACTIVSg500.m"); c1 = 2; c2 = 1
# system = SCOPF.System("data\\matpower\\ACTIVSg2000.m"); c1 = 1; c2 = 9
SCOPF.fix_generation_cost!(system);
voll = SCOPF.make_voll(system)
model, opf, pf, oplim, _, _, _ = SCOPF.opf(SCOPF.SC, system, HiGHS.Optimizer, voll=voll);
SCOPF.solve_model!(model)
idx = SCOPF.get_nodes_idx(opf.nodes)
bx = SCOPF.get_bus_idx.(opf.branches, [idx])
slack = SCOPF.find_slack(opf.nodes)[1]

# CALC MATRICES DIRECTLY
A = SCOPF.calc_A(opf.branches, idx)
D = SCOPF.calc_D(opf.branches)
DA = D*A
B = SCOPF.calc_B(A, DA)
X = SCOPF.calc_X(B, slack)
ϕ = SCOPF.calc_isf(DA, X)
Pᵢ = SCOPF.get_value(model, :p0)
θ = B \ Pᵢ
θ .-= θ[slack]
F = DA*θ
@test θ ≈ SCOPF.run_pf(pf.K, Pᵢ, pf.slack)

# CHECK DCPF-STRUCT CORRECTNESS
pf = SCOPF.DCPowerFlow(opf.nodes, opf.branches, idx, Pᵢ)
@test B ≈ pf.B
@test DA ≈ pf.DA
@test ϕ ≈ pf.ϕ
@test X ≈ pf.X
@test θ ≈ pf.θ
@test F ≈ SCOPF.calculate_line_flows!(pf, Pᵢ)
@test F ≈ SCOPF.calc_Pline!(pf)

K = SCOPF.get_klu(B, pf.slack)   
# CONTAINERS FOR POWER FLOW
flow1 = copy(pf.F)
flow2 = copy(pf.F)
flow3 = copy(pf.F)
flow4 = copy(pf.F)
flow5 = copy(pf.F)
flow6 = copy(pf.F)
flow7 = copy(pf.F)

# CONTINGENCY WITHOUT POWER INJECTION CHANGE OR ISLANDING
cont1 = bx[c1]
SCOPF.calculate_line_flows!(flow1, pf, cont1, c1) # IMML flow
SCOPF.calculate_line_flows!(flow2, θ, B, pf.DA, pf.B, Pᵢ, cont1, c1, pf.slack) # inverse with theta
@test flow1 ≈ flow2
SCOPF.get_isf!(ϕ, X, pf.X, pf.B, pf.DA, cont1, c1); LinearAlgebra.mul!(flow3, ϕ, Pᵢ); # IMML ptdf
@test flow2 ≈ flow3
SCOPF.get_isf!(ϕ, K, pf.DA, pf.B, cont1, c1, pf.slack); LinearAlgebra.mul!(flow4, ϕ, Pᵢ); # inverse with ptdf
@test flow3 ≈ flow4
SCOPF.calc_Pline!(flow5, θ, pf.X, pf.B, pf.DA, pf.θ, cont1, c1) # IMML theta
@test flow4 ≈ flow5 

# CONTINGENCY WITH POWER INJECTION CHANGE AND WITHOUT ISLANDING
ΔPc = rand(length(Pᵢ))
θ₂ = pf.X * (Pᵢ .+ ΔPc)
F₂ = pf.DA * θ₂
SCOPF.calculate_line_flows!(flow1, pf, cont1, c1, Pᵢ=(Pᵢ .+ ΔPc)) # IMML flow
SCOPF.calculate_line_flows!(flow2, θ, B, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont1, c1, pf.slack) # inverse with theta
@test flow1 ≈ flow2
SCOPF.get_isf!(ϕ, X, pf.X, pf.B, pf.DA, cont1, c1); LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc)); # IMML ptdf
@test flow2 ≈ flow3
SCOPF.get_isf!(ϕ, K, pf.DA, pf.B, cont1, c1, pf.slack); LinearAlgebra.mul!(flow4, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow3 ≈ flow4

# CONTINGENCY WITH POWER INJECTION CHANGE AND ISLANDING
cont2 = bx[c2]
islands, island, island_b = SCOPF.handle_islands(pf.B, pf.DA, cont2, c2, pf.slack)
SCOPF.calculate_line_flows!(flow1, pf, cont2, c2, Pᵢ=(Pᵢ .+ ΔPc), nodes=islands[island], branches=island_b) # IMML flow
SCOPF.calculate_line_flows!(flow2, θ, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont2, c2, pf.slack, islands[island], island_b) # inverse with theta
@test flow1 ≈ flow2
SCOPF.get_isf!(ϕ, pf.DA, pf.B, cont2, c2, pf.slack, islands[island], island_b); LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow2 ≈ flow3 
SCOPF.calc_Pline!(flow4, θ, pf.X, pf.B, pf.DA, θ₂, cont2, c2) # IMML theta
@test flow3 ≈ flow4
SCOPF.get_isf!(ϕ, X, pf.X, pf.B, pf.DA, cont2, c2); LinearAlgebra.mul!(flow5, ϕ, (Pᵢ .+ ΔPc)); # IMML ptdf
@test flow4 ≈ flow5 
SCOPF.get_isf!(ϕ, pf.ϕ, islands[island], island_b); LinearAlgebra.mul!(flow6, ϕ, (Pᵢ .+ ΔPc)); # hack
@test flow5 ≈ flow6 
SCOPF.calculate_line_flows!(flow7, ϕ, pf.ϕ, (Pᵢ .+ ΔPc), islands[island], island_b); # hack2
@test flow6 ≈ flow7 


println("Benchmarks")
immlF = @benchmark SCOPF.calculate_line_flows!($flow1, $pf, $cont1, $c1) # IMML flow
inv_theta = @benchmark SCOPF.calculate_line_flows!($flow2, $θ, $B, $pf.DA, $pf.B, $Pᵢ, $cont1, $c1, $pf.slack) # inverse with theta
imml_ptdf = @benchmark begin SCOPF.get_isf!($ϕ, $X, $pf.X, $pf.B, $pf.DA, $cont1, $c1); LinearAlgebra.mul!($flow3, $ϕ, $Pᵢ); end # IMML ptdf
inv_ptdf = @benchmark begin SCOPF.get_isf!($ϕ, $K, $pf.DA, $pf.B, $cont1, $c1, $pf.slack); LinearAlgebra.mul!($flow4, $ϕ, $Pᵢ); end # inverse with ptdf_pc
imml_theta = @benchmark SCOPF.calc_Pline!($flow5, $θ, $pf.X, $pf.B, $pf.DA, $pf.θ, $cont1, $c1)
println("        IMML flow; IMML thet; inv theta; imml ptdf;  inv ptdf")
@printf("Min:    %9.0f; %9.0f; %9.0f; %9.0f; %9.0f\n", minimum(immlF).time, minimum(imml_theta).time, minimum(inv_theta).time, minimum(imml_ptdf).time, minimum(inv_ptdf).time)
@printf("Median: %9.0f; %9.0f; %9.0f; %9.0f; %9.0f\n", median(immlF).time, median(imml_theta).time, median(inv_theta).time, median(imml_ptdf).time, median(inv_ptdf).time)

immlF = @benchmark SCOPF.calculate_line_flows!($flow1, $pf, $cont1, $c1, Pᵢ=($Pᵢ .+ $ΔPc)) # IMML flow
inv_theta = @benchmark SCOPF.calculate_line_flows!($flow2, $θ, $B, $pf.DA, $pf.B, ($Pᵢ .+ $ΔPc), $cont1, $c1, $pf.slack) # inverse with theta
imml_ptdf = @benchmark begin SCOPF.get_isf!($ϕ, $X, $pf.X, $pf.B, $pf.DA, $cont1, $c1); LinearAlgebra.mul!($flow3, $ϕ, ($Pᵢ .+ $ΔPc)); end # IMML ptdf
inv_ptdf = @benchmark begin SCOPF.get_isf!($ϕ, $K, $pf.DA, $pf.B, $cont1, $c1, $pf.slack); LinearAlgebra.mul!($flow4, $ϕ, ($Pᵢ .+ $ΔPc)); end # inverse with ptdf
println("         IMML flow; inv theta; imml ptdf;  inv ptdf")
@printf("Min:    %9.0f; %9.0f; %9.0f; %9.0f\n", minimum(immlF).time, minimum(inv_theta).time, minimum(imml_ptdf).time, minimum(inv_ptdf).time)
@printf("Median: %9.0f; %9.0f; %9.0f; %9.0f\n", median(immlF).time, median(inv_theta).time, median(imml_ptdf).time, median(inv_ptdf).time)

immlF = @benchmark SCOPF.calculate_line_flows!($flow1, $pf, $cont2, $c2, Pᵢ=($Pᵢ .+ $ΔPc), nodes=$islands[island], branches=$island_b) # IMML flow
inv_theta = @benchmark SCOPF.calculate_line_flows!($flow2, $θ, $pf.DA, $pf.B, ($Pᵢ .+ $ΔPc), $cont2, $c2, $pf.slack, $islands[island], $island_b) # inverse with theta
inv_ptdf = @benchmark begin SCOPF.get_isf!($ϕ, $pf.DA, $pf.B, $cont2, $c2, $pf.slack, $islands[island], $island_b); LinearAlgebra.mul!($flow3, ϕ, ($Pᵢ .+ $ΔPc)); end # inverse with ptdf
imml_theta = @benchmark begin θ₂ = SCOPF.run_pf($pf.K, ($Pᵢ .+ $ΔPc), $pf.slack); SCOPF.calc_Pline!($flow4, $θ, $pf.X, $pf.B, $pf.DA, θ₂, $cont2, $c2); end # IMML theta
imml_ptdf = @benchmark begin SCOPF.get_isf!($ϕ, $X, $pf.X, $pf.B, $pf.DA, $cont2, $c2); LinearAlgebra.mul!($flow5, ϕ, ($Pᵢ .+ $ΔPc)); end # IMML ptdf
hack = @benchmark begin SCOPF.get_isf!($ϕ, $pf.ϕ, $islands[island], $island_b); LinearAlgebra.mul!($flow6, ϕ, ($Pᵢ .+ $ΔPc)); end
hack2 = @benchmark SCOPF.calculate_line_flows!($flow7, $ϕ, $pf.ϕ, ($Pᵢ .+ $ΔPc), $islands[island], $island_b)
println("        IMML flow; IMML thet; imml ptdf;      hack;     hack2; inv theta;  inv ptdf")
@printf("Min:    %9.0f; %9.0f; %9.0f; %9.0f; %9.0f; %9.0f; %9.0f\n", minimum(immlF).time, minimum(imml_theta).time, minimum(imml_ptdf).time, minimum(hack).time, minimum(hack2).time, minimum(inv_theta).time, minimum(inv_ptdf).time)
@printf("Median: %9.0f; %9.0f; %9.0f; %9.0f; %9.0f; %9.0f; %9.0f\n", median(immlF).time, median(imml_theta).time, median(imml_ptdf).time, median(hack).time, median(hack2).time, median(inv_theta).time, median(inv_ptdf).time)
