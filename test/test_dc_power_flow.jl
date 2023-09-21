using Test
using Printf
using PowerSystems
import JuMP
import Gurobi
import LinearAlgebra

LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
system = System("data\\ELK14\\A5.m"); c1 = 1; c2 = 5
# system = SCOPF.System("data\\matpower\\IEEE_RTS.m"); c1 = 1; c2 = 11
# system = SCOPF.System("data\\matpower\\RTS_GMLC.m"); c1 = 1; c2 = 11
# system = SCOPF.System("data\\matpower\\ACTIVSg500.m"); c1 = 2; c2 = 1
# system = SCOPF.System("data\\matpower\\ACTIVSg2000.m"); c1 = 1; c2 = 9
SCOPF.fix_generation_cost!(system);
voll = SCOPF.make_voll(system)
model, opf, pf, oplim, _, _, _ = SCOPF.opf(SCOPF.SC, system, Gurobi.Optimizer, voll=voll);
SCOPF.solve_model!(model)
idx = SCOPF.get_nodes_idx(opf.nodes)
bx = SCOPF.get_bus_idx.(opf.branches, [idx])
slack = SCOPF.find_slack(opf.nodes)[1]

A = SCOPF.calc_A(opf.branches, idx)
D = SCOPF.calc_D(opf.branches)
DA = D*A
B = SCOPF.calc_B(A, DA)
X = SCOPF.calc_X(B, slack)
ϕ = SCOPF.calc_isf(D, A, X)
Pᵢ = SCOPF.get_value(model, :p0)
θ = SCOPF.run_pf(B, Pᵢ)
θ .-= θ[slack]

pf = SCOPF.DCPowerFlow(opf.nodes, opf.branches, Pᵢ, idx)
@test θ ≈ pf.θ
@test D*A*θ ≈ SCOPF.calculate_line_flows!(pf, Pᵢ)
@test D*A*θ ≈ SCOPF.calc_Pline!(pf)

flow1 = copy(pf.F)
flow2 = copy(pf.F)
flow3 = copy(pf.F)
flow4 = copy(pf.F)

cont = bx[c1]
SCOPF.calculate_line_flows!(flow1, pf, cont, c1) # IMML flow
SCOPF.calculate_line_flows!(flow2, θ, B, pf.DA, pf.B, Pᵢ, cont, c1, pf.slack) # inverse with theta
SCOPF.get_isf!(ϕ, X, pf.X, pf.B, pf.DA, cont, c1); LinearAlgebra.mul!(flow3, ϕ, Pᵢ); # IMML ptdf
SCOPF.get_isf!(ϕ, pf.DA, pf.B, cont, c1, pf.slack); LinearAlgebra.mul!(flow4, ϕ, Pᵢ); # inverse with ptdf
@test flow1 ≈ flow2
@test flow2 ≈ flow3
@test flow3 ≈ flow4

ΔPc = rand(length(Pᵢ))
SCOPF.calculate_line_flows!(flow2, θ, B, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont, c1, pf.slack) # inverse with theta
SCOPF.get_isf!(ϕ, X, pf.X, pf.B, pf.DA, cont, c1); LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc)); # IMML ptdf
SCOPF.get_isf!(ϕ, pf.DA, pf.B, cont, c1, pf.slack); LinearAlgebra.mul!(flow4, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow2 ≈ flow3
@test flow3 ≈ flow4

cont = bx[c2]
islands, island, island_b = SCOPF.handle_islands(pf.B, pf.DA, cont, c2, pf.slack)
SCOPF.calculate_line_flows!(flow2, θ, pf.DA, pf.B, (Pᵢ .+ ΔPc), cont, c2, pf.slack, islands[island], island_b) # inverse with theta
SCOPF.get_isf!(ϕ, pf.DA, pf.B, cont, c2, pf.slack, islands[island], island_b); LinearAlgebra.mul!(flow3, ϕ, (Pᵢ .+ ΔPc)); # inverse with ptdf
@test flow2 ≈ flow3

println("Benchmarks")
cont = bx[c1]
immlF = @benchmark SCOPF.calculate_line_flows!($flow1, $pf, $cont, $c1) # IMML flow
inv_theta = @benchmark SCOPF.calculate_line_flows!($flow2, $θ, $B, $pf.DA, $pf.B, $Pᵢ, $cont, $c1, $pf.slack) # inverse with theta
imml_ptdf = @benchmark begin SCOPF.get_isf!($ϕ, $X, $pf.X, $pf.B, $pf.DA, $cont, $c1); LinearAlgebra.mul!($flow3, $ϕ, $Pᵢ); end # IMML ptdf
inv_ptdf = @benchmark begin SCOPF.get_isf!($ϕ, $pf.DA, $pf.B, $cont, $c1, $pf.slack); LinearAlgebra.mul!($flow4, $ϕ, $Pᵢ); end # inverse with ptdf_pc
println("        IMML flow; inv theta; imml ptdf;  inv ptdf")
@printf("Min:    %9.0f; %9.0f; %9.0f; %9.0f\n", minimum(immlF).time, minimum(inv_theta).time, minimum(imml_ptdf).time, minimum(inv_ptdf).time)
@printf("Median: %9.0f; %9.0f; %9.0f; %9.0f\n", median(immlF).time, median(inv_theta).time, median(imml_ptdf).time, median(inv_ptdf).time)

ΔPc = rand(length(Pᵢ))
inv_theta = @benchmark SCOPF.calculate_line_flows!($flow2, $θ, $B, $pf.DA, $pf.B, ($Pᵢ .+ $ΔPc), $cont, $c1, $pf.slack) # inverse with theta
imml_ptdf = @benchmark begin SCOPF.get_isf!($ϕ, $X, $pf.X, $pf.B, $pf.DA, $cont, $c1); LinearAlgebra.mul!($flow3, $ϕ, ($Pᵢ .+ $ΔPc)); end # IMML ptdf
inv_ptdf = @benchmark begin SCOPF.get_isf!($ϕ, $pf.DA, $pf.B, $cont, $c1, $pf.slack); LinearAlgebra.mul!($flow4, $ϕ, ($Pᵢ .+ $ΔPc)); end # inverse with ptdf
println("        inv theta; imml ptdf;  inv ptdf")
@printf("Min:    %9.0f; %9.0f; %9.0f\n", minimum(inv_theta).time, minimum(imml_ptdf).time, minimum(inv_ptdf).time)
@printf("Median: %9.0f; %9.0f; %9.0f\n", median(inv_theta).time, median(imml_ptdf).time, median(inv_ptdf).time)

cont = bx[c2]
islands, island, island_b = SCOPF.handle_islands(pf.B, pf.DA, cont, c2, pf.slack)
inv_theta = @benchmark SCOPF.calculate_line_flows!($flow2, $θ, $pf.DA, $pf.B, ($Pᵢ .+ $ΔPc), $cont, $c2, $pf.slack, $islands[island], $island_b) # inverse with theta
inv_ptdf = @benchmark begin SCOPF.get_isf!($ϕ, $pf.DA, $pf.B, $cont, $c2, $pf.slack, $islands[island], $island_b); LinearAlgebra.mul!($flow3, ϕ, ($Pᵢ .+ $ΔPc)); end # inverse with ptdf
println("        inv theta;  inv ptdf")
@printf("Min:    %9.0f; %9.0f\n", minimum(inv_theta).time, minimum(inv_ptdf).time)
@printf("Median: %9.0f; %9.0f\n", median(inv_theta).time, median(inv_ptdf).time)
