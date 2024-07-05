# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2023

using BenchmarkTools
using Printf
using PowerSystems
import JuMP
import HiGHS
import LinearAlgebra
import Random
Random.seed!(42)

# SETUP
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
# system = System(joinpath("data","ELK14","A5.m")); c1 = 1; c2 = 5
# system = SCOPF.System(joinpath("data","matpower","IEEE_RTS.m")); c1 = 1; c2 = 11
system = SCOPF.System(joinpath("data","matpower","RTS_GMLC.m")); c1 = 1; c2 = 52
# system = SCOPF.System(joinpath("data","matpower","ACTIVSg500.m")); c1 = 2; c2 = 1
# system = SCOPF.System(joinpath("data","matpower","ACTIVSg2000.m")); c1 = 1; c2 = 9
# system = SCOPF.System(joinpath("data","matpower","case_ACTIVSg10k.m")); c1 = 1; c2 = 4
SCOPF.fix_generation_cost!(system);
voll = SCOPF.make_voll(system)
model, opf, pf, oplim, _, _, _ = SCOPF.opf_base(SCOPF.Base_SCOPF, system, HiGHS.Optimizer(), voll=voll);
SCOPF.solve_model!(model)
bx = SCOPF.get_bus_idx.(opf.branches, [opf.idx])
slack = SCOPF.find_slack(opf.nodes)[1]

Pᵢ = SCOPF.get_value(model, :p0)
pf = SCOPF.DCPowerFlow(opf.nodes, opf.branches, opf.idx, Pᵢ)
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
ΔPc = (rand(length(Pᵢ)) .- 0.5) ./ 20
cont2 = bx[c2]
islands, island, island_b = SCOPF.handle_islands(pf.B, pf.DA, cont2, c2, pf.slack)

println("Benchmarks")
# CONTINGENCY WITHOUT POWER INJECTION CHANGE OR ISLANDING
immlF = @benchmark SCOPF.calculate_line_flows!($flow1, $pf, $cont1, $c1) # IMML flow
inv_theta = @benchmark SCOPF.calculate_line_flows!($flow2, $θ, $B, $pf.DA, $pf.B, $Pᵢ, $cont1, $c1, $pf.slack) # inverse with theta
imml_ptdf = @benchmark begin SCOPF.calc_isf!($ϕ, $X, $pf.X, $pf.B, $pf.DA, $cont1, $c1); LinearAlgebra.mul!($flow3, $ϕ, $Pᵢ); end # IMML ptdf
inv_ptdf = @benchmark begin SCOPF.calc_isf!($ϕ, $K, $pf.DA, $pf.B, $cont1, $c1, $pf.slack); LinearAlgebra.mul!($flow4, $ϕ, $Pᵢ); end # inverse with ptdf_pc
imml_theta = @benchmark SCOPF.calc_Pline!($flow5, $θ, $pf.X, $pf.B, $pf.DA, $pf.θ, $cont1, $c1)
println("        IMML flow; IMML thet; inv theta; imml ptdf;  inv ptdf")
@printf("Min:    %9.0f; %9.0f; %9.0f; %9.0f; %9.0f\n", minimum(immlF).time, minimum(imml_theta).time, minimum(inv_theta).time, minimum(imml_ptdf).time, minimum(inv_ptdf).time)
@printf("Median: %9.0f; %9.0f; %9.0f; %9.0f; %9.0f\n", median(immlF).time, median(imml_theta).time, median(inv_theta).time, median(imml_ptdf).time, median(inv_ptdf).time)

# CONTINGENCY WITH POWER INJECTION CHANGE AND WITHOUT ISLANDING
immlF = @benchmark SCOPF.calculate_line_flows!($flow1, $pf, $cont1, $c1, Pᵢ=($Pᵢ .+ $ΔPc)) # IMML flow
inv_theta = @benchmark SCOPF.calculate_line_flows!($flow2, $θ, $B, $pf.DA, $pf.B, ($Pᵢ .+ $ΔPc), $cont1, $c1, $pf.slack) # inverse with theta
imml_theta = @benchmark SCOPF.calc_Pline!($flow4, $θ, $pf, $cont1, $c1, Pᵢ=($Pᵢ .+ $ΔPc)) # IMML theta
imml_ptdf = @benchmark begin SCOPF.calc_isf!($ϕ, $X, $pf.X, $pf.B, $pf.DA, $cont1, $c1); LinearAlgebra.mul!($flow3, $ϕ, ($Pᵢ .+ $ΔPc)); end # IMML ptdf
inv_ptdf = @benchmark begin SCOPF.calc_isf!($ϕ, $K, $pf.DA, $pf.B, $cont1, $c1, $pf.slack); LinearAlgebra.mul!($flow4, $ϕ, ($Pᵢ .+ $ΔPc)); end # inverse with ptdf
println("         IMML flow; IMML thet; inv theta; imml ptdf;  inv ptdf")
@printf("Min:    %9.0f; %9.0f; %9.0f; %9.0f; %9.0f\n", minimum(immlF).time, minimum(imml_theta).time, minimum(inv_theta).time, minimum(imml_ptdf).time, minimum(inv_ptdf).time)
@printf("Median: %9.0f; %9.0f; %9.0f; %9.0f; %9.0f\n", median(immlF).time, median(imml_theta).time, median(inv_theta).time, median(imml_ptdf).time, median(inv_ptdf).time)

# CONTINGENCY WITH POWER INJECTION CHANGE AND ISLANDING
immlF = @benchmark SCOPF.calculate_line_flows!($flow1, $pf, $cont2, $c2, Pᵢ=($Pᵢ .+ $ΔPc), nodes=$islands[island], branches=$island_b) # IMML flow
inv_theta = @benchmark SCOPF.calculate_line_flows!($flow2, $θ, $pf.DA, $pf.B, ($Pᵢ .+ $ΔPc), $cont2, $c2, $pf.slack, $islands[island], $island_b) # inverse with theta
inv_ptdf = @benchmark begin SCOPF.calc_isf!($ϕ, $pf.DA, $pf.B, $cont2, $c2, $pf.slack, $islands[island], $island_b); LinearAlgebra.mul!($flow3, ϕ, ($Pᵢ .+ $ΔPc)); end # inverse with ptdf
imml_theta = @benchmark SCOPF.calc_Pline!($flow4, $θ, $pf, $cont2, $c2, Pᵢ=($Pᵢ .+ $ΔPc), nodes=$islands[island], branches=$island_b) # IMML theta
imml_ptdf = @benchmark begin SCOPF.calc_isf!($ϕ, $X, $pf.X, $pf.B, $pf.DA, $cont2, $c2); LinearAlgebra.mul!($flow5, ϕ, ($Pᵢ .+ $ΔPc)); end # IMML ptdf
# hack = @benchmark begin SCOPF.calc_isf!($ϕ, $pf.ϕ, $islands[island], $island_b); LinearAlgebra.mul!($flow6, ϕ, ($Pᵢ .+ $ΔPc)); end
# hack2 = @benchmark SCOPF.calculate_line_flows!($flow7, $ϕ, $pf.ϕ, ($Pᵢ .+ $ΔPc), $islands[island], $island_b)
println("        IMML flow; IMML thet; imml ptdf;      hack;     hack2; inv theta;  inv ptdf")
@printf("Min:    %9.0f; %9.0f; %9.0f; %9.0f; %9.0f; %9.0f; %9.0f\n", minimum(immlF).time, minimum(imml_theta).time, minimum(imml_ptdf).time, minimum(inv_theta).time, minimum(inv_ptdf).time)
@printf("Median: %9.0f; %9.0f; %9.0f; %9.0f; %9.0f; %9.0f; %9.0f\n", median(immlF).time, median(imml_theta).time, median(imml_ptdf).time, median(inv_theta).time, median(inv_ptdf).time)

# CONTINGENCY WITH ISLANDING, ONLY PTDF-VECTOR
imml_ptdfvec = @benchmark SCOPF.calc_ptdf_vec!($pf.vn_tmp, $pf, $c1, $cont1[1], $cont1[2], 1, $bx[1][1], $bx[1][2]) # IMML ptdf-vec
inv_ptdfvec = @benchmark SCOPF.calc_isf_vec!($pf.vn_tmp, $K, $pf.DA, $pf.B, $cont1, $c1, $pf.slack, 1) # inverse with ptdf-vec
imml_ptdf = @benchmark SCOPF.calc_isf!($ϕ, $X, $pf.X, $pf.B, $pf.DA, $cont1, $c1) # IMML ptdf
inv_ptdf = @benchmark SCOPF.calc_isf!($ϕ, $K, $pf.DA, $pf.B, $cont1, $c1, $pf.slack) # inverse with ptdf_pc

inv_ptdfvec_i = @benchmark SCOPF.calc_isf_vec!($pf.vn_tmp, $pf.ϕ, $1, $islands[island], $island_b) # inverse with ptdf-vec
inv_ptdf_i = @benchmark SCOPF.calc_isf!($ϕ, $pf.DA, $pf.B, $cont2, $c2, $pf.slack, $islands[island], $island_b) # inverse with ptdf
imml_ptdf_i = @benchmark SCOPF.calc_isf!($ϕ, $X, $pf.X, $pf.B, $pf.DA, $cont2, $c2) # IMML ptdf

println("        imml vec; inv vec; imml; inv;   inv island vec; inv island; imml island")
@printf("Min:    %9.0f; %9.0f; %9.0f; %9.0f; %9.0f; %9.0f; %9.0f\n", minimum(imml_ptdfvec).time, minimum(inv_ptdfvec).time, minimum(imml_ptdf).time, minimum(inv_ptdf).time, minimum(inv_ptdfvec_i).time, minimum(inv_ptdf_i).time, minimum(imml_ptdf_i).time)
@printf("Median: %9.0f; %9.0f; %9.0f; %9.0f; %9.0f; %9.0f; %9.0f\n", median(imml_ptdfvec).time, median(inv_ptdfvec).time, median(imml_ptdf).time, median(inv_ptdf).time, median(inv_ptdfvec_i).time, median(inv_ptdf_i).time, median(imml_ptdf_i).time)
