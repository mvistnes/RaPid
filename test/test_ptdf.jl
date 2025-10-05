@testset "Test PTDF" begin

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
ϕ1 = Matrix{eltype(pf.DA)}(undef, size(pf.DA))
ϕ2 = Matrix{eltype(pf.DA)}(undef, size(pf.DA))

cont1 = bx[c1]

# CONTINGENCY WITHOUT ISLANDING
SCOPF.calc_ptdf!(ϕ1, X, pf.X, pf.B, pf.DA, cont1, c1); # IMML
SCOPF.calc_ptdf!(ϕ2, K, pf.DA, pf.B, cont1, c1, pf.slack); # inverse
@test ϕ1 ≈ ϕ2

# CONTINGENCY WITHOUT ISLANDING, ONLY PTDF-VECTOR
SCOPF.calc_ptdf_vec!(pf.vn_tmp, pf, c1, cont1[1], cont1[2], 3, bx[3][1], bx[3][2]) # IMML vec
@test ϕ1[3,:] ≈ pf.vn_tmp

SCOPF.calc_ptdf_vec!(pf.vn_tmp, K, pf.DA, pf.B, cont1, c1, pf.slack, 3) # inverse vec
@test ϕ1[3,:] ≈ pf.vn_tmp

# CONTINGENCY WITHOUT ISLANDING AND WITH DISTRIBUTED SLACK
SCOPF.set_dist_slack!(pf.ϕ, opf.mgx, oplim.pg_lim_max)
SCOPF.calc_ptdf!(ϕ1, X, pf.X, pf.B, pf.DA, cont1, c1); # IMML
SCOPF.calc_ptdf!(ϕ2, pf.K, pf.DA, pf.B, cont1, c1, pf.slack); # inverse 
@test ϕ1 ≈ ϕ2

# CONTINGENCY WITH ISLANDING
cont2 = bx[c2]
islands, island, island_b = SCOPF.handle_islands(pf.B, pf.DA, cont2, c2, pf.slack)
SCOPF.calc_ptdf!(ϕ1, pf.DA, pf.B, cont2, c2, pf.slack, islands[island], island_b); # inverse 
SCOPF.calc_ptdf!(ϕ2, X, pf.X, pf.B, pf.DA, cont2, c2); # IMML
@test ϕ1 ≈ ϕ2

# CONTINGENCY WITH ISLANDING, ONLY PTDF-VECTOR
SCOPF.calc_ptdf_vec!(pf.vn_tmp, pf.DA, pf.B, cont2, c2, pf.slack, islands[island], island_b, 3) # inverse vec
@test ϕ1[3,:] ≈ pf.vn_tmp

end