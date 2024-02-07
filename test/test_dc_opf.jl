using Test
using PowerSystems
import JuMP
import HiGHS
import LinearAlgebra
import Logging
import Random
Random.seed!(42)
Logging.disable_logging(Logging.Info)

# SETUP
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
# system = System(joinpath("data","ELK14","A5.m")); c1 = 1; c2 = 5
system = SCOPF.System(joinpath("cases","IEEE_RTS.m")); c1 = 1; c2 = 11
# system = SCOPF.System(joinpath("cases","RTS_GMLC.m")); c1 = 1; c2 = 52
# system = SCOPF.System(joinpath("cases","ACTIVSg500.m")); c1 = 2; c2 = 1
# system = SCOPF.System(joinpath("cases","ACTIVSg2000.m")); c1 = 1; c2 = 9
# system = SCOPF.System(joinpath("cases","case_ACTIVSg10k.m")); c1 = 1; c2 = 4
SCOPF.fix_generation_cost!(system);
voll = SCOPF.make_voll(system)
model, opf, pf, oplim, _, _, _ = SCOPF.opf_base(SCOPF.OPF(true, false, false, false, false), system, HiGHS.Optimizer(), voll=voll);
SCOPF.constrain_branches!(model, pf, oplim, 0.0)
model2, opf2, pf2, oplim2, _, _, _ = SCOPF.opf_base(SCOPF.OPF(true, false, false, false, false), system, HiGHS.Optimizer(), voll=voll);
SCOPF.add_branch_constraints!(model2, pf.ϕ, model2[:p0], oplim.branch_rating)
SCOPF.solve_model!(model2)

@test JuMP.objective_value(model) ≈ JuMP.objective_value(model2)
@test SCOPF.get_value(model, :p0) ≈ SCOPF.get_value(model2, :p0)
