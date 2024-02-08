import RaPidSCOPF as SCOPF
using Test
using PowerSystems
import JuMP
import LinearAlgebra
import Logging
import Random

# import Ipopt # LP, SOCP, NLP
# import Gurobi # LP, SOCP, NLP, MILP, MINLP
# import GLPK # LP, MILP
import HiGHS # LP, MILP
# import Tulip

Logging.disable_logging(Logging.Info)
# debug_logger = Logging.ConsoleLogger(stderr, Logging.Debug)
# Logging.global_logger(debug_logger); 
# Logging.disable_logging(Logging.Debug)
# Logging.disable_logging(Logging.Info)

## For timing functions (with allocations)
# using TimerOutputs
# const tmr = TimerOutput();
## SCOPF.tmr
## SCOPF.reset_timer!(SCOPF.tmr)

Random.seed!(42)
# const GUROBI_ENV = Gurobi.Env()
# optimizer = Gurobi.Optimizer(GUROBI_ENV)
# JuMP.set_optimizer_attribute(optimizer, "Threads", Threads.nthreads())
# optimizer = JuMP.optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.GLP_MSG_ON)


@testset "RaPidSCOPF" begin
    include("test_dc_opf.jl")

    include("test_dc_power_flow.jl")

    include("test_imml.jl")

    include("test_ieee_rts_benders.jl")
end