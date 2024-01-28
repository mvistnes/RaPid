using Test
using PowerSystems
import JuMP
# import Ipopt # LP, SOCP, NLP
import Gurobi # LP, SOCP, NLP, MILP, MINLP
# import GLPK # LP, MILP
import HiGHS # LP, MILP
# import Tulip
import Random
import Logging

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
const GUROBI_ENV = Gurobi.Env()
# optimizer = Gurobi.Optimizer(GUROBI_ENV)
# JuMP.set_optimizer_attribute(optimizer, "Threads", Threads.nthreads())
# optimizer = JuMP.optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.GLP_MSG_ON)

# open("output.txt", "w") do out
#     redirect_stdout(out) do

system = System(joinpath("data","matpower","IEEE_RTS.m"))
# SCOPF.fix_generation_cost!(system);
# nodes = SCOPF.sort_components!(SCOPF.get_nodes(system))
# idx = SCOPF.get_nodes_idx(nodes)
# SCOPF.set_ramp_limits!(system, 0.01)
SCOPF.set_rate!.(SCOPF.get_branches(system), SCOPF.get_rate.(SCOPF.get_branches(system)) * 0.8);
SCOPF.set_operation_cost!.(SCOPF.get_gens_h(system), [15.0, 16.0, 17.0, 18.0, 19.0, 20.0])

# voll, prob, contingencies = SCOPF.setup(system, 10, 40);
voll = [4304., 5098., 5245., 5419., 4834., 5585., 5785., 5192., 4575., 5244., 4478., 5698., 4465., 4859., 5032., 5256., 4598.]
branches = SCOPF.sort_components!(SCOPF.get_branches(system));
ctrl_generation = SCOPF.sort_components!(SCOPF.get_ctrl_generation(system));
# c = [7,12,13,21,22,23,27]
contingencies = branches
# contingencies = vcat(branches, ctrl_generation)
prob =
    [ # spesified for the RTS-96
        0.24, 0.51, 0.33, 0.39, 0.48, 0.38, 0.02, 0.36, 0.34, 0.33, 0.30, 0.44, 0.44,
        0.02, 0.02, 0.02, 0.02, 0.40, 0.39, 0.40, 0.52, 0.49, 0.38, 0.33, 0.41, 0.41,
        0.41, 0.35, 0.34, 0.32, 0.54, 0.35, 0.35, 0.38, 0.38, 0.34, 0.34, 0.45, # branches
        0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 
        0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 
        0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 # generators
    ]
# prob /= 8760
prob /= 100
# prob = prob[c]
short = 1.2
long = 1.0
ramp_minutes = 10.
max_shed = 0.1
ramp_mult = 2.

function test_benders(system, optimizer, voll, contingencies, prob, max_shed,ramp_mult, ramp_minutes, short, long)
    for case in [SCOPF.Base_SCOPF, SCOPF.P_SCOPF, SCOPF.OPF(true, false, true, false, false), SCOPF.PC_SCOPF, SCOPF.PCF_SCOPF, SCOPF.PC2_SCOPF, SCOPF.PC2F_SCOPF]
        p_failure = ifelse(case.C2F, 0.10, 0.00)
        println(case)
        println("Start PTDF")
        mod_ptdf, opf_ptdf, pf_ptdf, oplim_ptdf, Pc_ptdf, Pcc_ptdf, Pccx_ptdf = SCOPF.opf_base(case, system, HiGHS.Optimizer(), voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed,
            ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure);
        SCOPF.add_branch_constraints!(mod_ptdf, pf_ptdf.ϕ, mod_ptdf[:p0], oplim_ptdf.branch_rating)
        SCOPF.solve_model!(mod_ptdf);

        println("Start Contingency select")
        mod_cont, opf_cont, pf_cont, oplim_cont, Pc_cont, Pcc_cont, Pccx_cont, tot_t = SCOPF.run_contingency_select(case, system, HiGHS.Optimizer(), voll, prob, contingencies, 
            max_shed=max_shed, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure);

        println("Start Benders")
        model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = SCOPF.run_benders(case, system, HiGHS.Optimizer(), voll, prob, contingencies, max_shed=max_shed,
            ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure);

        @test JuMP.objective_value(mod_ptdf) ≈ JuMP.objective_value(mod_cont)
        @test JuMP.objective_value(mod_ptdf) ≈ JuMP.objective_value(model)
    end
    return
end

test_benders(system, HiGHS.Optimizer(), voll, contingencies, prob, max_shed,ramp_mult, ramp_minutes, short, long)

function benders_pc_scopf()
    println("\nPreventive-Corrective SCOPF")
    println("Start PTDF")
    @time mod_ptdf, opf_ptdf, pf_ptdf, oplim_ptdf, Pc_ptdf, Pcc_ptdf, Pccx_ptdf = SCOPF.opf_base(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GUROBI_ENV), voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    SCOPF.add_branch_constraints!(mod_ptdf, pf_ptdf.ϕ, mod_ptdf[:p0], oplim_ptdf.branch_rating)
    @time SCOPF.solve_model!(mod_ptdf);
    println("Objective value PTDF: ", JuMP.objective_value(mod_ptdf))
    SCOPF.print_contingency_results(opf_ptdf, Pc_ptdf, Pcc_ptdf)

    println("Start Contingency select")
    @time mod_cont, opf_cont, pf_cont, oplim_cont, Pc_cont, Pcc_cont, Pccx_cont, tot_t = SCOPF.run_contingency_select(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GUROBI_ENV), voll, prob, contingencies, 
        max_shed=max_shed, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00);
    println("Solver time: ", tot_t, " Objective value Contingency select: ", JuMP.objective_value(mod_cont))
    SCOPF.print_contingency_results(opf_cont, Pc_cont, Pcc_cont, Pccx_cont)

    println("Start Benders")
    @time model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = SCOPF.run_benders(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GUROBI_ENV), voll, prob, contingencies, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00);
    println("Solver time: ", tot_t, " Objective value Benders: ", JuMP.objective_value(model))
    SCOPF.print_contingency_results(opf, Pc, Pcc, Pccx)

    @test JuMP.objective_value(mod_ptdf) ≈ JuMP.objective_value(mod_cont)
    @test JuMP.objective_value(mod_ptdf) ≈ JuMP.objective_value(model)

    SCOPF.print_power_flow(opf, model)
    bx = SCOPF.typesort_component.(SCOPF.get_interarea(opf.branches), [opf])
    SCOPF.print_contingency_power_flow(opf, model, pf, Pc, Pcc, Pccx, short, long; subset=first.(bx))
end