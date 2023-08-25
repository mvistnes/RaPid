# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

using PowerSystems
import JuMP
import Gurobi
using BenchmarkTools

# open("output.txt", "w") do out
#     redirect_stdout(out) do

function test_ieee_rts()
    system = SCOPF.System("data\\matpower\\IEEE_RTS.m")
    SCOPF.set_rate!.(SCOPF.get_branches(system), SCOPF.get_rate.(SCOPF.get_branches(system)) * 0.8);
    benchmark_sys(system, 2, 11)

    SCOPF.set_operation_cost!.(SCOPF.get_gens_h(system), [15.0, 16.0, 17.0, 18.0, 19.0, 20.0])
    voll = [4304, 5098, 5245, 5419, 4834, 5585, 5785, 5192, 4575, 5244, 4478, 5698, 4465, 4859, 5032, 5256, 4598]
    branches = SCOPF.sort_components!(SCOPF.get_branches(system));
    contingencies = branches
    prob =
        [ # spesified for the RTS-96
            # 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
            # 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, # generators
            0.24, 0.51, 0.33, 0.39, 0.48, 0.38, 0.02, 0.36, 0.34, 0.33, 0.30, 0.44, 0.44,
            0.02, 0.02, 0.02, 0.02, 0.40, 0.39, 0.40, 0.52, 0.49, 0.38, 0.33, 0.41, 0.41,
            0.41, 0.35, 0.34, 0.32, 0.54, 0.35, 0.35, 0.38, 0.38, 0.34, 0.34, 0.45 # branches
        ]
    prob /= 100
    short = 1.2
    long = 1.0
    ramp_minutes = 10
    max_shed = 0.1
    ramp_mult = 10

    println("Start PTDF")
    @time opfm_ptdf, pf_ptdf, Pc_ptdf, Pcc_ptdf = SCOPF.opf(SCOPF.PCSC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_limit_multi=short, long_term_limit_multi=long);
    @time SCOPF.solve_model!(opfm_ptdf.mod);
    println("Objective value PTDF: ", JuMP.objective_value(opfm_ptdf.mod))
    SCOPF.print_contingency_P(opfm_ptdf, Pc_ptdf, Pcc_ptdf)
    # @time opfm_norm = SCOPF.scopf(SCOPF.PCSC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed, 
    #     ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_limit_multi=short);
    # @time SCOPF.solve_model!(opfm_norm.mod);

    println("Start Benders")
    @time opfm, pf, Pc, Pcc, Pccx = SCOPF.run_benders(SCOPF.PCSC, system, Gurobi.Optimizer, voll, prob, contingencies, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, branch_short_term_limit_multi=short, branch_long_term_limit_multi=long, p_failure=0.00);
    println("Objective value Benders: ", JuMP.objective_value(opfm.mod))
    SCOPF.print_contingency_P(opfm, Pc, Pcc, Pccx)

    println("Start Contingency select")
    @time opfm_cont, pf_cont, Pc_cont, Pcc_cont, Pccx_cont = SCOPF.run_contingency_select(SCOPF.PCSC, system, SCOPF.Gurobi.Optimizer, voll, prob, contingencies, 
        max_shed=max_shed, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, branch_short_term_limit_multi=short, branch_long_term_limit_multi=long, p_failure=0.00)
    println("Objective value Contingency select: ", JuMP.objective_value(opfm_cont.mod))
    SCOPF.print_contingency_P(opfm_cont, Pc_cont, Pcc_cont, Pccx_cont)
end

function test_big_sys()
    short = 1.5
    long = 1.25
    ramp_minutes = 10
    max_shed = 0.5
    ramp_mult = 10
    for (sysname, c1, c2) in zip(["data\\matpower\\ACTIVSg500.m", "data\\matpower\\ACTIVSg2000.m"], [2, 2], [1, 9])
        system = SCOPF.System(sysname)
        benchmark_sys(system, c1, c2)

        voll, prob, contingencies = SCOPF.setup(system, 100, 400);
        SCOPF.fix_generation_cost!(system);
        SCOPF.set_ramp_limits!(system, 0.01);
        prob = fill(0.01, length(contingencies))

        println("Start SC normal")
        @time opfm_norm = SCOPF.scopf(SCOPF.SC, system, SCOPF.Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=1.0,
            ramp_minutes=ramp_minutes, short_term_limit_multi=short, debug=true);
        @time SCOPF.solve_model!(opfm_norm.mod);
        println("Start SC PTDF")
        @time opfm_ptdf, pf_ptdf, Pc_ptdf, Pcc_ptdf = SCOPF.opf(SCOPF.SC, system, SCOPF.Gurobi.Optimizer, voll=voll, max_shed=max_shed);
        @time SCOPF.solve_model!(opfm_ptdf.mod);

        println("Start Benders")
        @time opfm, pf, Pc, Pcc, Pccx = SCOPF.run_benders(SCOPF.PCSC, system, Gurobi.Optimizer, voll, prob, contingencies, max_shed=max_shed,
            ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, branch_short_term_limit_multi=short, branch_long_term_limit_multi=long, p_failure=0.00);
            
        # println("Start Contingency select")
        # @time opfm_cont, pf_cont, Pc_cont, Pcc_cont, Pccx_cont = SCOPF.run_contingency_select(SCOPF.PCSC, system, SCOPF.Gurobi.Optimizer, voll, prob, contingencies, 
        #     max_shed=max_shed, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, branch_short_term_limit_multi=short, branch_long_term_limit_multi=long, p_failure=0.00)
        # println("Objective value Contingency select: ", JuMP.objective_value(opfm_cont.mod))
        # SCOPF.print_contingency_P(opfm_cont, Pc_cont, Pcc_cont, Pccx_cont)
    
        # println("Start PCSC PTDF")
        # @time opfm_ptdf, pf_ptdf, Pc_ptdf, Pcc_ptdf = SCOPF.opf(SCOPF.PCSC, system, Gurobi.Optimizer, voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed,
        #     ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_limit_multi=short, long_term_limit_multi=long);
        # @time SCOPF.solve_model!(opfm_ptdf.mod);
    end
end

function benchmark_sys(system::System, c1::Int, c2::Int)
    pf = SCOPF.DCPowerFlow(system)
    nodes = SCOPF.sort_components!(SCOPF.get_nodes(system))
    branches = SCOPF.sort_components!(SCOPF.get_branches(system))
    bx = SCOPF.get_bus_idx.(branches, [SCOPF.get_nodes_idx(nodes)])
    islands = SCOPF.island_detection(pf.B, bx[c2][1], bx[c2][2])
    island = SCOPF.find_ref_island(islands, pf.slack)
    island_b = SCOPF.find_island_branches(islands[island], pf.DA, c2)
    ptdf = similar(pf.ϕ)
    X = similar(pf.X)
    F = similar(pf.F)
    
    @time SCOPF.get_isf(pf.DA, pf.B, bx[c2], c2, pf.slack, islands[island], island_b)
    @time SCOPF.get_isf!(ptdf, X, pf.X, pf.B, pf.DA, bx[c1], c1)
    @time SCOPF.calculate_line_flows!(F, pf, bx[c1], c1)
end


test_ieee_rts()
test_big_sys()