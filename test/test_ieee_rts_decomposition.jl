@testset "Test SCOPF decomposition" begin
    
# open("output.txt", "w") do out
#     redirect_stdout(out) do

system = System(joinpath("cases","IEEE_RTS.m"))
# SCOPF.fix_generation_cost!(system);
# nodes = SCOPF.sort_components!(SCOPF.get_nodes(system))
# idx = SCOPF.get_nodes_idx(nodes)
# SCOPF.set_ramp_limits!(system, 0.01)
SCOPF.set_rating!.(SCOPF.get_branches(system), SCOPF.get_rating.(SCOPF.get_branches(system)) * 0.8);
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
max_shed = 1.
ramp_mult = 2.

function test_decomposition(system, optimizer, voll, contingencies, prob, max_shed,ramp_mult, ramp_minutes, short, long)
    for type in [SCOPF.Base_SCOPF, SCOPF.P_SCOPF, SCOPF.OPF(true, false, true, false, false), SCOPF.PC_SCOPF, SCOPF.PCF_SCOPF, SCOPF.PC2_SCOPF, SCOPF.PC2F_SCOPF]
        p_failure = ifelse(type.C2F, 0.10, 0.00)
        println(type)
        println("Start PTDF")
        mod, opf, pf, oplim, brc_up, brc_down, Pc, Pcc, Pccx = SCOPF.opf_base(type, system, HiGHS.Optimizer(), voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed,
            ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure);
        SCOPF.add_branch_constraints!(mod, pf.ϕ, mod[:p0], brc_up, brc_down, oplim.branch_rating)
        SCOPF.solve_model!(mod);

        # println("Start Contingency select")
        # case_cont, tot_t = SCOPF.run_contingency_select(type, system, HiGHS.Optimizer(), voll, prob, contingencies, 
        #     max_shed=max_shed, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure);

        println("Start SCOPF decomposition")
        case, tot_t = SCOPF.run_decomposed_optimization(type, system, HiGHS.Optimizer(), voll, prob, contingencies, max_shed=max_shed,
            ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure);

        # @test JuMP.objective_value(mod) ≈ JuMP.objective_value(case_cont.model)
        @test JuMP.objective_value(mod) ≈ JuMP.objective_value(case.model)
    end
    return
end

test_decomposition(system, HiGHS.Optimizer(), voll, contingencies, prob, max_shed,ramp_mult, ramp_minutes, short, long)

function decomposition_pc_scopf()
    println("\nPreventive-Corrective SCOPF")
    println("Start PTDF")
    @time case_ptdf = Case(SCOPF.opf_base(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GUROBI_ENV), voll=voll, contingencies=contingencies, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long)...)
    SCOPF.add_branch_constraints!(case_ptdf.model, case_ptdf.pf.ϕ, case_ptdf.model[:p0], case_ptdf.brc_up, case_ptdf.brc_down, case_ptdf.oplim.branch_rating)
    @time SCOPF.solve_model!(case_ptdf.model);
    println("Objective value PTDF: ", JuMP.objective_value(case_ptdf.model))
    SCOPF.print_contingency_results(case_ptdf)

    println("Start Contingency select")
    @time case_cont, tot_t = SCOPF.run_contingency_select(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GUROBI_ENV), voll, prob, contingencies, 
        max_shed=max_shed, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00);
    println("Solver time: ", tot_t, " Objective value Contingency select: ", JuMP.objective_value(case_cont.model))
    SCOPF.print_contingency_results(case_cont)

    println("Start SCOPF decomposition")
    @time case, tot_t = SCOPF.run_decomposed_optimization(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GUROBI_ENV), voll, prob, contingencies, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00);
    println("Solver time: ", tot_t, " Objective value SCOPF decomposition: ", JuMP.objective_value(case.model))
    SCOPF.print_contingency_results(case)

    @test JuMP.objective_value(case_ptdf.model) ≈ JuMP.objective_value(case_cont.model)
    @test JuMP.objective_value(case_ptdf.model) ≈ JuMP.objective_value(case.model)

    SCOPF.print_power_flow(case.opf, case.model)
    bx = SCOPF.typesort_component.(SCOPF.get_interarea(case.opf.branches), [case.opf])
    SCOPF.print_contingency_power_flow(case; subset=first.(bx))
end

end