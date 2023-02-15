# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
using PowerSystems
import JuMP
import Gurobi # LP, SOCP, Integer
import Test

module SCOPF

#=  
    Basic structure of the SCOPF problems
    Includes a basic formulation of a SCOPF, P-SCOPF and PC-SCOPF.
    Unit commitment and circuit breakers can be added to the problem
    making it a MILP.
    Variables can be set through the function call.
=#
include("N-1_SCOPF.jl")
export scopf, p_scopf, pc_scopf

#=
    A formulation of a PC-SCOPF where first a P-SCOPF is run,
    then a C-SCOPF is run using the results.
=#
include("short_long_SCOPF.jl")
export sl_scopf, c_scopf

#=
    A formulation of a PC_SCOPF using Benders' decomposition.
=#
include("benders.jl")
export run_benders, run_benders2, get_overload, get_islands, find_connected, calc_all_line_flows

#=
    Utilities for running and analyzing the SCOPFs.
=#
include("utils.jl")
export add_system_data_to_json, get_system, OPFmodel, opfmodel, beta, solve_model!, set_warm_start!, calc_severity, calc_line_severity


include("imml.jl")
export IMML, get_changed_angles, calculate_delta_line_flows, calculate_line_flows, get_overload, get_lodf

include("dc_power_flow.jl")
export get_net_Pᵢ, get_Pᵢ, calc_Pᵢ, build_adjacency, connectivitymatrix, calc_A, calc_X, calc_isf, calc_B, fast_calc_B, calc_D, get_isf, get_ptdf, find_overload, calculate_line_flows, calc_Pline, run_pf

include("post_process_opf.jl")
export scatterplot, make_save_plot, scatter_all, print_variabel, print_results, print_active_power, print_power_flow, print_contingency_power_flow, print_contingency_overflow

end