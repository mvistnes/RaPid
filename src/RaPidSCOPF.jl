# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2023

module RaPidSCOPF

using PowerSystems
using PowerModels
using DataFrames
using JuMP
# using Plots
# using StatsPlots
using Printf
# using PowerSystemCaseBuilder
# using Test
# using DelimitedFiles

import Statistics: mean
import Logging
import SparseArrays
# import StaticArrays
# import Octavian
# import SuiteSparseGraphBLAS
import MathOptInterface
const MOI = MathOptInterface
import KLU
import LinearAlgebra
import Dates

import Base: +, -

#=
    Definition of common types.
=#
include("types.jl")
export OPF, ExprC, ExprCC, ExprCCX

#=
    Utilities for running and analyzing the SCOPFs.
=#
include("utils.jl")
export add_system_data_to_json, get_system, OPFmodel, opfmodel, beta, solve_model!, set_warm_start!,
    calc_severity, calc_line_severity, get_branches, get_nodes, sort_components!,
    get_nodes_idx, get_bus_idx, find_slack, get_net_Pᵢ, get_Pᵢ

#=
    DC power flow functions.
=#
include("dc_power_flow.jl")
export DCPowerFlow, calc_Pᵢ, build_adjacency, connectivitymatrix, calc_A, calc_X, calc_isf,
    calc_B, calc_D, get_isf, get_ptdf, find_overload, filter_overload, calculate_line_flows,
    calculate_line_flows!, calc_Pline, run_pf

include("cont_dc_power_flow.jl")
export get_isf, get_ptdf, calculate_line_flows, calculate_line_flows!
    
#=
    AC power flow functions.
=#
include("ac_power_flow.jl")
export calc_Ybus, calc_B_mark, calc_B_doublemark

#=
    Functions using the Inverse Matrix Modification Lemma.
=#
include("imml.jl")
export IMML, get_changed_angles, calc_Pline, get_changed_X, get_isf, calculate_line_flows,
    calculate_line_flows!, get_overload, get_lodf

include("FDLF.jl")

#=
    Functions for handling islands in the system.
=#
include("island_handling.jl")
export handle_islands, island_detection_thread_safe, get_all_islands

#=  
    Basic structure of the SCOPF problems
    Includes a basic formulation of a SCOPF, P-SCOPF and PC-SCOPF.
    Unit commitment and circuit breakers can be added to the problem
    making it a MILP.
    Variables can be set through the function call.
=#
# include("N-1_SCOPF.jl")
# export scopf, p_scopf, pc_scopf

# include("AC_N-1_SCOPF.jl")
# export ac_scopf

include("OPF_PTDF.jl")
export opf

# #=
#     A formulation of a PC-SCOPF where first a P-SCOPF is run,
#     then a C-SCOPF is run using the results.
# =#
# include("short_long_SCOPF.jl")
# export sl_scopf, c_scopf

#=
    A formulation of a PC_SCOPF using Benders' decomposition.
=#
include("benders.jl")
export run_benders, find_system_state!, update_model!, find_overloads, print_benders_results

include("contingency_select.jl")
export run_contingency_select

#=
    Functions for running reliability calculations.
=#
include("reliability.jl")
export run_reliability_calculation

#=
    Functions for printing and plotting SCOPF results.
=#
include("post_process_opf.jl")
export print_variabels, print_corrective_values, print_results, print_contingency_results, 
    print_contingency_P, print_active_power, get_power_flow, print_power_flow, get_contingency_power_flow,
    print_contingency_power_flow, get_generation_results, print_generation_results

end
