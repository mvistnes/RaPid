# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

module RaPidSCOPF

using PowerSystems
using JuMP
using Plots
using StatsPlots
using Printf
# using PowerSystemCaseBuilder
using Test
using Statistics
using DelimitedFiles

import Logging
import SparseArrays
# import StaticArrays
# import Octavian
# import SuiteSparseGraphBLAS
import MathOptInterface
import KLU
const MOI = MathOptInterface
import LinearAlgebra
import Ipopt # LP, SOCP, Nonconvex
import Gurobi # LP, SOCP, Integer
import GLPK # LP, Integer
# import Tulip
import Test

## For timing functions (with allocations)
# using TimerOutputs
# const tmr = TimerOutput();
## SCOPF.tmr
## SCOPF.reset_timer!(SCOPF.tmr)

@enum OPF SC = 0 PSC = 1 PCSC = 2 PCFSC = 3 # SCOPF, P-SCOPF, PC-SCOPF

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
    Functions using the Inverse Matrix Modification Lemma.
=#
include("imml.jl")
export IMML, get_changed_angles, calc_Pline, get_changed_X, get_isf, calculate_line_flows,
    calculate_line_flows!, get_overload, get_lodf

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
export scatterplot, make_save_plot, scatter_all, print_variabel, print_results, print_active_power,
    print_power_flow, print_contingency_power_flow, print_contingency_overflow

end
