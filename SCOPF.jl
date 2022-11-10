# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

module SCOPF

#=  
    Basic structure of the SCOPF problems
    Includes a basic formulation of a SCOPF, P-SCOPF and PC-SCOPF.
    Unit commitment and circuit breakers can be added to the problem
    making it a MILP.
    Variables can be set through the function call.
=#
include("N-1_SCOPF.jl")
export OPFmodel, scopf, p_scopf, pc_scopf

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
export benders_scopf

#=
    Utilities for running and analyzing the SCOPFs.
=#
include("utils.jl")
export add_system_data_to_json, hot_start, comparison

end