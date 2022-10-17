using LinearAlgebra
using Printf
include("system_description.jl")
import .SystemDescription
include("dc_power_flow.jl")
import .DCPowerFlow

# NEED CONVERSION FROM PYTHON TO JULIA !!!!!!!
"""Fast Decoupled Power Flow

it_max is the maximum number of iterations allowed

lim is the convergance limit

method: Standard Fast Decoupled Power Flow = "std",
Primal Fast Decoupled Power Flow = "pri" and 
Dual Fast Decoupled Power Flow = "dual"

print_out = True -> state for each iteration"""
function fdlf(buses::Vector{Bus}, branches::Vector{Branch}, it_max::Int, lim::AbstractFloat; method = "std", print_out = false)
    numB = size(buses,1)
    sizeJ = sum(bus.type for bus in buses if bus.type != Int(ref::TypeB)) # size of jacobi matrix

    # Assuming 1 slack bus
    # building the two sub-matrices
    Ybus = makeYbus(buses, branches)
    B_id = makeB_id(buses)
    it = 1
    if method == "std"
        println("Standard Fast Decoupled Power Flow")
        dP = dQ = buildH(buses, branches)
    else if method == "pri"
        println("Primal Fast Decoupled Power Flow")
        dP = buildB(Ybus, B_id, numB, sizeJ)
        dQ = buildH(buses, branches)
    else if method == "dual"
        println("Dual Fast Decoupled Power Flow")
        dP = buildH(buses, branches)
        dQ = buildB(Ybus, B_id, numB, sizeJ)
    else
        return -1, buses, branches
    end

    while true
        buses, branches = PQ_update(buses, branches)
        if method == "std" || method == "pri"
            buses, branches, deltaT = solve_dP(buses, branches)
            buses, branches, deltaV = solve_dQ(buses, branches)
        else
            buses, branches, deltaV = solve_dQ(buses, branches)
            buses, branches, deltaT = solve_dP(buses, branches)
        end
        b = concatenate((deltaT,deltaV), axis=0)
        if print_out
            println("\n-------------------------  Iteration", it, " -------------------------")
            print_vectors(b, deltaPQ)
            print_buses(buses)
        end

        it += 1
        # check for convergance
        (convergeance(deltaPQ, lim) || it > it_max) && break
    end
    
    if it > it_max
        println("\nNo convergance in ", it, " iterations\n")
    else
        println("\nConvergance in ", it, " iterations\n")
    end
    print_jacobi_dpf(buses, deltaT, deltaV)
    return it, buses, branches
end

"""Builds a matrix with the Ybus-elements imaginary part"""
function buildB(Ybus, B_id, numB, sizeJ)
    B = zeros(Float64, sizeJ-numB-1, sizeJ-numB-1)
    for i in range(numB-1, sizeJ), j in range(numB-1, sizeJ)
        B[i-numB+1,j-numB+1] = - imag(Ybus[B_id[i],B_id[j]])
    end
    return B
end

"""Solves the dP/dVa-sub-matrix"""
function solve_dP(buses, branches, numB, sizeJ)
    b = deltaPQ[0:numB,:]
    deltaT = lu(dP) \ b
    b = concatenate((deltaT, zeros(sizeJ-numB,1)), axis=0)
    buses = TV_update(buses, b, B_id)
    return buses, branches, deltaT
end

"""Solves the dQ/dVm-sub-matrix"""
function solve_dQ(buses, branches, numB, sizeJ)
    deltaV = lu(dQ) \ deltaPQ[(numB-1):sizeJ]
    b = concatenate((zeros(numB-1,1), deltaV), axis=0)
    buses = TV_update(buses, b, B_id)
    return buses, branches, deltaV
end

"""Make a jacobi matrix out of the sub-matrices"""
function print_jacobi_dpf(dP, dQ)
    jac1 = concatenate((dP,zeros(size(dP))), axis=1)
    jac2 = concatenate((zeros(size(dQ)),dQ), axis=1)
    jac = concatenate((jac1,jac2), axis=0)
    print_jacobi(jac)
end  

