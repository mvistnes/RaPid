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
function fdlf(nodes::Vector{Bus}, branches::Vector{Branch}, it_max::Int, lim::AbstractFloat; method = "std", print_out = false)
    numB = size(nodes,1)
    sizeJ = sum(bus.type for bus in nodes if bus.type != Int(ref::TypeB)) # size of jacobi matrix

    # Assuming 1 slack bus
    # building the two sub-matrices
    Ybus = makeYbus(nodes, branches)
    B_id = makeB_id(nodes)
    it = 1
    if method == "std"
        println("Standard Fast Decoupled Power Flow")
        dP = dQ = buildH(nodes, branches)
    else if method == "pri"
        println("Primal Fast Decoupled Power Flow")
        dP = buildB(Ybus, B_id, numB, sizeJ)
        dQ = buildH(nodes, branches)
    else if method == "dual"
        println("Dual Fast Decoupled Power Flow")
        dP = buildH(nodes, branches)
        dQ = buildB(Ybus, B_id, numB, sizeJ)
    else
        return -1, nodes, branches
    end

    while true
        nodes, branches = PQ_update(nodes, branches)
        if method == "std" || method == "pri"
            nodes, branches, deltaT = solve_dP(nodes, branches)
            nodes, branches, deltaV = solve_dQ(nodes, branches)
        else
            nodes, branches, deltaV = solve_dQ(nodes, branches)
            nodes, branches, deltaT = solve_dP(nodes, branches)
        end
        b = concatenate((deltaT,deltaV), axis=0)
        if print_out
            println("\n-------------------------  Iteration", it, " -------------------------")
            print_vectors(b, deltaPQ)
            print_nodes(nodes)
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
    print_jacobi_dpf(nodes, deltaT, deltaV)
    return it, nodes, branches
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
function solve_dP(nodes, branches, numB, sizeJ)
    b = deltaPQ[0:numB,:]
    deltaT = lu(dP) \ b
    b = concatenate((deltaT, zeros(sizeJ-numB,1)), axis=0)
    nodes = TV_update(nodes, b, B_id)
    return nodes, branches, deltaT
end

"""Solves the dQ/dVm-sub-matrix"""
function solve_dQ(nodes, branches, numB, sizeJ)
    deltaV = lu(dQ) \ deltaPQ[(numB-1):sizeJ]
    b = concatenate((zeros(numB-1,1), deltaV), axis=0)
    nodes = TV_update(nodes, b, B_id)
    return nodes, branches, deltaV
end

"""Make a jacobi matrix out of the sub-matrices"""
function print_jacobi_dpf(dP, dQ)
    jac1 = concatenate((dP,zeros(size(dP))), axis=1)
    jac2 = concatenate((zeros(size(dQ)),dQ), axis=1)
    jac = concatenate((jac1,jac2), axis=0)
    print_jacobi(jac)
end  

