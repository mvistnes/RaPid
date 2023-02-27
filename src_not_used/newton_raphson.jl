# CC BY-NC-SA 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

module NewtonRaphson

using LinearAlgebra
using Printf
include("system_description.jl")
import .SystemDescription


"""
The Newton-Raphson method

it_max is the maximum number of iterations allowed

lim is the convergance limit

print_out = True -> state for each iteration
"""
function nr_method!(nodes::DataFrame, branches::DataFrame, it_max::UInt64, lim::Float64; print_out::Bool = false)
    nodes_vec, branches_vec, convert_bus = vectorize(nodes, branches)
    it, nodes_vec, branches_vec = nr_method!(nodes_vec, branches_vec, it_max, lim, print_out)
    nodes[!, :P] = zeros(size(nodes,1))
    nodes[!, :Q] = zeros(size(nodes,1))
    nodes[!, :Vm] = zeros(size(nodes,1))
    nodes[!, :Va] = zeros(size(nodes,1))
    for bus in eachrow(nodes)
        bus.P = P_bus[convert_bus[bus.ibus]]
        bus.Q = Q_bus[convert_bus[bus.ibus]]
        bus.Vm = Vm[convert_bus[bus.ibus]]
        bus.Va = Va[convert_bus[bus.ibus]]
    end
    return nodes, branches
end
function nr_method!(nodes::Vector{Bus}, branches::Vector{Branch}, it_max::UInt64, lim::Float64; print_out::Bool = false)
    numB = length(nodes) # number of nodes
    B_id, sizeJ = makeB_id(nodes)
    deltaPQ = zeros(Float64, sizeJ)
    PQ_init = [nodes[B_id[i]].P for i in range(numB-1); nodes[B_id[i]].Q for i in range(numB-1, sizeJ)]
    V = [nodes[B_id[i]].Va for i in range(numB-1); nodes[B_id[i]].Vm for i in range(numB-1, sizeJ)]

    it = 1
    jac = zeros(Float64, sizeJ,sizeJ)
    Ybus = makeYbus(nodes, branches)
    while true
        nodes, branches = PQ_update(Ybus, nodes, branches, B_id)
        jac = build_jacobi(Ybus, jac, nodes, branches, B_id)
        deltaV = linear_solve(deltaPQ, jac)
        nodes = V_update(nodes, deltaV, B_id)
        if print_out
            print("\n-------------------------  Iteration", it, " -------------------------")
            print_jacobi(jac)
            print_vectors(deltaV, deltaPQ)
            print_nodes(nodes)
        end

        # convergance and iteration check
        (convergeance(deltaPQ, lim) || it >= it_max) && break
        it += 1
    end

    if !print_out
        print_nodes(nodes)
    end
    if it >= it_max
        print("\nNo converge in ", it, " iterations\n")
        print_jacobi(jac)
        print_vectors(deltaV, deltaPQ)
    else
        print("\nConvergance in ", it, " iterations\n")
    end
    return it, nodes, branches
end

function print_nodes(nodes)
    println(" Bus \tVoltage mag \tVoltage ang \tActive power \tReactive power")
    println("-------------------------------------------------------------------")
    [println("% 4d", " \t% 8.4f"^4, nodes[i].id, nodes[i].Vm, nodes[i].Va, 
            nodes[i].P, nodes[i].Q) for i in range(1,length(nodes))];
end

function print_vectors(deltaV, deltaPQ)
    println("\nMismatch vector: deltaP and deltaQ")
    [print("% 8.6f    ", i) for i in deltaPQ];
    print("\n\nCorrection vector: deltaTheta and deltaV")
    [print("% 8.6f    ", i) for i in deltaV];
    println("\n")
end

function print_jacobi(jac)
    println("Jacobi matrix")
    for row in jac
        [print("% 8.4f    ", i) for i in row];
        println("")
    end
    println("")
end

#Intermediate calculations
Tij(Gij, Bij, Vai, Vaj) = return Gij*cos(Vai-Vaj) + Bij*sin(Vai-Vaj)
Uij(Gij, Bij, Vai, Vaj) = return Gij*sin(Vai-Vaj) - Bij*cos(Vai-Vaj)

"""Power equation at the bus"""
function power(Ybus, nodes, bus)
    p = 0
    q = 0
    for i in range(numB)
        y = Ybus[bus.id,i]
        g = real(y)
        b = imag(y)
        if g != 0 # If there is a line
            if i == bus.id
                p += nodes[i].Vm * g
                q -= nodes[i].Vm * b
            else
                p += nodes[i].Vm * Tij(g, b, bus.Va, nodes[i].Va)
                q += nodes[i].Vm * Uij(g, b, bus.Va, nodes[i].Va)
            end
        end
    end
    return p*bus.Vm, q*bus.Vm
end


"""Calculate the power and the power difference from initial at each bus"""
function PQ_update(Ybus, nodes, branches, B_id)
    #Calculating the power at each bus
    for i in range(numB)
        nodes[i].P, nodes[i].Q = power(Ybus, nodes, nodes[i])
    end
    
    #Find the differences from initial
    for i in range(numB-1)
        deltaPQ[i,0] = nodes[B_id[i]].P_init - nodes[B_id[i]].P
    end
    for i in range(numB-1, sizeJ)
        deltaPQ[i,0] = nodes[B_id[i]].Q_init - nodes[B_id[i]].Q
    end
    return nodes, branches
end
function PQ_update(Ybus, PQ, PQ_init)
    for i in 1:numB
        PQ[i], PQ[i+numB-1] = power(Ybus, PQ, i)
    end
    deltaPQ = PQ_init .- PQ
    return PQ, deltaPQ
end


"""Update the voltages and angles from deltaV"""
function V_update(nodes, deltaV, B_id)
    for i in range(numB-1)
        nodes[B_id[i]].Va += deltaV[i,0]
    end
    for i in range(numB-1, sizeJ)
        nodes[B_id[i]].Vm += deltaV[i,0]
    end
    return nodes
end
function V_update(V, deltaV)
    return V .+ deltaV
end

"""Jacobi matrix calculation"""
function build_jacobi(Ybus, jac, nodes, branches, B_id)
    # all sub-matrices have different equations for the diagonal
    # elements and the off-diagonal elements
    numB = length(nodes) # number of nodes

    # J1: dP/dTheta
    for i in range(numB-1)
        I = B_id[i]
        for j in range(numB-1)
            J = B_id[j]
            jac[i,j] = ifelse(i == j, 
                -nodes[I].Q - imag(Ybus[I,I]) * nodes[I].Vm^2, 
                nodes[I].Vm * nodes[J].Vm * Uij(real(Ybus[I,J]), 
                    imag(Ybus[I,J]), nodes[I].Va, nodes[J].Va)
            )
        end
    end
    
    # J2: dP/dV
    for i in range(numB-1)
        I = B_id[i]
        for j in range(numB-1, sizeJ)
            J = B_id[j]
            jac[i,j] = ifelse(I == J, 
                nodes[I].P / nodes[I].Vm + real(Ybus[I,I]) * nodes[I].Vm,
                nodes[I].Vm * Tij(real(Ybus[I,J]), imag(Ybus[I,J]), 
                    nodes[I].Va, nodes[J].Va)
            )
        end
    end

    # J3: dQ/dTheta
    for i in range(numB-1, sizeJ)
        I = B_id[i]
        for j in range(numB-1)
            J = B_id[j]
            jac[i,j] = ifelse(I == J, 
                nodes[I].P - real(Ybus[I,I]) * nodes[I].Vm^2,
                -nodes[I].Vm * nodes[J].Vm * Tij(real(Ybus[I,J]), 
                    imag(Ybus[I,J]), nodes[I].Va, nodes[J].Va)
            )
        end
    end

    # J4: dQ/dV
    for i in range(numB-1, sizeJ)
        I = B_id[i]
        for j in range(numB-1, sizeJ)
            J = B_id[j]
            jac[i,j] = ifelse(I == J, 
                nodes[I].Q / nodes[I].Vm - imag(Ybus[I,I]) * nodes[I].Vm,
                nodes[I].Vm * Uij(real(Ybus[I,J]), imag(Ybus[I,J]), 
                    nodes[I].Va, nodes[J].Va)
            )
        end
    end
end


"""Checks if the difference in power is small enough"""
function convergeance(deltaPQ, lim)
    for i in deltaPQ
        not(i > -lim && i < lim) && return false
    end
    return true
end


"""Solves the linear problem Ax=b """
function linear_solve(b, A)
    return lu(A) \ b
end

end