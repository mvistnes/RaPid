# CC BY-NC-SA 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

module NewtonRaphson

using LinearAlgebra
using Printf


"""
The Newton-Raphson method

it_max is the maximum number of iterations allowed

lim is the convergance limit

print_out = True -> state for each iteration
"""
function nr_method(system, it_max::Int64, lim::Float64, print_out::Bool = false)
    num_bus = length(system.buses) # number of buses
    size_J = 0 # size of jacobi matrix
    for bus in system.buses
        size_J += bus.type
    end
    deltaPQ = zeros(Float64, size_J)

    # needed when the slack bus is not the last bus number
    B_id = zeros(Int64, size_J) # bus number to matrix position
    i = 0 # counts bus numbers
    c = 0 # counts bus not slack
    c2 = 0 # counts PQ buses
    for bus in system.buses
        if bus.type == 1 || bus.type == 2
            B_id[c] = i
            c += 1
            if bus.type == 2
                B_id[c2+num_bus-1] = i
                c2 += 1
            end
        end
        i += 1
    end

    it = 1
    jac = zeros(Float64, size_J,size_J)
    while true
        system = PQ_update(system)
        build_jacobi(system)
        deltaTV = linear_solve(deltaPQ, jac)
        system.buses = TV_update(system.buses, deltaTV)
        if print_out
            print("\n-------------------------  Iteration", it, " -------------------------")
            print_jacobi(jac)
            print_vectors(deltaTV, deltaPQ)
            print_buses(system.buses)
        end

        # convergance and iteration check
        (convergeance(deltaPQ, lim) || it >= it_max) && break
        it += 1
    end

    if !print_out
        print_buses(system.buses)
    end
    if it >= it_max
        print("\nNo converge in ", it, " iterations\n")
        print_jacobi(jac)
        print_vectors(deltaTV, deltaPQ)
    else
        print("\nConvergance in ", it, " iterations\n")
    end
    return it, system
end

function print_buses(buses)
    println(" Bus \tVoltage mag \tVoltage ang \tActive power \tReactive power")
    println("-------------------------------------------------------------------")
    [println("% 4d", " \t% 8.4f"^4, buses[i].id, buses[i].vmag, buses[i].vang, 
            buses[i].P, buses[i].Q) for i in range(1,length(buses))];
end

function print_vectors(deltaTV, deltaPQ)
    println("\nMismatch vector: deltaP and deltaQ")
    [print("% 8.6f    ", i) for i in deltaPQ];
    print("\n\nCorrection vector: deltaTheta and deltaV")
    [print("% 8.6f    ", i) for i in deltaTV];
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
Tij(Gij, Bij, thetai, thetaj) = return Gij*cos(thetai-thetaj) + Bij*sin(thetai-thetaj)
Uij(Gij, Bij, thetai, thetaj) = return Gij*sin(thetai-thetaj) - Bij*cos(thetai-thetaj)


"""Power equation at the bus"""
function power(buses, bus)
    p = 0
    q = 0
    for i in range(num_bus)
        # g = real(Ybus[bus.id-1,i])
        # b = imag(Ybus[bus.id-1,i])
        y = yij(bus.id,i)
        g = real(y)
        b = imag(y)
        if g != 0 # If there is a line
            if i == bus.id
                p += buses[i].vmag * g
                q -= buses[i].vmag * b
            else
                p += buses[i].vmag * Tij(g, b, bus.vang, buses[i].vang)
                q += buses[i].vmag * Uij(g, b, bus.vang, buses[i].vang)
            end
        end
    end
    return p*bus.vmag, q*bus.vmag
end


"""Calculate the power and the power difference from initial at each bus"""
function PQ_update(system)
    #Calculating the power at each bus
    for i in range(num_bus)
        system.buses[i].P, system.buses[i].Q = power(system, system.buses[i])
    end
    
    #Find the differences from initial
    for i in range(num_bus-1)
        deltaPQ[i,0] = system.buses[B_id[i]].P_init - system.buses[B_id[i]].P
    end
    for i in range(num_bus-1, size_J)
        deltaPQ[i,0] = system.buses[B_id[i]].Q_init - system.buses[B_id[i]].Q
    end
    return system
end


"""Update the voltages and angles from deltaTV"""
function TV_update(buses, deltaTV)
    for i in range(num_bus-1)
        buses[B_id[i]].vang += deltaTV[i,0]
    end
    for i in range(num_bus-1, size_J)
        buses[B_id[i]].vmag += deltaTV[i,0]
    end
    return buses
end


"""Jacobi matrix calculation"""
function build_jacobi(system)
    sb = system.buses
    # all sub-matrices have different equations for the diagonal
    # elements and the off-diagonal elements

    # J1: dP/dTheta
    for i in range(num_bus-1)
        I = B_id[i]
        for j in range(num_bus-1)
            J = B_id[j]
            jac[i,j] = ifelse(I == J, 
                -sb[I].Q - np.imag(system.yij(I,I)) * sb[I].vmag^2, 
                sb[I].vmag * sb[J].vmag * Uij(np.real(system.yij(I,J)), 
                    np.imag(system.yij(I,J)), sb[I].vang, sb[J].vang)
            )
        end
    end
    
    # J2: dP/dV
    for i in range(num_bus-1)
        I = B_id[i]
        for j in range(num_bus-1, size_J)
            J = B_id[j]
            jac[i,j] = ifelse(I == J, 
                sb[I].P / sb[I].vmag + np.real(system.yij(I,I)) * sb[I].vmag,
                sb[I].vmag * Tij(np.real(system.yij(I,J)), np.imag(system.yij(I,J)), 
                    sb[I].vang, sb[J].vang)
            )
        end
    end

    # J3: dQ/dTheta
    for i in range(num_bus-1, size_J)
        I = B_id[i]
        for j in range(num_bus-1)
            J = B_id[j]
            jac[i,j] = ifelse(I == J, 
                sb[I].P - np.real(system.yij(I,I)) * sb[I].vmag^2,
                -sb[I].vmag * sb[J].vmag * Tij(np.real(system.yij(I,J)), 
                    np.imag(system.yij(I,J)), sb[I].vang, sb[J].vang)
            )
        end
    end

    # J4: dQ/dV
    for i in range(num_bus-1, size_J)
        I = B_id[i]
        for j in range(num_bus-1, size_J)
            J = B_id[j]
            jac[i,j] = ifelse(I == J, 
                sb[I].Q / sb[I].vmag - np.imag(system.yij(I,I)) * sb[I].vmag,
                sb[I].vmag * Uij(np.real(system.yij(I,J)), np.imag(system.yij(I,J)), 
                    sb[I].vang, sb[J].vang)
            )
        end
    end
end


"""Checks if the difference in power is small enough"""
function convergeance(deltaPQ, lim)
    for i in deltaPQ
        not(i > -lim && i < lim) && return False
    end
    return True
end


"""Solves the linear problem Ax=b """
function linear_solve(b, A)
    lu, piv = lu_factor(A)
    x = lu_solve((lu,piv), b)
    # A_inv = np.linalg.inv(A)
    # x = np.dot(A_inv, b)
    return x
end

end