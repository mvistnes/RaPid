# CC BY-NC-SA 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
module DCPowerFlow
using LinearAlgebra
using DataFrames
using Printf

"""
DC power flow calculation of a power system

Input: 
 - buses: a table with columns of ibus (bus number), type (PV=1, PQ=2, ref=3), and P.
 ! currently ref=0, change this!
 - branches: a table with columns of fbus (from bus number), tbus (to bus number), 
   and x (series reactance).

Output:
 - P: vector of real power into each bus (and with reactive power assemed 0 pu).
 - theta: vector of voltage angle at each bus (and with voltage magnutude assumed 1.0 pu).
"""
function dcopf(buses, branches)
    if buses.type[end] != 0 # Need the slack bus as the last bus
        for row in eachrow(buses)
            if row.type == 0
                push!(buses, row)
                deleteat!(buses, row.ibus)
                display(buses)
                break
            end
        end
    end

    H = buildH(buses, branches)
    P = copy.(buses.P[1:(end-1)])
    theta = calcangles(P, H)
    P = calcP(theta, branches)
    return P, theta
end


"""
Builds a matrix with the reactance of the lines

Input: as dcopf().

Output: an admittance matrix based on the line series reactances.
"""
function buildH(buses, branches)
    H = zeros(size(buses,1)-1, size(buses,1)-1)
    for branch in eachrow(branches)
        (f, t) = (branch.fbus, branch.tbus)
        y = 1 / branch.x
        if buses.type[f] != 0 # If the from bus is NOT a slack bus
            H[f,f] += y
            if buses.type[t] != 0 # If the to bus is NOT a slack bus
                H[t,t] += y
                H[f,t] = -y
                H[t,f] = -y
            end
        elseif buses.type[t] != 0
            H[t,t] += y
        end
    end
    return H
end


"""
Calculate voltage angles

Input:
 - P: vector of real power into each bus.
 - H: an admittance matrix based on the line series reactances.

Output: vector of voltage angle at each bus.
"""
function calcangles(P, H)
    luP = lu(H)
    return luP \ P
end


"""
Calculate active power

Input:
 - theta: vector of voltage angle at each bus.
 - branches: as in dcopf().

Output: vector of real power into each bus.
"""
function calcP(theta, branches)
    push!(theta, 0)
    P = zeros(size(theta,1))
    for branch in eachrow(branches)
        (f, t) = (branch.fbus, branch.tbus)
        p = (theta[f] - theta[t]) / branch.x
        # add/subtract the power to from/to bus
        P[f] += p
        P[t] -= p
    end
    return P
end


"""
Calculate and print the distribution factors of the branches

Input: as dcopf().

Output: a matrix of the distribution factors.
"""
function distr_factors(buses, branches)
    H = lu(buildH(buses, branches))
    a = zeros(size(branches,1), size(buses,1)-1) # Container for the distribution factors
    for branch in eachrow(branches)
        if branch.x != 0
            deltaP = zeros(size(buses,1)-1) # Container for the right side
            (f, t) = (branch.fbus, branch.tbus) 
            # (f)rom and (t)o bus at this branch
            if buses.type[f] != 0 # If the bus is NOT a slack bus
                deltaP[f] = 1 / branch.x
            end
            if buses.type[t] != 0 # If the bus is NOT a slack bus
                deltaP[t] = -1 / branch.x
            end
            a[rownumber(branch), :] = H \ deltaP # append factors to matrix
        end
    end

    # Nicely printed distribution factors
    print("Distribution factors\nBus  ")
    for e in eachrow(buses)
        if e.type != 0
            @printf("%12d", e.ibus)
        end
    end
    for branch in eachrow(branches)
        if branch.x != 0
            (f, t) = (branch.fbus, branch.tbus) 
            @printf("\nLine %d-%d: ", f, t)
            for e in a[rownumber(branch),:]
                @printf("%10.6f  ", e)
            end
        end
    end
    print("\n")
    return a
end


"""
Prints all buses with voltage magnitude and angle

Input: 
 - buses: as dcopf().
 - P: vector of real power into each bus.
 - theta: vector of voltage angle at each bus.
"""
function printsystem(buses, P, theta)
    println("\n Bus \tVoltage ang \tActive power")
    println("-------------------------------------")
    for i in 1:length(P)
        @printf("%4d \t%8.4f \t%8.4f\n", buses.ibus[i], theta[i], P[i])
    end
    println("")
end

end