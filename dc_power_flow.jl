# CC BY-NC-SA 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
module DCPowerFlow
using LinearAlgebra
using DataFrames
using Formatting

"""DC power flow calculation of a power system"""
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

"""Builds a matrix with the reactance of the lines"""
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

"""Calculate voltage angles"""
function calcangles(P, H)
    luP = lu(H)
    return luP \ P
end

"""Calculate active power"""
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

"""Calculate the distribution factors of the branches"""
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
            printfmt("{:12d}", e.ibus)
        end
    end
    for branch in eachrow(branches)
        if branch.x != 0
            (f, t) = (branch.fbus, branch.tbus) 
            printfmt("\nLine {:d}-{:d}: ", f, t)
            for e in a[rownumber(branch),:]
                printfmt("{:10.6f}  ", e)
            end
        end
    end
    print("\n")
end

"""Prints all buses with voltage magnitude and angle"""
function printsystem(buses, P, theta)
    println("\n Bus \tVoltage ang \tActive power")
    println("-------------------------------------")
    for i in 1:length(P)
        printfmt("{:4d} \t{:8.4f} \t{:8.4f}\n", buses.ibus[i], theta[i], P[i])
    end
    println("")
end

end