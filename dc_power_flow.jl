# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

module DCPowerFlow
using LinearAlgebra
using DataFrames
using Printf

include("system_description.jl")
import .SystemDescription

"""
DC power flow calculation of a power system.
No losses are accounted for.

Input: 
    - buses: a table with columns of ibus (bus number), type (PV=1, PQ=2, ref=3), and Pd.
    - branches: a table with columns of fbus (from bus number), tbus (to bus number), 
    and x (series reactance).

Output: 
    - buses: extended with
        - Pd: vector of real power into each bus (and with reactive power assemed 0 pu).
        - theta: vector of voltage angle at each bus (and with voltage magnitude assumed 1.0 pu).
"""
function dcpf!(buses::DataFrame, branches::DataFrame)
    buses_vec, branches_vec, convert_bus = vectorize(buses, branches)
    P_bus, theta, branches[!, :P] = dcpf!(buses_vec, branches_vec)
    buses[!, :P] = zeros(size(buses,1))
    buses[!, :theta] = zeros(size(buses,1))
    for bus in eachrow(buses)
        bus.P = P_bus[convert_bus[bus.ibus]]
        bus.theta = theta[convert_bus[bus.ibus]]
    end
    return buses, branches
end
function dcpf!(buses::Vector{Bus}, branches::Vector{Branch})
    H = buildH(buses, branches)
    P = [buses[i].Pd for i in range(1,length(buses)-1)]
    theta = calcangles(P, H)
    P_bus, P_branch = calcP(theta, branches)
    return P_bus, theta, P_branch
end


"""
Builds a matrix with the reactance of the lines

Input: as dcpf!().

Output: an admittance matrix based on the line series reactances.
"""
function buildH(buses::Vector{Bus}, branches::Vector{Branch})
    H = zeros(Float64, size(buses,1)-1, size(buses,1)-1)
    for branch in branches
        (f, t) = (branch.fbus, branch.tbus)
        y = 1 / branch.x
        if buses[f].type != ref::TypeB # If the from bus is NOT a slack bus
            H[f,f] += y
            if buses[t].type != ref::TypeB # If the to bus is NOT a slack bus
                H[t,t] += y
                H[f,t] = -y
                H[t,f] = -y
            end
        elseif buses[t].type != ref::TypeB
            H[t,t] += y
        end
    end
    return H
end


"""
Calculate voltage angles

Input:
 - Pd: vector of real power into each bus.
 - H: an admittance matrix based on the line series reactances.

Output: vector of voltage angle at each bus.
"""
function calcangles(P::Vector{AbstractFloat}, H::Matrix{AbstractFloat})
    luP = lu(H)
    return luP \ P
end


"""
Calculate active power

Input:
 - theta: vector of voltage angle at each bus.
 - branches: as in dcpf!().

Output: 
 - P_bus: vector of real power into each bus.
 - P_branch: vector of real power in a line.
"""
function calcP(theta::Vector{AbstractFloat}, branches::Vector{Branch})
    push!(theta, 0.0)
    P_bus = zeros(size(theta,1))
    P_branch = zeros(size(branches,1))
    for (i, branch) in enumerate(branches)
        (f, t) = (branch.fbus, branch.tbus)
        P_branch[i] = (theta[f] - theta[t]) / branch.x
        # add/subtract the power to from/to bus
        P_bus[f] += P_branch[i]
        P_bus[t] -= P_branch[i]
    end
    return P_bus, P_branch
end


"""
Calculate and print the distribution factors of the branches

Input: as dcpf!().

Output: a matrix of the distribution factors.
"""
function distr_factors(buses, branches)
    buses_vec, branches_vec, convert_bus = vectorize(buses, branches)
    H = lu(buildH(buses_vec, branches_vec))
    a = zeros(Float64, size(branches_vec,1), size(buses_vec,1)-1) # Container for the distribution factors
    for (i, branch) in enumerate(branches_vec)
        if branch.x != ref::TypeB
            deltaPd = zeros(Float64, size(buses_vec,1)-1) # Container for the right side
            (f, t) = (branch.fbus, branch.tbus) # (f)rom and (t)o bus at this branch
            if buses_vec[f].type != ref::TypeB # If the bus is NOT a slack bus
                deltaPd[f] = 1 / branch.x
            end
            if buses_vec[t].type != ref::TypeB # If the bus is NOT a slack bus
                deltaPd[t] = -1 / branch.x
            end
            a[i, :] = H \ deltaPd # append factors to matrix
        end
    end
    return a
end

function print_distr_factors(buses, branches, a)
    # Nicely printed distribution factors
    print("Distribution factors\nBus  ")
    for e in buses
        if e.type != ref::TypeB
            @printf("%12d", e.ibus)
        end
    end
    for (i, branch) in enumerate(branches)
        if branch.x != ref::TypeB
            @printf("\nLine %d-%d: ", branch.fbus, branch.tbus)
            for e in a[i,:]
                @printf("%10.6f  ", e)
            end
        end
    end
    print("\n")
end


"""
Prints all buses with voltage magnitude and angle

Input: 
 - buses: as dcpf!().
 - P: vector of real power into each bus.
 - theta: vector of voltage angle at each bus.
"""
function printsystem(buses, Pd::Vector{AbstractFloat}, theta::Vector{AbstractFloat})
    println("\n Bus \tVoltage ang \tActive power")
    println("-------------------------------------")
    for i in 1:length(Pd)
        @printf("%4d \t%8.4f \t%8.4f\n", buses.ibus[i], theta[i], P[i])
    end
    println("")
end

end