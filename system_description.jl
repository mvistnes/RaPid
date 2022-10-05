# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

module SystemDescription
import XLSX
using DataFrames
using SparseArrays

@enum TypeB PV=1 PQ=2 ref=3
struct Bus
    id::Int64 # Bus number
    type::TypeB
    Pd::Float64 # Active power inserted
    Qd::Float64 # Reactive power inserted
end
struct Branch
    fbus::Int64 # From bus number
    tbus::Int64 # To bus number
    r::Float64 # Line series resistance
    x::Float64 # Line series reactance
end

"""
Vectorize tables

Input:
    - buses: a table with columns of ibus (bus number), type (PV=1, PQ=2, ref=3), and Pd.
    ! currently ref=0, change this!
    - branches: a table with columns of fbus (from bus number), tbus (to bus number), 
    and x (series reactance).

Output: 
    - buses_vec, branches_vec: vectors of the input
    - convert_bus: conversion between the bus number and place in the vectors
"""
function vectorize(buses::DataFrame, branches::DataFrame)
    buses_vec = [Bus(i, TypeB(bus.type), bus.Pd, bus.Qd) for (i,bus) in enumerate(eachrow(buses))]
    slack = size(buses_vec,1)
    if buses.type[end] != Int(ref::TypeB) # Need the slack bus as the last bus
        for (i, row) in enumerate(eachrow(buses))
            if row.type == Int(ref::TypeB)
                slack = i
                buses_vec[i], buses_vec[end] = buses_vec[end], buses_vec[i]
                break
            end
        end
    end
    convert_bus = Dict()
    for (i,bus) in enumerate(eachrow(buses))
        convert_bus[bus.ibus] = ifelse(bus.type == Int(ref::TypeB), slack, i)
    end
    branches_vec = [Branch(convert_bus[branch.fbus], convert_bus[branch.tbus], branch.r, branch.x) 
                    for (i,branch) in enumerate(eachrow(branches))]
    return buses_vec, branches_vec, convert_bus
end
function move_ref_last!(buses::Vector{Bus}, branches::Vector{Bus})
    slack = size(buses,1)
    if buses.type[end] != Int(ref::TypeB) # Need the slack bus as the last bus
        for (i, row) in enumerate(eachrow(buses))
            if row.type == Int(ref::TypeB)
                slack = i
                buses[i], buses[end] = buses[end], buses[i]
                break
            end
        end
    end
    convert_bus = Dict()
    for (i,bus) in enumerate(eachrow(buses))
        convert_bus[bus.ibus] = ifelse(bus.type == Int(ref::TypeB), slack, i)
    end
    branches = [Branch(convert_bus[branch.fbus], convert_bus[branch.tbus], branch.r, branch.x) 
                    for (i,branch) in enumerate(branches)]
    return buses, branches, convert_bus
end

"""Bus to position"""
function makeB_id(buses)
    B_id = zeros(Int64, size_J) # bus number to matrix position
    c1 = 0 # counts bus not slack
    c2 = 0 # counts PQ buses
    for (i,bus) in enumerate(buses)
        if bus.type != ref::TypeB
            B_id[c1] = i
            c1 += 1
            if bus.type == ref::PQ
                B_id[c2+numB-1] = i
                c2 += 1
            end
        end
    end
    return B_id
end

"""Make the admittance bus matrix"""
function makeYbus(buses::DataFrame, branches::DataFrame)
    return makeYbus(collect(eachrow(buses)), collect(eachrow(branches)))
end
function makeYbus(buses::Vector{Bus}, branches::Vector{Branch})
    Ybus = spzeros(ComplexF32, size(buses,1), size(buses,1))
    
    # Calculating the diagonal elements
    for i in 1:size(buses,1), branch in branches
        if i == branch.fbus || i == branch.tbus
            Ybus[i,i] += 1/(branch.r + branch.x*im) + branch.b*im / 2
        end
    end

    # Calculating the off-diagonal elements
    for branch in branches
        (f, t) = (branch.fbus, branch.tbus)
        Ybus[f, t] = Ybus[t, f] = -1 / (branch.r + branch.x*im)
    end
    return Ybus
end

"""Return the yij element of the Ybus, 1-indexed"""
function yij(branches::DataFrame, i::Int, j::Int)
    return yij(collect(eachrow(branches)), i, j)
end
function yij(branches::Vector{branches}, i::Int, j::Int)
    if i == j
        # Calculating the diagonal element
        y = 0
        for branch in branches
            if i == branch.fbus || i == branch.tbus
                y += 1/(branch.r + branch.x*im) + branch.b*im / 2
            end
        end
        return y
    else
        # Calculating the off-diagonal element
        for branch in branches
            (f, t) = (branch.fbus, branch.tbus)
            if (f == i && t == j) || (f == j && t == i)
                return -1 / (branch.r + branch.x*im)
            end
        end
    end
end


"""
Import a Power System description to two DataFrames
(Bus and Branch) on the MatPower format
"""
function importXLSX(filename::String)
    buses = DataFrame(XLSX.readtable(filename, "BusData")...)
    rename!(buses, [:ibus, :type, :vmag, :vang, :Pd, :Qd, :Pmax, :Qmax])
    branches = DataFrame(XLSX.readtable(filename, "BranchData")...)
    rename!(branches, [:fbus, :tbus, :r, :x, :b])
    return buses, branches
end

# buses, branches = importXLSX("System_data_69bus.xlsx")
# @time Ybus = makeYbus(buses,branches)
# display(Ybus)

# using BenchmarkTools
# @benchmark DCPowerFlow.dcopf!(buses, branches) setup=(buses, branches = SystemDescription.importXLSX("System_data_69bus.xlsx"))
    
end