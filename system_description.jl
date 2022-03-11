# CC BY-NC-SA 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
module SystemDescription
import XLSX
using DataFrames
using SparseArrays

"""Make the admittance bus matrix"""
function makeYbus(buses, branches)
    Ybus = spzeros(ComplexF32, size(buses,1), size(buses,1))

    
    # Calculating the diagonal elements
    for i in 1:size(buses,1)
        for branch in eachrow(branches)
            if i == branch.fbus || i == branch.tbus
                Ybus[i,i] += 1/(branch.r + branch.x*im) + branch.b*im / 2
            end
        end
    end

    # Calculating the off-diagonal elements
    for branch in eachrow(branches)
        (f, t) = (branch.fbus, branch.tbus)
        Ybus[f, t] = Ybus[t, f] = -1 / (branch.r + branch.x*im)
    end
    return Ybus
end

"""Return the yij element of the Ybus, 1-indexed"""
function yij(branches, i, j)
    y = 0
    if i == j
        # Calculating the diagonal element
        for branch in eachrow(branches)
            if i == branch.fbus || i == branch.tbus
                y += 1/(branch.r + branch.x*im) + branch.b*im / 2
            end
        end
    else
        # Calculating the off-diagonal element
        for branch in eachrow(branches)
            (f, t) = (branch.fbus, branch.tbus)
            if (f == i && t == j) || (f == j && t == i)
                y = -1 / (branch.r + branch.x*im)
                break
                rem
            end
        end
    end
    return y
end


"""
Import a Power System description to two DataFrames
(Bus and Branch) on the MatPower format
"""
function importXLSX(filename)
    buses = DataFrame(XLSX.readtable(filename, "BusData")...)
    rename!(buses, [:ibus, :type, :vmag, :vang, :P, :Q, :Pmax, :Qmax])
    branches = DataFrame(XLSX.readtable(filename, "BranchData")...)
    rename!(branches, [:fbus, :tbus, :r, :x, :b])
    return buses, branches
end

# buses, branches = importXLSX("System_data_69bus.xlsx")
# @time Ybus = makeYbus(buses,branches)
# display(Ybus)
# @benchmark DCPowerFlow.dcopf(buses, branches) setup=(buses, branches = SystemDescription.importXLSX("System_data_69bus.xlsx"))
    
end