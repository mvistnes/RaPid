# CC BY-NC-SA 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

using DataFrames
import CSV

include("dc_power_flow.jl")
import .DCPowerFlow
include("linear_programming.jl")
import .LinearProgramming
include("system_description.jl")
import .SystemDescription

function import4area()
    buses = DataFrame(CSV.File(
        joinpath("\\\\home.ansatt.ntnu.no","matiaskv","Documents","03_Programming","4-area_test_network","4area_network_bus.csv")))
    rename!(buses, Dict(:bus_i => :ibus))
    branches = DataFrame(CSV.File(
        joinpath("\\\\home.ansatt.ntnu.no","matiaskv","Documents","03_Programming","4-area_test_network","4area_network_branch.csv"))) 
    gendata = DataFrame(CSV.File(
        joinpath("\\\\home.ansatt.ntnu.no","matiaskv","Documents","03_Programming","4-area_test_network","4area_network_gendata.csv"), 
        header=2, drop=[1])) 
    rename!(gendata, Dict(:Column2 => :ibus))
    loaddata = DataFrame(CSV.File(
        joinpath("\\\\home.ansatt.ntnu.no","matiaskv","Documents","03_Programming","4-area_test_network","4area_network_loaddata.csv"), 
        header=2, drop=[1])) 
    rename!(loaddata, Dict(:Column2 => :ibus, :Column3 => :type, :Column4 => :cost))

    buses[!, :Pd] = convert.(Float64, buses[!, :Pd])
    for gen in eachrow(gendata)
        for bus in eachrow(buses)
            if bus.ibus == gen.ibus
                bus.Pd += gen[2]
            end
        end
    end
    for load in eachrow(loaddata)
        for bus in eachrow(buses)
            if bus.ibus == load.ibus
                bus.Pd -= load[4]
            end
        end
    end
    return buses, branches
end

function test4area()
    buses, branches = import4area()
    buses, branches = DCPowerFlow.dcopf!(buses, branches)
    # display(buses)
    a = DCPowerFlow.distr_factors(buses, branches)
    limits = [1.5 for _ in 1:size(branches,1)]
    LinearProgramming.simplexmethod(a, DataFrames.getvector(buses.Pd), limits, false)
end

function testELK14A4()
    buses, branches = SystemDescription.importXLSX("System_data_ELK14_A4.xlsx")
    buses, branches = DCPowerFlow.dcopf!(buses, branches)
    a = DCPowerFlow.distr_factors(buses, branches)
    a = hcat(a, zeros(4))
    a = vcat(a, ones(1,4)) # part one of load balance equal constraint
    a = vcat(a, -ones(1,4)) # part two
    limits = [0.96, 0.64, -0.64, 2.1, 3.1, -3.1]
    a[2:4,:] .*= -1 # bigger than constraints
    cost = [4.0, 5.0, 3.0, 2.0]
    x, val = LinearProgramming.simplexmethod(a, limits, cost, false)
    display(x)
    println(val)
end

testELK14A4()