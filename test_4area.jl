# CC BY-NC-SA 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
module Test4area

using DataFrames
import CSV

include("dc_power_flow.jl")
import .DCPowerFlow
# include("linear_programming.jl")
# import .LinearProgramming

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

buses = DCPowerFlow.dcopf!(buses, branches)
display(buses)
end