using DataFrames
using PowerSystems
import CSV
using PowerSystemCaseBuilder

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

path = joinpath("\\\\home.ansatt.ntnu.no","matiaskv","Documents","03_Programming","4-area_test_network")
data = PowerSystemTableData(path, 100.0, joinpath(path, "./user_descriptors.yaml"))
system_data = System(data)
to_json(system_data, "4area.json", force=true)
