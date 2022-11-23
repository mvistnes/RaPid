# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
using PowerSystems
import JuMP
import Gurobi # LP, SOCP, Integer
import Test
include("utils.jl")
include("N-1_SCOPF.jl")
include("short_long_SCOPF.jl")
include("benders.jl")

# Initialize the system and run a base-case scopf
threebus = System(joinpath("data","threebus.json"))
voll = JuMP.Containers.DenseAxisArray(
    [rand(1000:3000, length(get_components(StaticLoad, threebus))); rand(1:30, length(get_components(RenewableGen, threebus)))], 
    [get_name.(get_components(StaticLoad, threebus)); get_name.(get_components(RenewableGen, threebus))]
)
opfm = scopf(SC, threebus, Gurobi.Optimizer, voll=voll)
solve_model!(opfm.mod)

prob = JuMP.Containers.DenseAxisArray([0.24, 0.51, 0.33],
    get_name.(get_components(ACBranch, threebus))
)
prob /= 8760
contingencies = get_name.(branches(threebus))

run_benders(opfm, voll, contingencies, prob)
print_results(opfm)