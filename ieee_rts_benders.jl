# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
using PowerSystems
import JuMP
import Gurobi # LP, SOCP, Integer
import Test
include("N-1_SCOPF.jl")
include("utils.jl")
include("post_process_opf.jl")
include("short_long_SCOPF.jl")
include("benders.jl")

# Initialize the system and run a base-case scopf
ieee_rts = System(joinpath("data","ieee_rts.json"))
voll = JuMP.Containers.DenseAxisArray(
    [rand(1000:3000, length(get_components(StaticLoad, ieee_rts))); rand(1:30, length(get_components(RenewableGen, ieee_rts)))], 
    [get_name.(get_components(StaticLoad, ieee_rts)); get_name.(get_components(RenewableGen, ieee_rts))]
)
prob = JuMP.Containers.DenseAxisArray(
    [ # spesified for the RTS-96
    # 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
    # 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, # generators
    0.24, 0.51, 0.33, 0.39, 0.48, 0.38, 0.02, 0.36, 0.34, 0.33, 0.30, 0.44, 0.44, 0.44, 
    0.02, 0.02, 0.02, 0.02, 0.40, 0.39, 0.40, 0.52, 0.49, 0.47, 0.38, 0.33, 0.41, 0.41, 
    0.41, 0.35, 0.34, 0.32, 0.54, 0.35, 0.35, 0.38, 0.38, 0.34, 0.34, 0.45, 0.46 # branches
    ],
    get_name.(get_components(ACBranch, ieee_rts))
)
prob /= 8760
contingencies = get_name.(get_branches(ieee_rts))


opfm = scopf(SC, ieee_rts, Gurobi.Optimizer, voll=voll)
solve_model!(opfm.mod)
run_benders(opfm, voll, contingencies, prob)
print_active_power(opfm)