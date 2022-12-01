# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
using PowerSystems
import JuMP
import Gurobi # LP, SOCP, Integer
import Test
include("utils.jl")
include("N-1_SCOPF.jl")
include("short_long_SCOPF.jl")
include("benders.jl")
include("post_process_opf.jl")

function run_benders(fname::String)
    system = System(joinpath("data",fname))
    voll = JuMP.Containers.DenseAxisArray(
        [rand(1000:3000, length(get_components(StaticLoad, system))); rand(1:30, length(get_components(RenewableGen, system)))], 
        [get_name.(get_components(StaticLoad, system)); get_name.(get_components(RenewableGen, system))]
    )

    prob = JuMP.Containers.DenseAxisArray(rand(0.1:0.4, length(get_components(ACBranch, system))),
        get_name.(get_components(ACBranch, system))
    )
    prob /= 8760
    contingencies = get_name.(get_branches(system))
    # contingencies = ["2-3-i_3"]

    opfm, contanal = run_benders(system, voll, contingencies, prob)
    print_active_power(opfm)
    print_power_flow(opfm)
    return opfm, contanal
end