# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
include("SCOPF.jl")
import .SCOPF

function setup(system::System)
    voll = JuMP.Containers.DenseAxisArray(
        [rand(1000:3000, length(get_components(StaticLoad, system))); rand(1:30, length(get_components(RenewableGen, system)))], 
        [get_name.(get_components(StaticLoad, system)); get_name.(get_components(RenewableGen, system))]
    )

    prob = JuMP.Containers.DenseAxisArray(rand(length(get_components(ACBranch, system))) .* 0.3 .+ 0.1,
        get_name.(get_components(ACBranch, system))
    )
    prob /= 8760
    contingencies = get_name.(SCOPF.get_branches(system))
    # contingencies = ["2-3-i_3"]
    return voll, prob, contingencies
end

function test()
    system = SCOPF.get_system("ELK14.json")
    voll, prob, contingencies = setup(system)
    opfm, imml = SCOPF.run_benders2(system, voll, prob)
    # print_active_power(opfm)
    # print_power_flow(opfm)
end