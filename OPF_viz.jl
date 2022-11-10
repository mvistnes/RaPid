# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022
using PowerSystems
using JuMP
# using Ipopt # LP, SOCP, Nonconvex
using Gurobi # LP, SOCP, Integer
# using GLPK # LP, Integer
using Plots
using StatsPlots
using Printf
# using PowerSystemCaseBuilder
using Test
include("utils.jl")
include("N-1_SCOPF.jl")
include("short_long_SCOPF.jl")
include("benders.jl")


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
opfm = scopf(SC, ieee_rts, Gurobi.Optimizer)
opfm = scopf(PCSC, ieee_rts, Gurobi.Optimizer, voll=voll, prob=prob)
sl_model = sl_scopf(ieee_rts, Gurobi.Optimizer, voll=voll, prob=prob)
scatterplot(model, ieee_rts, :pg0, Generator)
print_variabel(model, ieee_rts, :pg0, Generator)

get_active_power.(get_components(StaticLoad, ieee_rts)) - value.(pc_model[:lsd0]).data
[value.(model[:pgu]).data[1,:] , value.(model[:pgd]).data[1,:]]

for d in get_components(StaticLoad, threebus)
    set_active_power!(d, get_active_power(d)*2)
end

function print_cb()
    print("\t\t")
    for l in get_name.(get_components(Branch, ieee_rts))
        print(l[1:5], " ")
    end
    for c in get_name.(get_components(Branch, ieee_rts))[1:10]
        print("\n", c, "\t")
        for l in get_name.(get_components(Branch, ieee_rts))
            if c == l
                print("x     ")
            else
                print(value(model[:cbc][l,c]) > 0.5 ? 1 : 0, "     ")
            end
        end
    end
end

plt = scatter(
    [get_name.(get_components(Branch, ieee_rts))], 
    [value.(model[:pf0]).data ./ get_rate.(get_components(Branch, ieee_rts))], 
    dpi=100, 
    size=(600,600), 
    label = false, 
    rotation=90, 
    title = "pf0/rate"
)
path = mkpath(joinpath("results","N-1_SCOPF_Gurobi"))
savefig(plt, joinpath(path,"pf0_rate.pdf"))

objective_value(model)
objective_value(sl_model[2])

# for g in gens
#     @printf("%s = %.4f <= %.2f \n", get_name(g), value(pg0[get_name(g)]), get_active_power_limits(g).max)
# end
# for l in get_name.(demands)
#     value(ls0[l]) > 0.00001 ? @printf("%s = %.4f \n", l,value(ls0[l])) : print()
# end
# for l in branches
#     @printf("%s = %.4f <= %.2f \n", get_name(l), value(pf0[get_name(l)]), get_rate(l))
# end

solve_time(model)

