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
using Statistics
# include("N-1_SCOPF.jl")
# include("short_long_SCOPF.jl")

scatterplot(model, system, name, type) = 
scatter(
    [get_name.(get_components(type, system))], 
    [value.(model[name]).data], 
    dpi=100, 
    size=(600,600), 
    label = false, 
    rotation=90, 
    title = name
)

function make_save_plot(model, system, sys_name, name)
plt = scatterplot(model, system, name[1], name[2])
display(plt)
path = mkpath(joinpath("results",sys_name))
savefig(plt, joinpath(path,"$(name[1]).pdf"))
end

function scatter_all(model, system; sys_name = "")
names = [
    (:pg0, Generator), (:pgu, Generator), (:pgd, Generator), (:u, ThermalGen),
    (:pf0, Branch), (:pfc, Branch), (:pfcc, Branch), 
    (:ls0, StaticLoad), (:lsc, StaticLoad), (:lscc, StaticLoad), 
    (:qg0, Generator), (:qgu, Generator), (:qgd, Generator), 
    (:qf0, Branch), (:qfc, Branch), (:qfcc, Branch), 
    (:va0, Bus), (:vac, Bus), (:vacc, Bus), 
    (:cbc, Bus), (:cbcc, Bus)
]
for name in names
    try
        plt = scatterplot(model, system, name[1], name[2])
        display(plt)
        path = mkpath(joinpath("results",sys_name))
        savefig(plt, joinpath(path,"$(name[1]).pdf"))
    catch e
        print("No $(name[1]). ")
    end
end
end

function print_variabel(model, system, name, type)
    for i in sort(get_name.(get_components(type, system)))
        @printf("%s: %.3f\n", i, value(model[name][i]))
    end
end

function print_results(opfm::OPFmodel)
    print_variabel(opfm.mod, opfm.sys, :pg0, Generator)
    print_variabel(opfm.mod, opfm.sys, :ls0, StaticLoad)
    print_variabel(opfm.mod, opfm.sys, :va0, Bus)
end

function print_active_power(opfm::OPFmodel)
    nodes = get_sorted_nodes(opfm.sys)
    list_gen = make_list(opfm, get_ctrl_generation, nodes)
    list_d = make_list(opfm, get_demands, nodes)
    list_r = make_list(opfm, get_renewables, nodes)
    xprint(x, val) = @printf("%7s:%6.2f, ", split(x.name, "_", limit=2)[end], val)
    sort_x!(list) = sort!(list, by = x -> parse(Int64, split(x.name, "_")[end]))
    sort_alt!(list) = sort!(list, by = x -> x.name)
    get_injection(g::Union{ThermalGen, HydroGen}) = value(opfm.mod[:pg0][get_name(g)])
    get_injection(r::RenewableGen) = get_active_power(r) - value(opfm.mod[:ls0][get_name(r)])
    get_injection(d::StaticLoad) =  -get_active_power(d) + value(opfm.mod[:ls0][get_name(d)])
    print(" Bus   Total  Injections...")
    tot = 0.0
    for (n, g, d, r) in zip(nodes, list_gen, list_d, list_r)
        try 
            sort_x!(g)
            sort_x!(d)
            sort_x!(r)
        catch
            sort_alt!(g)
            sort_alt!(d)
            sort_alt!(r)
        end
        val_g = get_injection.(g)
        val_d = get_injection.(d)
        val_r = get_injection.(r)
        b_tot = sum(val_g, init=0.0) + sum(val_d, init=0.0) + sum(val_r, init=0.0)
        tot += b_tot
        @printf "\n%4s  %6.2f  " n.name b_tot
        xprint.(g, val_g)
        xprint.(d, val_d)
        xprint.(r, val_r)
    end
    @printf("\n Sum  %6.2f\n", tot)
end

function print_power_flow(opfm::OPFmodel)
    branches = get_sorted_branches(opfm.sys)
    println("    Branch    Flow  Rating")
    for branch in branches
        @printf("%10s  %6.2f  %6.2f\n", 
                get_name(branch), 
                value(opfm.mod[:pf0][get_name(branch)]), 
                get_rate(branch)
            )
    end
end

# # corrective control failure probability
# phi(p, n) = sum((-1)^k * p^k * binomial(n,k) for k in 1:n)

# # severity function
# @expression(opf_m, severity, sum(voll[d] * lse[d] for d in get_name.(demands)))

# add_system_data_to_json()
# system_data = System("system_data.json")
# results = opf_model(system_data, Ipopt.Optimizer)
# value.(results[:pl])
