using PowerSystems
const PSY = PowerSystems
using JuMP
using Ipopt # LP, SOCP, Nonconvex
using Gurobi # LP, SOCP, Integer
using GLPK # LP, Integer
using Plots
using StatsPlots
using Printf
using PowerSystemCaseBuilder
using Test

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

function scatter_all(model, system; sys_name = "")
    names = [
        (:pg0, Generator), (:pgu, Generator), (:pgd, Generator), (:pf0, Branch), 
        (:pfc, Branch), (:pfcc, Branch), (:ls0, StaticLoad), (:lsc, StaticLoad), (:qg0, Generator), 
        (:qgu, Generator), (:qgd, Generator), (:qf0, Branch), (:qfc, Branch), (:qfcc, Branch), 
        (:va0, Bus), (:vac, Bus), (:vacc, Bus)
    ]
    for name in names
        try
            plt = scatterplot(model, system, name[1], name[2])
            display(plt)
            savefig(plt, @sprintf("%s%s.pdf", sys_name, name[1]))
        catch e
            @printf("No %s. ", name[1])
        end
    end
end

function print_variabel(model, system, name, type)
    for (i,x) in zip(get_name.(get_components(type, system)), value.(model[name]).data)
        @printf("%s: %.3f\n", i, x)
    end
end

ieee_rts = System("ieee_rts.json")
scopf_m = scopf(ieee_rts, Gurobi.Optimizer, corrective=true)
scopf_m = sl_scopf(ieee_rts, Gurobi.Optimizer)
scatterplot(scopf_m, ieee_rts, :pg0, Generator)
print_variabel(scopf_m, ieee_rts, :pg0, Generator)

get_active_power.(get_components(StaticLoad, ieee_rts)) - value.(pc_scopf_m[:lsd0]).data
[value.(scopf_m[:pgu]).data[1,:] , value.(scopf_m[:pgd]).data[1,:]]

for d in get_components(StaticLoad, ieee_rts)
    set_active_power!(d, get_active_power(d)*1.5) 
end

plt = scatter(
    [get_name.(get_components(Branch, ieee_rts))], 
    [value.(scopf_m[:pf0]).data ./ get_rate.(get_components(Branch, ieee_rts))], 
    dpi=100, 
    size=(600,600), 
    label = false, 
    rotation=90, 
    title = "pf0/rate"
)
savefig(plt, "pf0_rate_Ipopt.pdf")

objective_value(p_opf_m) + objective_value(pc_opf_m)
objective_value(scopf_m)


A = JuMP.Containers.DenseAxisArray((rand(length(gens_h)).*(0.5-0.02).+0.02)./8760, get_name.(gens_h))