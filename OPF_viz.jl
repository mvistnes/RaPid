using PowerSystems
const PSY = PowerSystems
using JuMP
using Ipopt # LP, SOCP, Nonconvex
using Gurobi # LP, SOCP, Integer
using Plots
using StatsPlots
using Printf
using PowerSystemCaseBuilder

scatterplot(model, name, type) = 
    scatter(
        [get_name.(get_components(type, ieee_rts))], 
        [value.(model[name]).data], 
        dpi=100, 
        size=(600,600), 
        label = false, 
        rotation=90, 
        title = name
    )

function scatter_all(model)
    names = [
        (:pg0, Generator), (:pgu, Generator), (:pgd, Generator), (:pf0, Branch), 
        (:pfc, Branch), (:pfcc, Branch), (:ls0, StaticLoad), (:lsc, StaticLoad), (:qg0, Generator), 
        (:qgu, Generator), (:qgd, Generator), (:qf0, Branch), (:qfc, Branch), (:qfcc, Branch), 
        (:va0, Bus), (:vac, Bus), (:vacc, Bus)
    ]
    for name in names
        try
            plt = scatterplot(model, name[1], name[2])
            display(plt)
        catch e
            @printf("No %s. ", name[1])
        end
    end
end

ieee_rts = System("ieee_rts.json")
scopf_m = scopf(ieee_rts, Ipopt.Optimizer)

scatterplot(scopf_m, :pg0, Generator)
scatterplot(scopf_m, :pgu, Generator)
scatterplot(scopf_m, :pgd, Generator)
scatterplot(scopf_m, :pf0, Branch)
scatterplot(scopf_m, :pfc, Branch)
scatterplot(scopf_m, :pfcc, Branch)
scatterplot(scopf_m, :ls0, StaticLoad)
scatterplot(scopf_m, :lsc, StaticLoad)

scatterplot(scopf_m, :qg0, Generator)
scatterplot(scopf_m, :qgu, Generator)
scatterplot(scopf_m, :qgd, Generator)
scatterplot(scopf_m, :qf0, Branch)
scatterplot(scopf_m, :qfc, Branch)
scatterplot(scopf_m, :qfcc, Branch)
scatterplot(scopf_m, :va0, Bus)
scatterplot(scopf_m, :vac, Bus)
scatterplot(scopf_m, :vacc, Bus)

get_active_power.(get_components(StaticLoad, ieee_rts)) - value.(pc_scopf_m[:lsd0]).data

scatter(
    [get_name.(get_components(Branch, ieee_rts))], 
    [value.(scopf_m[:pfcc]).data ./ get_rate.(get_components(Branch, ieee_rts))], 
    dpi=100, 
    size=(600,600), 
    label = false, 
    rotation=90, 
    title = "pfcc/rate"
)