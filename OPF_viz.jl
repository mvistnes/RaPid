using PowerSystems
const PSY = PowerSystems
using JuMP
using Ipopt
using Plots
using StatsPlots
using Printf
using PowerSystemCaseBuilder

scatterplot(model, name, type) = scatter([get_name.(get_components(type, ieee_rts))], [value.(model[name]).data], dpi=100, size=(600,600), label = false, rotation=90)

ieee_rts = System("ieee_rts.json")
p_scopf_m = p_scopf(ieee_rts, Ipopt.Optimizer)
scatterplot(p_scopf_m, :pg0, ThermalStandard)
scatterplot(p_scopf_m, :qg0, ThermalStandard)
scatterplot(p_scopf_m, :pf0, Line)
scatterplot(p_scopf_m, :pfc, Line)
scatterplot(p_scopf_m, :vm0, Bus)
scatterplot(p_scopf_m, :vmc, Bus)
scatterplot(p_scopf_m, :va0, Bus)
scatterplot(p_scopf_m, :vac, Bus)
scatterplot(p_scopf_m, :ls0, StaticLoad)

pc_scopf_m = pc_scopf(ieee_rts, Ipopt.Optimizer)
scatterplot(pc_scopf_m, :pg0, ThermalStandard)
scatterplot(pc_scopf_m, :pgu, ThermalStandard)
scatterplot(pc_scopf_m, :pgd, ThermalStandard)
scatterplot(pc_scopf_m, :pf0, Line)
scatterplot(pc_scopf_m, :pfc, Line)
scatterplot(pc_scopf_m, :pfcc, Line)
scatterplot(pc_scopf_m, :angle0, Bus)
scatterplot(pc_scopf_m, :anglec, Bus)
scatterplot(pc_scopf_m, :anglecc, Bus)
scatterplot(pc_scopf_m, :ls0, StaticLoad)
scatterplot(pc_scopf_m, :lsc, StaticLoad)

get_active_power.(get_components(StaticLoad, ieee_rts)) - value.(pc_scopf_m[:lsd0]).data