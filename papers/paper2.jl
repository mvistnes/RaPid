import Statistics
import Plots
import Measures
import CSV
using LaTeXStrings

include("../test/runtests.jl");

system = include("../data/RTS_GMLC/config.jl")
n = "g500"
fname = joinpath("data", "matpower", "ACTIVSg500.m")
system, voll, contingencies, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec = setup_system(joinpath("data", "matpower", "ACTIVSg500.m"))


Plots.scatter(df[!,:max_shed]./10, df[!,:pg0]./10, leg=:none, xlabel=L"\Gamma ~ [MWh]", ylabel="Base cost [k\$/h]", size=(500,300))
Plots.savefig("../results/500bus_pg0_vs_max_shed.pdf")

a = Plots.scatter(df[!,:pg0]./10, df[!,:p], ylabel=L"\pi ~ [-]", leg=:none)
b = Plots.scatter(df[!,:pg0]./10, df[!,:VOLL], ylabel=L"VOLL ~ [-]", leg=:none)
c = Plots.scatter(df[!,:pg0]./10, df[!,:VOLL].*df[!,:p], ylabel=L"VOLL \cdot \pi ~ [-]", xlabel="Base costs [k\$/h]", leg=:none)
Plots.plot(a,b,c,layout=Plots.grid(3,1), size=(500,900), bottom_margin=7Measures.mm, left_margin=7Measures.mm)
Plots.savefig("../results/500bus_pg0_vs_p_VOLL.pdf")

df_names = [L"VOLL" L"\pi" L"\Gamma" L"u_k" L"b_{k1}" L"b_{k2}" L"EOC" L"\sum^{\mathcal{G}}_g c_g p_g" L"\sum^\mathcal{K} \sum^\mathcal{D} u_{dk1}" L"\sum^\mathcal{K} \sum^\mathcal{D} u_{dk2}" L"\sum^\mathcal{K} \pi_k \sum^\mathcal{D} u_{dk1}" L"\sum^\mathcal{K} \pi_k \sum^\mathcal{D} u_{dk2}"]
filter(:obj_val => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df)
itr = [1,2,3,4,5,6,8,9,10,7];
Plots.heatmap(Statistics.cor(Matrix(df)[:,itr]), xticks=(1:length(df_names[itr]), df_names[itr]), yticks=(1:length(df_names[itr]), df_names[itr]), xrotation=90, right_margin=7Measures.mm, bottom_margin=7Measures.mm, aspect_ratio=1, size=(500,500))
df = DataFrames.DataFrame(CSV.File("../results/500_bus_random.csv"))
Plots.heatmap(Statistics.cor(Matrix(df)[:,itr]), xticks=(1:length(df_names[itr]), df_names[itr]), yticks=(1:length(df_names[itr]), df_names[itr]), xrotation=90, right_margin=7Measures.mm, bottom_margin=7Measures.mm, aspect_ratio=1, size=(500,500))
df = DataFrames.DataFrame(CSV.File("../results/500bus_random_no_split.csv")) 
Plots.heatmap(Statistics.cor(Matrix(df)[:,itr]), xticks=(1:length(df_names[itr]), df_names[itr]), yticks=(1:length(df_names[itr]), df_names[itr]), xrotation=90, right_margin=7Measures.mm, bottom_margin=7Measures.mm, aspect_ratio=1, size=(500,500))

Plots.scatter(df[!,:max_shed]./10, df[!,:pg0]./10, leg=:none, xlabel=L"\Gamma ~ [MWh]", ylabel="Base cost [k\$/h]", size=(500,300))


df_names = [L"VOLL" L"\pi" L"\Gamma" L"c_k" L"\bar h_b^{S,k}" L"\bar h_b^{L,k}" L"EOC" L"\sum^{\mathcal{G}}_g c_g p_g" L"\sum^\mathcal{K} \sum^\mathcal{D} u_{d}^{S,k}" L"\sum^\mathcal{K} \sum^\mathcal{D} u_{d}^{L,k}"]
Plots.heatmap(Statistics.cor(Matrix(df)[:,itr]), xticks=(1:length(df_names[itr]), df_names[itr]), yticks=(1:length(df_names[itr]), df_names[itr]), xrotation=90, right_margin=7Measures.mm, bottom_margin=7Measures.mm, aspect_ratio=1, size=(500,500))
Plots.savefig("../results/500bus_heatmap.pdf")

df = DataFrames.DataFrame(CSV.File("../results/ACTIVSg2000_bus_random.csv")) 
df = filter(:obj_val => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df)
Plots.heatmap(Statistics.cor(Matrix(df)), xticks=(1:length(names(df)), names(df)), yticks=(1:length(names(df)), names(df)), xrotation=90, right_margin=7Measures.mm, bottom_margin=7Measures.mm, aspect_ratio=1, size=(500,500))

scatter_plot(df, :long, :pg0)
Plots.scatter(df[!,:VOLL].*100, df[!,:pg0]./10, label="C-SCOPF", xlabel="Branch overload limit [%]", ylabel="Base cost [k\$/h]")
Plots.scatter(df[!,:VOLL].*df[!,:p], df[!,:pg0]./10, xlabel="VOLL*p", ylabel="Base cost [k\$/h]")
Plots.scatter(df[!,:long].*df[!,:short], df[!,:pg0]./10, xlabel="long*short", ylabel="Base cost [k\$/h]")
Plots.scatter(df[!,:long].*df[!,:short], df[!,:obj_val]./10, xlabel="long*short", ylabel="Obj val [k\$/h]")
Plots.scatter(df[!,:long].*df[!,:p], df[!,:obj_val]./10, xlabel="long*p", ylabel="Obj val [k\$/h]")