import Dates, DataFrames, CSV, Plots, StatsPlots, Measures, Gurobi, Printf, FileIO, JLD2, Random, Logging, JuMP
using LaTeXStrings
import RaPidSCOPF as SCOPF
const GRB_ENV = Gurobi.Env()
Random.seed!(53715)

## Run this line in REPL
system = include(joinpath(pwd(), "data", "RTS_GMLC", "config.jl"))
include(joinpath(pwd(), "src", "scenario.jl"))
include(joinpath(pwd(), "src", "scenario_utils.jl"))

df = DataFrames.DataFrame(CSV.File("data/RTS_GMLC/branch.csv"))
df_n = DataFrames.DataFrame(CSV.File("data/RTS_GMLC/bus.csv"))
branches = SCOPF.sort_components!(SCOPF.get_branches(system));
ia_br = SCOPF.get_interarea(branches)
brc = reduce(vcat, [ia_br, branches[findall(x->x.name ∈ ["A30", "A21", "B30", "B21", "C30", "C21"], branches)]])

max_lng = maximum(df_n[!,:lng])
delay = Dict{String, Float64}()
for b in brc
    fbus = df[df.UID .== b.name, "From Bus"][1]
    tbus = df[df.UID .== b.name, "To Bus"][1]
    delay[b.name] = round.((max_lng - df_n[df_n."Bus ID" .== fbus,:lng][1] + max_lng - df_n[df_n."Bus ID" .== tbus,:lng][1]) / 2)
end

prob = reduce(vcat, [[df[!,"Perm OutRate"][i] for i in 1:length(branches) if x == branches[i].name] for x in df.UID])
w = [branches[i] ∈ brc ? Storm(x/100/8760, 0.02*(0.5+rand()), rand(20:25), rand(1:3), rand(10:14), rand(10:14)) : storm() for (i,x) in enumerate(prob)]
# w = [branches[i] ∈ brc ? Storm(x/100/8760, 0.02*(0.5+rand()), 20+delay[branches[i].name], rand(1:3), rand(10:14), rand(10:14)) : storm() for (i,x) in enumerate(prob)]
w2 = [branches[i] ∈ brc ? Storm(x/100/8760, 0.02*(0.5+rand()), rand(20:25), rand(1:3), rand(10:14), rand(10:14)) : storm() for (i,x) in enumerate(prob)]
# w2 = [branches[i] ∈ brc ? Storm(x/100/8760, 0.02*(0.5+rand()), rand(25:30), rand(1:3), rand(10:14), rand(10:14)) : storm() for (i,x) in enumerate(prob)]
p = generate_p_from_weather(prob/8760, 60, w)
p2 = generate_p_from_weather(prob/8760, 60, w2)
scenarioes = generate_contingency_probabilities(p, 5, 10000)
scenarioes2 = generate_contingency_probabilities(p2, 5, 10000)

style = Plots.supported_styles()[2:end]
style = reshape(style, 1, length(style))
bnames = first.(sort(get_name.(branches) .=> w, by=x->x[2].max)[1:11])
Plots.plot([[last(sort(x[2])) for x in scenarioes] [last(sort(x[2])) for x in scenarioes2]], labels=[])
Plots.plot(reduce(hcat, generate_weather.([60], sort(w, by=x->x.max)[1:11])).*100, label=permutedims(bnames), 
    xlabel="Time [h]", ylabel="Outage probability [%]", palette = Plots.palette(:seaborn_colorblind6), size=(500,250))
Plots.boxplot(reduce(hcat, generate_weather.([60], sort(w2, by=x->x.max)[1:11]))', 
    xlabel="Time", ylabel="Probability", palette = Plots.palette([:lightgrey],1), leg=:none)

vals = [collect(values(sort(SCOPF.sumvals(length.(last.(x[1])), x[2]), by=i->i[1]))) for x in scenarioes]
c_size = maximum(length.(vals))
outs = zeros(60,c_size);
for (i,x) in enumerate(vals)
    for (j,y) in enumerate(x)
        outs[i,j] = y
    end
end
StatsPlots.groupedbar([sum(outs[:,1:3], dims=2) outs[:,4:end]], bar_position = :stack, xlabel="Time [h]", ylabel="Probability", 
    palette = Plots.palette(:seaborn_colorblind6, rev=true), labels=["N-k > 3" "N-3" "N-2" "N-1" "N-0"], leg=:bottomright, grid=:none)
Plots.plot(cumsum([outs[:,1:4] sum(outs[:,5:end], dims=2)], dims = 2)[:,end:-1:1].*100, fill=0, lc=:black, legend=:none, grid=:none,
    xlabel="Time [h]", ylabel="Probability [%]", palette = Plots.palette(:seaborn_colorblind6, rev=true), size=(500,300),
    annotation=[(27, 20, ("N-0", 8)), (27, 60, ("N-1", 8)), (27, 85, ("N-2", 8)), (27, 97, ("N-3", 8)), (27, 103, ("N-k > 3", 8))])
Plots.plot([collect(values(sort(SCOPF.countmemb(length.(last.(x[1]))), by=i->i[1]))) for x in scenarioes])
Plots.plot([collect(values(sort(SCOPF.sumvals(length.(last.(x[1])), x[2]), by=i->i[1]))) for x in scenarioes])

scenarioes_n1 = scenario_plus_n1(scenarioes, length(branches))
scenarioes2_n1 = scenario_plus_n1(scenarioes2, length(branches))

contingencies = ["branch" => [i] for i in 1:length(SCOPF.get_branches(system))]
demands = SCOPF.sort_components!(SCOPF.get_demands(system))
base_voll = [6.20 4.89 5.30 5.62 6.11 5.50 5.41 5.40 2.30 4.14 5.39 3.41 3.01 3.54 3.75 2.29 3.64]
voll = vec([base_voll... base_voll... base_voll...] * 100)
SCOPF.set_ramp_limits!(system)
demands = SCOPF.sort_components!((SCOPF.get_demands(system)))
renewables = SCOPF.sort_components!((SCOPF.get_renewables(system)))
gen = SCOPF.sort_components!(SCOPF.get_generation(system))
SCOPF.set_operation_cost!.(SCOPF.get_gens_h(system), rand(20) .* (55 - 45) .+ 45)
# SCOPF.set_rate!.(SCOPF.get_branches(system), SCOPF.get_rate.(SCOPF.get_branches(system)) * 0.8);
short = 1.2
long = 1.0
# n = "long"
prev_lim = 1.0
ramp_minutes = 10.0
max_shed = 1.5
reserve = 0.0
ramp_mult = 2.0
demand_mult = 1.0
renew_mult = 1.0
renew_cost = 0.00
renew_ramp = 0
i = 9
t = Dates.DateTime("2020-04-11T13:00:00")
# t = Dates.DateTime("2020-08-25T14:00:00")
time_limit_sec = length(contingencies)^2 + 10
labels = ["NkC" "NkCz" "NkCz2" "NkP" "N1C" "N1P" "N0"]
# xguidefont=fsize, yguidefont=fsize, xtickfont=fsize, ytickfont=fsize, legendfont=fsize
name = Printf.@sprintf "prev%d_ls%d_pd%d_pr%d_ramp%d_r%d_%s" prev_lim*100 max_shed*100 demand_mult*100 renew_mult*100 renew_ramp reserve*100 string(t)[1:end-6]
println("\n\n"*name*"\n")
results = run_plotsave_weather(name, system, labels, scenarioes_n1, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim, demands, demand_mult, renewables, renew_mult, t);
results_w2 = run_plotsave_weather2(results, "s2"*name, system, labels[1:6], scenarioes2_n1, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim);
results2 = run_plotsave_weather("w2_"*name, system, labels, scenarioes2_n1, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim, demands, demand_mult, renewables, renew_mult, t);

sx = repeat(["W1", "W2"], inner = length(results))
df_results = table_results(results, system, scenarioes_n1)
costs = distributed_costs(results, case.opf, scenarioes_n1)
StatsPlots.groupedbar!(Matrix(costs)[:,[9,8,7,6,4,3,2,1]]./10,
        bar_position = :stack, bar_width=0.5, xticks=(1:length(results), getindex.(results, 1)), leg=:outertopright,
        label=reverse([L"G^0" L"\ell^0" L"G^{r+}" L"G^{r-}" L"\ell^{S,\mathcal{K}}" L"G^{L+,\mathcal{K}}" L"G^{L-,\mathcal{K}}" L"\ell^{L-,\mathcal{K}}"]),
        # label=[L"p_\mathcal{G}^0" L"\ell_\mathcal{D}^0" L"r_\mathcal{G}^+" L"r_\mathcal{G}^-" L"\ell_\mathcal{D}^{S,\mathcal{K}}" L"p_\mathcal{G}^{L+,\mathcal{K}}" L"p_\mathcal{G}^{L-,\mathcal{K}}" L"\ell_\mathcal{D}^0"],
        ylabel="Cost [k\$]", palette = Plots.palette(:seaborn_colorblind, rev=true)[[10,1,2,3,4,5,6,7,8,9]], size=(500,200), xrotation = 45, bottommargin=5*Plots.Measures.mm)
StatsPlots.groupedbar(sum(Matrix(costs2)[:,[9,8,7,6,4,3,2,1]]./10, dims=2), label=:none, bar_width=0.8, 
        xticks=(1:length(results), getindex.(results, 1)), palette = Plots.palette(:Greys_3, rev=true))
forced_ls = [SCOPF.calc_forced_ls(case, x[2], x[1]) for x in scenarioes_n1];
pl, xax, data = plot_probability_contingency_data(results, :Pc, :lsc, labels, forced_ls, scenarioes_n1)

# for max_shed in [1.5, 3.0], demand_mult in [1.0, 1.5], renew_mult in [0.5, 1.0], renew_ramp in [0.0, 10.0]
labels = ["N-k-C" "N-k-zC" "N-k-zVOLL"]
for demand_mult in 0.5:0.25:1.5, renew_mult in 0.0:0.25:1.0
    name = Printf.@sprintf "prev%d_ls%d_pd%d_pr%d_ramp%d_r%d_%s" prev_lim*100 max_shed*100 demand_mult*100 renew_mult*100 renew_ramp reserve*100 string(t)[1:end-6]
    println("\n\n"*name*"\n")
    # results = run_plotsave_weather(name, system, labels, scenarioes_n1, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim, demands, demand_mult, renewables, renew_mult, t)
    results = run_plotsave_weather_zc("zc_"*name, system, labels, scenarioes_n1, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, demands, demand_mult, renewables, renew_mult, t)
end
for max_shed in 0.0:5.0
    name = Printf.@sprintf "prev%d_ls%d_pd%d_pr%d_ramp%d_r%d_%s" prev_lim*100 max_shed*100 demand_mult*100 renew_mult*100 renew_ramp reserve*100 string(t)[1:end-6]
    println("\n\n"*name*"\n")
    results = run_plotsave_weather(name, system, labels, scenarioes_n1, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim, demands, demand_mult, renewables, renew_mult, t)
end
for renew_mult in 0.0:0.25:1.5
    name = Printf.@sprintf "prev%d_ls%d_pd%d_pr%d_ramp%d_r%d_%s" prev_lim*100 max_shed*100 demand_mult*100 renew_mult*100 renew_ramp reserve*100 string(t)[1:end-6]
    println("\n\n"*name*"\n")
    results = run_plotsave_weather(name, system, labels, scenarioes_n1, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim, demands, demand_mult, renewables, renew_mult, t)
end
for prev_lim in 0.6:0.2:1.4
    name = Printf.@sprintf "prev%d_ls%d_pd%d_pr%d_ramp%d_r%d_%s" prev_lim*100 max_shed*100 demand_mult*100 renew_mult*100 renew_ramp reserve*100 string(t)[1:end-6]
    println("\n\n"*name*"\n")
    results = run_plotsave_weather(name, system, labels, scenarioes_n1, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim, demands, demand_mult, renewables, renew_mult, t)
end
for renew_ramp in [0.0, 1.0, 10.0, 100.0, 1000.0]
    name = Printf.@sprintf "prev%d_ls%d_pd%d_pr%d_ramp%d_r%d_%s" prev_lim*100 max_shed*100 demand_mult*100 renew_mult*100 renew_ramp reserve*100 string(t)[1:end-6]
    println("\n\n"*name*"\n")
    results = run_plotsave_weather(name, system, labels, scenarioes_n1, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim, demands, demand_mult, renewables, renew_mult, t)
end
for reserve in 0.0:0.01:0.1
    name = Printf.@sprintf "prev%d_ls%d_pd%d_pr%d_ramp%d_r%d_%s" prev_lim*100 max_shed*100 demand_mult*100 renew_mult*100 renew_ramp reserve*100 string(t)[1:end-6]
    println("\n\n"*name*"\n")
    results = run_plotsave_weather(name, system, labels, scenarioes_n1, voll, prob, contingencies, max_shed, reserve, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim, demands, demand_mult, renewables, renew_mult, t)
end

sensitivity = []
itr = 0.5:5.0 # max_shed
itr = 0.25:0.25:1.5 # renew_mult
itr = 0.6:0.2:1.4 # prev_lim
itr = [0.0, 1.0, 10.0, 100.0, 1000.0] # renew_ramp
for renew_ramp in itr
    name = Printf.@sprintf "prev%d_ls%d_pd%d_pr%d_ramp%d_r%d_%s" prev_lim*100 max_shed*100 demand_mult*100 renew_mult*100 renew_ramp reserve*100 string(t)[1:end-6]
    results = FileIO.load("results/"*name*".jld2", "results");
    base_costs, corrective_costs = calc_costs(results, scenarioes_n1)
    push!(sensitivity, prev_lim => sum.(base_costs) .+ sum.(corrective_costs))
    # push!(sensitivity, prev_lim => sum.(reduce(vcat, [base_costs, base_costs2])) .+ sum.(reduce(vcat, [corrective_costs, corrective_costs2])))
    # push!(sensitivity, prev_lim => Dict(:base_costs => base_costs, :corrective_costs => corrective_costs))
end
Plots.plot(reduce(hcat, [reduce(vcat, x[2][:base_costs]) for x in sensitivity])',labels=labels, xticks=(1:length(itr), itr), xlabel="Renewable corrective cost", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6), thickness_scaling=1.2)
Plots.plot(reduce(hcat, [reduce(vcat, x[2][:corrective_costs]) for x in sensitivity])',labels=labels, xticks=(1:length(itr), itr), xlabel="Renewable corrective cost", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6), thickness_scaling=1.2)
Plots.plot(reduce(hcat, [reduce(vcat, x[2][:base_costs] .+ x[2][:corrective_costs]) for x in sensitivity])',labels=labels, xticks=(1:length(itr), Int.(itr.*100)), xlabel="Renewable availability [%]", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6), thickness_scaling=1.2, leg=:topright)
Plots.plot(reduce(hcat, [reduce(vcat, x[2]) for x in sensitivity])'[:,1:5], line=:dash, alpha=0.5, labels=labels, xticks=(1:length(itr), Int.(itr.*100)), xlabel="Renewable availability [%]", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6)[1:5], thickness_scaling=1.2, leg=:topright)
Plots.plot!(reduce(hcat, [reduce(vcat, x[2]) for x in sensitivity])'[:,6:end],labels=labels .* " 2S", xticks=(1:length(itr), Int.(itr.*100)), xlabel="Load shed limit [MW]", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6)[1:5], thickness_scaling=1.2, leg=:topright)

# Linear sensitivity plot
function load_costs(prev_lim, max_shed, demand_mult, renew_mult, renew_ramp, reserve, t, scenarioes_n1)
    name = Printf.@sprintf "prev%d_ls%d_pd%d_pr%d_ramp%d_r%d_%s" prev_lim*100 max_shed*100 demand_mult*100 renew_mult*100 renew_ramp reserve*100 string(t)[1:end-6]
    results = FileIO.load("results/"*name*".jld2", "results");
    base_costs, corrective_costs = calc_costs(results, scenarioes_n1)
    return getindex.(last.(base_costs),26) .+ getindex.(last.(corrective_costs),26)
    # return (sum.(last.(base_costs)) .+ sum.(last.(corrective_costs))) ./1000
end

rm_vals = []; pl_vals = []; ms_vals = []; rr_vals = []
rm_itr = 0.25:0.25:1.5; pl_itr = 0.6:0.2:1.4; ms_itr = 0.0:5.0; rr_itr = [0.0, 1.0, 10.0, 100.0]
for renew_mult in rm_itr
    push!(rm_vals, load_costs(prev_lim, max_shed, demand_mult, renew_mult, renew_ramp, reserve, t, scenarioes_n1))
end
rm = reduce(hcat, rm_vals)'
a = Plots.plot(rm, labels=:none, xlabel="Availability [%]", ylabel="Cost [\$]", xticks=(1:length(rm_itr), Int.(rm_itr.*100)))

# for prev_lim in pl_itr
#     push!(pl_vals, load_costs(prev_lim, max_shed, demand_mult, renew_mult, renew_ramp, reserve, t, scenarioes_n1))
# end
# pl = reduce(hcat, pl_vals)'
# b = Plots.plot(pl, labels=:none, xlabel="Preventive limit [%]", xticks=(1:length(pl_itr), Int.(pl_itr.*100)))

# for max_shed in ms_itr
#     push!(ms_vals, load_costs(prev_lim, max_shed, demand_mult, renew_mult, renew_ramp, reserve, t, scenarioes_n1))
# end
# ms = reduce(hcat, ms_vals)'
# c = Plots.plot(ms, labels=:none, xlabel="Load shed limit [MW]", xticks=(1:length(ms_itr), Int.(ms_itr.*100)))

for renew_ramp in rr_itr
    push!(rr_vals, load_costs(prev_lim, max_shed, demand_mult, renew_mult, renew_ramp, reserve, t, scenarioes_n1))
end
rr = reduce(hcat, rr_vals)'
d = Plots.plot(rr, labels=labels, leg=:outertopright, xlabel="Corr. cost [\$/MWh]", xticks=(1:length(rr_itr), Int.(rr_itr)))
# Plots.plot(a,b,c,d, layout=Plots.grid(1,4, widths=(5/17,4/17,4/17,4/17)), size=(1200,300), yaxis=[0,500], 
Plots.plot(a,d, layout=Plots.grid(1,2, widths=(4/8,4/8)), size=(500,300), yaxis=[0,500], 
    bottom_margin=7Measures.mm, left_margin=7Measures.mm, palette = Plots.palette(:seaborn_colorblind6), thickness_scaling=1.0)

ty = [["short", "long"], [0.5, 5.0], [1.5, 1.0], [0.0, 1000.0]]
choose(ty, i, c) = c == i ? ty[i][2] : ty[i][1]
for c in 1:length(ty)
    name = Printf.@sprintf "prev-%s_ls%d_pd150_pr%dramp%d_time%d" choose(ty,1,c) choose(ty,2,c)*100 choose(ty,3,c)*100 choose(ty,4,c) i
    results = FileIO.load("results/"*name*".jld2", "results");
    base_costs = [[isnan(x[:obj_val]) ? NaN : x[:costs][1,:Base] for x in r] for (_,r) in results]
    corrective_costs = [[isnan(x[:obj_val]) ? NaN : (x[:costs][1,:r] + sum(p[2] .* (x[:costs][:,:Pc] + x[:costs][:,:Pcc]))) for (x,p) in zip(r,scenarioes_n1)] for (_,r) in results]
    Plots.plot(results_base_costs + results_corrective_costs, labels=labels, line=(:steppre, :dash), alpha=0.5, xlabel="Time", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6)[1:5])
    Plots.plot!(base_costs + corrective_costs, labels=labels.*" 2", line=:steppre, xlabel="Time", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6)[1:5], leg=:topright)
    Plots.savefig("results/"*name*"_objective_2.pdf")
end

# HEATMAP
vals = DataFrames.DataFrame(prev_lim=Float64[], max_shed=Float64[], demand_mult=Float64[], renew_mult=Float64[], renew_ramp=Float64[], i=Int64[], 
    Nk_C=Float64[], Nk_P=Float64[], N1_C=Float64[], N1_P=Float64[], N0=Float64[])
for max_shed in [1.5, 3.0], demand_mult in [1.0, 1.5], renew_mult in [0.5, 1.0], renew_ramp in [0.0, 10.0]
    name = Printf.@sprintf "prev%d_ls%d_pd%d_pr%d_ramp%d_r%d_%s" prev_lim*100 max_shed*100 demand_mult*100 renew_mult*100 renew_ramp reserve*100 string(t)[1:end-6]
    valb = CSV.read("results/"*name*"_base_cost.csv", DataFrames.DataFrame)
    valc = CSV.read("results/"*name*"_corrective_cost.csv", DataFrames.DataFrame)
    # push!(vals, (prev_lim, max_shed, demand_mult, renew_mult, renew_ramp, i, valb[1,1:5]..., valc[1,1:5]...), promote=true)
    push!(vals, (prev_lim, max_shed, demand_mult, renew_mult, renew_ramp, i, sum(Matrix(valb .+ valc), dims=1)[1:5]...), promote=true)
end
function make_heatmap(vals)
    short_labels = ["NkC", "NkP", "N1C", "N1P", "N0"]
    short_types = [L"1.0\Gamma, ~ c_g=0", L"1.0\Gamma, ~ c_g=10", L"2.0\Gamma, ~ c_g=0", L"2.0\Gamma, ~ c_g=10"]
    data = Matrix(vals[!,r"N"]) ./ 1000
    pl = Plots.heatmap(data', yticks=(1:5,short_labels), xticks=(1:size(data,1),short_types), xrotation = 45, c=Plots.cgrad(:island, rev=true))
    Plots.annotate!([(j,i,Plots.text(Int(round(data[j,i], digits=0)), 8, :black, :center)) for i in 1:size(data,2) for j in 1:size(data,1)], linecolor=:white)
    return pl
end
dm1r05 = vals[vals.demand_mult .== 1.0 .&& vals.renew_mult .== 0.5 .&& vals.prev_lim .== 1.0,:]
a = make_heatmap(dm1r05)
Plots.plot!(a, ylabel="Renew 50%, 2348 MW")
dm15r05 = vals[vals.demand_mult .== 1.5 .&& vals.renew_mult .== 0.5 .&& vals.prev_lim .== 1.0,:]
b = make_heatmap(dm15r05)
dm1r1 = vals[vals.demand_mult .== 1.0 .&& vals.renew_mult .== 1.0 .&& vals.prev_lim .== 1.0,:]
c = make_heatmap(dm1r1)
Plots.plot!(c, ylabel="Renew 100%, 4697 MW", xlabel="Load 100%, 3616 MW")
dm15r1 = vals[vals.demand_mult .== 1.5 .&& vals.renew_mult .== 1.0 .&& vals.prev_lim .== 1.0,:]
d = make_heatmap(dm15r1)
Plots.plot!(d, xlabel="Load 150%, 5424 MW")
Plots.plot(a,b,c,d, size=(700,750), layout=Plots.grid(2,2), thickness_scaling=1.2)

# Single heatmap
vals = DataFrames.DataFrame(demand_mult=Float64[], renew_mult=Float64[], i=Int64[], Nk_C=Float64[], Nk_zC=Float64[], Nk_zVOLL=Float64[])
for demand_mult in 0.5:0.25:1.5, renew_mult in 0.0:0.25:1.0
    name = Printf.@sprintf "zc_prev%d_ls%d_pd%d_pr%d_ramp%d_r%d_%s" prev_lim*100 max_shed*100 demand_mult*100 renew_mult*100 renew_ramp reserve*100 string(t)[1:end-6]
    valb = CSV.read("results/"*name*"_base_cost.csv", DataFrames.DataFrame)
    valc = CSV.read("results/"*name*"_corrective_cost.csv", DataFrames.DataFrame)
    # push!(vals, (prev_lim, max_shed, demand_mult, renew_mult, renew_ramp, i, valb[1,1:5]..., valc[1,1:5]...), promote=true)
    push!(vals, (demand_mult, renew_mult, i, Matrix(valb .+ valc)[27,:][1:3]...), promote=true)
    # push!(vals, (demand_mult, renew_mult, i, sum(Matrix(valb .+ valc), dims=1)[1:3]...), promote=true)
end
max_pd = sum(get_active_power.(demands))
max_renewable = sum(get_active_power(g) for g in SCOPF.sort_components!(SCOPF.get_generation(system)) if typeof(g) <: RenewableGen)
sort!(vals, :demand_mult)
data = reshape(Matrix(vals[!,r"Nk_zVOLL"]), (5,5)) ./ 10
pl = Plots.heatmap(data', 
    yticks=(1:5,Int.(round.(vals[!,:demand_mult][1:5:end].*max_pd.*100, digits=0))), 
    xticks=(1:size(data,1),Int.(round.(vals[!,:renew_mult][1:5].*max_renewable.*100, digits=0))), 
    c=Plots.cgrad(:vik), colorbar_scale=:log10, colorbar_title="Total cost [k\$]",
    xlabel="Renewable penetration [MW]", ylabel="Load level [MW]", size=(430,400), rightmargin=0*Plots.Measures.mm, aspect_ratio=:equal)
Plots.annotate!([(j,i,Plots.text(Int(round(data[j,i], digits=0)), 8, :black, :center)) for j in 1:size(data,1) for i in 1:size(data,2)], linecolor=:white)


bx = SCOPF.get_bus_idx.(SCOPF.sort_components!(SCOPF.get_branches(system)), [SCOPF.get_nodes_idx(SCOPF.sort_components!(SCOPF.get_nodes(system)))]);
forced_ls = [SCOPF.calc_forced_ls(case, x[2], x[1]) for x in scenarioes_n1]
Plots.plot([sum(p * SCOPF.is_islanded(case.pf, bx[c[2]], c[2]) for (c,p) in zip(s...) if length(c[2]) > 0) for s in scenarioes_n1].*100, xlabel="Time [h]", ylabel="Post-contingency islanding [%]", labels="Left", leg=:topleft)
Plots.plot!(Plots.twinx(), c=:red, [sum(s[2][i]*x for (i,x) in l) for (s,l) in zip(scenarioes_n1,forced_ls)].*100, labels="Right", leg=:topright, size=(500,250), ylabel="Forced load shed [MW]")
# Plots.plot([count(SCOPF.is_islanded(case.pf, bx[c[2]], c[2]) for c in s[1] if length(c[2]) > 0) for s in scenarioes_n1], xlabel="Time [h]", ylabel="Count of islands", labels="Left", leg=:topleft)
# Plots.plot!(Plots.twinx(), c=:red, sum.(x for (i,x) in forced_ls), labels="Right", leg=:topright, size=(500,300), ylabel="Forced load shed [MW]")

res_new = []
for (i,(r, s)) in enumerate(zip(res, scenarioes_n1))
    println("\n",i)
    case, tot_t = SCOPF.run_benders(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, s[2], s[1], max_shed=max_shed, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10)
    push!(res_new, run_corrective_with_base(case, r))
end
base_costs_new = [isnan(x[:obj_val]) ? NaN : x[:costs][1,:Base] for x in res_new]
corrective_costs_new = [isnan(x[:obj_val]) ? NaN : (x[:costs][1,:r] + sum(p[2] .* (x[:costs][:,:Pc] + x[:costs][:,:Pcc]))) for (x,p) in zip(res_new,scenarioes_n1)]

Plots.plot(base_costs[2] + corrective_costs[2], labels=labels, line=(:steppre, :dash), alpha=0.5, xlabel="Time", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6)[1:5])
Plots.plot!(base_costs_new + corrective_costs_new, labels=labels.*" no corr cost", line=:steppre, xlabel="Time", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6)[1:5], leg=:topright)
Plots.savefig("results/"*name*"_objective_no-corr-cost.pdf")



consequences = Dict()
for (s, case) in zip(scenarioes_n1, res_sn1)
    for (c,v) in zip(s[1], case[:costs][:,:Pcc]) 
        if v != 0.0
            for i in c[2]
                if get(consequences, i, 0) == 0
                    consequences[i] = [v] 
                else
                    push!(consequences[i], v) 
                end
            end
        end
    end
end
for (_,c) in consequences
    sort!(unique!(c), rev=true)
end
sort(consequences, by=x->x[2])




vals3 = []
for i in 0:23
    t = Dates.DateTime("2020-01-01T"*string(i)*":00:00")
    demands = SCOPF.sort_components!((SCOPF.get_demands(system)))
    renewables = SCOPF.sort_components!((SCOPF.get_renewables(system)))
    set_time_series_value!(demands, t)
    set_active_power!.(demands, get_active_power.(demands)*1.5)
    set_time_series_value!(renewables, t)
    set_active_power!.(renewables, get_active_power.(renewables)*2.0)
    @time case, tot_t = SCOPF.run_benders(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(GRB_ENV), voll, prob/8760, contingencies, max_shed=max_shed, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10, debug=true);
    max_demand = sum(get_active_power(g) for g in demands)
    max_renewable = sum(get_active_power(g) for g in gen if typeof(g) <: RenewableGen)
    prod = SCOPF.get_type_prod(system, JuMP.value.(case.model[:pg0]))
    push!(vals3, (prod, max_renewable, max_demand))
end


import FileIO, JLD2
a = 1
FileIO.save("myfile.jld2","a",a)
b = FileIO.load("myfile.jld2","a")