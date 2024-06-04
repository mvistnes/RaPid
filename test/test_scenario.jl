import Dates, DataFrames, CSV, Plots, StatsPlots, Gurobi, Printf, FileIO, JLD2, Random, Logging, JuMP
import RaPidSCOPF as SCOPF
Random.seed!(53715)

## Run this line in REPL
# system = include("data/RTS_GMLC/config.jl")

df = DataFrames.DataFrame(CSV.File("data/RTS_GMLC/branch.csv"))
branches = SCOPF.sort_components!(SCOPF.get_branches(system));
ia_br = SCOPF.get_interarea(branches)
brc = reduce(vcat, [ia_br, branches[findall(x->x.name ∈ ["A30", "A21", "B30", "B21", "C30", "C21"], branches)]])

prob = reduce(vcat, [[df[!,"Perm OutRate"][i] for i in 1:length(branches) if x == branches[i].name] for x in df.UID])
w = [branches[i] ∈ brc ? Storm(x/100/8760, 0.02*(0.5+rand()), rand(20:25), rand(1:3), rand(10:14), rand(10:14)) : storm() for (i,x) in enumerate(prob)]
w2 = [branches[i] ∈ brc ? Storm(x/100/8760, 0.02*(0.5+rand()), rand(25:30), rand(1:3), rand(10:14), rand(10:14)) : storm() for (i,x) in enumerate(prob)]
p = generate_p_from_weather(prob/8760, 60, w)
p2 = generate_p_from_weather(prob/8760, 60, w2)
scenarioes = generate_contingency_probabilities(p, 5, 10000)
scenarioes2 = generate_contingency_probabilities(p2, 5, 10000)

style = Plots.supported_styles()[2:end]
style = reshape(style, 1, length(style))
bnames = first.(sort(get_name.(branches) .=> w, by=x->x[2].max)[1:11])
Plots.plot([[last(sort(x[2])) for x in scenarioes] [last(sort(x[2])) for x in scenarioes2]], labels=[])
Plots.plot(reduce(hcat, generate_weather.([60], sort(w, by=x->x.max)[1:11])), label=permutedims(bnames), 
    xlabel="Time", ylabel="Probability", palette = Plots.palette(:seaborn_colorblind6))
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
StatsPlots.groupedbar([sum(outs[:,1:3], dims=2) outs[:,4:end]], bar_position = :stack, xlabel="Time", ylabel="Probability", 
    palette = Plots.palette(:seaborn_colorblind6, rev=true), labels=["N-k > 3" "N-3" "N-2" "N-1" "N-0"], leg=:bottomright)
Plots.plot(cumsum([outs[:,1:4] sum(outs[:,5:end], dims=2)], dims = 2)[:,end:-1:1], fill=0, lc=:black, legend=:none, 
    xlabel="Time", ylabel="Probability", palette = Plots.palette(:seaborn_colorblind6, rev=true), leg=:bottomright, 
    annotation=[(27, 0.2, ("N-0", 10)), (27, 0.6, ("N-1", 10)), (27, 0.83, ("N-2", 10)), (27, 0.96, ("N-3", 10)), (27, 1.03, ("N-k > 3", 10))])
Plots.plot([collect(values(sort(SCOPF.countmemb(length.(last.(x[1]))), by=i->i[1]))) for x in scenarioes])
Plots.plot([collect(values(sort(SCOPF.sumvals(length.(last.(x[1])), x[2]), by=i->i[1]))) for x in scenarioes])

function scenario_plus_n1(scenarioes)
    scenarioes_n1 = deepcopy(scenarioes)
    for (i,x) in enumerate(scenarioes)
        diff = setdiff(1:120, reduce(vcat, [x for x in getindex.(x[1], 2) if length(x) == 1]))
        for c in diff
            push!(scenarioes_n1[i][1], "branch" => [c])
            push!(scenarioes_n1[i][2], prob[c]/8760)
        end
    end
    return scenarioes_n1
end
scenarioes_n1 = scenario_plus_n1(scenarioes)
scenarioes2_n1 = scenario_plus_n1(scenarioes2)

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
ramp_mult = 2.0
demand_mult = 1.0
renew_mult = 1.0
renew_cost = 0.00
renew_ramp = 0
i = 9
t = Dates.DateTime("2020-04-11T13:00:00")
# t = Dates.DateTime("2020-08-25T14:00:00")
time_limit_sec = length(contingencies)^2 + 10
labels = ["S-P" "S-C" "N-1-P" "N-1-C" "N-0"]
# xguidefont=fsize, yguidefont=fsize, xtickfont=fsize, ytickfont=fsize, legendfont=fsize
# for #=prev_lim in [short, long],=# max_shed in [1.5, 3.0], demand_mult in [1.0, 1.5], renew_mult in [0.5, 1.0], renew_ramp in [0.0, 10.0]# , i in [2,9,17]
for max_shed in 2.5:5.0
    name = Printf.@sprintf "prev-%d_ls%d_pd%d_pr%dramp%d_%s" prev_lim*100 max_shed*100 demand_mult*100 renew_mult*100 renew_ramp string(t)[1:end-6]
    println("\n\n"*name*"\n")
    # t = Dates.DateTime("2020-01-01T"*string(i)*":00:00")
    set_time_series_value!(demands, t)
    set_active_power!.(demands, get_active_power.(demands)*demand_mult)
    set_time_series_value!(renewables, t)
    set_active_power!.(renewables, get_active_power.(renewables)*renew_mult)
    # res, res_prev, res_prev_n1, res_n1, res_base = run_cases(system, scenarioes_n1, voll, prob, contingencies, max_shed, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim);
    res_prev, res, res_prev_n1, res_n1, res_base = FileIO.load("results/"*name*".jld2", "Scen_Prev", "Scen_Prob","N-1_Prev", "N-1_Prob", "N-0");

    base_costs = [[isnan(x[:obj_val]) ? NaN : x[:costs][1,:Base] for x in r] for r in [res_prev, res, res_prev_n1, res_n1, res_base]]
    corrective_costs = [[isnan(x[:obj_val]) ? NaN : (x[:costs][1,:r] + sum(p[2] .* (x[:costs][:,:Pc] + x[:costs][:,:Pcc]))) for (x,p) in zip(r,scenarioes_n1)] for r in [res_prev,res,res_prev_n1,res_n1,res_base]]
    
    CSV.write("results/"*name*"_base_cost.csv", DataFrames.DataFrame((base_costs, 
        ["Scen_Prev", "Scen_Prob", "N-1_Prev", "N-1_Prob", "N-0"])...))
    CSV.write("results/"*name*"_corrective_cost.csv", DataFrames.DataFrame((corrective_costs, 
        ["Scen_Prev", "Scen_Prob", "N-1_Prev", "N-1_Prob", "N-0"])...)) 
    FileIO.save("results/"*name*".jld2", "Scen_Prev", res_prev, "Scen_Prob", res, "N-1_Prev", res_prev_n1, "N-1_Prob", res_n1, "N-0", res_base)

    Plots.savefig(Plots.plot(base_costs, labels=labels, 
        line=:steppre, xlabel="Time", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6), thickness_scaling=1.2),
        "results/"*name*"_base_cost.pdf")

    Plots.savefig(Plots.plot(corrective_costs,labels=labels, 
        line=:steppre, xlabel="Time", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6), thickness_scaling=1.2),
        "results/"*name*"_corrective_cost.pdf")

    Plots.savefig(Plots.plot(base_costs + corrective_costs,  labels=labels, 
        line=:steppre, xlabel="Time", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6), thickness_scaling=1.2),
        "results/"*name*"_objective.pdf")

    max_renewable = sum(get_active_power(g) for g in gen if typeof(g) <: RenewableGen)
    Plots.savefig(Plots.plot([[SCOPF.get_all_type_prod(system, r)[3,:] for r in [res_prev, res, res_prev_n1, res_n1, res_base]]... fill(max_renewable, size(scenarioes_n1,1))], 
        label=[labels... "Max renewable"], 
        line=:steppre, xlabel="Time", ylabel="Amount", palette = Plots.palette(:seaborn_colorblind6), leg=:topright, thickness_scaling=1.2),
        "results/"*name*"_renewable.pdf")
    
    res2 = []
    for (i,(r, s)) in enumerate(zip(res, scenarioes2_n1))
        println("\nScen.2 ", i, " prob")
        case, tot_t = SCOPF.run_benders(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(), voll, s[2], s[1], max_shed=max_shed, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10)
        push!(res2, run_corrective_with_base(case, r))
    end
    res2_prev = []
    for (i,(r, s)) in enumerate(zip(res_prev, scenarioes2_n1))
        println("\nScen.2 ", i, " preventive")
        case, tot_t = SCOPF.run_benders(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(), voll, s[2], s[1], max_shed=max_shed, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10)
        push!(res2_prev, run_corrective_with_base(case, r))
    end
    _, _, res2_prev_n1, res2_n1, res2_base = run_cases(system, scenarioes2_n1, voll, prob, contingencies, max_shed, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim);
    base_costs2 = [[isnan(x[:obj_val]) ? NaN : x[:costs][1,:Base] for x in r] for r in [res2_prev, res2, res2_prev_n1, res2_n1, res2_base]]
    corrective_costs2 = [[isnan(x[:obj_val]) ? NaN : (x[:costs][1,:r] + sum(p[2] .* (x[:costs][:,:Pc] + x[:costs][:,:Pcc]))) for (x,p) in zip(r,scenarioes2_n1)] for r in [res2_prev,res2,res2_prev_n1,res2_n1,res2_base]]    

    Plots.plot(corrective_costs, labels=labels, line=(:steppre, :dash), alpha=0.5, xlabel="Time", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6)[1:5])
    Plots.plot!(corrective_costs2, labels=labels.*" S2", line=:steppre, xlabel="Time", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6)[1:5], leg=:topright, thickness_scaling=1.2)
    Plots.savefig("results/"*name*"_corrective_cost_s2.pdf")

    Plots.plot(base_costs + corrective_costs, labels=labels, line=(:steppre, :dash), alpha=0.5, xlabel="Time", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6)[1:5])
    Plots.plot!(base_costs2 + corrective_costs2, labels=labels.*" S2", line=:steppre, xlabel="Time", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6)[1:5], leg=:topright, thickness_scaling=1.2)
    Plots.savefig("results/"*name*"_objective_s2.pdf")
    
    CSV.write("results/"*name*"_base_cost.csv", DataFrames.DataFrame(([base_costs... base_costs2...], 
        ["Scen_Prev", "Scen_Prob", "N-1_Prev", "N-1_Prob", "N-0", "Scen_Prev_S2", "Scen_Prob_S2", "N-1_Prev_S2", "N-1_Prob_S2", "N-0_S2"])...))
    CSV.write("results/"*name*"_corrective_cost.csv", DataFrames.DataFrame(([corrective_costs... corrective_costs2...], 
        ["Scen_Prev", "Scen_Prob", "N-1_Prev", "N-1_Prob", "N-0", "Scen_Prev_S2", "Scen_Prob_S2", "N-1_Prev_S2", "N-1_Prob_S2", "N-0_S2"])...)) 
    FileIO.save("results/"*name*".jld2", "Scen_Prev", res_prev, "Scen_Prob", res, "N-1_Prev", res_prev_n1, "N-1_Prob", res_n1, "N-0", res_base, "Scen_Prev_S2", res2_prev_n1, "Scen_Prob_S2", res2, "N-1_Prev_S2", res2_prev, "N-1_Prob_S2", res2_n1, "N-0_S2", res2_base)
    
end

sensitivity = []
itr = 0.5:5.0 # max_shed
itr = 0.25:0.25:1.5 # renew_mult
itr = 0.6:0.2:1.4 # prev_lim
itr = [0.0, 1.0, 10.0, 100.0, 1000.0] # renew_ramp
for renew_ramp in itr
    name = Printf.@sprintf "prev-%d_ls%d_pd%d_pr%dramp%d_%s" prev_lim*100 max_shed*100 demand_mult*100 renew_mult*100 renew_ramp string(t)[1:end-6]
    # println("\n\n"*name*"\n")
    # t = Dates.DateTime("2020-01-01T"*string(i)*":00:00")
    # set_time_series_value!(demands, t)
    # set_active_power!.(demands, get_active_power.(demands)*1.5)
    # set_time_series_value!(renewables, t)
    # set_active_power!.(renewables, get_active_power.(renewables)*renew_mult)
    # res, res_prev, res_prev_n1, res_n1, res_base = run_cases(system, scenarioes_n1[[30]], voll, prob, contingencies, max_shed, ramp_mult, renew_cost, renew_ramp, ramp_minutes, short, long, prev_lim);
    # res_prev, res, res_prev_n1, res_n1, res_base, res2_prev_n1, res2, res2_prev, res2_n1, res2_base = FileIO.load("results/"*name*".jld2", "Scen_Prev", "Scen_Prob","N-1_Prev", "N-1_Prob", "N-0", "Scen_Prev_S2", "Scen_Prob_S2", "N-1_Prev_S2", "N-1_Prob_S2", "N-0_S2");
    res_prev, res, res_prev_n1, res_n1, res_base = FileIO.load("results/"*name*".jld2", "Scen_Prev", "Scen_Prob","N-1_Prev", "N-1_Prob", "N-0");

    base_costs = [[isnan(x[:obj_val]) ? NaN : x[:costs][1,:Base] for x in r] for r in [res_prev, res, res_prev_n1, res_n1, res_base]]
    corrective_costs = [[isnan(x[:obj_val]) ? NaN : (x[:costs][1,:r] + sum(p[2] .* (x[:costs][:,:Pc] + x[:costs][:,:Pcc]))) for (x,p) in zip(r,scenarioes_n1)] for r in [res_prev,res,res_prev_n1,res_n1,res_base]]
    # base_costs2 = [[isnan(x[:obj_val]) ? NaN : x[:costs][1,:Base] for x in r] for r in [res2_prev, res2, res2_prev_n1, res2_n1, res2_base]]
    # corrective_costs2 = [[isnan(x[:obj_val]) ? NaN : (x[:costs][1,:r] + sum(p[2] .* (x[:costs][:,:Pc] + x[:costs][:,:Pcc]))) for (x,p) in zip(r,scenarioes2_n1)] for r in [res2_prev,res2,res2_prev_n1,res2_n1,res2_base]]    

    push!(sensitivity, prev_lim => sum.(base_costs) .+ sum.(corrective_costs))
    # push!(sensitivity, prev_lim => sum.(reduce(vcat, [base_costs, base_costs2])) .+ sum.(reduce(vcat, [corrective_costs, corrective_costs2])))
    # push!(sensitivity, prev_lim => Dict(:base_costs => base_costs, :corrective_costs => corrective_costs))
end
Plots.plot(reduce(hcat, [reduce(vcat, x[2][:base_costs]) for x in sensitivity])',labels=labels, xticks=(1:length(itr), itr), xlabel="Renewable corrective cost", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6), thickness_scaling=1.2)
Plots.plot(reduce(hcat, [reduce(vcat, x[2][:corrective_costs]) for x in sensitivity])',labels=labels, xticks=(1:length(itr), itr), xlabel="Renewable corrective cost", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6), thickness_scaling=1.2)
Plots.plot(reduce(hcat, [reduce(vcat, x[2][:base_costs] .+ x[2][:corrective_costs]) for x in sensitivity])',labels=labels, xticks=(1:length(itr), Int.(itr.*100)), xlabel="Renewable availability [%]", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6), thickness_scaling=1.2, leg=:topright)
Plots.plot(reduce(hcat, [reduce(vcat, x[2]) for x in sensitivity])'[:,1:5], line=:dash, alpha=0.5, labels=labels, xticks=(1:length(itr), Int.(itr.*100)), xlabel="Renewable availability [%]", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6)[1:5], thickness_scaling=1.2, leg=:topright)
Plots.plot!(reduce(hcat, [reduce(vcat, x[2]) for x in sensitivity])'[:,6:end],labels=labels .* " 2S", xticks=(1:length(itr), Int.(itr.*100)), xlabel="Load shed limit [MW]", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6)[1:5], thickness_scaling=1.2, leg=:topright)

# save all with => reduce(hcat, [reduce(vcat, x[2][:base_costs] .+ x[2][:corrective_costs]) for x in sensitivity])'
# then 
a = Plots.plot(rm ./1000, labels=:none, xlabel="Renewable availability [%]", ylabel="Cost", xticks=(1:6, Int.((0.25:0.25:1.5).*100)))
b = Plots.plot(pl ./1000, labels=:none, xlabel="Preventive limit [%]", xticks=(1:5, Int.((0.6:0.2:1.4).*100)))
c = Plots.plot(ms ./1000, labels=:none, xlabel="Load shed limit [MW]", xticks=(1:5, Int.((0.5:5.0).*100)))
d = Plots.plot(rr[1:end-1,:] ./1000, labels=labels, leg=:outertopright, xlabel="Renewable corrective cost", xticks=(1:4, Int.([0.0, 1.0, 10.0, 100.0])))
Plots.plot(a,b,c,d, layout=Plots.grid(1,4, widths=(5/17,4/17,4/17,4/17)), size=(1200,300), yaxis=[0,15], bottom_margin=7Measures.mm, left_margin=7Measures.mm, palette = Plots.palette(:seaborn_colorblind6), thickness_scaling=1.0)

ty = [["short", "long"], [0.5, 5.0], [1.5, 1.0], [0.0, 1000.0]]
choose(ty, i, c) = c == i ? ty[i][2] : ty[i][1]
for c in 1:length(ty)
    name = Printf.@sprintf "prev-%s_ls%d_pd150_pr%dramp%d_time%d" choose(ty,1,c) choose(ty,2,c)*100 choose(ty,3,c)*100 choose(ty,4,c) i
    res_prev, res, res_prev_n1, res_n1, res_base = FileIO.load("results/"*name*".jld2", "Scen_Prev", "Scen_Prob","N-1_Prev", "N-1_Prob", "N-0");
    base_costs = [[isnan(x[:obj_val]) ? NaN : x[:costs][1,:Base] for x in r] for r in [res_prev, res, res_prev_n1, res_n1, res_base]]
    corrective_costs = [[isnan(x[:obj_val]) ? NaN : (x[:costs][1,:r] + sum(p[2] .* (x[:costs][:,:Pc] + x[:costs][:,:Pcc]))) for (x,p) in zip(r,scenarioes_n1)] for r in [res_prev,res,res_prev_n1,res_n1,res_base]]
    Plots.plot(results_base_costs + results_corrective_costs, labels=labels, line=(:steppre, :dash), alpha=0.5, xlabel="Time", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6)[1:5])
    Plots.plot!(base_costs + corrective_costs, labels=labels.*" 2", line=:steppre, xlabel="Time", ylabel="Cost", palette = Plots.palette(:seaborn_colorblind6)[1:5], leg=:topright)
    Plots.savefig("results/"*name*"_objective_2.pdf")
end

# HEATMAP
vals30 = DataFrames.DataFrame(prev_lim=Float64[], max_shed=Float64[], demand_mult=Float64[], renew_mult=Float64[], renew_ramp=Float64[], i=Int64[], Scen_Prev_B=Float64[], Scen_Prob_B=Float64[], N_1_Prev_B=Float64[], N_1_Prob_B=Float64[], N_0_B=Float64[], Scen_Prev_C=Float64[], Scen_Prob_C=Float64[], N_1_Prev_C=Float64[], N_1_Prob_C=Float64[], N_0_C=Float64[])
for prev_lim in [short, long], max_shed in [1.5, 3.0], demand_mult in [1.0, 1.5], renew_mult in [0.5, 1.0], renew_ramp in [0.0, 10.0]# , i in [2,9,17]
    name = Printf.@sprintf "prev-%d_ls%d_pd%d_pr%dramp%d_%s" prev_lim*100 max_shed*100 demand_mult*100 renew_mult*100 renew_ramp string(t)[1:end-6]    
    valb = CSV.read("results/"*name*"_base_cost.csv", DataFrames.DataFrame)
    valc = CSV.read("results/"*name*"_corrective_cost.csv", DataFrames.DataFrame)
    # push!(vals30, (prev_lim, max_shed, demand_mult, renew_mult, renew_ramp, i, valb[1,1:5]..., valc[1,1:5]...), promote=true)
    push!(vals30, (prev_lim, max_shed, demand_mult, renew_mult, renew_ramp, i, sum(Matrix(valb), dims=1)[1:5]..., sum(Matrix(valc), dims=1)[1:5]...), promote=true)
end
# dataB = vals30[!,r"B"]
# dataC = vals30[!,r"C"]
# data = Matrix(dataB) .+ Matrix(dataC);
dm1r05 = vals30[vals30.demand_mult .== 1.0 .&& vals30.renew_mult .== 0.5 .&& vals30.prev_lim .== 1.0,:]
data = Matrix(dm1r05[!,r"B"]) .+ Matrix(dm1r05[!,r"C"])
data .= data ./ 1000
a = Plots.heatmap(data, xticks=(1:5,["SP", "SC", "N1P", "N1C", "N0"]), yticks=(1:size(data,1)), c=Plots.cgrad(:seaborn_icefire_gradient), ylabel="Renew 50%")
Plots.annotate!([(j,i,Plots.text(Int(round(data[i,j], digits=0)), 8, :black, :center)) for i in 1:size(data,1) for j in 1:size(data,2)], linecolor=:white)
dm15r05 = vals30[vals30.demand_mult .== 1.5 .&& vals30.renew_mult .== 0.5 .&& vals30.prev_lim .== 1.0,:]
data = Matrix(dm15r05[!,r"B"]) .+ Matrix(dm15r05[!,r"C"])
data .= data ./ 1000
b = Plots.heatmap(data, xticks=(1:5,["SP", "SC", "N1P", "N1C", "N0"]), yticks=(1:size(data,1)), c=Plots.cgrad(:seaborn_icefire_gradient))
Plots.annotate!([(j,i,Plots.text(Int(round(data[i,j], digits=0)), 8, :black, :center)) for i in 1:size(data,1) for j in 1:size(data,2)], linecolor=:white)
dm1r1 = vals30[vals30.demand_mult .== 1.0 .&& vals30.renew_mult .== 1.0 .&& vals30.prev_lim .== 1.0,:]
data = Matrix(dm1r1[!,r"B"]) .+ Matrix(dm1r1[!,r"C"])
data .= data ./ 1000
c = Plots.heatmap(data, xticks=(1:5,["SP", "SC", "N1P", "N1C", "N0"]), yticks=(1:size(data,1)), c=Plots.cgrad(:seaborn_icefire_gradient), ylabel="Renew 100%", xlabel="Load 100%")
Plots.annotate!([(j,i,Plots.text(Int(round(data[i,j], digits=0)), 8, :black, :center)) for i in 1:size(data,1) for j in 1:size(data,2)], linecolor=:white)
dm15r1 = vals30[vals30.demand_mult .== 1.5 .&& vals30.renew_mult .== 1.0 .&& vals30.prev_lim .== 1.0,:]
data = Matrix(dm15r1[!,r"B"]) .+ Matrix(dm15r1[!,r"C"])
data .= data ./ 1000
d = Plots.heatmap(data, xticks=(1:5,["SP", "SC", "N1P", "N1C", "N0"]), yticks=(1:size(data,1)), c=Plots.cgrad(:seaborn_icefire_gradient), xlabel="Load 150%")
Plots.annotate!([(j,i,Plots.text(Int(round(data[i,j], digits=0)), 8, :black, :center)) for i in 1:size(data,1) for j in 1:size(data,2)], linecolor=:white)
Plots.plot(a,b,c,d, size=(800,600), layout=Plots.grid(2,2), thickness_scaling=1.2)


bx = SCOPF.get_bus_idx.(SCOPF.sort_components!(SCOPF.get_branches(system)), [SCOPF.get_nodes_idx(SCOPF.sort_components!(SCOPF.get_nodes(system)))]);
forced_ls = [SCOPF.calc_forced_ls(case, x[2], x[1]) for x in scenarioes_n1]
Plots.plot([count(SCOPF.is_islanded(case.pf, bx[c[2]], c[2]) for c in s[1] if length(c[2]) > 0) for s in scenarioes_n1], xlabel="Time", ylabel="Count of islands", leg=:none)
Plots.plot!(Plots.twinx(), c=:red, sum.(forced_ls), leg=:none, ylabel="Forced load shed by islands")

res_new = []
for (i,(r, s)) in enumerate(zip(res, scenarioes_n1))
    println("\n",i)
    case, tot_t = SCOPF.run_benders(SCOPF.Base_SCOPF, system, Gurobi.Optimizer(), voll, s[2], s[1], max_shed=max_shed, ramp_mult=ramp_mult, renew_cost=renew_cost, renew_ramp=renew_ramp, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10)
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
    @time case, tot_t = SCOPF.run_benders(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(), voll, prob/8760, contingencies, max_shed=max_shed, ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10, debug=true);
    max_demand = sum(get_active_power(g) for g in demands)
    max_renewable = sum(get_active_power(g) for g in gen if typeof(g) <: RenewableGen)
    prod = SCOPF.get_type_prod(system, JuMP.value.(case.model[:pg0]))
    push!(vals3, (prod, max_renewable, max_demand))
end


import FileIO, JLD2
a = 1
FileIO.save("myfile.jld2","a",a)
b = FileIO.load("myfile.jld2","a")