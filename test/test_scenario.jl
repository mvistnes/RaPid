import Dates, DataFrames, CSV, Plots, StatsPlots
Random.seed!(53715)

# system = include("data/RTS_GMLC/config.jl")
# system, voll, contingencies, prob, short, long, ramp_minutes, ramp_mult, max_shed, time_limit_sec = setup_system(joinpath("data","matpower","RTS_GMLC.m"));

df = DataFrames.DataFrame(CSV.File("data/RTS_GMLC/branch.csv"))
voll, _, contingencies = SCOPF.setup(system, 100.0, 400.0)
SCOPF.set_ramp_limits!(system)
SCOPF.set_renewable_prod!(system, 0.5)
SCOPF.set_operation_cost!.(SCOPF.get_gens_h(system), rand(20) .* (55 - 45) .+ 45)
SCOPF.set_rate!.(SCOPF.get_branches(system), SCOPF.get_rate.(SCOPF.get_branches(system)) * 0.8);
short = 1.2
long = 1.0
ramp_minutes = 10.0
max_shed = 1.5
ramp_mult = 2.0
renew_cost = 0.00
time_limit_sec = length(contingencies)^2 + 10

branches = SCOPF.sort_components!(SCOPF.get_branches(system));
ia_br = SCOPF.get_interarea(branches)
brc = reduce(vcat, [ia_br, branches[findall(x->x.name ∈ ["A30", "A21", "B30", "B21", "C30", "C21"], branches)]])

prob = reduce(vcat, [[df[!,"Perm OutRate"][i] for i in 1:length(branches) if x == branches[i].name] for x in df.UID])
w = [branches[i] ∈ ia_br ? Storm(x/8760, 0.03*(0.8+rand()*0.4), rand(20:25), rand(1:3), rand(10:14), rand(10:14)) : storm() for (i,x) in enumerate(prob)];
p = generate_p_from_weather(prob/8760, 60, w);
scenarioes = generate_contingency_probabilities(p, 5, 100000);
Plots.plot([last(sort(x[2])) for x in scenarioes])

vals = [collect(values(sort(SCOPF.sumvals(length.(last.(x[1])), x[2]), by=i->i[1]))) for x in scenarioes]
c_size = maximum(length.(vals))
outs = zeros(60,c_size);
for (i,x) in enumerate(vals)
    for (j,y) in enumerate(x)
        outs[i,c_size-j+1] = y
    end
end
StatsPlots.groupedbar(outs, bar_position = :stack, palette = Plots.palette(:rainbow))
Plots.plot([collect(values(sort(SCOPF.countmemb(length.(last.(x[1]))), by=i->i[1]))) for x in scenarioes])
Plots.plot([collect(values(sort(SCOPF.sumvals(length.(last.(x[1])), x[2]), by=i->i[1]))) for x in scenarioes])

function extract_results(case::SCOPF.Case)
    costs = DataFrames.DataFrame(SCOPF.get_costs(case)...)
    pg0 = SCOPF.get_value(case.model, :pg0)
    ls0 = SCOPF.get_value(case.model, :ls0)
    Pc = [SCOPF.get_value(case.model, x) for (i,x) in case.Pc]
    Pcc = [SCOPF.get_value(case.model, x) for (i,x) in case.Pcc]
    obj_val = JuMP.objective_value(case.model)
    return costs, pg0, ls0, Pc, Pcc, obj_val
end

res = [@time case, tot_t = SCOPF.run_benders(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(), voll, s[2], s[1], max_shed=max_shed, ramp_mult=ramp_mult, renew_cost=renew_cost, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10) 
    for s in scenarioes];
res_rob = [@time case, tot_t = SCOPF.run_benders(SCOPF.P_SCOPF, system, Gurobi.Optimizer(), voll, s[2], s[1], max_shed=max_shed, ramp_mult=ramp_mult, renew_cost=renew_cost, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10) 
    for s in scenarioes];
res_rob_n1 = @time case, tot_t = SCOPF.run_benders(SCOPF.P_SCOPF, system, Gurobi.Optimizer(), voll, prob/8760, contingencies, max_shed=max_shed, ramp_mult=ramp_mult, renew_cost=renew_cost, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10);
Plots.plot(JuMP.objective_value.(getproperty.(first.(res),:model))) 

scenarioes_n1 = deepcopy(scenarioes);
for (i,x) in enumerate(scenarioes)
    diff = setdiff(1:120, reduce(vcat, [x for x in getindex.(x[1], 2) if length(x) == 1]))
    for c in diff
        push!(scenarioes_n1[i][1], "branch" => [c])
        push!(scenarioes_n1[i][2], prob[c]/8760)
    end
end
res_sn1 = [@time case, tot_t = SCOPF.run_benders(SCOPF.PC2_SCOPF, system, Gurobi.Optimizer(), voll, s[2], s[1], max_shed=max_shed, ramp_mult=ramp_mult, renew_cost=renew_cost, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10) 
    for s in scenarioes_n1];
res_rob_sn1 = [@time case, tot_t = SCOPF.run_benders(SCOPF.P_SCOPF, system, Gurobi.Optimizer(), voll, s[2], s[1], max_shed=max_shed, ramp_mult=ramp_mult, renew_cost=renew_cost, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=0.00, max_itr=10) 
    for s in scenarioes_n1];

df = [DataFrames.DataFrame(SCOPF.get_costs(x[1])...) for x in res];
df_rob = [DataFrames.DataFrame(SCOPF.get_costs(x[1])...) for x in res_rob];
df_sn1 = [DataFrames.DataFrame(SCOPF.get_costs(x[1])...) for x in res_sn1];
df_rob_sn1 = [DataFrames.DataFrame(SCOPF.get_costs(x[1])...) for x in res_rob_sn1];
style = Plots.supported_styles()[2:end]
style = reshape(style, 1, length(style))
Plots.plot([[x[1,:Base] for x in df_rob] [x[1,:Base] for x in df]], label=["Scenario Robust" "Scenario Prob"], line=style)
Plots.plot([[x[1,:Base] for x in df_sn1] [x[1,:Base] for x in df_rob_sn1]], label=["Scenario Robust N-1" "Scenario Prob N-1"], line=style)

function get_type_prod(system::System, case::SCOPF.Case)
    gens = SCOPF.sort_components!(SCOPF.get_generation(system))
    type_prod = zeros(3)
    for (i,g) in enumerate(gens)
        type = typeof(g)
        val = JuMP.value(case.model[:pg0][i])
        if type <: ThermalGen 
            type_prod[1] += val
        elseif type <: HydroGen 
            type_prod[2] += val
        elseif type <: RenewableGen
            type_prod[3] += val
        else
            @error "Unknown generation type"
        end
    end
    return [ThermalGen, HydroGen, RenewableGen] .=> type_prod
end

function get_all_type_prod(system::System, cases::Vector{<:SCOPF.Case})
    return reduce(hcat, [getindex.(get_type_prod(system, x), 2) for x in cases])
end

slack = SCOPF.find_slack(islands[2], pf.slack, case.oplim.pg_lim_max, case.opf.mgx)
ptdf = SCOPF.calc_isf(case.pf.DA, case.pf.B, bx[c_i], c_i, slack, islands[2], islands_b[2])
flow = ptdf * JuMP.value.(res[40][1].model[:p0])
(SCOPF.calc_Pᵢ_from_flow(branches, flow, length(nodes), case.opf.idx) .- p0)[islands[2]]

function generate_states(branches::Integer, probability::Real)
    vals = Dict()
    for _ in 1:10000
        x = branches - count([probability < rand() for _ in 1:branches])
        vals[x] = get(vals, x, 0) + 1
    end
    return [i => v/10000 for (i,v) in sort(vals, by=x->x[1])]
end

