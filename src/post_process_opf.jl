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
    sort_x!(list) = sort!(list, by = x -> x.name)
    get_injection(g::Union{ThermalGen, HydroGen}) = value(opfm.mod[:pg0][get_name(g)])
    get_injection(r::RenewableGen) = get_active_power(r) - value(opfm.mod[:ls0][get_name(r)])
    get_injection(d::StaticLoad) =  -get_active_power(d) + value(opfm.mod[:ls0][get_name(d)])
    print("       Bus   Total  Injections...")
    tot = 0.0
    for (n, g, d, r) in zip(nodes, list_gen, list_d, list_r)
        sort_x!(g)
        sort_x!(d)
        sort_x!(r)
        val_g = get_injection.(g)
        val_d = get_injection.(d)
        val_r = get_injection.(r)
        b_tot = sum(val_g, init=0.0) + sum(val_d, init=0.0) + sum(val_r, init=0.0)
        tot += b_tot
        @printf "\n%10s  %6.2f  " n.name b_tot
        xprint.(g, val_g)
        xprint.(d, val_d)
        xprint.(r, val_r)
    end
    @printf("\n       Sum  %6.2f\n", tot)
end

function print_power_flow(opfm::OPFmodel)
    branches = get_sorted_branches(opfm.sys)
    print_power_flow(
            get_name.(branches), 
            value.(opfm.mod[:pf0][get_name.(branches)]).data, 
            get_rate.(branches)
        )
end
function print_power_flow(names::AbstractVector{String}, flow::AbstractVector{<:Real}, rate::AbstractVector{<:Real})
    println("         Branch    Flow  Rating")
    string_line(n, f, r) = @sprintf("%15s  %6.2f  %6.2f\n", n, f, r)
    for (n,f,r) in zip(names, flow, rate)
        if abs(f) > r + 0.0001
            printstyled(string_line(n, f, r), color = :red)
        elseif abs(f) > r * 0.9
            printstyled(string_line(n, f, r), color = :yellow)
        else
            print(string_line(n, f, r))
        end
    end
end

function print_contingency_power_flow(opfm::OPFmodel, pf, ΔP, short_term_limit_multi::Float64 = 1.5)
    branches = get_sorted_branches(opfm.sys)
    nodes = get_sorted_nodes(opfm.sys)
    idx = get_nodes_idx(nodes)
    b_names = get_name.(branches)
    linerates = get_rate.(branches) * short_term_limit_multi
    println("Base case")
    print_power_flow(b_names, pf.F, linerates)
    flow = Vector{Float64}()
    for (c,cont) in enumerate(branches)
        println("Contingency ", cont.name)
        (f, t) = get_bus_idx(cont, idx)
        δP = get_ΔP(opfm, length(nodes), idx, ΔP, c)
        try
            flow = calculate_line_flows(pf, (f, t), c, 
                        (pf.Pᵢ .+ get_ΔP(opfm, length(nodes), idx, ΔP, c)))
        catch DivideError
            island, island_b = handle_islands(pf, get_bus_idx(cont, idx), c)
            ptdf = zeros(eltype(pf.ϕ), size(pf.ϕ))
            ptdf[island_b, island] = get_isf(pf.X, pf.B, pf.DA, (f, t), c, island,)
            flow = ptdf * (pf.Pᵢ .+ δP)
        end
        if isempty(flow)
            println("Contingency ", cont.name, " resulted in a disconnected reference bus.")
        else
            print_power_flow(b_names, flow, linerates)
        end
    end
end

function print_contingency_power_flow(opfm::OPFmodel, short_term_limit_multi::Float64 = 1.5)
    branches = get_sorted_branches(opfm.sys)
    nodes = get_sorted_nodes(opfm.sys)
    idx = get_nodes_idx(nodes)
    slack = find_slack(nodes)
    b_names = get_name.(branches)
    linerates = get_rate.(branches) * short_term_limit_multi
    P = get_net_Pᵢ(opfm.mod, opfm.sys, nodes, idx)
    println("Base case")
    print_power_flow(b_names, value.(opfm.mod[:pf0][get_name.(branches)]).data, linerates)
    for (c,cont) in enumerate(branches)
        println("Contingency ", cont.name)
        flow = get_ptdf(branches[1:end .!= c], length(nodes), idx, slack[1]) * P
        print_power_flow(b_names[1:end .!= c], flow, linerates[1:end .!= c])
    end
end

function print_contingency_overflow(opfm::OPFmodel, pf, ΔP, short_term_limit_multi::Float64 = 1.5)
    string_line(c, n, f, r) = @sprintf("%15s %15s  %6.2f  %6.2f\n", c, n, f, r)
    println("    Contingency        Branch    Flow  Rating")
    branches = get_sorted_branches(opfm.sys)
    nodes = get_sorted_nodes(opfm.sys)
    idx = get_nodes_idx(nodes)
    b_names = get_name.(branches)
    linerates = get_rate.(branches) * short_term_limit_multi
    flow = Vector{Float64}()
    for (c,cont) in enumerate(branches)
        (f, t) = get_bus_idx(cont, idx)
        δP = get_ΔP(opfm, length(nodes), idx, ΔP, c)
        try
            flow = calculate_line_flows(pf, (f, t), c, (pf.Pᵢ .+ δP))
        catch DivideError
            island, island_b = handle_islands(pf, get_bus_idx(cont, idx), c)
            ptdf = zeros(eltype(pf.ϕ), size(pf.ϕ))
            try
                ptdf[island_b, island] = get_isf(pf.X, pf.B, pf.DA, (f, t), c, island)
                flow = ptdf * (pf.Pᵢ .+ δP)
                catch DivideError
                    @warn "Islands forming due to contingency on line $(f)-$(t)-i_$c."
                    continue
                end
        end
        for (n,f,r) in zip(b_names, flow, linerates)
            if abs(f) > r + 0.0001
                print(string_line(cont.name, n, f, r))
            end
        end
    end
end

function print_contingency_overflow(opfm::OPFmodel, short_term_limit_multi::Float64 = 1.5)
    string_line(c, n, f, r) = @sprintf("%15s %15s  %6.2f  %6.2f\n", c, n, f, r)
    println("    Contingency          Branch    Flow  Rating")
    branches = get_sorted_branches(opfm.sys)
    nodes = get_sorted_nodes(opfm.sys)
    idx = get_nodes_idx(nodes)
    slack = find_slack(nodes)
    b_names = get_name.(branches)
    linerates = get_rate.(branches) * short_term_limit_multi
    P = get_net_Pᵢ(opfm.mod, opfm.sys, nodes, idx)
    for (c,cont) in enumerate(branches)
        try
            flow = get_ptdf(branches[1:end .!= c], length(nodes), idx, slack[1]) * P
            for (n,f,r) in zip(b_names[1:end .!= c], flow, linerates[1:end .!= c])
                if abs(f) > r + 0.0001
                    print(string_line(cont.name, n, f, r))
                end
            end
        catch SingularException
            continue
        end
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
