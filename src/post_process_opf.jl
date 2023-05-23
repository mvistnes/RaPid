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

function print_variabel(model::Model, list::AbstractVector, name::Symbol)
    for (i,g) in enumerate(list)
        @printf("%12s: %.3f\n", g.name, JuMP.value(model[name][i]))
    end
end

function print_active_power(opfm::OPFmodel)
    idx = get_nodes_idx(opfm.nodes)
    list = make_list(opfm, idx, opfm.nodes)
    xprint(x, val) = @printf("%7s:%6.2f, ", split(x.name, "_", limit=2)[end], val)
    sort_x!(list) = sort!(list, by = x -> x.name)
    get_injection(g::Union{ThermalGen, HydroGen}) = value(opfm.mod[:pg0][get_name(g)])
    get_injection(r::RenewableGen) = get_active_power(r) - value(opfm.mod[:ls0][get_name(r)])
    get_injection(d::StaticLoad) =  -get_active_power(d) + value(opfm.mod[:ls0][get_name(d)])
    print("       Bus   Total  Injections...")
    tot = 0.0
    for (n, g, d, r) in zip(opfm.nodes, list_gen, list_d, list_r)
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

function print_power_flow(opfm::OPFmodel, sep::String = " ")
    try
        print_power_flow(get_name.(opfm.branches), value.(opfm.mod[:pf0]), 
            get_rate.(opfm.branches), sep)
    catch KeyError
        print_power_flow(get_name.(opfm.branches), 
            calculate_line_flows(get_isf(opfm.branches, opfm.nodes), get_net_Pᵢ(opfm)), 
            get_rate.(opfm.branches), sep)
    end
end
function print_power_flow(names::AbstractVector{String}, flow::AbstractVector, 
        rate::AbstractVector{<:Real}, sep::String = " ")
    println("         Branch    Flow  Rating")
    string_line(n, f, r, sep=" ") = @sprintf("%15s%s %6.2f%s %6.2f\n", n, sep, f, sep, r)
    for (n,f,r) in zip(names, flow, rate)
        if abs(f) > r + 0.0001
            printstyled(string_line(n, f, r, sep), color = :red)
        elseif abs(f) > r * 0.9
            printstyled(string_line(n, f, r, sep), color = :yellow)
        else
            print(string_line(n, f, r, sep))
        end
    end
end

function print_contingency_power_flow(opfm::OPFmodel, pf, Pcc, short_term_limit_multi::Float64 = 1.5)
    idx = get_nodes_idx(opfm.nodes)
    list = make_list(opfm, idx, opfm.nodes)
    b_names = get_name.(opfm.branches)
    linerates = get_rate.(opfm.branches)
    Pᵢ = calc_Pᵢ(pf)
    println("Base case")
    print_power_flow(b_names, pf.F, linerates)
    flow = Vector{Float64}()
    for (c,cont) in enumerate(opfm.branches)
        println("Contingency ", cont.name)
        (f, t) = get_bus_idx(cont, idx)
        P = Pᵢ .+ get_ΔPcc(opfm, length(opfm.nodes), list, Pcc, c)
        try
            flow = calculate_line_flows(pf, (f, t), c, P)
        catch DivideError
            island, island_b = handle_islands(pf.B, pf.DA, (f,t), c, pf.slack)
            ptdf = zeros(eltype(pf.ϕ), size(pf.ϕ))
            ptdf[island_b, island] = get_isf(pf.DA, pf.B, (f,t), c, pf.slack, island, island_b)
            flow = ptdf * P
        end
        if isempty(flow)
            println("Contingency ", cont.name, " resulted in a disconnected reference bus.")
        else
            print_power_flow(b_names, flow, linerates)
        end
    end
end

function print_contingency_power_flow(opfm::OPFmodel, short_term_limit_multi::Float64 = 1.5)
    idx = get_nodes_idx(opfm.nodes)
    b_names = get_name.(opfm.branches)
    linerates = get_rate.(opfm.branches)
    P = calc_Pᵢ(pf)
    println("Base case")
    print_power_flow(b_names, value.(opfm.mod[:pf0]).data, linerates)
    for (c,cont) in enumerate(opfm.branches)
        println("Contingency ", cont.name)
        flow = get_ptdf(opfm.branches[1:end .!= c], length(opfm.nodes), idx, slack[1]) * P
        print_power_flow(b_names[1:end .!= c], flow, linerates[1:end .!= c])
    end
end

function print_contingency_overflow(opfm::OPFmodel, pf, Pcc, short_term_limit_multi::Float64 = 1.5)
    string_line(c, n, f, r) = @sprintf("%15s %15s  %6.2f  %6.2f\n", c, n, f, r)
    println("    Contingency        Branch    Flow  Rating")
    idx = get_nodes_idx(opfm.nodes)
    list = make_list(opfm, idx, opfm.nodes)
    b_names = get_name.(opfm.branches)
    linerates = get_rate.(opfm.branches)
    Pᵢ = calc_Pᵢ(pf)
    flow = Vector{Float64}()
    for (c,cont) in enumerate(opfm.branches)
        (f, t) = get_bus_idx(cont, idx)
        P = Pᵢ .+ get_ΔPcc(opfm, length(opfm.nodes), list, Pcc, c)
        try
            flow = calculate_line_flows(pf, (f, t), c, P)
        catch DivideError
            island, island_b = handle_islands(pf.B, pf.DA, (f,t), c, pf.slack)
            ptdf = zeros(eltype(pf.ϕ), size(pf.ϕ))
            ptdf[island_b, island] = get_isf(pf.DA, pf.B, (f,t), c, pf.slack, island, island_b)
            flow = ptdf * P
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
    idx = get_nodes_idx(opfm.nodes)
    slack = find_slack(opfm.nodes)
    b_names = get_name.(opfm.branches)
    linerates = get_rate.(opfm.branches)
    θ = SCOPF.get_sorted_angles(opfm.mod)
    P = calc_Pᵢ(calc_B(opfm.branches, length(opfm.nodes), idx), θ)
    for (c,cont) in enumerate(opfm.branches)
        try
            flow = get_ptdf(opfm.branches[1:end .!= c], length(opfm.nodes), idx, slack[1]) * P
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

function print_results(opfm::OPFmodel)
    print_variabel(opfm.mod, opfm.ctrl_generation, :pg0)
    print_variabel(opfm.mod, opfm.dc_branches, :pfdc0)
    print_variabel(opfm.mod, opfm.renewables, :pr0)
    print_variabel(opfm.mod, opfm.demands, :ls0)
end

function print_corrective_results(opfm::OPFmodel)
    for (i_g,g) in enumerate(opfm.ctrl_generation)
        @printf("%12s: %.3f\n", g.name, JuMP.value(opfm.mod[:pg0][i_g]))
        for (i,c) in enumerate(JuMP.value.(opfm.mod[:pgc][i_g,:]))
            if c > 0.0001
                @printf("           %4dc: pgc: %.3f\n", i, c)
            end
        end
        for (i,c) in enumerate(zip(JuMP.value.(opfm.mod[:pgu][i_g,:]), JuMP.value.(opfm.mod[:pgd][i_g,:])))
            if c[1] > 0.0001
                @printf("           %4dc: pgu: %.3f\n", i, c[1])
            end
            if c[2] > 0.0001
                @printf("           %4dc: pgd: %.3f\n", i, c[2])
            end
        end
    end
    for (i_g,g) in enumerate(opfm.dc_branches)
        @printf("%12s: %.3f\n", g.name, JuMP.value(opfm.mod[:pfdc0][i_g]))
        for (i,c) in enumerate(JuMP.value.(opfm.mod[:pfdccc][i_g,:]))
            if c > 0.0001
                @printf("           %4dc: pfdccc: %.3f\n", i, c)
            end
        end
    end
    for (i_g,g) in enumerate(opfm.renewables)
        @printf("%12s: %.3f\n", g.name, JuMP.value(opfm.mod[:pr0][i_g]))
        for (i,c) in enumerate(JuMP.value.(opfm.mod[:prc][i_g,:]))
            if c > 0.0001
                @printf("           %4dc: prc: %.3f\n", i, c)
            end
        end
        for (i,c) in enumerate(JuMP.value.(opfm.mod[:prcc][i_g,:]))
            if c > 0.0001
                @printf("           %4dc: prcc: %.3f\n", i, c)
            end
        end
    end
    for (i_g,g) in enumerate(opfm.demands)
        @printf("%12s: %.3f\n", g.name, JuMP.value(opfm.mod[:ls0][i_g]))
        for (i,c) in enumerate(JuMP.value.(opfm.mod[:lsc][i_g,:]))
            if c > 0.0001
                @printf("           %4dc: lsc: %.3f\n", i, c)
            end
        end
        for (i,c) in enumerate(JuMP.value.(opfm.mod[:lscc][i_g,:]))
            if c > 0.0001
                @printf("           %4dc: lscc: %.3f\n", i, c)
            end
        end
    end
end

function print_sorted_corrective_results(opfm::OPFmodel)
    for (i_g,g) in enumerate(opfm.ctrl_generation)
        @printf("%12s: %.3f\n", g.name, JuMP.value(opfm.mod[:pg0][i_g]))
    end
    for (i_g,g) in enumerate(opfm.dc_branches)
        @printf("%12s: %.3f\n", g.name, JuMP.value(opfm.mod[:pfdc0][i_g]))
    end
    for (i_g,g) in enumerate(opfm.renewables)
        @printf("%12s: %.3f\n", g.name, JuMP.value(opfm.mod[:pr0][i_g]))
    end
    for (i_g,g) in enumerate(opfm.demands)
        @printf("%12s: %.3f\n", g.name, JuMP.value(opfm.mod[:ls0][i_g]))
    end

    for (i,c) in enumerate(opfm.contingencies)
        for (i_g,g) in enumerate(opfm.ctrl_generation)
            if JuMP.value(opfm.mod[:pgc][i_g,i]) > 0.0001
                @printf("%4dc: %12s:  pgc: %.3f\n", i, g.name, JuMP.value(opfm.mod[:pgc][i_g,i]))
            end
        end
        for (i_g,g) in enumerate(opfm.renewables)
            if JuMP.value(opfm.mod[:prc][i_g,i]) > 0.0001
                @printf("%4dc: %12s:  prc: %.3f\n", i, g.name, JuMP.value(opfm.mod[:prc][i_g,i]))
            end
        end
        for (i_g,g) in enumerate(opfm.demands)
            if JuMP.value(opfm.mod[:lsc][i_g,i]) > 0.0001
                @printf("%4dc: %12s:  lsc: %.3f\n", i, g.name, JuMP.value(opfm.mod[:lsc][i_g,i]))
            end
        end
    end

    for (i,c) in enumerate(opfm.contingencies)
        for (i_g,g) in enumerate(opfm.ctrl_generation)
            if JuMP.value(opfm.mod[:pgu][i_g,i]) > 0.0001
                @printf("%4dc: %12s:  pgu: %.3f\n", i, g.name, JuMP.value(opfm.mod[:pgu][i_g,i]))
            end
            if JuMP.value(opfm.mod[:pgd][i_g,i]) > 0.0001
                @printf("%4dc: %12s:  pgd: %.3f\n", i, g.name, JuMP.value(opfm.mod[:pgd][i_g,i]))
            end
        end
        for (i_g,g) in enumerate(opfm.dc_branches)
            if JuMP.value(opfm.mod[:pfdccc][i_g,i]) > 0.0001
                @printf("%4dc: %12s: pfdccc: %.3f\n", i, g.name, JuMP.value(opfm.mod[:pfdccc][i_g,i]))
            end
        end
        for (i_g,g) in enumerate(opfm.renewables)
            if JuMP.value(opfm.mod[:prcc][i_g,i]) > 0.0001
                @printf("%4dc: %12s: prcc: %.3f\n", i, g.name, JuMP.value(opfm.mod[:prcc][i_g,i]))
            end
        end
        for (i_g,g) in enumerate(opfm.demands)
            if JuMP.value(opfm.mod[:lscc][i_g,i]) > 0.0001
                @printf("%4dc: %12s: lscc: %.3f\n", i, g.name, JuMP.value(opfm.mod[:lscc][i_g,i]))
            end
        end
    end
end

function print_contingency_P(opfm, idx)
    println("Pc")
    for (c,cont) in enumerate(opfm.contingencies)
        P = zeros(length(opfm.nodes))
        for (i_g,g) in enumerate(opfm.ctrl_generation)
            if JuMP.value(opfm.mod[:pgc][i_g,c]) > 0.0001
                P[idx[g.bus.number]] -= JuMP.value(opfm.mod[:pgc][i_g,c])
            end
        end
        for (i_g,g) in enumerate(opfm.renewables)
            if JuMP.value(opfm.mod[:prc][i_g,c]) > 0.0001
                P[idx[g.bus.number]] -= JuMP.value(opfm.mod[:prc][i_g,c])
            end
        end
        for (i_g,g) in enumerate(opfm.demands)
            if JuMP.value(opfm.mod[:lsc][i_g,c]) > 0.0001
                P[idx[g.bus.number]] += JuMP.value(opfm.mod[:lsc][i_g,c])
            end
        end
        for (i,x) in enumerate(P)
            if abs(x) > 0.00001
                Printf.@printf "%s,%2d,%5.2f\n" cont i x
            end
        end
    end

    
    println("Pcc")
    for (c,cont) in enumerate(opfm.contingencies)
        P = zeros(length(opfm.nodes))
        for (i_g,g) in enumerate(opfm.ctrl_generation)
            if JuMP.value(opfm.mod[:pgu][i_g,c]) > 0.0001
                P[idx[g.bus.number]] += JuMP.value(opfm.mod[:pgu][i_g,c])
            end
            if JuMP.value(opfm.mod[:pgd][i_g,c]) > 0.0001
                P[idx[g.bus.number]] -= JuMP.value(opfm.mod[:pgd][i_g,c])
            end
        end
        for (i_g,g) in enumerate(opfm.dc_branches)
            if JuMP.value(opfm.mod[:pfdccc][i_g,c]) > 0.0001
                P[idx[g.bus.number]] -= beta(g.bus, g) * JuMP.value(opfm.mod[:pfdccc][i_g,c])
            end
        end
        for (i_g,g) in enumerate(opfm.renewables)
            if JuMP.value(opfm.mod[:prcc][i_g,c]) > 0.0001
                P[idx[g.bus.number]] -= JuMP.value(opfm.mod[:prcc][i_g,c])
            end
        end
        for (i_g,g) in enumerate(opfm.demands)
            if JuMP.value(opfm.mod[:lscc][i_g,c]) > 0.0001
                P[idx[g.bus.number]] += JuMP.value(opfm.mod[:lscc][i_g,c])
            end
        end
        for (i,x) in enumerate(P)
            if abs(x) > 0.00001
                Printf.@printf "%s,%2d,%5.2f\n" cont i x
            end
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
