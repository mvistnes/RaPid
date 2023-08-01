# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

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

function get_power_flow(opfm::OPFmodel, sep::String = " ")
    try
        return value.(opfm.mod[:pf0])
    catch KeyError
        return calculate_line_flows(get_isf(opfm.branches, opfm.nodes), get_net_Pᵢ(opfm))
    end
end
function print_power_flow(opfm::OPFmodel, sep::String = " ")
    print_power_flow(get_name.(opfm.branches), get_power_flow(opfm, sep), 
            get_rate.(opfm.branches), sep=sep)
end
function print_power_flow(names::AbstractVector{String}, flow::AbstractVector, 
        rate::AbstractVector{<:Real}; sep::String = " ", risky_flow = 0.9)
    println("         Branch     Flow  Rating")
    string_line(n, f, r, sep=" ") = @sprintf("%15s%s %7.3f%s %5.1f\n", n, sep, f, sep, r)
    for (n,f,r) in zip(names, flow, rate)
        if abs(f) > r + 0.0001
            printstyled(string_line(n, f, r, sep), color = :red)
        elseif abs(f) > r * risky_flow
            printstyled(string_line(n, f, r, sep), color = :yellow)
        else
            print(string_line(n, f, r, sep))
        end
    end
end

function get_contingency_power_flow(opfm::OPFmodel)
    idx = get_nodes_idx(opfm.nodes)
    slack = find_slack(opfm.nodes)
    θ = SCOPF.get_sorted_angles(opfm.mod)
    P = calc_Pᵢ(calc_B(opfm.branches, length(opfm.nodes), idx), θ)
    flow = Vector{Vector{Float64}}(undef, length(opfm.contingencies))
    for c in 1:length(opfm.contingencies)
        flow[c] = get_isf(opfm.branches[1:end .!= c], opfm.nodes, idx, slack[1]) * P
    end
    return flow
end

function print_contingency_power_flow(opfm::OPFmodel, b_names::Vector{String}, flow::Vector{Vector{Float64}}, linerates::Vector{<:Float64})
    for c in 1:length(opfm.contingencies)
        if isempty(flow)
            println("Contingency ", opfm.contingencies[c].name, " resulted in a disconnected reference bus.")
        else
            println("Contingency ", opfm.contingencies[c].name)
            print_power_flow(b_names, flow[c], linerates)
        end
    end
end

function print_contingency_power_flow(opfm::OPFmodel, pf::DCPowerFlow, Pc, Pcc, Pccx = Dict{Int, NTuple{2, Any}}(), short_term_limit::Float64 = 1.5, long_term_limit::Float64 = 1.0)
    b_names = get_name.(opfm.branches)
    linerates = get_rate.(opfm.branches)
    ptdf = get_contingency_ptdf(opfm, pf)
    idx = get_nodes_idx(opfm.nodes)
    list = make_list(opfm, idx, opfm.nodes)
    Pᵢ = calc_Pᵢ(pf)

    println("Base case")
    print_power_flow(b_names, pf.F, linerates)

    c_flow = [ptdf[i] * (Pᵢ .+ get_ΔPc(list, zeros(length(opfm.nodes)), Pc, c)) for (i,(c,cont)) in enumerate(get_branch_bus_idx(opfm.branches, opfm.contingencies, idx))]
    println("Short-term post-contingency")
    print_contingency_power_flow(opfm, b_names, c_flow, linerates * short_term_limit)
    
    cc_flow = [ptdf[i] * (Pᵢ .+ get_ΔPcc(opfm, list, zeros(length(opfm.nodes)), Pcc, c)) for (i,(c,cont)) in enumerate(get_branch_bus_idx(opfm.branches, opfm.contingencies, idx))]
    println("Long-term post-contingency")
    print_contingency_power_flow(opfm, b_names, cc_flow, linerates * long_term_limit)
    
    ccx_flow = [ptdf[i] * (Pᵢ .+ get_ΔPccx(list, zeros(length(opfm.nodes)), Pccx, c)) for (i,(c,cont)) in enumerate(get_branch_bus_idx(opfm.branches, opfm.contingencies, idx))]
    println("Long-term post-contingency corrective failed")
    print_contingency_power_flow(opfm, b_names, ccx_flow, linerates * long_term_limit)
end

function print_contingency_power_flow(opfm::OPFmodel, rate_limit_multi::Float64 = 1.5)
    b_names = get_name.(opfm.branches)
    linerates = get_rate.(opfm.branches)
    idx = get_nodes_idx(opfm.nodes)
    println("Base case")
    print_power_flow(b_names, JuMP.value.(opfm.mod[:pf0]), linerates)
    println("Short-term post-contingency")
    for (c,cont) in get_branch_bus_idx(opfm.branches, opfm.contingencies, idx)
        println("Contingency ", opfm.branches[c].name)
        print_power_flow(b_names, JuMP.value.(opfm.mod[:pfc][:,c]), linerates * rate_limit_multi)
    end
    println("Long-term post-contingency")
    for (c,cont) in get_branch_bus_idx(opfm.branches, opfm.contingencies, idx)
        println("Contingency ", opfm.branches[c].name)
        print_power_flow(b_names, JuMP.value.(opfm.mod[:pfcc][:,c]), linerates)
    end
end

function print_contingency_overflow(opfm::OPFmodel, pf::DCPowerFlow, Pc, Pcc, Pccx = Dict{Int, NTuple{2, Any}}(), short_term_limit::Float64 = 1.5, long_term_limit::Float64 = 1.0)
    string_line(c, n, f, r) = @sprintf("%15s %15s  %6.2f  %6.2f\n", c, n, f, r)
    print_string_line(c, n, f, r) = abs(f) > r + 0.0001 && print(string_line(c, n, f, r))
    print_all_string_line(c, n, f, r) = print_string_line.(c, n, f, r)
    b_names = get_name.(opfm.branches)
    linerates = get_rate.(opfm.branches)
    ptdf = get_contingency_ptdf(opfm, pf)
    idx = get_nodes_idx(opfm.nodes)
    list = make_list(opfm, idx, opfm.nodes)
    Pᵢ = calc_Pᵢ(pf)

    c_flow = [ptdf[i] * (Pᵢ .+ get_ΔPc(list, zeros(length(opfm.nodes)), Pc, c)) for (i,(c,cont)) in enumerate(get_branch_bus_idx(opfm.branches, opfm.contingencies, idx))]
    println("Short-term post-contingency")
    println("    Contingency          Branch    Flow  Rating")
    print_all_string_line.(get_name.(opfm.contingencies), [b_names], c_flow, [linerates * short_term_limit])
    
    cc_flow = [ptdf[i] * (Pᵢ .+ get_ΔPcc(opfm, list, zeros(length(opfm.nodes)), Pcc, c)) for (i,(c,cont)) in enumerate(get_branch_bus_idx(opfm.branches, opfm.contingencies, idx))]
    println("Long-term post-contingency")
    println("    Contingency          Branch    Flow  Rating")
    print_all_string_line.(get_name.(opfm.contingencies), [b_names], cc_flow, [linerates * long_term_limit])
    
    ccx_flow = [ptdf[i] * (Pᵢ .+ get_ΔPccx(list, zeros(length(opfm.nodes)), Pccx, c)) for (i,(c,cont)) in enumerate(get_branch_bus_idx(opfm.branches, opfm.contingencies, idx))]
    println("Long-term post-contingency corrective failed")
    println("    Contingency          Branch    Flow  Rating")
    print_all_string_line.(get_name.(opfm.contingencies), [b_names], ccx_flow, [linerates * long_term_limit])
    return
end

function print_contingency_overflow(opfm::OPFmodel, rate_limit_multi::Float64 = 1.5)
    string_line(c, n, f, r) = @sprintf("%15s %15s  %6.2f  %6.2f\n", c, n, f, r)
    print_string_line(c, n, f, r) = abs(f) > r + 0.0001 && print(string_line(c, n, f, r))
    println("    Contingency          Branch    Flow  Rating")
    idx = get_nodes_idx(opfm.nodes)
    slack = find_slack(opfm.nodes)
    b_names = get_name.(opfm.branches)
    linerates = get_rate.(opfm.branches)
    θ = SCOPF.get_sorted_angles(opfm.mod)
    P = calc_Pᵢ(calc_B(opfm.branches, length(opfm.nodes), idx), θ)
    print_string_line(opfm.branches[c].name, n, f, r)
end

function print_results(opfm::OPFmodel)
    print_variabel(opfm.mod, opfm.ctrl_generation, :pg0)
    print_variabel(opfm.mod, opfm.dc_branches, :pfdc0)
    print_variabel(opfm.mod, opfm.renewables, :pr0)
    print_variabel(opfm.mod, opfm.demands, :ls0)
end

function print_preventive_results(opfm::OPFmodel)
    for (i_g,g) in enumerate(opfm.ctrl_generation)
        @printf("%12s: %.3f\n", g.name, JuMP.value(opfm.mod[:pg0][i_g]))
        for (i,c) in enumerate(JuMP.value.(opfm.mod[:pgc][i_g,:]))
            if c > 0.0001
                @printf("           %4dc: pgc: %.3f\n", i, c)
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
    end
    for (i_g,g) in enumerate(opfm.demands)
        @printf("%12s: %.3f\n", g.name, JuMP.value(opfm.mod[:ls0][i_g]))
        for (i,c) in enumerate(JuMP.value.(opfm.mod[:lsc][i_g,:]))
            if c > 0.0001
                @printf("           %4dc: lsc: %.3f\n", i, c)
            end
        end
    end
end

function print_corrective_results(opfm::OPFmodel)
    for (i_g,g) in enumerate(opfm.ctrl_generation)
        @printf("%12s: %.3f (%.3f)\n", g.name, JuMP.value(opfm.mod[:pg0][i_g]), get_active_power_limits(g).max)
        for (i,c) in enumerate(JuMP.value.(opfm.mod[:pgc][i_g,:]))
            if c > 0.0001
                @printf("          c %12s: pgc: %.3f\n",opfm.contingencies[i].name, c)
            end
        end
        for (i,c) in enumerate(zip(JuMP.value.(opfm.mod[:pgu][i_g,:]), JuMP.value.(opfm.mod[:pgd][i_g,:])))
            if c[1] > 0.0001
                @printf("          c %12s: pgu: %.3f\n",opfm.contingencies[i].name, c[1])
            end
            if c[2] > 0.0001
                @printf("          c %12s: pgd: %.3f\n",opfm.contingencies[i].name, c[2])
            end
        end
    end
    for (i_g,g) in enumerate(opfm.dc_branches)
        @printf("%12s: %.3f (%.3f)\n", g.name, JuMP.value(opfm.mod[:pfdc0][i_g]), get_active_power_limits(g).max)
        for (i,c) in enumerate(JuMP.value.(opfm.mod[:pfdccc][i_g,:]))
            if c > 0.0001
                @printf("          c %12s: pfdccc: %.3f\n",opfm.contingencies[i].name, c)
            end
        end
    end
    for (i_g,g) in enumerate(opfm.renewables)
        @printf("%12s: %.3f (%.3f)\n", g.name, JuMP.value(opfm.mod[:pr0][i_g]), get_active_power_limits(g).max)
        for (i,c) in enumerate(JuMP.value.(opfm.mod[:prc][i_g,:]))
            if c > 0.0001
                @printf("          c %12s: prc: %.3f\n",opfm.contingencies[i].name, c)
            end
        end
        for (i,c) in enumerate(JuMP.value.(opfm.mod[:prcc][i_g,:]))
            if c > 0.0001
                @printf("          c %12s: prcc: %.3f\n",opfm.contingencies[i].name, c)
            end
        end
    end
    for (i_g,g) in enumerate(opfm.demands)
        @printf("%12s: %.3f (%.3f)\n", g.name, JuMP.value(opfm.mod[:ls0][i_g]), get_active_power(g))
        for (i,c) in enumerate(JuMP.value.(opfm.mod[:lsc][i_g,:]))
            if c > 0.0001
                @printf("          c %12s: lsc: %.3f\n",opfm.contingencies[i].name, c)
            end
        end
        for (i,c) in enumerate(JuMP.value.(opfm.mod[:lscc][i_g,:]))
            if c > 0.0001
                @printf("          c %12s: lscc: %.3f\n",opfm.contingencies[i].name, c)
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

    try
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
    catch
    end
end

function print_contingency_P(opfm, idx)
    println("Pc")
    for (c,cont) in enumerate(opfm.contingencies)
        P = zeros(length(opfm.nodes))
        for (i_g,g) in enumerate(opfm.ctrl_generation)
            P[idx[g.bus.number]] -= JuMP.value(opfm.mod[:pgc][i_g,c])
        end
        for (i_g,g) in enumerate(opfm.renewables)
            P[idx[g.bus.number]] -= JuMP.value(opfm.mod[:prc][i_g,c])
        end
        for (i_g,g) in enumerate(opfm.demands)
            P[idx[g.bus.number]] += JuMP.value(opfm.mod[:lsc][i_g,c])
        end
        for (i,x) in enumerate(P)
            if abs(x) > 0.00001
                Printf.@printf "%s,%2d,%7.4f\n" cont.name i x
            end
        end
    end

    
    println("Pcc")
    for (c,cont) in enumerate(opfm.contingencies)
        P = zeros(length(opfm.nodes))
        for (i_g,g) in enumerate(opfm.ctrl_generation)
            P[idx[g.bus.number]] += JuMP.value(opfm.mod[:pgu][i_g,c])
            P[idx[g.bus.number]] -= JuMP.value(opfm.mod[:pgd][i_g,c])
        end
        for (i_g,g) in enumerate(opfm.dc_branches)
            P[idx[g.bus.number]] -= beta(g.bus, g) * JuMP.value(opfm.mod[:pfdccc][i_g,c])
        end
        for (i_g,g) in enumerate(opfm.renewables)
            P[idx[g.bus.number]] -= JuMP.value(opfm.mod[:prcc][i_g,c])
        end
        for (i_g,g) in enumerate(opfm.demands)
            P[idx[g.bus.number]] += JuMP.value(opfm.mod[:lscc][i_g,c])
        end
        for (i,x) in enumerate(P)
            if abs(x) > 0.00001
                Printf.@printf "%s,%2d,%7.4f\n" cont.name i x
            end
        end
    end
end

function print_contingency_P(opfm, Pc, Pcc, Pccx, idx)
    println("Pc")
    for (c,cont) in enumerate(opfm.contingencies)
        P = zeros(length(opfm.nodes))
        x = get(Pc, c, 0)
        if x != 0
            for (i_g,g) in enumerate(opfm.ctrl_generation)
                P[idx[g.bus.number]] -= JuMP.value(x[1][i_g])
            end
            for (i_g,g) in enumerate(opfm.renewables)
                P[idx[g.bus.number]] -= JuMP.value(x[2][i_g])
            end
            for (i_g,g) in enumerate(opfm.demands)
                P[idx[g.bus.number]] += JuMP.value(x[3][i_g])
            end
            for (i,val) in enumerate(P)
                if abs(val) > 0.00001
                    Printf.@printf "%s,%2d,%7.4f\n" cont.name i val
                end
            end
        end
    end

    
    println("Pcc")
    for (c,cont) in enumerate(opfm.contingencies)
        P = zeros(length(opfm.nodes))
        x = get(Pcc, c, 0)
        if x != 0
            for (i_g,g) in enumerate(opfm.ctrl_generation)
                P[idx[g.bus.number]] += JuMP.value(x[1][i_g])
                P[idx[g.bus.number]] -= JuMP.value(x[2][i_g])
            end
            for (i_g,g) in enumerate(opfm.dc_branches)
                P[idx[g.bus.number]] -= beta(g.bus, g) * JuMP.value(x[3][i_g])
            end
            for (i_g,g) in enumerate(opfm.renewables)
                P[idx[g.bus.number]] -= JuMP.value(x[4][i_g])
            end
            for (i_g,g) in enumerate(opfm.demands)
                P[idx[g.bus.number]] += JuMP.value(x[5][i_g])
            end
            for (i,val) in enumerate(P)
                if abs(val) > 0.00001
                    Printf.@printf "%s,%2d,%7.4f\n" cont.name i val
                end
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
