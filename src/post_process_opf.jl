# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

# scatterplot(model, system, name, type) =
#     scatter(
#         [get_name.(get_components(type, system))],
#         [value.(model[name]).data],
#         dpi=100,
#         size=(600, 600),
#         label=false,
#         rotation=90,
#         title=name
#     )

# function make_save_plot(model, system, sys_name, name)
#     plt = scatterplot(model, system, name[1], name[2])
#     display(plt)
#     path = mkpath(joinpath("results", sys_name))
#     savefig(plt, joinpath(path, "$(name[1]).pdf"))
# end

# function scatter_all(model, system; sys_name="")
#     names = [
#         (:pg0, Generator), (:pgu, Generator), (:pgd, Generator), (:u, ThermalGen),
#         (:pf0, Branch), (:pfc, Branch), (:pfcc, Branch),
#         (:ls0, StaticLoad), (:lsc, StaticLoad), (:lscc, StaticLoad),
#         (:qg0, Generator), (:qgu, Generator), (:qgd, Generator),
#         (:qf0, Branch), (:qfc, Branch), (:qfcc, Branch),
#         (:va0, ACBus), (:vac, ACBus), (:vacc, ACBus),
#         (:cbc, ACBus), (:cbcc, ACBus)
#     ]
#     for name in names
#         try
#             plt = scatterplot(model, system, name[1], name[2])
#             display(plt)
#             path = mkpath(joinpath("results", sys_name))
#             savefig(plt, joinpath(path, "$(name[1]).pdf"))
#         catch e
#             print("No $(name[1]). ")
#         end
#     end
# end

function print_variabel(model::Model, list::AbstractVector, name::Symbol)
    for (i, g) in enumerate(list)
        @printf("%12s: %.3f\n", g.name, JuMP.value(model[name][i]))
    end
end

function print_active_power(opf::OPFsystem)
    idx = get_nodes_idx(opf.nodes)
    list = make_list(opf, idx, opf.nodes)
    xprint(x, val) = @printf("%7s:%6.2f, ", split(x.name, "_", limit=2)[end], val)
    sort_x!(list) = sort!(list, by=x -> x.name)
    get_injection(g::Union{ThermalGen,HydroGen}) = value(mod[:pg0][get_name(g)])
    get_injection(r::RenewableGen) = get_active_power(r) - value(mod[:ls0][get_name(r)])
    get_injection(d::StaticLoad) = -get_active_power(d) + value(mod[:ls0][get_name(d)])
    print("       Bus   Total  Injections...")
    tot = 0.0
    for (n, g, d, r) in zip(opf.nodes, list_gen, list_d, list_r)
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

function get_power_flow(opf::OPFsystem, mod::Model, sep::String=" ")
    try
        return value.(mod[:pf0])
    catch KeyError
        return calculate_line_flows(get_isf(opf.branches, opf.nodes), get_net_Pᵢ(opf, mod))
    end
end
function print_power_flow(opf::OPFsystem, mod::Model, sep::String=" ")
    print_power_flow(get_name.(opf.branches), get_power_flow(opf, mod, sep),
        get_rate.(opf.branches), sep=sep)
end
function print_power_flow(names::AbstractVector{String}, flow::AbstractVector,
    rate::AbstractVector{<:Real}; sep::String=" ", risky_flow=0.9)
    println("         Branch     Flow  Rating")
    string_line(n, f, r, sep=" ") = @sprintf("%15s%s %7.3f%s %5.1f\n", n, sep, f, sep, r)
    for (n, f, r) in zip(names, flow, rate)
        if abs(f) > r + 0.0001
            printstyled(string_line(n, f, r, sep), color=:red)
        elseif abs(f) > r * risky_flow
            printstyled(string_line(n, f, r, sep), color=:yellow)
        else
            print(string_line(n, f, r, sep))
        end
    end
end

function get_contingency_power_flow(opf::OPFsystem, mod::Model)
    idx = get_nodes_idx(opf.nodes)
    slack = find_slack(opf.nodes)
    θ = SCOPF.get_sorted_angles(mod)
    P = calc_Pᵢ(calc_B(opf.branches, idx), θ)
    flow = Vector{Vector{Float64}}(undef, length(opf.contingencies))
    for c in 1:length(opf.contingencies)
        flow[c] = get_isf(opf.branches[1:end.!=c], opf.nodes, idx, slack[1]) * P
    end
    return flow
end

function print_contingency_power_flow(opf::OPFsystem, b_names::Vector{String}, flow::Vector{Vector{Float64}}, linerates::Vector{<:Float64})
    for c in 1:length(opf.contingencies)
        if isempty(flow)
            println("Contingency ", opf.contingencies[c].name, " resulted in a disconnected reference bus.")
        else
            println("Contingency ", opf.contingencies[c].name)
            print_power_flow(b_names, flow[c], linerates)
        end
    end
end

function print_contingency_power_flow(opf::OPFsystem, mod::Model, pf::DCPowerFlow, Pc, Pcc, Pccx=Dict{Int,NTuple{2,Any}}(), short_term_limit::Float64=1.5, long_term_limit::Float64=1.0)
    b_names = get_name.(opf.branches)
    linerates = get_rate.(opf.branches)
    ptdf = get_contingency_ptdf(opf, pf)
    idx = get_nodes_idx(opf.nodes)
    list = make_list(opf, idx, opf.nodes)
    Pᵢ = calc_Pᵢ(pf)

    println("Base case")
    print_power_flow(b_names, pf.F, linerates)

    c_flow = [ptdf[i] * (Pᵢ .+ get_ΔPc(mod, list, zeros(length(opf.nodes)), Pc, c)) for (i, (c, cont)) in enumerate(get_branch_bus_idx(opf.branches, opf.contingencies, idx))]
    println("Short-term post-contingency")
    print_contingency_power_flow(opf, b_names, c_flow, linerates * short_term_limit)

    cc_flow = [ptdf[i] * (Pᵢ .+ get_ΔPcc(mod, opf.nodes, opf.dc_branches, list, zeros(length(opf.nodes)), Pcc, c)) for (i, (c, cont)) in enumerate(get_branch_bus_idx(opf.branches, opf.contingencies, idx))]
    println("Long-term post-contingency")
    print_contingency_power_flow(opf, b_names, cc_flow, linerates * long_term_limit)

    ccx_flow = [ptdf[i] * (Pᵢ .+ get_ΔPccx(mod, list, zeros(length(opf.nodes)), Pccx, c)) for (i, (c, cont)) in enumerate(get_branch_bus_idx(opf.branches, opf.contingencies, idx))]
    println("Long-term post-contingency corrective failed")
    print_contingency_power_flow(opf, b_names, ccx_flow, linerates * long_term_limit)
end

function print_contingency_power_flow(opf::OPFsystem, mod::Model, rate_limit_multi::Float64=1.5)
    b_names = get_name.(opf.branches)
    linerates = get_rate.(opf.branches)
    idx = get_nodes_idx(opf.nodes)
    println("Base case")
    print_power_flow(b_names, JuMP.value.(mod[:pf0]), linerates)
    println("Short-term post-contingency")
    for (c, cont) in get_branch_bus_idx(opf.branches, opf.contingencies, idx)
        println("Contingency ", opf.branches[c].name)
        print_power_flow(b_names, JuMP.value.(mod[:pfc][:, c]), linerates * rate_limit_multi)
    end
    println("Long-term post-contingency")
    for (c, cont) in get_branch_bus_idx(opf.branches, opf.contingencies, idx)
        println("Contingency ", opf.branches[c].name)
        print_power_flow(b_names, JuMP.value.(mod[:pfcc][:, c]), linerates)
    end
end

function print_contingency_overflow(opf::OPFsystem, mod::Model, pf::DCPowerFlow, Pc, Pcc, Pccx=Dict{Int,NTuple{2,Any}}(), short_term_limit::Float64=1.5, long_term_limit::Float64=1.0)
    string_line(c, n, f, r) = @sprintf("%15s %15s  %6.2f  %6.2f\n", c, n, f, r)
    print_string_line(c, n, f, r) = abs(f) > r + 0.0001 && print(string_line(c, n, f, r))
    print_all_string_line(c, n, f, r) = print_string_line.(c, n, f, r)
    b_names = get_name.(opf.branches)
    linerates = get_rate.(opf.branches)
    ptdf = get_contingency_ptdf(opf, pf)
    idx = get_nodes_idx(opf.nodes)
    list = make_list(opf, idx, opf.nodes)
    Pᵢ = calc_Pᵢ(pf)

    c_flow = [ptdf[i] * (Pᵢ .+ get_ΔPc(mod, list, zeros(length(opf.nodes)), Pc, c)) for (i, (c, cont)) in enumerate(get_branch_bus_idx(opf.branches, opf.contingencies, idx))]
    println("Short-term post-contingency")
    println("    Contingency          Branch    Flow  Rating")
    print_all_string_line.(get_name.(opf.contingencies), [b_names], c_flow, [linerates * short_term_limit])

    cc_flow = [ptdf[i] * (Pᵢ .+ get_ΔPcc(mod, opf.nodes, opf.dc_branches, list, zeros(length(opf.nodes)), Pcc, c)) for (i, (c, cont)) in enumerate(get_branch_bus_idx(opf.branches, opf.contingencies, idx))]
    println("Long-term post-contingency")
    println("    Contingency          Branch    Flow  Rating")
    print_all_string_line.(get_name.(opf.contingencies), [b_names], cc_flow, [linerates * long_term_limit])

    ccx_flow = [ptdf[i] * (Pᵢ .+ get_ΔPccx(mod, list, zeros(length(opf.nodes)), Pccx, c)) for (i, (c, cont)) in enumerate(get_branch_bus_idx(opf.branches, opf.contingencies, idx))]
    println("Long-term post-contingency corrective failed")
    println("    Contingency          Branch    Flow  Rating")
    print_all_string_line.(get_name.(opf.contingencies), [b_names], ccx_flow, [linerates * long_term_limit])
    return
end

function print_contingency_overflow(opf::OPFsystem, mod::Model, rate_limit_multi::Float64=1.5)
    string_line(c, n, f, r) = @sprintf("%15s %15s  %6.2f  %6.2f\n", c, n, f, r)
    print_string_line(c, n, f, r) = abs(f) > r + 0.0001 && print(string_line(c, n, f, r))
    println("    Contingency          Branch    Flow  Rating")
    idx = get_nodes_idx(opf.nodes)
    slack = find_slack(opf.nodes)
    b_names = get_name.(opf.branches)
    linerates = get_rate.(opf.branches)
    θ = SCOPF.get_sorted_angles(mod)
    P = calc_Pᵢ(calc_B(opf.branches, idx), θ)
    print_string_line(opf.branches[c].name, n, f, r)
end

function print_results(opf::OPFsystem, mod::Model)
    print_variabel(mod, opf.ctrl_generation, :pg0)
    print_variabel(mod, opf.dc_branches, :pfdc0)
    print_variabel(mod, opf.renewables, :pr0)
    print_variabel(mod, opf.demands, :ls0)
end

function print_generation_results(opf::OPFsystem, mod::Model)
    (pg_lim_min, pg_lim_max) = split_pair(get_active_power_limits.(opf.ctrl_generation))
    Pg = get_value(mod, :pg0)
    gen_bus = opf.ctrl_generation .|> get_bus .|> get_number
    println(" Gen   Bus        Active Power Limits \n" *
            "  #     #       Pmin       Pg       Pmax  \n" *
            "----  -----   --------  --------  --------")
    for (i, (n, min, pg, max)) in enumerate(zip(gen_bus, pg_lim_min, Pg, pg_lim_max))
        @printf("%4d  %4d  %8.2f  %8.2f  %8.2f\n", i, n, min, pg, max)
    end
    
    pr_lim = get_active_power.(opf.renewables)
end

function print_preventive_results(opf::OPFsystem, mod::Model)
    for (i_g, g) in enumerate(opf.ctrl_generation)
        @printf("%12s: %.3f\n", g.name, JuMP.value(mod[:pg0][i_g]))
        for (i, c) in enumerate(JuMP.value.(mod[:pgc][i_g, :]))
            if c > 0.0001
                @printf("           %4dc: pgc: %.3f\n", i, c)
            end
        end
    end
    for (i_g, g) in enumerate(opf.renewables)
        @printf("%12s: %.3f\n", g.name, JuMP.value(mod[:pr0][i_g]))
        for (i, c) in enumerate(JuMP.value.(mod[:prc][i_g, :]))
            if c > 0.0001
                @printf("           %4dc: prc: %.3f\n", i, c)
            end
        end
    end
    for (i_g, g) in enumerate(opf.demands)
        @printf("%12s: %.3f\n", g.name, JuMP.value(mod[:ls0][i_g]))
        for (i, c) in enumerate(JuMP.value.(mod[:lsc][i_g, :]))
            if c > 0.0001
                @printf("           %4dc: lsc: %.3f\n", i, c)
            end
        end
    end
end

function print_corrective_results(opf::OPFsystem, mod::Model)
    for (i_g, g) in enumerate(opf.ctrl_generation)
        @printf("%12s: %.3f (%.3f)\n", g.name, JuMP.value(mod[:pg0][i_g]), get_active_power_limits(g).max)
        for (i, c) in enumerate(JuMP.value.(mod[:pgc][i_g, :]))
            if c > 0.0001
                @printf("          c %12s: pgc: %.3f\n", opf.contingencies[i].name, c)
            end
        end
        for (i, c) in enumerate(zip(JuMP.value.(mod[:pgu][i_g, :]), JuMP.value.(mod[:pgd][i_g, :])))
            if c[1] > 0.0001
                @printf("          c %12s: pgu: %.3f\n", opf.contingencies[i].name, c[1])
            end
            if c[2] > 0.0001
                @printf("          c %12s: pgd: %.3f\n", opf.contingencies[i].name, c[2])
            end
        end
    end
    for (i_g, g) in enumerate(opf.dc_branches)
        @printf("%12s: %.3f (%.3f)\n", g.name, JuMP.value(mod[:pfdc0][i_g]), get_active_power_limits(g).max)
        for (i, c) in enumerate(JuMP.value.(mod[:pfdccc][i_g, :]))
            if c > 0.0001
                @printf("          c %12s: pfdccc: %.3f\n", opf.contingencies[i].name, c)
            end
        end
    end
    for (i_g, g) in enumerate(opf.renewables)
        @printf("%12s: %.3f (%.3f)\n", g.name, JuMP.value(mod[:pr0][i_g]), get_active_power_limits(g).max)
        for (i, c) in enumerate(JuMP.value.(mod[:prc][i_g, :]))
            if c > 0.0001
                @printf("          c %12s: prc: %.3f\n", opf.contingencies[i].name, c)
            end
        end
        for (i, c) in enumerate(JuMP.value.(mod[:prcc][i_g, :]))
            if c > 0.0001
                @printf("          c %12s: prcc: %.3f\n", opf.contingencies[i].name, c)
            end
        end
    end
    for (i_g, g) in enumerate(opf.demands)
        @printf("%12s: %.3f (%.3f)\n", g.name, JuMP.value(mod[:ls0][i_g]), get_active_power(g))
        for (i, c) in enumerate(JuMP.value.(mod[:lsc][i_g, :]))
            if c > 0.0001
                @printf("          c %12s: lsc: %.3f\n", opf.contingencies[i].name, c)
            end
        end
        for (i, c) in enumerate(JuMP.value.(mod[:lscc][i_g, :]))
            if c > 0.0001
                @printf("          c %12s: lscc: %.3f\n", opf.contingencies[i].name, c)
            end
        end
    end
end

function print_sorted_corrective_results(opf::OPFsystem, mod::Model)
    for (i_g, g) in enumerate(opf.ctrl_generation)
        @printf("%12s: %.3f\n", g.name, JuMP.value(mod[:pg0][i_g]))
    end
    for (i_g, g) in enumerate(opf.dc_branches)
        @printf("%12s: %.3f\n", g.name, JuMP.value(mod[:pfdc0][i_g]))
    end
    for (i_g, g) in enumerate(opf.renewables)
        @printf("%12s: %.3f\n", g.name, JuMP.value(mod[:pr0][i_g]))
    end
    for (i_g, g) in enumerate(opf.demands)
        @printf("%12s: %.3f\n", g.name, JuMP.value(mod[:ls0][i_g]))
    end

    for (i, c) in enumerate(opf.contingencies)
        for (i_g, g) in enumerate(opf.ctrl_generation)
            if JuMP.value(mod[:pgc][i_g, i]) > 0.0001
                @printf("%4dc: %12s:  pgc: %.3f\n", i, g.name, JuMP.value(mod[:pgc][i_g, i]))
            end
        end
        for (i_g, g) in enumerate(opf.renewables)
            if JuMP.value(mod[:prc][i_g, i]) > 0.0001
                @printf("%4dc: %12s:  prc: %.3f\n", i, g.name, JuMP.value(mod[:prc][i_g, i]))
            end
        end
        for (i_g, g) in enumerate(opf.demands)
            if JuMP.value(mod[:lsc][i_g, i]) > 0.0001
                @printf("%4dc: %12s:  lsc: %.3f\n", i, g.name, JuMP.value(mod[:lsc][i_g, i]))
            end
        end
    end

    try
        for (i, c) in enumerate(opf.contingencies)
            for (i_g, g) in enumerate(opf.ctrl_generation)
                if JuMP.value(mod[:pgu][i_g, i]) > 0.0001
                    @printf("%4dc: %12s:  pgu: %.3f\n", i, g.name, JuMP.value(mod[:pgu][i_g, i]))
                end
                if JuMP.value(mod[:pgd][i_g, i]) > 0.0001
                    @printf("%4dc: %12s:  pgd: %.3f\n", i, g.name, JuMP.value(mod[:pgd][i_g, i]))
                end
            end
            for (i_g, g) in enumerate(opf.dc_branches)
                if JuMP.value(mod[:pfdccc][i_g, i]) > 0.0001
                    @printf("%4dc: %12s: pfdccc: %.3f\n", i, g.name, JuMP.value(mod[:pfdccc][i_g, i]))
                end
            end
            for (i_g, g) in enumerate(opf.renewables)
                if JuMP.value(mod[:prcc][i_g, i]) > 0.0001
                    @printf("%4dc: %12s: prcc: %.3f\n", i, g.name, JuMP.value(mod[:prcc][i_g, i]))
                end
            end
            for (i_g, g) in enumerate(opf.demands)
                if JuMP.value(mod[:lscc][i_g, i]) > 0.0001
                    @printf("%4dc: %12s: lscc: %.3f\n", i, g.name, JuMP.value(mod[:lscc][i_g, i]))
                end
            end
        end
    catch
    end
end

function print_contingency_P(opf::OPFsystem, mod::Model, idx)
    println("Pc")
    for (c, cont) in enumerate(opf.contingencies)
        P = zeros(length(opf.nodes))
        for (i_g, g) in enumerate(opf.ctrl_generation)
            P[idx[g.bus.number]] -= JuMP.value(mod[:pgc][i_g, c])
        end
        for (i_g, g) in enumerate(opf.renewables)
            P[idx[g.bus.number]] -= JuMP.value(mod[:prc][i_g, c])
        end
        for (i_g, g) in enumerate(opf.demands)
            P[idx[g.bus.number]] += JuMP.value(mod[:lsc][i_g, c])
        end
        for (i, x) in enumerate(P)
            if abs(x) > 0.00001
                Printf.@printf "%s,%2d,%7.4f\n" cont.name i x
            end
        end
    end


    println("Pcc")
    for (c, cont) in enumerate(opf.contingencies)
        P = zeros(length(opf.nodes))
        for (i_g, g) in enumerate(opf.ctrl_generation)
            P[idx[g.bus.number]] += JuMP.value(mod[:pgu][i_g, c])
            P[idx[g.bus.number]] -= JuMP.value(mod[:pgd][i_g, c])
        end
        for (i_g, g) in enumerate(opf.dc_branches)
            P[idx[g.bus.number]] -= beta(g.bus, g) * JuMP.value(mod[:pfdccc][i_g, c])
        end
        for (i_g, g) in enumerate(opf.renewables)
            P[idx[g.bus.number]] -= JuMP.value(mod[:prcc][i_g, c])
        end
        for (i_g, g) in enumerate(opf.demands)
            P[idx[g.bus.number]] += JuMP.value(mod[:lscc][i_g, c])
        end
        for (i, x) in enumerate(P)
            if abs(x) > 0.00001
                Printf.@printf "%s,%2d,%7.4f\n" cont.name i x
            end
        end
    end
end

function print_contingency_P(opf::OPFsystem, Pc, Pcc, Pccx=Dict())
    idx = SCOPF.get_nodes_idx(opf.nodes)
    println("Pc")
    for (c, cont) in enumerate(opf.contingencies)
        P = zeros(length(opf.nodes))
        x = get(Pc, c, 0)
        if x != 0
            for (i_g, g) in enumerate(opf.ctrl_generation)
                P[idx[g.bus.number]] -= JuMP.value(x.pgc[i_g])
            end
            for (i_g, g) in enumerate(opf.renewables)
                P[idx[g.bus.number]] -= JuMP.value(x.prc[i_g])
            end
            for (i_g, g) in enumerate(opf.demands)
                P[idx[g.bus.number]] += JuMP.value(x.lsc[i_g])
            end
            for (i, val) in enumerate(P)
                if abs(val) > 0.00001
                    Printf.@printf "%s,%2d,%7.4f\n" cont.name i val
                end
            end
        end
    end

    println("Pcc")
    for (c, cont) in enumerate(opf.contingencies)
        P = zeros(length(opf.nodes))
        x = get(Pcc, c, 0)
        if x != 0
            for (i_g, g) in enumerate(opf.ctrl_generation)
                P[idx[g.bus.number]] += JuMP.value(x.pgu[i_g])
                P[idx[g.bus.number]] -= JuMP.value(x.pgd[i_g])
            end
            for (i_g, g) in enumerate(opf.dc_branches)
                P[idx[g.arc.from.number]] += JuMP.value(x.pfdccc[i_g])
                P[idx[g.arc.to.number]] -= JuMP.value(x.pfdccc[i_g])
            end
            for (i_g, g) in enumerate(opf.renewables)
                P[idx[g.bus.number]] -= JuMP.value(x.prcc[i_g])
            end
            for (i_g, g) in enumerate(opf.demands)
                P[idx[g.bus.number]] += JuMP.value(x.lscc[i_g])
            end
            for (i, val) in enumerate(P)
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
