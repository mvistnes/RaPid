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

function add_to!(vec::AbstractVector, varref::AbstractVector{VariableRef}, idx::Dict, 
    components::AbstractVector, func::Function
)
    for (i, x) in enumerate(components)
        ix = idx[x.bus.number]
        vec[ix] = func(vec[ix], JuMP.value(varref[i]))
    end
end

function print_variabels(components::AbstractVector, varref::AbstractVector{VariableRef})
    for (i, x) in enumerate(components)
        @printf(" %.3f, %s\n", JuMP.value(varref[i]), x.name)
    end
end

function print_corrective_values(components::AbstractVector{<:Component}, varref::AbstractVector{VariableRef}, 
    symb::Symbol
)
    for (i, x) in enumerate(components)
        if JuMP.value(varref[i]) > 0.0001
            @printf("   %4s: %.3f, %s\n", string(symb), JuMP.value(varref[i]), x.name)
        end
    end
end

function print_results(opf::OPFsystem, m::Model)
    print_variabels(opf.generation, m[:pg0])
    print_variabels(opf.dc_branches, m[:pfdc0])
    print_variabels(opf.demands, m[:ls0])
end

function print_contingency_results(opf::OPFsystem, Pc::Dict{<:Integer, ExprC}, 
    c::Integer
)
    x = get(Pc, c, 0)
    if x != 0
        println("Contingency ", c)
        print_corrective_values(opf.generation, x.pgu, :pguc)
        print_corrective_values(opf.generation, x.pgd, :pgdc)
        print_corrective_values(opf.demands, x.lsc, :lsc)
    end
end

function print_contingency_results(opf::OPFsystem, Pcc::Dict{<:Integer, ExprCC}, 
    c::Integer
)
    x = get(Pcc, c, 0)
    if x != 0
        println("Contingency ", c)
        print_corrective_values(opf.generation, x.pgu, :pgucc)
        print_corrective_values(opf.generation, x.pgd, :pgdcc)
        print_corrective_values(opf.dc_branches, x.pfdccc, :pfdccc)
        print_corrective_values(opf.demands, x.lscc, :lscc)
    end
end

function print_contingency_results(opf::OPFsystem, Pccx::Dict{<:Integer, ExprCCX}, 
    c::Integer
)
    x = get(Pccx, c, 0)
    if x != 0
        println("Contingency ", c)
        print_corrective_values(opf.generation, x.pgdx, :pgdx)
        print_corrective_values(opf.demands, x.lsccx, :lsccx)
    end
end

function print_contingency_results(case::Case)
    print_contingency_results.([case.opf], [case.Pc], keys(case.Pc))
    print_contingency_results.([case.opf], [case.Pcc], keys(case.Pcc))
    print_contingency_results.([case.opf], [case.Pcc], keys(case.Pccx))
    return
end

function print_contingency_P(case::Case)
    opf = case.opf
    println("Pc; cont, node, val")
    for (c, cont) in enumerate(opf.contingencies)
        P = zeros(length(opf.nodes))
        x = get(case.Pc, c, 0)
        if x != 0
            add_to!(P, x.pgu, opf.idx, opf.generation, +)
            add_to!(P, x.pgd, opf.idx, opf.generation, -)
            add_to!(P, x.lsc, opf.idx, opf.demands, +)
            for (i, val) in enumerate(P)
                if abs(val) > 0.00001
                    Printf.@printf "%s,%2d,%7.4f\n" cont.name i val
                end
            end
        end
    end

    println("Pcc; cont, node, val")
    for (c, cont) in enumerate(opf.contingencies)
        P = zeros(length(opf.nodes))
        x = get(case.Pcc, c, 0)
        if x != 0
            add_to!(P, x.pgu, opf.idx, opf.generation, +)
            add_to!(P, x.pgd, opf.idx, opf.generation, -)
            add_to!(P, x.lscc, opf.idx, opf.demands, +)
            for (i_g, g) in enumerate(opf.dc_branches)
                P[opf.idx[g.arc.from.number]] += JuMP.value(x.pfdccc[i_g])
                P[opf.idx[g.arc.to.number]] -= JuMP.value(x.pfdccc[i_g])
            end
            for (i, val) in enumerate(P)
                if abs(val) > 0.00001
                    Printf.@printf "%s,%2d,%7.4f\n" cont.name i val
                end
            end
        end
    end
end

function print_active_power(opf::OPFsystem, m::Model)
    list = make_list(opf, opf.nodes)
    xprint(x, val) = @printf("%6.3f:%s,\t", val, x.name)
    sort_x!(list) = sort!(list, by=x -> x.name)
    print("       Bus   Total  Injections...")
    tot = 0.0
    for x in list
        val_g = value.(m[:pg0][x.generation])
        val_d = -value.(m[:pd][x.demands]) + value.(m[:ls0][x.demands])
        b_tot = sum(val_g, init=0.0) + sum(val_d, init=0.0) 
        tot += b_tot
        @printf "\n%10s  %6.3f  " x.node.name b_tot
        xprint.(opf.generation[x.generation], val_g)
        xprint.(opf.demands[x.demands], val_d)
    end
    @printf("\n       Sum %10.6f\n", tot)
end

function get_power_flow(opf::OPFsystem, pf::DCPowerFlow, oplim::Oplimits, m::Model; subset::AbstractVector{<:Integer}=Int64[])
    br_names = [Tuple(opf.mbx[i,:].nzind) for i in axes(opf.mbx,1)]
    flow = calculate_line_flows(pf.ϕ, get_net_Pᵢ(m))
    rate = oplim.branch_rating
    if !isempty(subset)
        br_names = br_names[subset]
        flow = flow[subset]
        rate = rate[subset]
    end
    return Dict("name" => br_names, "flow" => flow, "rate" => rate)
end

function print_power_flow(case::Case; subset::AbstractVector{<:Integer}=Int64[],
    all=true, sep::String=" ", risky_flow=0.9, atol=1e-14
)
    vals = get_power_flow(case.opf, case.pf, case.oplim, case.model, subset=subset)
    if !isempty(subset)
        all=true
    end
    print_power_flow(
        vals["name"], vals["flow"], vals["rate"]; all=all, sep=sep, risky_flow=risky_flow, atol=atol
    )
end

function print_power_flow(names::AbstractVector, flow::AbstractVector,
    rate::AbstractVector{<:Real}; all=true, sep::String=" ", risky_flow=0.9, atol=1e-14
)
    all && println("     Branch     Flow  Rating")
    string_line(n, f, r, sep) = @sprintf("%5d-%5d%s %7.3f%s %7.3f\n", n[1], n[2], sep, f, sep, r)
    for (n, f, r) in zip(names, flow, rate)
        if abs(f) > r + atol
            printstyled(string_line(n, f, r, sep), color=:red)
        elseif abs(f) > r * risky_flow
            printstyled(string_line(n, f, r, sep), color=:yellow)
        elseif all
            print(string_line(n, f, r, sep))
        end
    end
end

function get_contingency_power_flow(opf::OPFsystem, m::Model, pf::DCPowerFlow, 
    Pᵢ::AbstractVector{<:Real}, Pc::Dict, name::String; subset::AbstractVector{<:Integer}=Int64[]
)
    islands = Vector{Vector{Int}}()
    island = 0
    island_b = Int[]
    ΔP = similar(Pᵢ)
    vals = Dict{String, Vector}()
    for (i, c_obj) in enumerate(opf.contingencies)
        cont = typesort_component(c_obj, opf)
        if is_islanded(pf, cont[2], cont[1])
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont[2], cont[1], pf.slack)
            inodes = islands[island]
        else
            empty!(island_b)
            inodes = Int[]
        end
        flow = calculate_contingency_line_flows!(ΔP, Pc, opf, m, pf, Pᵢ, cont, i, c_obj, inodes, island_b)
        vals[name * c_obj.name] = isempty(subset) ? flow : flow[subset] 
    end
    return vals
end

function print_contingency_power_flow(opf::OPFsystem, m::Model, pf::DCPowerFlow,
    Pᵢ::AbstractVector{<:Real}, b_names::Vector{String}, linerates::Vector{<:Float64}, Pc::Dict; 
    subset::AbstractVector{<:Integer}=Int64[], all=true, sep::String=" ", risky_flow=0.9, atol=1e-14
)
    islands = Vector{Vector{Int}}()
    island = 0
    island_b = Int[]
    ΔP = similar(Pᵢ)
    !all && println("         Branch     Flow  Rating")
    for (i, c_obj) in enumerate(opf.contingencies)
        cont = typesort_component(c_obj, opf)
        if is_islanded(pf, cont[2], cont[1])
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont[2], cont[1], pf.slack)
            inodes = islands[island]
        else
            empty!(islands)
            island = 0
            empty!(island_b)
            inodes = Int[]
        end
        flow = calculate_contingency_line_flows!(ΔP, Pc, opf, m, pf, Pᵢ, cont, i, c_obj, inodes, island_b)
        if all || any(abs.(flow) .> linerates * risky_flow)
            println("Contingency ", c_obj.name)
            if isempty(subset)
                print_power_flow(b_names, flow, linerates; all=all, sep=sep, risky_flow=risky_flow, atol=atol)
            else
                print_power_flow(b_names[subset], flow[subset], linerates[subset]; all=all, sep=sep, risky_flow=risky_flow, atol=atol)
            end
        end
    end
end

function get_contingency_power_flow(case::Case; 
    subset::AbstractVector{<:Integer}=Int64[], all=true, sep::String=" ", risky_flow=0.9, atol=1e-14
)
    b_names = PowerSystems.get_name.(case.opf.branches)
    Pᵢ = get_value(case.model, :p0)
    brst = case.oplim.branch_rating * case.oplim.short_term_multi
    brlt = case.oplim.branch_rating * case.oplim.long_term_multi
    vals = Dict("name" => b_names, "st_rate" => brst, "lt_rate" => brlt)

    if !isempty(case.Pc) 
        merge!(vals, calc_contingency_power_flow(case.opf, case.model, case.pf, Pᵢ, case.Pc, "short_"; subset=subset))
    end

    if !isempty(case.Pcc)
        merge!(vals, calc_contingency_power_flow(case.opf, case.model, case.pf, Pᵢ, case.Pcc, "long_"; subset=subset))
    end

    if !isempty(case.Pccx)
        merge!(vals, calc_contingency_power_flow(case.opf, case.model, case.pf, Pᵢ, case.Pccx, "long_x_"; subset=subset))
    end
    return vals
end

function print_contingency_power_flow(case::Case; 
    subset::AbstractVector{<:Integer}=Int64[], all=true, sep::String=" ", risky_flow=0.9, atol=1e-14
)
    b_names = PowerSystems.get_name.(case.opf.branches)
    Pᵢ = get_value(case.model, :p0)
    brst = case.oplim.branch_rating * case.oplim.short_term_multi
    brlt = case.oplim.branch_rating * case.oplim.long_term_multi

    if !isempty(case.Pc) 
        println("Short-term post-contingency")
        print_contingency_power_flow(case.opf, case.model, case.pf, Pᵢ, b_names, brst, case.Pc; 
                                     subset=subset, all=all, sep=sep, risky_flow=risky_flow, atol=atol)
    end

    if !isempty(case.Pcc)
        println("Long-term post-contingency")
        print_contingency_power_flow(case.opf, case.model, case.pf, Pᵢ, b_names, brlt, case.Pcc; 
                                     subset=subset, all=all, sep=sep, risky_flow=risky_flow, atol=atol)
    end

    if !isempty(case.Pccx)
        println("Long-term post-contingency corrective failed")
        print_contingency_power_flow(case.opf, case.model, case.pf, Pᵢ, b_names, brlt, case.Pccx; 
                                     subset=subset, all=all, sep=sep, risky_flow=risky_flow, atol=atol)
    end
    return
end

function get_generation_results(case::Case)
    Pg = get_value(m, :pg0)
    gen_bus = opf.generation .|> get_bus .|> get_number
    return Dict("Bus_num" => gen_bus, "P" => Pg, "Min_P" => case.oplim.pg_lim_min, "Max_P" => case.oplim.pg_lim_max)
end

function print_generation_results(case::Case)
    vals = get_generation_results(case)
    println(" Gen   Bus        Active Power Limits \n" *
            "  #     #       Pmin       Pg       Pmax  \n" *
            "----  -----   --------  --------  --------")
    for (i, (n, min, pg, max)) in enumerate(zip(vals["Bus_num"], vals["Min_P"], vals["P"], vals["Max_P"]))
        @printf("%4d  %4d  %8.2f %s%8.2f %s%8.2f\n", i, n, min, min ≈ pg ? "=" : " ", pg, max ≈ pg ? "=" : " ", max)
    end
end


function print_benders_results(case::Case, system::System, lim::Real=1e-14)
    function print_c(itr, symb::String, x::Int)
        for (i, c_obj) in enumerate(case.opf.contingencies)
            c = get(itr, i, 0)
            if c != 0 && JuMP.value(getfield(c, Symbol(symb))[x]) > lim
                @printf("          c %5d: %s: %.3f\n", i, symb, JuMP.value(getfield(c, Symbol(symb))[x]))
            end
        end
    end

    for (x, g) in enumerate(sort_components!(get_generation(system)))
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(case.model[:pg0][x]), case.oplim.pg_lim_max[x])
        print_c(case.Pc, "pgu", x)
        print_c(case.Pc, "pgd", x)
        print_c(case.Pcc, "pgu", x)
        print_c(case.Pcc, "pgd", x)
        print_c(case.Pccx, "pgdx", x)
    end
    for (x, g) in enumerate(sort_components!(get_dc_branches(system)))
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(case.model[:pfdc0][x]), case.oplim.dc_lim_max[x])
        print_c(case.Pcc, "pfdccc", x)
    end
    for (x, g) in enumerate(sort_components!(get_demands(system)))
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(case.model[:ls0][x]), case.oplim.pd_lim[x])
        print_c(case.Pc, "lsc", x)
        print_c(case.Pcc, "lscc", x)
        print_c(case.Pccx, "lsccx", x)
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
