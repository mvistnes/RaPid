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
# include("N-1_SCOPF.jl")
# include("short_long_SCOPF.jl")

function add_system_data_to_json(;
        file_name="system_data.json",
        data_name=joinpath("matpower","IEEE_RTS.m"),
        # time_series=joinpath("forecasts","5bus_ts","timeseries_pointers_da.json"),
        DATA_DIR="data"
    )
    # if !isdir(DATA_DIR)
    #     download(PowerSystems.UtilsData.TestData, folder = DATA_DIR)
    # end
    system_data = System(joinpath(DATA_DIR, data_name))
    # add_time_series!(system_data, joinpath(DATA_DIR,time_series))
    to_json(system_data, file_name, force=true)
end

function hot_start(model::JuMP.Model)
    x = JuMP.all_variables(model)
    x_solution = JuMP.value.(x)
    JuMP.set_start_value.(x, x_solution)
    JuMP.optimize!(model)
    return model
end

function comparison(range, fname, optimizer; 
        sys_name = "",
        voll = nothing, 
        contingencies = nothing, 
        prob = nothing,
        time_limit_sec = 600,
        unit_commit::Bool=false,
        max_shed = 0.1,
        max_curtail::Float64 = 1.0,
        ratio = 0.5, 
        short_term_limit_multi = 1.5,
        ramp_minutes::Int64 = 10,
        repair_time::Float64 = 1.0)

    @info "----------- Start of simulation -----------"
    system = System(fname)
    p_names = [(:pg0, Generator), (:pf0, Branch), (:pfc, Branch), (:ls0, StaticLoad), (:lsc, StaticLoad), (:va0, Bus), (:vac, Bus)]
    c_names = [(:pgu, Generator), (:pgd, Generator), (:pfcc, Branch), (:lscc, StaticLoad), (:vacc, Bus)]
    result = Vector()
    pg = Vector()
    x = 1
    for i in range
        try
            for d in get_components(StaticLoad, system)
                set_active_power!(d, get_active_power(d)*i/x)
            end
            x = i
            opfm = scopf(PCSC, system, optimizer, voll=voll, prob=prob, contingencies=contingencies, time_limit_sec=time_limit_sec,
                unit_commit=unit_commit, max_shed=max_shed, max_curtail=max_curtail, ratio=ratio, short_term_limit_multi=short_term_limit_multi, 
                ramp_minutes=ramp_minutes, repair_time=repair_time)
            solve_model!(opfm.mod)
            sl_model = sl_scopf(system, optimizer, voll=voll, prob=prob, contingencies=contingencies, time_limit_sec=time_limit_sec,
                unit_commit=unit_commit, max_shed=max_shed, max_curtail=max_curtail, ratio=ratio, short_term_limit_multi=short_term_limit_multi, 
                ramp_minutes=ramp_minutes, repair_time=repair_time)
            push!(result, (i,
                solve_time(opfm.mod), objective_value(opfm.mod),
                solve_time(sl_model[1].mod)+solve_time(sl_model[2]), objective_value(sl_model[2])
            ))
            push!(pg, (i, value.(opfm.mod[:pg0]).data, value.(sl_model[1].mod[:pg0]).data))
            for name in p_names
                make_save_plot(opfm.mod, system, string("norm",sys_name,string(i)), name)
                make_save_plot(sl_model[1].mod, system, string("sl",sys_name,string(i)), name)
            end
            for name in c_names
                make_save_plot(opfm.mod, system, string("norm",sys_name,string(i)), name)
                make_save_plot(sl_model[2], system, string("sl",sys_name,string(i)), name)
            end
        catch e
            @warn "$(e) when load is x$(i)"
            break
        end
    end
    return result, pg
end

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
    for (i,x) in zip(get_name.(get_components(type, system)), value.(model[name]).data)
        @printf("%s: %.3f\n", i, x)
    end
end

# # corrective control failure probability
# phi(p, n) = sum((-1)^k * p^k * binomial(n,k) for k in 1:n)

# # severity function
# @expression(opf_m, severity, sum(voll[d] * lse[d] for d in get_name.(demands)))


# for l in get_name.(get_components(ACBranch, ieee_rts))
#     for g in get_name.(get_components(Generator, ieee_rts))
#         if value.(results[:pgu])[g,l] > 0.00001 
#             @printf("%s: \t%s \t= %.5f \n", l, g, value.(results[:pgu])[g,l])
#         end
#     end
# end
        
# add_system_data_to_json()
# system_data = System("system_data.json")
# results = opf_model(system_data, Ipopt.Optimizer)
# value.(results[:pl])
