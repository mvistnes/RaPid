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
    system = System(joinpath("data", fname))
    p_names = [(:pg0, Generator), (:pf0, Branch), (:pfc, Branch), (:ls0, StaticLoad), (:lsc, StaticLoad), (:va0, Bus), (:vac, Bus)]
    c_names = [(:pgu, Generator), (:pgd, Generator), (:pfcc, Branch), (:lscc, StaticLoad), (:vacc, Bus)]
    result = Vector()
    pg = Vector()
    ls = zeros(length(demands(system))+length(renewables(system)))
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
                value.(opfm.mod[:pf0]).data ./ get_rate.(get_components(Branch, system)).|> abs |> mean,
                solve_time(sl_model[1].mod)+solve_time(sl_model[2]), objective_value(sl_model[2]), 
                value.(sl_model[1].mod[:pf0]).data ./ get_rate.(get_components(Branch, system)).|> abs |> mean
            ))
            push!(pg, (i, value.(opfm.mod[:pg0]).data, value.(sl_model[1].mod[:pg0]).data))
            ls .+= (value.(opfm.mod[:ls0]) .+ sum(value.(opfm.mod[:lsc]).data, dims=2) .+ sum(value.(opfm.mod[:lscc]).data, dims=2)).data
            # for name in p_names
            #     make_save_plot(opfm.mod, system, string("norm",sys_name,string(i)), name)
            #     make_save_plot(sl_model[1].mod, system, string("sl",sys_name,string(i)), name)
            # end
            # for name in c_names
            #     make_save_plot(opfm.mod, system, string("norm",sys_name,string(i)), name)
            #     make_save_plot(sl_model[2], system, string("sl",sys_name,string(i)), name)
            # end
        catch e
            @warn "$(e) when load is x$(i)"
            break
        end
    end
    return result, pg, ls
end
