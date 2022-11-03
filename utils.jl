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
include("utils.jl")
include("N-1_SCOPF.jl")
include("short_long_SCOPF.jl")

function add_system_data_to_json(;
        file_name="system_data.json",
        data_name=joinpath("matpower","IEEE_RTS.m"),
        # time_series=joinpath("forecasts","5bus_ts","timeseries_pointers_da.json"),
        DATA_DIR="data"
    )
    # if !isdir(DATA_DIR)
    #     download(PowerSystems.UtilsData.TestData, folder = DATA_DIR)
    # end
    system_data = PSY.System(joinpath(DATA_DIR, data_name))
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

function comparison(range, fname)
    @info "----------- Start of simulation -----------"
    system = System(fname)
    result = Vector()
    pg = Vector()
    x = 1
    for i in range
        try
            for d in get_components(StaticLoad, system)
                set_active_power!(d, get_active_power(d)*i/x)
            end
            x = i
            opfm = scopf(PCSC, system, Gurobi.Optimizer, voll=voll, prob=prob, max_shed=0.5)
            solve_model!(opfm.mod)
            sl_model = sl_scopf(system, Gurobi.Optimizer, voll=voll, prob=prob, max_shed=0.5)
            push!(result, (i,
                solve_time(opfm.mod), objective_value(opfm.mod),
                solve_time(sl_model[2]), objective_value(sl_model[2])
            ))
            push!(pg, (i, value.(opfm.mod[:pg0]).data, value.(sl_model[1].mod[:pg0]).data))
        catch e
            @warn "Error when load is x$(i)"
            break
        end
    end
    return result, pg
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
