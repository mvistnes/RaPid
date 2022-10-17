using PowerSystems
using JuMP

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
