# Example from PowerSystems docs https://nrel-siip.github.io/PowerSystems.jl/stable/modeler_guide/generatopf_modeling_with_JuMP/

using PowerSystems
const PSY = PowerSystems
using JuMP
using Ipopt
import Plots
using PowerSystemCaseBuilder

function add_system_data_to_json(;
        file_name="system_data.json",
        data_name=joinpath("matpower","case5_re_uc_pwl.m"),
        time_series=joinpath("forecasts","5bus_ts","timeseries_pointers_da.json"),
        DATA_DIR="data"
    )
    # if !isdir(DATA_DIR)
    #     download(PowerSystems.UtilsData.TestData, folder = DATA_DIR)
    # end
    system_data = System(joinpath(DATA_DIR, data_name))
    add_time_series!(system_data, joinpath(DATA_DIR,time_series))
    to_json(system_data, file_name, force=true)
end

function installed_capacity(system::System; technology::Type{T} = Generator) where T <: Generator
    capacity = 0.0
    for g in get_components(T, system)
        capacity += get_max_active_power(g)
    end
    return capacity
end

function opf_model(system::System, optimizer)
    opf_m = Model(optimizer)
    set_time_limit_sec(opf_m, 60.0)

    thermal_gens_names = get_name.(get_components(ThermalStandard, system))
    line_names = get_name.(get_components(Line, system))
    bus_names = get_name.(get_components(Bus, system))
    # ybus = Ybus(system)
    ptdf = PTDF(system)

    # make active power variables for the generators
    @variable(opf_m, pg[g in thermal_gens_names]) # TODO: add load curtailment variable
    # make voltage angle and active power variables for the nodes
    @variable(opf_m, angle[b in bus_names])
    @variable(opf_m, pb[b in bus_names])

    # restrict active power generation to min and max values
    for g in get_components(ThermalStandard, system)
        set_lower_bound(pg[get_name(g)], get_active_power_limits(g).min) 
        set_upper_bound(pg[get_name(g)], get_active_power_limits(g).max)
    end

    # incerted power at each bus
    for b in get_components(Bus, system)
        @constraint(opf_m, pb[get_name(b)] == sum(get_bus(g) == b ? pg[get_name(g)] : 0 for g in get_components(ThermalStandard, system))
                                             - sum(get_bus(l) == b ? get_max_active_power(l) : 0 for l in get_components(StaticLoad, system)))
    end

    # sum up the load and sum up renewable active power generation as negative load
    net_load = sum(get_max_active_power(g) for g in get_components(StaticLoad, system)) # - sum(get_max_active_power(g) for g in get_components(RenewableGen, system))
    # power balance
    @constraint(opf_m, power_balance, sum(pg[g] for g in thermal_gens_names) == net_load)

    # power flow on line and line limit
    # for l in get_components(Line, system)
    #     bus_from = get_name(get_arc(l).from)
    #     bus_to = get_name(get_arc(l).to)
    #     @constraint(opf_m, pb[bus_from] - pb[bus_to] == (angle[bus_from] - angle[bus_to]) / get_x(l))
    #     @constraint(opf_m, -get_rate(l) <= pb[bus_from] - pb[bus_to] <= get_rate(l))
    # end
    @constraint(opf_m, get_rate.(get_components(Line, system)) <= ptdf * pb <= get_rate.(get_components(Line, system)))

    # minimize cost of generation, quadratic costs
    @objective(opf_m, Min, sum(
        pg[get_name(g)]^2 * get_cost(get_variable(get_operation_cost(g)))[1] +
        pg[get_name(g)] * get_cost(get_variable(get_operation_cost(g)))[2]
        for g in get_components(ThermalGen, system)
        )
    )

    optimize!(opf_m)
    for g in thermal_gens_names
        println(g, " = ", value(pg[g]))
    end

    return opf_m
end

# add_system_data_to_json()
# system_data = System("system_data.json")
# results = opf_model(system_data, Ipopt.Optimizer)
# value.(results[:pl])
# write_to_file(results, "model.mps")
# latex_formulation(results)
# model = read_from_file("model.mps") # the names in the model will not be registered, use variable_by_name or constraint_by_name.