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
    to_json(system_data, joinpath(DATA_DIR, file_name), force=true)
end

@enum OPF SC=0 PSC=1 PCSC=2 # SCOPF, P-SCOPF, SC-SCOPF

""" OPF model type """
mutable struct OPFmodel
    sys::PowerSystems.System
    mod::JuMP.Model
    # beta::Union{Nothing, JuMP.Containers.DenseAxisArray} # Should convert to sparse some time in the future
    voll::JuMP.Containers.DenseAxisArray
    cost::JuMP.Containers.DenseAxisArray
    contingencies::Union{Nothing, Vector{String}}
    prob::Union{Nothing, JuMP.Containers.DenseAxisArray}
end

""" Constructor for OPFmodel """
function opfmodel(sys::System, optimizer, time_limit_sec, voll=nothing, contingencies=nothing, prob=nothing)
    mod = Model(optimizer)
    # if GLPK.Optimizer == optimizer 
    #     set_optimizer_attribute(mod, "msg_lev", GLPK.GLP_MSG_ON)
    # end
    c = [get_generator_cost(g)[2] for g in get_gens_t(sys)]
    cost = JuMP.Containers.DenseAxisArray(
            [c; [Statistics.mean(c) for _ in get_gens_h(sys)]],
            get_name.(get_ctrl_generation(sys))
        )
    # β = JuMP.Containers.DenseAxisArray{Int64}(undef, get_name.(get_nodes(sys)), get_name.(get_branches(sys)))
    # for l in get_branches(sys)
    #     β[get_name(l.arc.from), get_name(l)] = 1
    #     β[get_name(l.arc.to), get_name(l)] = -1
    # end

    # @assert length(voll) == length(get_nonctrl_generation(sys)) # should be implemented
    # @assert isempty(prob) || length(contingencies) == length(prob)

    return OPFmodel(sys, mod, voll, cost, contingencies, prob)
end

"""Find value of β, β = 1 if from-bus is the bus, β = -1 if to-bus is the bus, 0 else"""
beta(bus::Bus, branch::Branch) = bus == branch.arc.from ? 1 : bus == branch.arc.to ? -1 : 0
beta(bus::String, branch::Branch) = bus == branch.arc.from.name ? 1 : bus == branch.arc.to.name ? -1 : 0
function beta(sys::System, branch::String) 
    l = get_component(ACBranch, sys, branch)
    return [beta(b, l) for b in get_nodes(sys)]
end

a(line::String,contingency::String) = line != contingency ? 1 : 0

get_generator_cost(gen) = get_operation_cost(gen) |> get_variable |> get_cost |> get_value
get_value(x::Vector{Tuple{Float64, Float64}}) = x[1]
get_value(x::Tuple{Float64, Float64}) = x

""" An iterator to a type of power system component """
get_gens_t(sys::System) = get_components(ThermalGen, sys)
get_gens_h(sys::System) = get_components(HydroGen, sys)
get_branches(sys::System) = get_components(ACBranch, sys) # Includes both Line and Phase Shifting Transformer
get_nodes(sys::System) = get_components(Bus, sys)
get_demands(sys::System) = get_components(StaticLoad, sys)
get_renewables(sys::System) = get_components(RenewableGen, sys) # Renewable modelled as negative demand
get_ctrl_generation(sys::System) = Iterators.flatten((get_gens_t(sys), get_gens_h(sys))) # An iterator of all controllable generators
get_nonctrl_generation(sys::System) = Iterators.flatten((get_demands(sys), get_renewables(sys))) # An iterator of all non-controllable load and generation

""" A sorted vector to a type of power system component """
get_sorted_nodes(sys::System) = sort_nodes!(collect(get_components(Bus, sys)))
get_sorted_branches(sys::System) = sort_branches!(collect(get_components(ACBranch, sys)))
sort_nodes!(nodes::Vector{Bus}) = sort!(nodes,by = x -> x.number)
sort_branches!(branches::Vector{<: Branch}) = sort!(branches,
        by = x -> (get_number(get_arc(x).from), get_number(get_arc(x).to))
    )

# it_name(::Type{T}, mos::Model) where {T <: Component} = get_name.(get_components(T,mod))

""" Make a DenseAxisArray using the list and function for the value of each element """
make_named_array(value_func, list) = JuMP.Containers.DenseAxisArray(
    [value_func(x) for x in list], get_name.(list) 
)

""" An array of the Value Of Lost Load for the demand and get_renewables """
make_voll(sys::System) = JuMP.Containers.DenseAxisArray(
        [rand(1000:3000, length(get_demands(sys))); rand(1:30, length(get_renewables(sys)))], 
        [get_name.(get_demands(sys)); get_name.(get_renewables(sys))]
    )

""" An array of the outage probability of the contingencies """
make_prob(contingencies::Vector{String}) = JuMP.Containers.DenseAxisArray(
        (rand(length(contingencies)).*(0.5-0.02).+0.02)./8760, 
        contingencies
    )

""" A list with type_func componentes distributed on their node """
function make_list(opfm::OPFmodel, type_func, nodes = get_nodes(opfm.sys)) 
    list = JuMP.Containers.DenseAxisArray(
        [[] for _ in 1:length(nodes)], get_name.(nodes)
    )
    for g in type_func(opfm.sys)
        push!(list[g.bus.name], g)
    end
    return list
end

""" Return the (first) slack bus in the system. """
function find_slack(nodes)
    for (i,x) in enumerate(nodes)
        x.bustype == BusTypes.REF && return i,x
    end
    @warn "No slack bus found!"
end
find_slack(sys::System) = find_slack(get_sorted_nodes(sys))

""" Run optimizer to solve the model and check for optimality """
function solve_model!(model::Model)
    optimize!(model)
    if termination_status(model) != MOI.OPTIMAL 
        @warn "Model not optimally solved with status $(termination_status(model))!"
    else
        @info "Model solved in $(solve_time(model)) seconds with an objective value of $(objective_value(model))"
    end
    return model
end

" Sets the start values of symbols to the previous solution "
function set_warm_start!(model::JuMP.Model, symb::Symbol) 
    x_sol = JuMP.value.(model[symb])
    JuMP.set_start_value.(model[symb], x_sol)
    return model
end
function set_warm_start!(model::JuMP.Model, symbs::Vector{Symbol}) 
    x_sol = [JuMP.value.(model[symb]) for symb in symbs]
    for symb in symbs
        JuMP.set_start_value.(model[symb], x_sol)
    end
    return model
end

" Get a dict of bus-number-keys and vector-position-values "
get_nodes_idx(nodes::Vector{Bus}) = PowerSystems._make_ax_ref(nodes)

" Get the number of the from-bus and to-bus from a branch"
get_bus_id(branch::ACBranch) = (branch.arc.from.number, branch.arc.to.number)

" Return the overload of a line, else return 0.0 "
find_overload(flow::Float64, lim::Float64) = abs(flow)-lim > 0 ? sign(flow)*(abs(flow)-lim) : 0.0

" DC line flow calculation"
calculate_line_flows(isf::Array{Float64,2}, Pᵢ::Array{Float64}) = isf*Pᵢ

" DC power flow calculation "
run_pf(X::Array{Float64, 2}, Pᵢ::Array{Float64}) = X*Pᵢ 

""" Return the net power injected at each node sorted by bus number. """
function get_net_Pᵢ(opfm::OPFmodel, nodes::Vector{Bus}, Pᵢ = copy(get_Pᵢ(opfm, nodes)))
    p = JuMP.value.(opfm.mod[:ls0])
    for r in get_components(RenewableGen, opfm.sys)
        Pᵢ[r.bus.number] += get_active_power(r) - p[get_name(r)]
    end
    for d in get_components(StaticLoad, opfm.sys)
        Pᵢ[d.bus.number] -= get_active_power(d) + p[get_name(d)]
    end
    # @assert abs(sum(Pᵢ)) < 0.001
    return Pᵢ
end

""" Return the power injected at each node sorted by bus number. """
function get_Pᵢ(opfm::OPFmodel, nodes::Vector{Bus})
    Pᵢ = zeros(length(nodes))
    p = JuMP.value.(opfm.mod[:pg0])
    for g in get_components(Generator, opfm.sys)
        Pᵢ[g.bus.number] += p[get_name(g)]
    end
    return Pᵢ
end

