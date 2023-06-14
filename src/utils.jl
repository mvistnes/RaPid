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

function get_system(
        data_dir = "data\\RTS_GMLC",
        base_power = 100.0,
        descriptors = "data\\RTS_GMLC\\user_descriptors.yaml",
        timeseries_metadata_file = "data\\RTS_GMLC\\timeseries_pointers.json",
        generator_mapping_file = "data\\RTS_GMLC\\generator_mapping.yaml"
        )
    data = PowerSystemTableData(
        data_dir,
        base_power,
        descriptors;
        timeseries_metadata_file = timeseries_metadata_file,
        generator_mapping_file = generator_mapping_file,
    )
    return System(data)
end

function setup(system::System, prob_min = 0.1, prob_max = 0.4)
    voll = make_voll(system)
    contingencies = sort_components!(get_branches(system))
    prob = make_prob(contingencies, prob_min, prob_max)
    # contingencies = ["2-3-i_3"]
    return voll, prob, contingencies
end

function fix_generation_cost!(sys::System)
    # set_operation_cost!.(get_renewables(sys), make_cost_renewables(sys))
    cost = [get_generator_cost(g)[2] for g in get_gens_t(sys)]
    c_mean = Statistics.mean(cost)
    set_operation_cost!.(get_ctrl_generation(sys), 
        [[iszero(c) ? c_mean * (rand() + 0.5) : c for c in cost]; 
        [c_mean * (rand() + 0.1) for _ in get_gens_h(sys)]])
end

set_operation_cost!(gen::Generator, val::Real) = 
    PowerSystems.set_operation_cost!(gen, ThreePartCost(VariableCost((0.0, val)), 0.0, 0.0, 0.0))
set_operation_cost!(gen::RenewableGen, val::Real) = 
    PowerSystems.set_operation_cost!(gen, TwoPartCost(VariableCost((0.0, val)), 0.0))

""" Set the renewable production to a ratio of maximum active power """
function set_renewable_prod!(system::System, ratio::Real=0.5)
    for g in get_renewables(system)
        set_active_power!(g, get_max_active_power(g)*ratio)
    end
end

""" Set the active power demand """
function set_active_power_demand!(system::System, demands = get_max_active_power.(get_demands(system)))
    set_active_power!.(get_demands(system), demands)
end

PowerSystems.set_active_power!(dem::StandardLoad, val::Real) = 
    PowerSystems.set_current_active_power!(dem, val)

function set_ramp_limits!(system::System, ramp_mult::Real = 0.01)
    gens = get_ctrl_generation(system)
    ramp_lim = [(up=(x.max*ramp_mult), down=(x.max*ramp_mult)) for x in get_active_power_limits.(gens)]
    PowerSystems.set_ramp_limits!.(gens,ramp_lim)
end

""" An array of the Value Of Lost Load for the demand and renewables """
make_voll(sys::System, x_range::UnitRange=1000:3000) = rand(x_range, length(get_demands(sys)))
make_cost_renewables(sys::System, x_range::UnitRange=1:30) = rand(x_range, length(get_renewables(sys)))

""" An array of the outage probability of the contingencies """
make_prob(contingencies::AbstractVector, prob_min = 0.1, prob_max = 0.4) = 
    (rand(length(contingencies)).*(prob_max - prob_min) .+ prob_min)./8760

get_system(fname::String) = System(joinpath("data",fname))

""" OPF model type """
mutable struct OPFmodel
    sys::PowerSystems.System
    mod::JuMP.Model

    cost_ctrl_gen::Vector{<:Real}
    cost_renewables::Vector{<:Real}
    voll::Vector{<:Real}
    prob::Vector{<:Real}

    ctrl_generation::Vector{Generator}
    branches::Vector{<:ACBranch}
    dc_branches::Vector{<:DCBranch}
    nodes::Vector{Bus}
    demands::Vector{StaticLoad}
    renewables::Vector{RenewableGen}
    contingencies::Vector{<:Component}
end

""" Constructor for OPFmodel """
function opfmodel(sys::System, optimizer, time_limit_sec, voll::Vector{<:Real}, 
        contingencies::Vector{<:Branch}=Component[], prob::Vector{<:Real}=[])
    mod = Model(optimizer)
    # mod = Model(optimizer; add_bridges = false)
    set_string_names_on_creation(mod, false)
    # @debug begin
    #     set_string_names_on_creation(mod, true)
    #     "Variable names is on."
    # end

    # if GLPK.Optimizer == optimizer 
    #     set_optimizer_attribute(mod, "msg_lev", GLPK.GLP_MSG_ON)
    # end

    ctrl_generation = collect(get_ctrl_generation(sys))
    branches = sort_components!(get_branches(sys))
    dc_branches = sort_components!(get_dc_branches(sys))
    nodes = sort_components!(get_nodes(sys))
    demands = collect(get_demands(sys))
    renewables = collect(get_renewables(sys))

    cost_ctrl_gen = Vector{Float64}([get_generator_cost(g)[2] for g in ctrl_generation])
    cost_renewables = Vector{Float64}() # Vector{Float64}([get_generator_cost(g)[2] for g in renewables])

    check_values(cost_ctrl_gen, ctrl_generation, "cost")
    check_values(voll, demands, "voll")
    check_values(getindex.(get_active_power_limits.(ctrl_generation), :max), ctrl_generation, "max_active_power")
    check_values(get_active_power.(demands), demands, "max_active_power")
    check_values(get_active_power.(renewables), renewables, "max_active_power")
    (rampup, rampdown) = split_pair(get_ramp_limits.(ctrl_generation))
    check_values(rampup, ctrl_generation, "rampup")
    check_values(rampdown, ctrl_generation, "rampdown")
    check_values(get_r.(branches), branches, "r")
    check_values(get_x.(branches), branches, "x")
    check_values(get_rate.(branches), branches, "rate")
    limits = split_pair(get_active_power_limits_from.(dc_branches))
    check_values(limits[1], dc_branches, "rate_dc_branches_min")
    check_values(limits[2], dc_branches, "rate_dc_branches_max")

    islands = island_detection(nodes, branches)
    if length(islands) > 1
        @info "The system is separated into islands" islands
    end

    return OPFmodel(sys, mod, cost_ctrl_gen, cost_renewables, voll, prob, 
        ctrl_generation, branches, dc_branches, nodes, demands, renewables, contingencies)
end

""" Automatic constructor for OPFmodel where voll, prob, and contingencies are automatically computed. """
function opfmodel(sys::System, optimizer, time_limit_sec)
    voll, prob, contingencies = setup(sys, 1, 4)
    fix_generation_cost!(sys)
    return opfmodel(sys, optimizer, time_limit_sec, voll, contingencies, prob)
end

function check_values(val::Real, var::Component, name::String; atol=1e-5)
    if isapprox(val, zero(typeof(val)); atol=atol)
        @info "$(get_name(var)) has $val value in $name."
    end
end
check_values(val, var, name; atol=1e-5) = check_values.(val, var, [name], atol=atol)

"""Find value of β, β = 1 if from-bus is the bus, β = -1 if to-bus is the bus, 0 else"""
beta(bus::Bus, branch::Branch) = bus == branch.arc.from ? 1 : bus == branch.arc.to ? -1 : 0
beta(bus::String, branch::Branch) = bus == branch.arc.from.name ? 1 : bus == branch.arc.to.name ? -1 : 0
function beta(sys::System, branch::String) 
    l = get_component(ACBranch, sys, branch)
    return [beta(b, l) for b in get_nodes(sys)]
end

a(line::String, contingency::String) = ifelse(line != contingency, 1, 0)

get_generator_cost(gen) = get_operation_cost(gen) |> get_variable |> get_cost |> _get_g_value
_get_g_value(x::AbstractVector{<:Tuple{Real, Real}}) = x[1]
_get_g_value(x::Tuple{<:Real, <:Real}) = x

""" An iterator to a type of power system component """
get_gens_t(sys::System) = get_components(ThermalGen, sys)
get_gens_h(sys::System) = get_components(HydroGen, sys)
get_branches(sys::System) = get_components(ACBranch, sys) # Includes both Line, Transformer, and Phase Shifting Transformer
get_dc_branches(sys::System) = get_components(DCBranch, sys) 
get_nodes(sys::System) = get_components(Bus, sys)
get_demands(sys::System) = get_components(StaticLoad, sys)
get_renewables(sys::System) = get_components(RenewableGen, sys) # Renewable modelled as negative demand
get_ctrl_generation(sys::System) = Iterators.flatten((get_gens_t(sys), get_gens_h(sys))) # An iterator of all controllable generators
get_generation(sys::System) = Iterators.flatten((get_ctrl_generation(sys), get_renewables(sys))) # An iterator of all generation

""" A sorted vector to a type of power system component """
sort_components!(list::PowerSystems.FlattenIteratorWrapper) = sort_components!(collect(list))
sort_components!(nodes::AbstractVector{Bus}) = sort!(nodes, by = x -> x.number)
sort_components!(components::AbstractVector{<:Component}) = sort!(components, by = x -> x.bus.number)
sort_components!(branches::AbstractVector{<:Branch}) = sort!(branches,
        by = x -> (get_number(get_arc(x).from), get_number(get_arc(x).to))
    )

get_sorted_angles(model) = JuMP.value.(model[:va0])

# it_name(::Type{T}, mos::Model) where {T <:Component} = get_name.(get_components(T,mod))


""" Make a DenseAxisArray using the list and function for the value of each element """
make_named_array(value_func, list) = JuMP.Containers.DenseAxisArray(
    [value_func(x) for x in list], get_name.(list) 
)

struct CTypes
    ctrl_generation
    renewables
    demands
    branches
    dc_branches
end

""" Make a list where all components distributed on their node """
function make_list(opfm::OPFmodel, idx::Dict{<:Any, <:Int}, nodes::AbstractVector{Bus}) 
    list = [CTypes([Int[] for _ in 1:5]...) for _ in 1:length(nodes)]
    for (i,r) in enumerate(opfm.ctrl_generation)
        push!(list[idx[r.bus.number]].ctrl_generation, i)
    end
    for (i,r) in enumerate(opfm.renewables)
        push!(list[idx[r.bus.number]].renewables, i)
    end
    for (i,d) in enumerate(opfm.demands)
        push!(list[idx[d.bus.number]].demands, i)
    end
    for (i,d) in enumerate(opfm.branches)
        push!(list[idx[d.arc.from.number]].branches, i)
        push!(list[idx[d.arc.to.number]].branches, i)
    end
    for (i,d) in enumerate(opfm.dc_branches)
        push!(list[idx[d.arc.from.number]].dc_branches, i)
        push!(list[idx[d.arc.to.number]].dc_branches, i)
    end
    return list
end

""" A list with type_func componentes distributed on their node """
function make_list(system::System, type_func, nodes = get_nodes(system)) 
    list = JuMP.Containers.DenseAxisArray(
        [[] for _ in 1:length(nodes)], get_name.(nodes)
    )
    for g in type_func(system)
        push!(list[g.bus.name], g)
    end
    return list
end

find_in_model(mod::Model, ctrl::ThermalGen, name::String) = mod[:pg0][name]
find_in_model(mod::Model, ctrl::HydroGen, name::String) = mod[:pg0][name]
find_in_model(mod::Model, r::RenewableGen, name::String) = mod[:pr0][name]
find_in_model(mod::Model, d::StaticLoad, name::String) = mod[:ls0][name]
find_in_model(mod::Model, dc::DCBranch, name::String) = mod[:pfdc0][name]

""" Return the (first) slack bus in the system. 
Return: slack bus number in nodes. The slack bus.
"""
function find_slack(nodes::AbstractVector{Bus})
    for (i,x) in enumerate(nodes)
        x.bustype == BusTypes.REF && return i,x
    end
    @error "No slack bus found!"
end
find_slack(sys::System) = find_slack(sort_components!(get_nodes(sys)))

""" Run optimizer to solve the model and check for optimality """
function solve_model!(model::Model)
    optimize!(model)
    if termination_status(model) != MOI.OPTIMAL 
        @warn "Model not optimally solved with status $(termination_status(model))!"
    else
        @info @sprintf "Model solved in %.6fs with objective value %.4f" solve_time(model) objective_value(model)
    end
    return model
end

" Sets the start values of symbols to the previous solution "
function set_warm_start!(model::JuMP.Model, symb::Symbol) 
    x_sol = JuMP.value.(model[symb])
    JuMP.set_start_value.(model[symb], x_sol)
    return model
end
function set_warm_start!(model::JuMP.Model, symbs::AbstractVector{Symbol}) 
    x_sol = [JuMP.value.(model[symb]) for symb in symbs]
    for symb in symbs
        JuMP.set_start_value.(model[symb], x_sol)
    end
    return model
end

" Get a dict of bus-number-keys and vector-position-values "
get_nodes_idx(nodes::AbstractVector{Bus}) = _make_ax_ref(get_number.(nodes))

""" Return a Dict of bus name number to index. Deprecated in PowerSystems """
function _make_ax_ref(ax::AbstractVector)
    ref = Dict{eltype(ax), Int}()
    for (ix, el) in enumerate(ax)
        if haskey(ref, el)
            @error("Repeated index element $el. Index sets must have unique elements.")
        end
        ref[el] = ix
    end
    return ref
end

# TODO: remove the need for idx in most functions, use array

" Get the bus number of the from-bus and to-bus from a branch "
get_bus_id(branch::Branch) = (branch.arc.from.number, branch.arc.to.number)

" Get the bus number index of the from-bus and to-bus from a branch "
get_bus_idx(branch::Branch, idx::Dict{<:Any, <:Int}) = 
    (idx[branch.arc.from.number], idx[branch.arc.to.number])

" Get the bus number index of the from-bus and to-bus from all branches "
get_bus_idx(branches::AbstractVector{<:Branch}, idx::Dict{<:Any, <:Int}) =
    split_pair(get_bus_idx.(branches, [idx]))

get_branch_bus_idx(branches::AbstractVector{<:Branch}, contingencies::AbstractVector{<:Branch}, idx::Dict{<:Any, <:Int}) = 
    [(findfirst(x -> x == b, branches), get_bus_idx(b, idx)) for b in contingencies]

function extract_bus_idx(branch::String, idx::Dict{<:Any, <:Int})
    (f,t) = extract_bus(branch)
    return (idx[f], idx[t])
end

function extract_bus(branch::String)
    x = split(branch, "-")
    return (parse(Int,x[1]), parse(Int,x[2]))
end

" Get the bus number index of the from-bus and to-bus from all branches "
function get_bus_idx(A::AbstractMatrix)
    m, n = size(A)
    ix = [[0,0] for _ in 1:m]
    rows = SparseArrays.rowvals(A)
    for j = 1:n
        for i in SparseArrays.nzrange(A, j)
            row = rows[i]
            if iszero(ix[row][1])
                ix[row][1] = j
            else
                ix[row][2] = j
            end
        end
    end
    return ix
end

" Split a Vector{Pair} into a Pair{Vector}"
split_pair(val::AbstractVector) = map(first, val), map(last, val)

" zip a Pair{Vector, Vector} into a Vector{Pair}"
zip_pair(val::Pair{AbstractVector, AbstractVector}) = [(a,b) for (a,b) in zip(first(val),last(val))]
zip_pair(val1::AbstractVector, val2::AbstractVector) = zip_pair((val1, val2))

" Get dual value (upper or lower bound) from model reference "
get_low_dual(varref::VariableRef) = dual(LowerBoundRef(varref))
get_high_dual(varref::VariableRef) = dual(UpperBoundRef(varref))

" Fix for name difference in StandardLoad "
PowerSystems.get_active_power(value::StandardLoad) = value.current_active_power

""" Return the net power injected at each node. """
function get_net_Pᵢ(opfm::OPFmodel, idx::Dict{<:Any, <:Int} = get_nodes_idx(opfm.nodes)
    )
    P = zeros(length(opfm.nodes))
    for (i,r) in enumerate(opfm.renewables)
        P[idx[r.bus.number]] += get_active_power(r) - JuMP.value(opfm.mod[:pr0][i])
    end
    for (i,d) in enumerate(opfm.demands)
        P[idx[d.bus.number]] -= get_active_power(d) + JuMP.value(opfm.mod[:ls0][i])
    end
    for (i,r) in enumerate(opfm.ctrl_generation)
        P[idx[r.bus.number]] += JuMP.value(opfm.mod[:pg0][i])
    end
    for (i,d) in enumerate(opfm.dc_branches)
        P[idx[d.arc.from.number]] += JuMP.value(opfm.mod[:pfdc0][i])
        P[idx[d.arc.to.number]] -= JuMP.value(opfm.mod[:pfdc0][i])
    end
    # @assert abs(sum(Pᵢ)) < 0.001
    return P
end

function get_Pd!(P::Vector{<:Real}, opfm::OPFmodel, idx::Dict{<:Any, <:Int})
    for r in opfm.renewables
        P[idx[r.bus.number]] += get_active_power(r)
    end
    for d in opfm.demands
        P[idx[d.bus.number]] -= get_active_power(d)
    end
end
function get_Pd(opfm::OPFmodel, idx::Dict{<:Any, <:Int} = get_nodes_idx(opfm.nodes))
    P = zeros(length(opfm.nodes))
    get_Pd!(P, opfm, idx)
    return P
end

""" Return power shed at each node. """
function get_Pshed!(P::Vector{<:Real}, opfm::OPFmodel, idx::Dict{<:Any, <:Int})
    for (i,r) in enumerate(opfm.renewables)
        P[idx[r.bus.number]] -= JuMP.value(opfm.mod[:pr0][i])
    end
    for (i,d) in enumerate(opfm.demands)
        P[idx[d.bus.number]] += JuMP.value(opfm.mod[:ls0][i])
    end
end
function get_Pshed(opfm::OPFmodel, idx::Dict{<:Any, <:Int} = get_nodes_idx(opfm.nodes))
    P = zeros(length(opfm.nodes))
    get_Pshed!(P, opfm, idx)
    return P
end

""" Return the power injected by controlled generation at each node. """
function get_Pgen!(P::Vector{<:Real}, opfm::OPFmodel, idx::Dict{<:Any, <:Int})
    for (i,r) in enumerate(opfm.ctrl_generation)
        P[idx[r.bus.number]] += JuMP.value(opfm.mod[:pg0][i])
    end
    for (i,d) in enumerate(opfm.dc_branches)
        P[idx[d.arc.from.number]] += JuMP.value(opfm.mod[:pfdc0][i])
        P[idx[d.arc.to.number]] -= JuMP.value(opfm.mod[:pfdc0][i])
    end
    return P
end
function get_Pgen(opfm::OPFmodel, idx::Dict{<:Any, <:Int} = get_nodes_idx(opfm.nodes))
    P = zeros(length(opfm.nodes))
    get_Pgen!(P, opfm, idx)
    return P
end

function get_controllable(opfm::OPFmodel, idx::Dict{<:Any, <:Int} = get_nodes_idx(opfm.nodes))
    P = get_Pgen(opfm, idx)
    get_Pshed!(P, opfm, idx)
    return P
end

" Calculate the severity index for the system based on line loading "
function calc_severity(opfm::OPFmodel, lim::Real = 0.9)
    rate = make_named_array(get_rate, get_branches(opfm.sys))
    sev = 0
    for c in 1:length(opfm.contingencies)
        for l in get_name.(get_branches(opfm.sys))
            sev += calc_line_severity(value(opfm.mod[:pfc][l,c]), rate[l], lim)
        end
    end 
    return sev
end

calc_severity(
        rate::AbstractVector{<:Real}, 
        contingencies::AbstractVector, 
        branches::AbstractVector, 
        pfc, 
        lim::Real = 0.9
        ) = 
    sum(calc_line_severity(pfc[l,c], rate[l], lim) for c in 1:length(opfm.contingencies), l in branches)

calc_severity(values::AbstractVector{<:Real}, rate::AbstractVector{<:Real}, lim::Real = 0.9) = 
    sum(calc_line_severity.(values, rate, [lim]))

" Calculate the severity index for a component based on the rating "
calc_line_severity(value::Real, rate::Real, lim::Real = 0.9) = 
    abs(value) > lim * rate ? (1/(1-lim)) * (abs(value) / rate - lim) : 0
