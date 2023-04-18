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

function setup(system::System, prob_min = 0.1, prob_max = 0.4)
    voll = make_voll(system)
    cost_renewables = make_cost_renewables(system)
    contingencies = get_name.(SCOPF.get_branches(system))
    prob = make_prob(contingencies, prob_min, prob_max)
    # contingencies = ["2-3-i_3"]
    return voll, cost_renewables, prob, contingencies
end

""" An array of the Value Of Lost Load for the demand and renewables """
make_voll(sys::System) = JuMP.Containers.DenseAxisArray(
        rand(1000:3000, length(get_demands(sys))), 
        get_name.(get_demands(sys))
    )
make_cost_renewables(sys::System) = JuMP.Containers.DenseAxisArray(
        rand(1:30, length(get_renewables(sys))), 
        get_name.(get_renewables(sys))
    )

""" An array of the outage probability of the contingencies """
make_prob(contingencies::AbstractVector{String}, prob_min = 0.1, prob_max = 0.4) = JuMP.Containers.DenseAxisArray(
        (rand(length(contingencies)).*(prob_max - prob_min) .+ prob_min)./8760, 
        contingencies
    )

get_system(fname::String) = System(joinpath("data",fname))

""" OPF model type """
mutable struct OPFmodel
    sys::PowerSystems.System
    mod::JuMP.Model
    voll::JuMP.Containers.DenseAxisArray
    cost_gen::JuMP.Containers.DenseAxisArray
    cost_renewables::JuMP.Containers.DenseAxisArray
    contingencies::Union{Nothing, Vector{String}}
    prob::Union{Nothing, JuMP.Containers.DenseAxisArray}
end

""" Constructor for OPFmodel """
function opfmodel(sys::System, optimizer, time_limit_sec, voll=nothing, cost_renewables=nothing, contingencies=nothing, prob=nothing)
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
    c = [get_generator_cost(g)[2] for g in get_gens_t(sys)]
    c_mean = Statistics.mean(c)
    cost_gen = JuMP.Containers.DenseAxisArray(
            [c; [c_mean for _ in get_gens_h(sys)]],
            get_name.(get_ctrl_generation(sys))
        )
    # β = JuMP.Containers.DenseAxisArray{Int64}(undef, get_name.(get_nodes(sys)), get_name.(get_branches(sys)))
    # for l in get_branches(sys)
    #     β[get_name(l.arc.from), get_name(l)] = 1
    #     β[get_name(l.arc.to), get_name(l)] = -1
    # end

    return OPFmodel(sys, mod, voll, cost_gen, cost_renewables, contingencies, prob)
end

"""Find value of β, β = 1 if from-bus is the bus, β = -1 if to-bus is the bus, 0 else"""
beta(bus::Bus, branch::Branch) = bus == branch.arc.from ? 1 : bus == branch.arc.to ? -1 : 0
beta(bus::String, branch::Branch) = bus == branch.arc.from.name ? 1 : bus == branch.arc.to.name ? -1 : 0
function beta(sys::System, branch::String) 
    l = get_component(ACBranch, sys, branch)
    return [beta(b, l) for b in get_nodes(sys)]
end

a(line::String,contingency::String) = ifelse(line != contingency, 1, 0)

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
get_sorted_nodes(sys::System) = sort_nodes!(collect(get_nodes(sys)))
get_sorted_branches(sys::System) = sort_branches!(collect(get_branches(sys)))
sort_nodes!(nodes::AbstractVector{Bus}) = sort!(nodes,by = x -> x.number)
sort_branches!(branches::AbstractVector{<:Branch}) = sort!(branches,
        by = x -> (get_number(get_arc(x).from), get_number(get_arc(x).to))
    )

# it_name(::Type{T}, mos::Model) where {T <:Component} = get_name.(get_components(T,mod))


""" Make a DenseAxisArray using the list and function for the value of each element """
make_named_array(value_func, list) = JuMP.Containers.DenseAxisArray(
    [value_func(x) for x in list], get_name.(list) 
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

find_in_model(mod::Model, ctrl::ThermalGen, name::String) = mod[:pg0][name]
find_in_model(mod::Model, ctrl::HydroGen, name::String) = mod[:pg0][name]
find_in_model(mod::Model, r::RenewableGen, name::String) = mod[:pr0][name]
find_in_model(mod::Model, d::StaticLoad, name::String) = -mod[:ls0][name]
find_in_model(mod::Model, dc::DCBranch, name::String) = -mod[:pfdc0][name]

""" Return the (first) slack bus in the system. 
Return: slack bus number in nodes. The slack bus.
"""
function find_slack(nodes::AbstractVector{Bus})
    for (i,x) in enumerate(nodes)
        x.bustype == BusTypes.REF && return i,x
    end
    @error "No slack bus found!"
end
find_slack(sys::System) = find_slack(get_sorted_nodes(sys))

""" Set the renewable production to a ratio of maximum active power """
function set_renewable_prod!(system::System, ratio::Real=0.5)
    for g in get_components(RenewableGen, system)
        set_active_power!(g, get_max_active_power(g)*ratio)
    end
    return system
end

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
get_nodes_idx(nodes::AbstractVector{Bus}) = _make_ax_ref(nodes)

""" Return a Dict of bus name number to index. Deprecated in PowerSystems """
_make_ax_ref(buses::AbstractVector{Bus}) = _make_ax_ref(get_number.(buses))

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
split_pair(val) = map(first, val), map(last, val)

" Get dual value (upper or lower bound) from model reference "
get_low_dual(varref::VariableRef) = dual(LowerBoundRef(varref))
get_high_dual(varref::VariableRef) = dual(UpperBoundRef(varref))

" Fix for name difference in StandardLoad "
PowerSystems.get_active_power(value::StandardLoad) = value.current_active_power

""" Return the net power injected at each node. """
function get_net_Pᵢ(opfm::OPFmodel, nodes::AbstractVector{Bus}, 
        idx::Dict{<:Any, <:Int} = get_nodes_idx(nodes), P = get_Pgen(opfm, nodes)
    )
    Pᵢ = copy(P)
    p = JuMP.value.(opfm.mod[:pr0])
    for r in get_renewables(opfm.sys)
        Pᵢ[idx[r.bus.number]] += get_active_power(r) - p[get_name(r)]
    end
    p = JuMP.value.(opfm.mod[:ls0])
    for d in get_demands(opfm.sys)
        Pᵢ[idx[d.bus.number]] -= get_active_power(d) + p[get_name(d)]
    end
    # @assert abs(sum(Pᵢ)) < 0.001
    return Pᵢ
end

function get_Pd(opfm::OPFmodel, nodes::AbstractVector{Bus}, 
        idx::Dict{<:Any, <:Int} = get_nodes_idx(nodes))
    P = zeros(length(nodes))
    for r in get_renewables(opfm.sys)
        add_device_value_to_node!(P, get_number(get_bus(r)), get_active_power(r), idx)
    end
    for d in get_demands(opfm.sys)
        add_device_value_to_node!(P, get_number(get_bus(d)), -get_active_power(d), idx)
    end
    # @assert abs(sum(P)) < 0.001
    return P
end

""" Return the power injected by controlled generation at each node. """
function get_Pgen(opfm::OPFmodel, nodes::AbstractVector{Bus}, 
        idx::Dict{<:Any, <:Int} = get_nodes_idx(nodes)
    )
    P = zeros(length(nodes))
    p = JuMP.value.(opfm.mod[:pg0])
    pn = get_number.(get_bus.(get_ctrl_generation(opfm.sys)))
    add_device_value_to_node!.([P], pn, p, [idx])

    p = JuMP.value.(opfm.mod[:pfdc0])
    a = get_arc.(get_dc_branches(opfm.sys))
    from = get_number.(get_from.(a))
    to = get_number.(get_to.(a))
    add_device_value_to_node!.([P], from, p, [idx])
    add_device_value_to_node!.([P], to, -p, [idx])
    return P
end

add_device_value_to_node!(arr::AbstractVector, device::Integer, val::Real, idx::Dict{<:Any, <:Int}) = arr[idx[device]] += val


" Calculate the severity index for the system based on line loading "
function calc_severity(opfm::OPFmodel, lim::Real = 0.9)
    rate = make_named_array(get_rate, get_branches(opfm.sys))
    sev = 0
    for c in opfm.contingencies
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
    sum(calc_line_severity(pfc[l,c], rate[l], lim) for c in contingencies, l in branches)

calc_severity(values::AbstractVector{<:Real}, rate::AbstractVector{<:Real}, lim::Real = 0.9) = 
    sum(calc_line_severity.(values, rate, [lim]))

" Calculate the severity index for a component based on the rating "
calc_line_severity(value::Real, rate::Real, lim::Real = 0.9) = 
    abs(value) > lim * rate ? (1/(1-lim)) * (abs(value) / rate - lim) : 0
