function add_system_data_to_json(;
    file_name="system_data.json",
    data_name=joinpath("matpower", "IEEE_RTS.m"),
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
    data_dir="data\\RTS_GMLC",
    base_power=100.0,
    descriptors="data\\RTS_GMLC\\user_descriptors.yaml",
    timeseries_metadata_file="data\\RTS_GMLC\\timeseries_pointers.json",
    generator_mapping_file="data\\RTS_GMLC\\generator_mapping.yaml"
)
    data = PowerSystemTableData(
        data_dir,
        base_power,
        descriptors;
        timeseries_metadata_file=timeseries_metadata_file,
        generator_mapping_file=generator_mapping_file
    )
    return System(data)
end

function setup(system::System, prob_min=0.1, prob_max=0.4)
    voll = make_voll(system)
    contingencies = sort_components!(get_branches(system))
    prob = make_prob(contingencies, prob_min, prob_max)
    return voll, prob, contingencies
end

function fix_generation_cost!(sys::System)
    # set_operation_cost!.(get_renewables(sys), make_cost_renewables(sys))
    cost = [get_generator_cost(g)[2] for g in get_gens_t(sys)]
    c_mean = mean(cost)
    set_operation_cost!.(get_ctrl_generation(sys),
        [[iszero(c) ? c_mean * (rand() + 0.5) : c for c in cost]
            [c_mean * (rand() + 0.1) for _ in get_gens_h(sys)]])
end

set_operation_cost!(gen::Generator, val::Real) =
    PowerSystems.set_operation_cost!(gen, ThreePartCost(VariableCost((0.0, val)), 0.0, 0.0, 0.0))
set_operation_cost!(gen::RenewableGen, val::Real) =
    PowerSystems.set_operation_cost!(gen, TwoPartCost(VariableCost((0.0, val)), 0.0))

""" Set the renewable production to a ratio of maximum active power """
function set_renewable_prod!(system::System, ratio::Real=0.5)
    for g in get_renewables(system)
        set_active_power!(g, get_max_active_power(g) * ratio)
    end
end

""" Set the active power demand """
function set_active_power_demand!(system::System, demands=get_max_active_power.(get_demands(system)))
    set_active_power!.(get_demands(system), demands)
end

PowerSystems.set_active_power!(dem::StandardLoad, val::Real) =
    PowerSystems.set_current_active_power!(dem, val)

function set_ramp_limits!(system::System, ramp_mult::Real=0.01)
    set_ramp_limit!.(get_ctrl_generation(system), ramp_mult)
end
function set_ramp_limit!(gen::Generator, ramp_mult::Real=0.01, min_ramp=1e-6)
    ramp_lim = PowerSystems.get_ramp_limits(gen)
    if ramp_lim < (up=min_ramp, down=min_ramp)
        p_lim = get_active_power_limits(gen)
        PowerSystems.set_ramp_limits!(gen, (up=(p_lim.max * ramp_mult), down=(p_lim.max * ramp_mult)))
    end
end

get_generator_cost(gen::Generator) = get_operation_cost(gen) |> get_variable |> get_cost |> _get_g_value
_get_g_value(x::AbstractVector{<:Tuple{Real,Real}}) = x[1]
_get_g_value(x::Tuple{<:Real,<:Real}) = x

function get_generator_cost(ctrl_generation::AbstractVector{<:Generator}, ramp_cost::Float64)
    cost = Vector{NamedTuple{(:fix, :var, :ramp)}}(undef, length(ctrl_generation))
    for (i,g) in enumerate(ctrl_generation)
        c = get_generator_cost(g)
        cost[i] = (fix=c[1], var=c[2], ramp=c[2]*ramp_cost)
    end
    return cost
end

_get_g_value(x::T) where {T<:Real} = (zero(T), x)

get_y(value::Branch) = 1 / (get_r(value) + get_x(value) * im)

function get_specs(branch::Line)
    y = get_y(branch)
    return real(y), imag(y), get_b(branch), 1.0, 1.0, 0.0
end
function get_specs(branch::Transformer2W)
    y = get_y(branch)
    B = (from=PowerSystems.get_primary_shunt(branch), to=0.0)
    return real(y), imag(y), B, 1.0, 1.0, 0.0
end
function get_specs(branch::TapTransformer)
    y = get_y(branch)
    B = (from=PowerSystems.get_primary_shunt(branch), to=0.0)
    tap = get_tap(branch)
    tr = tap
    return real(y), imag(y), B, tap, tr, 0.0
end
function get_specs(branch::PhaseShiftingTransformer)
    y = get_y(branch)
    B = (from=PowerSystems.get_primary_shunt(branch), to=0.0)
    tap = get_tap(branch)
    angle = get_α(branch)
    tr = tap * cos(angle)
    ti = tap * sin(angle)
    return real(y), imag(y), B, tap, tr, ti
end

""" An array of the Value Of Lost Load for the demand and renewables """
make_voll(sys::System, x_range::UnitRange=UnitRange(1000.0,3000.0)) = rand(x_range, length(get_demands(sys)))
make_cost_renewables(sys::System, x_range::UnitRange=UnitRange(1.0,30.0)) = rand(x_range, length(get_renewables(sys)))

""" An array of the outage probability of the contingencies """
make_prob(contingencies::AbstractVector, prob_min=0.1, prob_max=0.4) =
    (rand(length(contingencies)) .* (prob_max - prob_min) .+ prob_min) ./ 8760

get_system(fname::String) = System(joinpath("data", fname))

function create_model(optimizer; time_limit_sec::Integer=10000, silent::Bool=true, debug::Bool=false)
    mod = Model(optimizer; add_bridges=false)
    set_string_names_on_creation(mod, debug)
    MOI.set(mod, MOI.Silent(), silent) # supress output from the solver
    set_time_limit_sec(mod, time_limit_sec)
    return mod
end

function check_values(val::Real, var::Component, typename::String; comp=<, atol=1e-5)
    if comp(val, atol)
    # if isapprox(val, zero(typeof(val)); atol=atol)
        @info "$(get_name(var)) has $val value in $typename."
    end
end
check_values(val::AbstractVector{<:Real}, var::AbstractVector{<:Component}, typename::String; atol=1e-5) =
    check_values.(val, var, [typename], atol=atol)

function check_values(val::Generator)
    lmt = get_active_power_limits(val)
    check_values(lmt[:min], val, "min_active_power")
    check_values(lmt[:max], val, "max_active_power")
    ramp = get_ramp_limits(val)
    check_values(ramp[1], val, "rampup")
    check_values(ramp[2], val, "rampdown")
end

function check_values(val::ACBranch)
    check_values(get_r(val), val, "r")
    check_values(get_x(val), val, "x")
    check_values(get_rate(val), val, "rate")
end

function check_values(val::TwoTerminalHVDCLine)
    lmt = get_active_power_limits_from(val)
    check_values(lmt[1], val, "rate_dc_branches_min", comp=≈)
    check_values(lmt[2], val, "rate_dc_branches_max")
end

"""Find value of β, β = 1 if from-bus is the bus, β = -1 if to-bus is the bus, 0 else"""
beta(bus::ACBus, branch::Branch) = bus == branch.arc.from ? 1 : bus == branch.arc.to ? -1 : 0
beta(bus::String, branch::Branch) = bus == branch.arc.from.name ? 1 : bus == branch.arc.to.name ? -1 : 0
function beta(sys::System, branch::String)
    l = get_component(ACBranch, sys, branch)
    return [beta(b, l) for b in get_nodes(sys)]
end

a(line::String, contingency::String) = ifelse(line != contingency, 1, 0)

""" An iterator to a type of power system component """
get_gens_t(sys::System) = get_components(ThermalGen, sys) |> collect
get_gens_h(sys::System) = get_components(HydroGen, sys) |> collect
get_lines(sys::System) = get_components(Line, sys) |> collect 
get_transformers2w(sys::System) = get_components(Transformer2W, sys) |> collect
get_taptransformers(sys::System) = get_components(TapTransformer, sys) |> collect
get_phaseshiftingtransformers(sys::System) = get_components(PhaseShiftingTransformer, sys) |> collect
get_branches(sys::System) = vcat(get_lines(sys), get_transformers2w(sys), get_taptransformers(sys), get_phaseshiftingtransformers(sys))
get_dc_branches(sys::System) = get_components(TwoTerminalHVDCLine, sys) |> collect
get_nodes(sys::System) = get_components(ACBus, sys) |> collect
get_demands(sys::System) = get_components(StaticLoad, sys) |> collect
get_renewables(sys::System) = get_components(RenewableGen, sys) |> collect # Renewable modelled as negative demand
get_ctrl_generation(sys::System) = vcat(get_gens_t(sys), get_gens_h(sys)) # An iterator of all controllable generators
get_generation(sys::System) = vcat(get_ctrl_generation(sys), get_renewables(sys)) # An iterator of all generation

""" A sorted vector to a type of power system component """
sort_components!(list::PowerSystems.FlattenIteratorWrapper{T}) where {T} = sort_components!(collect(T, list))
sort_components!(nodes::Vector{ACBus}) = sort!(nodes, by=x -> x.number)
sort_components!(components::Vector{<:Component}) = sort!(components, by=x -> x.bus.number)
sort_components!(branches::Vector{<:Branch}) = sort!(branches,
    by=x -> (get_number(get_arc(x).from), get_number(get_arc(x).to))
)

""" OPF system type """
mutable struct OPFsystem{T<:Real}
    cost_ctrl_gen::Vector{NamedTuple{(:fix, :var, :ramp), Tuple{T, T, T}}}
    cost_renewables::Vector{T}
    voll::Vector{T}
    prob::Vector{T}

    ctrl_generation::Vector{Generator}
    branches::Vector{ACBranch}
    dc_branches::Vector{TwoTerminalHVDCLine}
    nodes::Vector{ACBus}
    demands::Vector{StaticLoad}
    renewables::Vector{RenewableGen}
    contingencies::Vector{Component}
end

""" Constructor for OPFsystem """
function opfsystem(sys::System, voll::Vector{T}, contingencies::Vector{<:Component}=Component[], prob::Vector{T}=[],
    ramp_mult::Real=1.0; check=false
) where {T<:Real}

    ctrl_generation = sort_components!(get_ctrl_generation(sys))
    branches = sort_components!(get_branches(sys))
    dc_branches = sort_components!(get_dc_branches(sys))
    nodes = sort_components!(get_nodes(sys))
    demands = sort_components!(get_demands(sys))
    renewables = sort_components!(get_renewables(sys))

    cost_ctrl_gen = get_generator_cost(ctrl_generation, ramp_mult)
    cost_renewables = Vector{Float64}() # Vector{Float64}([get_generator_cost(g)[2] for g in renewables])

    if check
        check_values(getfield.(cost_ctrl_gen, :fix), ctrl_generation, "fixedcost")
        check_values(getfield.(cost_ctrl_gen, :var), ctrl_generation, "varcost")
        check_values(getfield.(cost_ctrl_gen, :ramp), ctrl_generation, "rampcost")
        check_values(voll, demands, "voll")
        check_values.(ctrl_generation)
        check_values.(branches)
        check_values.(dc_branches)
    end

    A = calc_B(branches, nodes)
    islands = island_detection(A'A)
    if length(islands) > 1
        @warn "The system is separated into islands" islands
    end

    return OPFsystem{T}(cost_ctrl_gen, cost_renewables, voll, prob,
        ctrl_generation, branches, dc_branches, nodes, demands, renewables, contingencies)
end

""" Automatic constructor for OPFsystem where voll, prob, and contingencies are automatically computed. """
function opfsystem(sys::System)
    voll, prob, contingencies = setup(sys, 1, 4)
    fix_generation_cost!(sys)
    return opfsystem(sys, voll, contingencies, prob)
end

get_cost_ctrl_gen(opf::OPFsystem) = opf.cost_ctrl_gen
get_cost_renewables(opf::OPFsystem) = opf.cost_renewables
get_voll(opf::OPFsystem) = opf.voll
get_prob(opf::OPFsystem) = opf.prob

get_ctrl_generation(opf::OPFsystem) = opf.ctrl_generation
get_branches(opf::OPFsystem) = opf.branches
get_dc_branches(opf::OPFsystem) = opf.dc_branches
get_nodes(opf::OPFsystem) = opf.nodes
get_demands(opf::OPFsystem) = opf.demands
get_renewables(opf::OPFsystem) = opf.renewables
get_contingencies(opf::OPFsystem) = opf.contingencies

function get_interarea(branches::AbstractVector{T}) where {T<:Branch}
    vals = Vector{T}()
    for br in branches
        from = get_area(get_from(get_arc(br)))
        to = get_area(get_to(get_arc(br)))
        if from != to
            push!(vals, br)
        end
    end
    return vals
end

function typesort_components(list::AbstractVector{<:Component})
    gen = Generator[]
    demand = StaticLoad[]
    bus = ACBus[]
    branch = ACBranch[]
    for comp in list
        t = typeof(comp)
        if t <: Generator
            push!(gen, comp)
        elseif t <: StaticLoad
            push!(demand, comp)
        elseif t <: ACBus
            push!(bus, comp)
        elseif t <: ACBranch
            push!(branch, comp)
        else
            error(string(typeof(comp)) + " is not supported")
        end
    end
    return (gen=gen, demand=demand, bus=bus, branch=branch)
end

function find_component(val::Component, list::AbstractVector{<:Component})
    i = findfirst(in([val]), list)
    i === nothing && error(string(get_name(val)) * " is not found")
    return i
end

typesort_component(val::Generator, opf::OPFsystem, idx::Dict{<:Int,<:Int}) =
    (:ctrl_generation, find_component(val, get_ctrl_generation(opf)), get_bus_idx(val, idx))

typesort_component(val::StaticLoad, opf::OPFsystem, idx::Dict{<:Int,<:Int}) =
    (:demands, find_component(val, get_demands(opf)), get_bus_idx(val, idx))

typesort_component(val::ACBus, opf::OPFsystem, idx::Dict{<:Int,<:Int}) =
    (:nodes, find_component(val, get_nodes(opf)), get_bus_idx(val, idx))

typesort_component(val::ACBranch, opf::OPFsystem, idx::Dict{<:Int,<:Int}) =
    (:branches, find_component(val, get_branches(opf)), get_bus_idx(val, idx))

""" Make a DenseAxisArray using the list and function for the value of each element """
make_named_array(value_func, list) = JuMP.Containers.DenseAxisArray(
    [value_func(x) for x in list], get_name.(list)
)

struct CTypes{T<:Integer}
    node::ACBus
    ctrl_generation::Vector{T}
    renewables::Vector{T}
    demands::Vector{T}
    branches::Vector{T}
    dc_branches::Vector{T}
end

""" Make a list where all component numbers are distributed on their node """
function make_list(opf::OPFsystem, idx::Dict{<:Int,<:Int}, nodes::AbstractVector{ACBus})
    list = [CTypes(n, [Int[] for _ in 1:5]...) for n in nodes]
    for (i, r) in enumerate(opf.ctrl_generation)
        push!(list[idx[r.bus.number]].ctrl_generation, i)
    end
    for (i, r) in enumerate(opf.renewables)
        push!(list[idx[r.bus.number]].renewables, i)
    end
    for (i, d) in enumerate(opf.demands)
        push!(list[idx[d.bus.number]].demands, i)
    end
    for (i, d) in enumerate(opf.branches)
        push!(list[idx[d.arc.from.number]].branches, i)
        push!(list[idx[d.arc.to.number]].branches, i)
    end
    for (i, d) in enumerate(opf.dc_branches)
        push!(list[idx[d.arc.from.number]].dc_branches, i)
        push!(list[idx[d.arc.to.number]].dc_branches, i)
    end
    return list
end

""" A list with type_func componentes distributed on their node """
function make_list(system::System, type_func, nodes=get_nodes(system))
    list = JuMP.Containers.DenseAxisArray(
        [[] for _ in 1:length(nodes)], get_name.(nodes)
    )
    for g in type_func(system)
        push!(list[g.bus.name], g)
    end
    return list
end

function get_in_list(type::Symbol, nodes::Vector{Int}, list::Vector)
    res = Int[]
    for n in nodes 
        push!(res, getproperty(getindex(list, n), type)...)
    end
    return res
end

find_in_model(mod::Model, ::ThermalGen, name::String) = mod[:pg0][name]
find_in_model(mod::Model, ::HydroGen, name::String) = mod[:pg0][name]
find_in_model(mod::Model, ::RenewableGen, name::String) = mod[:pr0][name]
find_in_model(mod::Model, ::StaticLoad, name::String) = mod[:ls0][name]
find_in_model(mod::Model, ::DCBranch, name::String) = mod[:pfdc0][name]

""" Return the (first) slack bus in the system. 
Return: slack bus number in nodes. The slack bus.
"""
function find_slack(nodes::AbstractVector{ACBus})
    for (i, x) in enumerate(nodes)
        x.bustype == ACBusTypes.REF && return i, x
    end
    @error "No slack bus found!"
    return 0, ACBus()
end
find_slack(sys::System) = find_slack(sort_components!(get_nodes(sys)))

""" Run optimizer to solve the model and check for optimality """
function solve_model!(model::Model)
    optimize!(model)
    if !has_values(model)
        @warn "Model not optimally solved with status $(termination_status(model))!"
    else
        @info @sprintf "Model solved in %.6fs with objective value %.10f" MOI.get(model, MOI.SolveTimeSec()) objective_value(model)
    end
    return model
end

" Sets the start values of symbols to the previous solution "
function set_warm_start!(model::JuMP.Model, symb::Symbol)
    x_sol = get_value(model, symb)
    JuMP.set_start_value.(model[symb], x_sol)
    return model
end

" Get a dict of bus-number-keys and vector-position-values "
get_nodes_idx(nodes::AbstractVector{ACBus}) = _make_ax_ref(get_number.(nodes))

""" Return a Dict of bus name number to index. Deprecated in PowerSystems """
function _make_ax_ref(ax::AbstractVector{T}) where {T}
    ref = Dict{T,Int}()
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
get_bus_idx(branch::Branch, idx::Dict{<:Int,T}) where {T<:Int} =
    (idx[branch.arc.from.number]::T, idx[branch.arc.to.number]::T)

get_bus_idx(val::StaticInjection, idx::Dict{<:Int,T}) where {T<:Int} =
    idx[val.bus.number]::T

" Get the bus number index of the from-bus and to-bus from all branches "
function get_bus_idx(branches::AbstractVector{<:Branch}, idx::Dict{<:Int,<:Int})
    return get_bus_idx.(branches, [idx])
end

function get_component_bus_idx(comp::Component, contingencies::AbstractVector{<:Branch}, idx::Dict{<:Int,<:Int})
    ix = indexin(contingencies, comp)
    return (ix, get_bus_idx(comp, idx))
end
get_component_bus_idx(comps::AbstractVector{<:Component}, contingencies::AbstractVector{<:Component}, idx::Dict{<:Int,<:Int}) =
    [(x, get_bus_idx(comps[x], idx)) for x in indexin(contingencies, comps)]

" Get the bus number index of the from-bus and to-bus from all branches "
function get_bus_idx(A::SparseArrays.SparseMatrixCSC{T}) where {T}
    m, n = size(A)
    ix = [(0, 0) for _ in 1:m]
    rows = SparseArrays.rowvals(A)
    for j = 1:n
        for i in SparseArrays.nzrange(A, j)
            row = rows[i]
            if iszero(ix[row][1])
                ix[row] = (j, 0)
            else
                ix[row] = (ix[row][1], j)
            end
        end
    end
    return ix
end

function sorted_missing(vals::AbstractVector{T}, maxval::Int) where {T<:Int}
    res = T[]
    i = 1
    for x in 1:maxval
        if get(vals, i, 0) == x
            i += 1
        else
            push!(res, x)
        end
    end
    return res
end

function not_insorted_nodes(from_bus::Integer, to_bus::Integer, nodes::AbstractVector{<:Integer})
    if insorted(from_bus, nodes)
        return 2
    elseif insorted(to_bus, nodes)
        return 1
    end
    throw(DivideError())
end

function zero_not_in_array!(array::AbstractArray{T}, subarray::AbstractVector{<:Integer}, n::Integer = ndims(array)
) where {T<:Real}
    i = firstindex(array, n)
    x = iterate(subarray)
    while lastindex(array, n) >= i
        if x === nothing 
            set_zero!.([array], i:lastindex(array, n), [Val(n)])
            break
        end
        val, state = x
        if i != val
            set_zero!(array, i, Val(n))
        elseif i == val
            x = iterate(subarray, state)
        else
            array .= zero(T)
            break
        end
        i += 1
    end
    return array
end

set_zero!(mx::AbstractMatrix{T}, i::Integer, ::Val{2}) where {T<:Real} = mx[:,i] .= zero(T)
set_zero!(mx::AbstractMatrix{T}, i::Integer, ::Val{1}) where {T<:Real} = mx[i,:] .= zero(T)
set_zero!(vx::AbstractVector{T}, i::Integer, ::Val{1}) where {T<:Real} = vx[i] = zero(T)

" Split a Vector{Pair} into a Pair{Vector}"
split_pair(val::AbstractVector) = map(first, val), map(last, val)

" zip a Pair{Vector, Vector} into a Vector{Pair}"
zip_pair(val::Pair{AbstractVector,AbstractVector}) = [(a, b) for (a, b) in zip(first(val), last(val))]
zip_pair(val1::AbstractVector, val2::AbstractVector) = zip_pair((val1, val2))

""" Return the overload of a line, else return 0.0 """
find_overload(flow::T, rate::Real, atol::Real=1e-6) where {T<:Real} =
    abs(flow) - rate > atol ? sign(flow) * (abs(flow) - rate) : zero(T)

filter_overload(flow::AbstractVector{<:Real}, linerating::AbstractVector{<:Real}, atol::Real=1e-6) =
    [(i, ol) for (i, ol) in enumerate(find_overload.(flow, linerating)) if abs(ol) > atol]

filter_overload(Δflow::AbstractVector{<:Tuple}, linerating::AbstractVector{<:Real}, atol::Real=1e-6) =
    [(i, find_overload(ol, linerating[i])) for (i, ol) in Δflow if abs(find_overload(ol, linerating[i])) > atol]

" Get dual value (upper or lower bound) from model reference "
get_low_dual(varref::VariableRef) = dual(LowerBoundRef(varref))
get_high_dual(varref::VariableRef) = dual(UpperBoundRef(varref))

""" This is a faster version of JuMP.value """
get_value(mod::Model, symb::Symbol) =
    MOI.get.([JuMP.backend(mod)], [MOI.VariablePrimal()], JuMP.index.(mod[symb]))
get_value(mod::Model, var::JuMP.VariableRef) =
    MOI.get(JuMP.backend(mod), MOI.VariablePrimal(), JuMP.index(var))
get_value(mod::Model, var::Vector{JuMP.VariableRef})::Vector{Float64} =
    MOI.get(JuMP.backend(mod), MOI.VariablePrimal(), JuMP.index.(var))

function get_variable_values(mod::Model)
    vals = Dict{Symbol, Vector{Float64}}()
    for obj in mod.obj_dict
        if typeof(last(obj)) <: Vector{VariableRef}
            vals[first(obj)] = get_value(mod, first(obj))
        end
    end
    return vals
end

function get_variable_values(mod::Model, Pc::Dict{<:Int,T}) where {T<:Union{ExprC, ExprCC, ExprCCX}}
    vals = Dict(x => Dict{Int, Vector{Float64}}() for x in fieldnames(T))
    for (c, cont) in Pc
        for symb in fieldnames(T)
            vals[symb][c] = get_value(mod, getfield(cont, symb))
        end
    end
    return vals
end

" Fix for name difference in StandardLoad "
PowerSystems.get_active_power(value::StandardLoad) = value.constant_active_power

""" Return the net power injected at each node. """
function get_net_Pᵢ(opf::OPFsystem, mod::Model, idx::Dict{<:Int,<:Int}=get_nodes_idx(opf.nodes))
    P = zeros(length(opf.nodes))
    vals = get_value(mod, :pr0)
    for (i, r) in enumerate(opf.renewables)
        P[idx[r.bus.number]] += get_active_power(r) - vals[i]
    end
    vals = get_value(mod, :ls0)
    for (i, d) in enumerate(opf.demands)
        P[idx[d.bus.number]] -= get_active_power(d) + vals[i]
    end
    vals = get_value(mod, :pg0)
    for (i, r) in enumerate(opf.ctrl_generation)
        P[idx[r.bus.number]] += vals[i]
    end
    vals = get_value(mod, :pfdc0)
    for (i, d) in enumerate(opf.dc_branches)
        P[idx[d.arc.from.number]] += vals[i]
        P[idx[d.arc.to.number]] -= vals[i]
    end
    # @assert abs(sum(Pᵢ)) < 0.001
    return P
end
function get_net(power_func::Function, opf::OPFsystem, idx::Dict{<:Int,<:Int})
    vals = zeros(length(opf.nodes))
    for r in opf.renewables
        vals[idx[r.bus.number]] += power_func(r)
    end
    for d in opf.demands
        vals[idx[d.bus.number]] -= power_func(d)
    end
    for r in opf.ctrl_generation
        vals[idx[r.bus.number]] += power_func(r)
    end
    for d in opf.dc_branches
        vals[idx[d.arc.from.number]] += power_func(d)
        vals[idx[d.arc.to.number]] -= power_func(d)
    end
    return vals
end
get_net_P(opf::OPFsystem, idx::Dict{<:Int,<:Int}=get_nodes_idx(opf.nodes)) = 
    get_net(get_active_power, opf, idx)
get_net_Q(opf::OPFsystem, idx::Dict{<:Int,<:Int}=get_nodes_idx(opf.nodes)) = 
    get_net(get_reactive_power, opf, idx)


function get_Pd!(P::Vector{<:Real}, opf::OPFsystem, idx::Dict{<:Int,<:Int})
    for r in opf.renewables
        P[idx[r.bus.number]] += get_active_power(r)
    end
    for d in opf.demands
        P[idx[d.bus.number]] -= get_active_power(d)
    end
    return P
end
function get_Pd(opf::OPFsystem, idx::Dict{<:Int,<:Int}=get_nodes_idx(opf.nodes))
    get_Pd!(zeros(length(opf.nodes)), opf, idx)
end

""" Return power shed at each node. """
function get_Pshed!(P::Vector{<:Real}, opf::OPFsystem, mod::Model, idx::Dict{<:Int,<:Int})
    vals = get_value(mod, :pr0)
    for (i, r) in enumerate(opf.renewables)
        P[idx[r.bus.number]] -= vals[i]
    end
    vals = get_value(mod, :ls0)
    for (i, d) in enumerate(opf.demands)
        P[idx[d.bus.number]] += vals[i]
    end
    return P
end
function get_Pshed(opf::OPFsystem, mod::Model, idx::Dict{<:Int,<:Int}=get_nodes_idx(opf.nodes))
    get_Pshed!(zeros(length(opf.nodes)), opf, mod, idx)
end

""" Return the power injected by controlled generation at each node. """
function get_Pgen!(P::Vector{<:Real}, opf::OPFsystem, mod::Model, idx::Dict{<:Int,<:Int})
    vals = get_value(mod, :pg0)
    for (i, r) in enumerate(opf.ctrl_generation)
        P[idx[r.bus.number]] += vals[i]
    end
    vals = get_value(mod, :pfdc0)
    for (i, d) in enumerate(opf.dc_branches)
        P[idx[d.arc.from.number]] += vals[i]
        P[idx[d.arc.to.number]] -= vals[i]
    end
    return P
end
function get_Pgen(opf::OPFsystem, mod::Model, idx::Dict{<:Int,<:Int}=get_nodes_idx(opf.nodes))
    get_Pgen!(zeros(length(opf.nodes)), opf, mod, idx)
end

function get_controllable(opf::OPFsystem, mod::Model, idx::Dict{<:Int,<:Int}=get_nodes_idx(opf.nodes))
    P = get_Pgen(opf, mod, idx)
    get_Pshed!(P, opf, mod, idx)
    return P
end

" Calculate the severity index for the system based on line loading "
function calc_severity(opf::OPFsystem, mod::Model, lim::Real=0.9)
    rate = make_named_array(get_rate, get_branches(opf.sys))
    sev = 0
    for c in 1:length(opf.contingencies)
        for l in get_name.(get_branches(opf.sys))
            sev += calc_line_severity(value(mod[:pfc][l, c]), rate[l], lim)
        end
    end
    return sev
end

calc_severity(
    rate::AbstractVector{<:Real},
    contingencies::AbstractVector,
    branches::AbstractVector,
    pfc,
    lim::Real=0.9
) =
    sum(calc_line_severity(pfc[l, c], rate[l], lim) for c in 1:length(contingencies), l in branches)

calc_severity(values::AbstractVector{<:Real}, rate::AbstractVector{<:Real}, lim::Real=0.9) =
    sum(calc_line_severity.(values, rate, [lim]))

" Calculate the severity index for a component based on the rating "
calc_line_severity(value::Real, rate::Real, lim::Real=0.9) =
    abs(value) > lim * rate ? (1 / (1 - lim)) * (abs(value) / rate - lim) : 0
