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
    contingencies = ["branch" => [i] for i in 1:length(get_branches(system))]
    prob = make_prob(contingencies, prob_min, prob_max)
    return voll, prob, contingencies
end

function fix_generation_cost!(sys::System)
    # set_operation_cost!.(get_renewables(sys), make_cost_renewables(sys))
    cost = [get_generator_cost(g)[2] for g in get_gens_t(sys)]
    c_mean = mean(cost)
    set_operation_cost!.(get_generation(sys),
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

function set_renewable_prod!(system::System, t::Dates.DateTime)
    for g in get_renewables(system)
        v = get_time_series_values(SingleTimeSeries, g, "max_active_power", 
            start_time = t, len=1)
        set_active_power!(g, get_max_active_power(g) * v)
    end
end

""" Set the active power demand """
function set_active_power_demand!(system::System, demands=get_max_active_power.(get_demands(system)))
    set_active_power!.(get_demands(system), demands)
end

PowerSystems.set_active_power!(dem::StandardLoad, val::Real) =
    PowerSystems.set_constant_active_power!(dem, val)

" Fix for name difference in StandardLoad "
PowerSystems.get_active_power(value::StandardLoad) = value.constant_active_power

const fueldict = Dict(
    ThermalFuels.COAL => 0.02, 
    ThermalFuels.DISTILLATE_FUEL_OIL => 0.03, 
    ThermalFuels.NATURAL_GAS => 0.04, 
    ThermalFuels.NUCLEAR => 0.005, 
    ThermalFuels.OTHER => 0.01
)

""" Set the generator ramping ability p.u./min based on ramp (or fuel type) and maximal active power """
function set_ramp_limit!(gen::Generator, ramp::Real)
    p_lim = get_active_power_limits(gen)
    PowerSystems.set_ramp_limits!(gen, (up=(p_lim.max * ramp), down=(p_lim.max * ramp)))
end

set_ramp_limit!(gen::ThermalGen) = set_ramp_limit!(gen, fueldict[get_fuel(gen)])
set_ramp_limit!(gen::HydroDispatch) = set_ramp_limit!(gen, 0.2)
set_ramp_limit!(gen::RenewableGen) = nothing # set_ramp_limit!(gen, 1.0)

set_ramp_limits!(system::System) = set_ramp_limit!.(get_generation(system))

PowerSystems.get_active_power_limits(gen::RenewableGen) = (min=get_active_power(gen), max=get_active_power(gen))
PowerSystems.get_ramp_limits(gen::RenewableGen) = get_active_power_limits(gen)

get_generator_cost(gen::Generator) = get_operation_cost(gen) |> get_variable |> get_cost |> _get_g_value
_get_g_value(x::AbstractVector{<:Tuple{Real,Real}}) = mean([y[1] for y in x if y[1] != 0.0]), mean([y[2] for y in x if y[2] != 0.0])
_get_g_value(x::Tuple{<:Real,<:Real}) = x

function get_generator_cost(generation::AbstractVector{<:Generator}, ramp_mult::Real, renew_cost::Real, renew_ramp::Real)
    cost = Vector{NamedTuple{(:fix, :var, :ramp)}}(undef, length(generation))
    for (i,g) in enumerate(generation)
        if typeof(g) <: RenewableGen
            cost[i] = (fix=0.0, var=renew_cost, ramp=renew_ramp)
        else
            c = get_generator_cost(g)
            cost[i] = (fix=c[1], var=c[2], ramp=c[2]*ramp_mult)
        end
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
    m = direct_model(optimizer)
    set_string_names_on_creation(m, debug)
    MOI.set(m, MOI.Silent(), silent) # supress output from the solver
    set_time_limit_sec(m, time_limit_sec)
    return m
end

function check_values(val::Real, var::Component, typename::String; comp=<, atol=1e-5)
    if comp(val, atol)
    # if isapprox(val, zero(typeof(val)); atol=atol)
        @info "$(PowerSystems.get_name(var)) has $val value in $typename."
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

function nan_to_zero!(vals::AbstractArray{T}) where {T<:Real}
    @inbounds for i in eachindex(vals)
        if isequal(vals[i], NaN)
            vals[i] = zero(T)
        end
    end
    return vals
end

"""Find value of β, β = 1 if from-bus is the bus, β = -1 if to-bus is the bus, 0 else"""
beta(bus::ACBus, branch::Branch) = bus == branch.arc.from ? 1 : bus == branch.arc.to ? -1 : 0
beta(bus::String, branch::Branch) = bus == branch.arc.from.name ? 1 : bus == branch.arc.to.name ? -1 : 0
function beta(sys::System, branch::String)
    l = get_component(ACBranch, sys, branch)
    return [beta(b, l) for b in get_nodes(sys)]
end

beta(val::Number) = ifelse(val < 0.0, -1.0, 1.0)

a(line::String, contingency::String) = ifelse(line != contingency, 1, 0)

""" An iterator to a type of power system component """
get_gens_t(sys::System) = get_components(ThermalGen, sys) |> collect
get_gens_h(sys::System) = get_components(HydroGen, sys) |> collect
get_lines(sys::System) = get_components(Line, sys) |> collect 
get_transformers2w(sys::System) = get_components(Transformer2W, sys) |> collect
get_taptransformers(sys::System) = get_components(TapTransformer, sys) |> collect
get_arcs(sys::System) = get_components(Arc, sys) |> collect
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
    by=x -> (x.arc.from.number, x.arc.to.number)
)
sort_components!(arcs::Vector{<:Arc}) = sort!(arcs, by=x -> (x.from.number, x.to.number))

filter_active!(list::AbstractVector{<:Component}) = filter(comp -> comp.available, list)

""" OPF system type """
mutable struct OPFsystem{TR<:Real,TI<:Integer}
    cost_gen::Vector{NamedTuple{(:fix, :var, :ramp), Tuple{TR, TR, TR}}}
    cost_renewables::Vector{TR}
    voll::Vector{TR}
    prob::Vector{TR}

    idx::Dict{TI,TI}
    mgx::SparseArrays.SparseMatrixCSC{Int8, TI}
    mdx::SparseArrays.SparseMatrixCSC{Int8, TI}
    mbx::SparseArrays.SparseMatrixCSC{Int8, TI}
    mdcx::SparseArrays.SparseMatrixCSC{Int8, TI}

    contingencies::Vector{Pair{String, Vector{TI}}}
end

""" Constructor for OPFsystem """
function opfsystem(sys::System, voll::Vector{TR}, contingencies::Vector{Pair{String, Vector{Int64}}}=Pair{String, Int64}[], 
    prob::Vector{TR}=[]; ramp_mult::Real=1.0, renew_cost::Real=0.0, renew_ramp::Real=0.0, check=false
) where {TR<:Real}

    generation = sort_components!(get_generation(sys))
    branches = sort_components!(get_branches(sys))
    # arcs = sort_components!(get_arcs(sys))
    dc_branches = sort_components!(get_dc_branches(sys))
    nodes = sort_components!(get_nodes(sys))
    demands = sort_components!(get_demands(sys))

    idx = get_nodes_idx(nodes)
    mgx = calc_connectivity(generation, length(nodes), idx)
    mbx = calc_A(branches, length(nodes), idx)
    mdcx = calc_A(dc_branches, length(nodes), idx)
    mdx = calc_connectivity(demands, length(nodes), idx)

    cost_gen = get_generator_cost(generation, ramp_mult, renew_cost, renew_ramp)
    cost_renewables = Vector{TR}() # Vector{Float64}([get_generator_cost(g)[2] for g in renewables])

    if check
        check_values(getfield.(cost_gen, :fix), generation, "fixedcost")
        check_values(getfield.(cost_gen, :var), generation, "varcost")
        check_values(voll, demands, "voll")
        check_values.(generation)
        check_values.(branches)
        check_values.(dc_branches)
    end

    A = calc_B(branches, length(nodes), idx)
    islands = island_detection(A'A)
    if length(islands) > 1
        @error "The system is separated into islands" islands
    end

    return OPFsystem{TR, Int64}(cost_gen, cost_renewables, voll, prob, idx, mgx, mdx, mbx, mdcx, contingencies)
end

""" Automatic constructor for OPFsystem where voll, prob, and contingencies are automatically computed. """
function opfsystem(sys::System)
    voll, prob, contingencies = setup(sys, 1, 4)
    fix_generation_cost!(sys)
    return opfsystem(sys, voll, contingencies, prob)
end

get_cost_gen(opf::OPFsystem) = opf.cost_gen
get_cost_renewables(opf::OPFsystem) = opf.cost_renewables
get_voll(opf::OPFsystem) = opf.voll
get_prob(opf::OPFsystem) = opf.prob

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

function get_area_branches(branches::AbstractVector{T}, area::Area) where {T<:Branch}
    vals = Vector{T}()
    for br in branches
        from = get_area(get_from(get_arc(br)))
        to = get_area(get_to(get_arc(br)))
        if from == area || to == area
            push!(vals, br)
        end
    end
    return vals
end

function get_area_nodes(nodes::AbstractVector{T}) where {T<:Bus}
    vals = Dict{Area, Vector{T}}()
    for n in nodes
        area = get_area(n)
        if get(vals, area, 0) == 0
            vals[area] = T[]
        end
        push!(vals[area], n)
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
    i === nothing && error(string(PowerSystems.get_name(val)) * " is not found")
    return i
end

typesort_component(val::Generator, opf::OPFsystem) =
    (find_component(val, get_generation(opf)), get_bus_idx(val, opf.idx))

typesort_component(val::StaticLoad, opf::OPFsystem) =
    (find_component(val, get_demands(opf)), get_bus_idx(val, opf.idx))

typesort_component(val::ACBus, opf::OPFsystem) =
    (find_component(val, get_nodes(opf)), get_bus_idx(val, opf.idx))

typesort_component(val::ACBranch, opf::OPFsystem) =
    (find_component(val, get_branches(opf)), get_bus_idx(val, opf.idx))

function typesort_component(val::Pair{String, <:Any}, opf::OPFsystem)
    if first(val) == "branch"
        return (last(val), opf.mbx[last(val),:].nzval)
    elseif first(val) == "gen"
        return (last(val), opf.mbx[last(val),:].nzval)
    else
        return (-1,-1)
    end
end

""" Make a DenseAxisArray using the list and function for the value of each element """
make_named_array(value_func, list) = JuMP.Containers.DenseAxisArray(
    [value_func(x) for x in list], PowerSystems.get_name.(list)
)

find_in_model(m::Model, ::ThermalGen, name::String) = m[:pg0][name]
find_in_model(m::Model, ::HydroGen, name::String) = m[:pg0][name]
find_in_model(m::Model, ::StaticLoad, name::String) = m[:ls0][name]
find_in_model(m::Model, ::DCBranch, name::String) = m[:pfdc0][name]

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

function find_slack(nodes::AbstractVector{<:Integer}, slack::Integer, pg_lim_max::AbstractVector{<:Real}, mgx::AbstractMatrix{<:Real})
    s = searchsorted(nodes, slack)
    if length(s) == 0
        gens = reduce(vcat, [mgx[:,i].nzind for i in nodes]) # generators in island
        _, sg = findmax(pg_lim_max[gens]) # slack generator (biggest in island)
        si = mgx[gens[sg], :].nzind[1] # slack node
        s = searchsorted(nodes, si)
    end
    return first(s)
end

""" Run optimizer to solve the model and check for optimality """
function solve_model!(model::Model)
    optimize!(model)
    if !is_solved_and_feasible(model)
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

uniqueidx(v::AbstractVector{<:Branch}) = unique(i -> v[i].arc, eachindex(v))
uniqueidx(v) = unique(i -> v[i], eachindex(v))

" Get the bus number of the from-bus and to-bus from a branch "
get_bus_id(branch::Branch) = (branch.arc.from.number, branch.arc.to.number)

" Get the bus number index of the from-bus and to-bus from a branch "
get_bus_idx(branch::Branch, idx::Dict{<:Int,<:Int}) = get_bus_idx(branch.arc, idx)
get_bus_idx(arc::Arc, idx::Dict{<:Int,T}) where {T<:Int} =
    (idx[arc.from.number]::T, idx[arc.to.number]::T)
get_bus_idx(branch::SparseArrays.SparseVector{<:Integer, T}) where {T<:Integer} = 
    first(branch.nzval) < 0 ? (last(branch.nzind), first(branch.nzind)) : (first(branch.nzind), last(branch.nzind))

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

function zero_not_in_array!(array::AbstractArray{T}, subarray::AbstractVector{<:Integer}, nval::Val
) where {T<:Real}
    n = ndims(array)
    i = firstindex(array, n)
    x = iterate(subarray)
    while lastindex(array, n) >= i
        if x === nothing 
            set_zero!.([array], i:lastindex(array, n), [nval])
            break
        end
        val, state = x
        if i != val
            set_zero!(array, i, nval)
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
zero_not_in_array!(array::AbstractArray{<:Real}, subarray::AbstractVector{<:Integer}) = 
    zero_not_in_array!(array, subarray, Val(ndims(array)))

set_zero!(mx::AbstractMatrix{T}, i::Integer, ::Val{2}) where {T<:Real} = mx[:,i] .= zero(T)
set_zero!(mx::AbstractMatrix{T}, i::Integer, ::Val{1}) where {T<:Real} = mx[i,:] .= zero(T)
set_zero!(vx::AbstractVector{T}, i::Integer, ::Val{1}) where {T<:Real} = vx[i] = zero(T)

set_tol_zero!(vals::AbstractArray{T}, atol::Float64 = 1e-14) where {T<:Real} = 
    @. vals[abs(vals) < atol && vals != zero(T)] = zero(T)

" Split a Vector{Pair} into a Pair{Vector}"
split_pair(val::AbstractVector) = map(first, val), map(last, val)

" zip a Pair{Vector, Vector} into a Vector{Pair}"
zip_pair(val::Tuple{AbstractVector,AbstractVector}) = [(a, b) for (a, b) in zip(first(val), last(val))]
zip_pair(val1::AbstractVector, val2::AbstractVector) = zip_pair((val1, val2))

""" Return the overload of a line, else return 0.0 """
find_overload(flow::T, rate::Real, atol::Real=1e-6) where {T<:Real} =
    abs(flow) - rate > atol ? sign(flow) * (abs(flow) - rate) : zero(T)

filter_overload(flow::AbstractVector{<:Real}, linerating::AbstractVector{<:Real}, atol::Real=1e-6) =
    [(i, ol) for (i, ol) in enumerate(find_overload.(flow, linerating)) if abs(ol) > atol]

filter_overload(Δflow::AbstractVector{<:Tuple}, linerating::AbstractVector{<:Real}, atol::Real=1e-6) =
    [(i, find_overload(ol, linerating[i])) for (i, ol) in Δflow if abs(find_overload(ol, linerating[i])) > atol]

find_overloaded_branches(flow::AbstractVector{<:Real}, linerating::AbstractVector{<:Real}, atol::Real=1e-6) =
    findall(abs.(flow) .- linerating .> atol)

" Get dual value (upper or lower bound) from model reference "
get_low_dual(varref::VariableRef) = dual(LowerBoundRef(varref))
get_high_dual(varref::VariableRef) = dual(UpperBoundRef(varref))

""" This is a faster version of JuMP.value """
get_value(m::Model, symb::Symbol) =
    MOI.get.([JuMP.backend(m)], [MOI.VariablePrimal()], JuMP.index.(m[symb]))
get_value(m::Model, var::JuMP.VariableRef) =
    MOI.get(JuMP.backend(m), MOI.VariablePrimal(), JuMP.index(var))
get_value(m::Model, var::Vector{JuMP.VariableRef})::Vector{Float64} =
    MOI.get(JuMP.backend(m), MOI.VariablePrimal(), JuMP.index.(var))

function get_variable_values(m::Model)
    vals = Dict{Symbol, Vector{Float64}}()
    for obj in m.obj_dict
        if typeof(last(obj)) <: Vector{VariableRef}
            vals[first(obj)] = get_value(m, first(obj))
        end
    end
    return vals
end

function get_variable_values(m::Model, Pc::Dict{<:Int,T}) where {T<:Union{ExprC, ExprCC, ExprCCX}}
    vals = Dict(x => Dict{Int, Vector{Float64}}() for x in fieldnames(T))
    for (c, cont) in Pc
        for symb in fieldnames(T)
            vals[symb][c] = get_value(m, getfield(cont, symb))
        end
    end
    return vals
end

""" Return the net power injected at each node. """
get_net_Pᵢ(m::Model) = get_value(m, :p0)

""" Return active power from demands at each node. """
function get_Pd!(P::Vector{<:Real}, opf::OPFsystem, m::Model)
    P .-= opf.mdx' * get_value(m, :pd)
    return P
end
get_Pd(opf::OPFsystem, m::Model) = get_Pd!(zeros(size(opf.mgx, 2)), opf, m)

""" Return power shed from demands at each node. """
function get_Pshed!(P::Vector{<:Real}, opf::OPFsystem, m::Model)
    P .+= opf.mdx' * get_value(m, :ls0)
    return P
end
get_Pshed(opf::OPFsystem, m::Model) = get_Pshed!(zeros(size(opf.mgx, 2)), opf, m)

""" Return the power injected by controlled generation at each node. """
function get_Pgen!(P::Vector{<:Real}, opf::OPFsystem, m::Model)
    P .+= opf.mgx' * get_value(m, :pg0) + opf.mdcx' * get_value(m, :pfdc0)
    return P
end
get_Pgen(opf::OPFsystem, m::Model) = get_Pgen!(zeros(size(opf.mgx, 2)), opf, m)

" Return the controlled generation and power shedding at each node. "
function get_controllable(opf::OPFsystem, m::Model)
    P = get_Pgen(opf, m)
    get_Pshed!(P, opf, m)
    return P
end

" Calculate the severity index for the system based on line loading "
function calc_severity(opf::OPFsystem, m::Model, lim::Real=0.9)
    rate = make_named_array(get_rate, get_branches(opf.sys))
    sev = 0
    for c in 1:length(opf.contingencies)
        for l in PowerSystems.get_name.(get_branches(opf.sys))
            sev += calc_line_severity(value(m[:pfc][l, c]), rate[l], lim)
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

" An AbstractJuMPScalar nicely formatted to a string "
sprint_expr(expr::AbstractJuMPScalar, lim=1e-14) =
    join(Printf.@sprintf("%s%5.2f %s ", (x[2] > 0 ? "+" : "-"), abs(x[2]), x[1])
            for x in expr.terms if abs(x[2]) > lim) *
    Printf.@sprintf("<= %s%.2f", (expr.constant > 0 ? "-" : " "), abs(expr.constant)
    )

function print_cuts(iterations::Integer, type::OPF, pre::Integer, corr1::Integer, corr2::Integer, corr2f::Integer)
    print(iterations, " iterations with cuts added:")
    type.P && print(" pre=", pre)
    type.C1 && print(" corr1=", corr1)
    type.C2 && print(" corr2=", corr2)
    type.C2F && print(" corr2f=", corr2f)
    println("")
end

function countmemb(itr::AbstractVector{T}) where {T<:Real}
    d = Dict{T, Int}()
    for val in itr
        d[val] = get(d, val, 0) + 1
    end
    return d
end

function sumvals(itr::AbstractVector{T1}, vals::AbstractVector{T2}) where {T1<:Real, T2<:Real}
    d = Dict{T1, T2}()
    for (i,val) in enumerate(itr)
        d[val] = get(d, val, 0) + vals[i]
    end
    return d
end