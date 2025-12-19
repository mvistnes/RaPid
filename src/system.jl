""" 
    OPF system type 
    Contains `all` data needed for the SCOPF model
"""
mutable struct OPFsystem{TR<:Real,TI<:Integer}
    cost_ctrl_gen::Vector{NamedTuple{(:fix, :var, :ramp), Tuple{TR, TR, TR}}}
    cost_renewables::Vector{TR}
    voll::Vector{TR}
    prob::Vector{TR}

    idx::Dict{TI,TI}
    mgx::SparseArrays.SparseMatrixCSC{Int8, TI}
    mdx::SparseArrays.SparseMatrixCSC{Int8, TI}
    mbx::SparseArrays.SparseMatrixCSC{Int8, TI}
    mdcx::SparseArrays.SparseMatrixCSC{Int8, TI}
    mrx::SparseArrays.SparseMatrixCSC{Int8, TI}

    ctrl_generation::Vector{Generator}
    branches::Vector{ACBranch}
    arcs::Vector{Arc}
    dc_branches::Vector{TwoTerminalHVDCLine}
    nodes::Vector{ACBus}
    demands::Vector{StaticLoad}
    renewables::Vector{RenewableGen}
    contingencies::Vector{Component}
end

""" Constructor for OPFsystem """
function opfsystem(sys::System, voll::Vector{TR}, contingencies::Vector{<:Component}=Component[], prob::Vector{TR}=[],
    ramp_mult::Real=1.0; check=false
) where {TR<:Real}

    ctrl_generation = sort_components!(get_ctrl_generation(sys))
    branches = sort_components!(get_branches(sys))
    arcs = sort_components!(get_arcs(sys))
    dc_branches = sort_components!(get_dc_branches(sys))
    nodes = sort_components!(get_nodes(sys))
    demands = sort_components!(get_demands(sys))
    renewables = sort_components!(get_renewables(sys))

    idx = get_nodes_idx(nodes)
    mgx = calc_connectivity(ctrl_generation, length(nodes), idx)
    mbx = calc_A(branches, length(nodes), idx)
    mdcx = calc_A(dc_branches, length(nodes), idx)
    mdx = calc_connectivity(demands, length(nodes), idx)
    mrx = calc_connectivity(renewables, length(nodes), idx)

    cost_ctrl_gen = get_generator_cost(ctrl_generation, ramp_mult)
    cost_renewables = Vector{TR}() # Vector{Float64}([get_generator_cost(g)[2] for g in renewables])

    if check
        check_values(getfield.(cost_ctrl_gen, :fix), ctrl_generation, "fixedcost")
        check_values(getfield.(cost_ctrl_gen, :var), ctrl_generation, "varcost")
        check_values(getfield.(cost_ctrl_gen, :ramp), ctrl_generation, "rampcost")
        check_values(voll, demands, "voll")
        check_values.(ctrl_generation)
        check_values.(branches)
        check_values.(dc_branches)
    end

    A = calc_B(branches, length(nodes), idx)
    islands = island_detection(A'A)
    if length(islands) > 1
        @error "The system is separated into islands" islands
    end

    # sump = sum(prob)
    return OPFsystem{TR, Int}(cost_ctrl_gen, cost_renewables, voll, prob, idx, mgx, mdx, mbx, mdcx, mrx,
        ctrl_generation, branches, arcs, dc_branches, nodes, demands, renewables, contingencies)
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

typesort_component(val::Generator, opf::OPFsystem) =
    (find_component(val, get_ctrl_generation(opf)), get_bus_idx(val, opf.idx))

typesort_component(val::StaticLoad, opf::OPFsystem) =
    (find_component(val, get_demands(opf)), get_bus_idx(val, opf.idx))

typesort_component(val::ACBus, opf::OPFsystem) =
    (find_component(val, get_nodes(opf)), get_bus_idx(val, opf.idx))

typesort_component(val::ACBranch, opf::OPFsystem) =
    (find_component(val, get_branches(opf)), get_bus_idx(val, opf.idx))

""" Make a DenseAxisArray using the list and function for the value of each element """
make_named_array(value_func, list) = JuMP.Containers.DenseAxisArray(
    [value_func(x) for x in list], PowerSystems.get_name.(list)
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
function make_list(opf::OPFsystem, nodes::AbstractVector{ACBus})
    list = [CTypes(n, [Int[] for _ in 1:5]...) for n in nodes]
    for (i, r) in enumerate(opf.ctrl_generation)
        push!(list[opf.idx[r.bus.number]].ctrl_generation, i)
    end
    for (i, r) in enumerate(opf.renewables)
        push!(list[opf.idx[r.bus.number]].renewables, i)
    end
    for (i, d) in enumerate(opf.demands)
        push!(list[opf.idx[d.bus.number]].demands, i)
    end
    for (i, d) in enumerate(opf.branches)
        push!(list[opf.idx[d.arc.from.number]].branches, i)
        push!(list[opf.idx[d.arc.to.number]].branches, i)
    end
    for (i, d) in enumerate(opf.dc_branches)
        push!(list[opf.idx[d.arc.from.number]].dc_branches, i)
        push!(list[opf.idx[d.arc.to.number]].dc_branches, i)
    end
    return list
end

function get_net(power_func::Function, opf::OPFsystem)
    vals = zeros(length(opf.nodes))
    for r in opf.renewables
        vals[opf.idx[r.bus.number]] += power_func(r)
    end
    for d in opf.demands
        vals[opf.idx[d.bus.number]] -= power_func(d)
    end
    for r in opf.ctrl_generation
        vals[opf.idx[r.bus.number]] += power_func(r)
    end
    for d in opf.dc_branches
        vals[opf.idx[d.arc.from.number]] += power_func(d)
        vals[opf.idx[d.arc.to.number]] -= power_func(d)
    end
    return vals
end
get_net_P(opf::OPFsystem) = get_net(get_active_power, opf)
get_net_Q(opf::OPFsystem) = get_net(get_reactive_power, opf)

""" Return active power from renweables and demands at each node. """
function get_Pd!(P::Vector{<:Real}, opf::OPFsystem, m::Model)
    P .+= opf.mrx' * get_value(m, :pr) - opf.mdx' * get_value(m, :pd)
    return P
end
get_Pd(opf::OPFsystem, m::Model) = get_Pd!(zeros(length(opf.nodes)), opf, m)

""" Return power shed from renewables and demands at each node. """
function get_Pshed!(P::Vector{<:Real}, opf::OPFsystem, m::Model)
    P .+= opf.mdx' * get_value(m, :ls0) - opf.mrx' * get_value(m, :pr0)
    return P
end
get_Pshed(opf::OPFsystem, m::Model) = get_Pshed!(zeros(length(opf.nodes)), opf, m)

""" Return the power injected by controlled generation at each node. """
function get_Pgen!(P::Vector{<:Real}, opf::OPFsystem, m::Model)
    P .+= opf.mgx' * get_value(m, :pg0) + opf.mdcx' * get_value(m, :pfdc0)
    return P
end
get_Pgen(opf::OPFsystem, m::Model) = get_Pgen!(zeros(length(opf.nodes)), opf, m)

" Return the controlled generation and power shedding at each node. "
function get_controllable(opf::OPFsystem, m::Model)
    P = get_Pgen(opf, m)
    get_Pshed!(P, opf, m)
    return P
end

" Calculate the severity index for the system based on line loading "
function calc_severity(opf::OPFsystem, m::Model, lim::Real=0.9)
    rate = make_named_array(get_rating, get_branches(opf.sys))
    sev = 0
    for c in 1:length(opf.contingencies)
        for l in PowerSystems.get_name.(get_branches(opf.sys))
            sev += calc_line_severity(value(m[:pfc][l, c]), rate[l], lim)
        end
    end
    return sev
end
