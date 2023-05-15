""" System """
mutable struct System
    name::String
    baseMVA::float

    nodes::Vector{Node}
    branches::Vector{Branch}
    generators::Vector{Generator}
    demands::Vector{Demand}
    renewables::Vector{Renewable}
    dc_branches::Vector{DCBranch}
end

function System()
    return System(name, baseMVA, nodes, branches, generators, demands, renewables, dc_branches)
end

mutable struct Arc
    from::Node
    to::Node
    R::Float
    X::Float
    B::Float
end

mutable struct Cost
    fixed::Float
    variable::Float
end

mutable struct MinMax
    min::Float
    max::Float
end

mutable struct Phasor
    r::Real
    ϕ::Real
end

Phasor(val::Complex) = Phasor(abs(val), angle(val))

rec(r::Real, ϕ::Real) = r*cos(ϕ) + 1im*r*sin(ϕ)

∠(ϕ::Real) = cos(ϕ) + 1im*sin(ϕ)

""" Base component type """
mutable struct Component
    number::Integer
    id::string
    mttf::Float
    mttr::Float
end

mutable struct Node <: Component
    type::NodeType
    BaseKv::Float
    voltage_setpoint::Phasor
    G::Float
    B::Float
end

mutable struct Branch <: Component
    arc::Arc
    continous_rating::Float
    long_term_rating::Float
    short_term_rating::Float
end

mutable struct ACBranch <: Branch
    
end

mutable struct Line <: ACBranch
    
end

mutable struct Transformer <: ACBranch
    
end

mutable struct DCBranch <: Branch
    ramp_rate::MinMax
    cost::Cost
end

mutable struct PowerDevice <: Component
    node::Node
    voltage_magnitude_setpoint::Float
    active_power_limits::MinMax
    reactive_power_limits::MinMax
    ramp_rate::MinMax
    cost::Cost
end

mutable struct Generator <: PowerDevice
    
end

mutable struct Demand <: Component
    
end

mutable struct Renewable <: Component
    
end
