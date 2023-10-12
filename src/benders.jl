# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022 

""" Benders type """
mutable struct Benders{TR<:Real,TI<:Integer}
    idx::Dict{TI,TI}
    list::Vector{CTypes{TI}}

    obj::AbstractJuMPScalar

    Pᵢ::Vector{TR}
    Pg::Vector{TR}
    Pd::Vector{TR}
end

""" Constructor for Benders type """
function benders(opf::OPFsystem, mod::Model)
    idx = get_nodes_idx(opf.nodes)
    list = make_list(opf, idx, opf.nodes)
    # cgen = connectivitymatrix(system, length(nodes), idx)

    obj = objective_function(mod)
    # Pᵢ = get_net_Pᵢ(opf, idx)
    Pᵢ = get_value(mod, :p0)
    Pg = get_controllable(opf, mod, idx)
    Pd = get_Pd(opf, idx) # Fixed injection
    @assert isapprox(Pg, (Pᵢ - Pd); atol=1e-8) string(Pg - (Pᵢ - Pd))
    return Benders(idx, list, obj, Pᵢ, Pg, Pd)
end

function run_benders(
    type::OPF, 
    system::System, 
    optimizer, 
    voll::Vector{<:Number}, 
    prob::Vector{Float64}, 
    contingencies::Vector{ACBranch};
    time_limit_sec::Int64=600, 
    ramp_minutes=10, 
    ramp_mult=10, 
    max_shed=1.0, 
    max_curtail=1.0, 
    lim=1e-6,
    short_term_multi::Real=1.5, 
    long_term_multi::Real=1.2, 
    max_itr=length(contingencies),
    p_failure=0.0, 
    branch_c=nothing, 
    rate_c=0.0, 
    silent=true,
    debug=false
)
    # LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
    # opf = scopf(SC, system, optimizer, voll=voll, prob=prob, contingencies=contingencies, ramp_minutes=ramp_minutes, max_shed=max_shed, 
    #     renewable_prod = 1.0, debug = debug)
    mod, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SC, system, optimizer, voll=voll, contingencies=contingencies, prob=prob,
        time_limit_sec=time_limit_sec, ramp_minutes=ramp_minutes, ramp_mult=ramp_mult, max_shed=max_shed, max_curtail=max_curtail,
        short_term_multi=short_term_multi, long_term_multi=long_term_multi, p_failure=p_failure, silent=silent, debug=debug)
    return run_benders!(SC, type, mod, opf, pf, oplim, Pc, Pcc, Pccx, lim, max_itr, branch_c, rate_c, debug)
end

""" 
Solve the optimization model using Benders decomposition.

Function creates extra production constraints on generators in the system
based on the power transfer distribution factors of the generators in the
system. 
"""
function run_benders!(
    basetype::OPF, 
    type::OPF, 
    mod::Model, 
    opf::OPFsystem, 
    pf::DCPowerFlow, 
    oplim::Oplimits, 
    Pc::Dict{<:Integer, Main.SCOPF.ExprC}, 
    Pcc::Dict{<:Integer, Main.SCOPF.ExprCC}, 
    Pccx::Dict{<:Integer, Main.SCOPF.ExprCCX},
    lim=1e-6,
    max_itr=length(opf.contingencies),
    branch_c=nothing, 
    rate_c=0.0, 
    debug=false
)
    @assert isempty(Pc) || basetype == PSC::OPF
    @assert (isempty(Pcc) && isempty(Pccx)) || type < PCSC::OPF
    # LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
    # if !isnothing(branch_c) 
    #     @info "Flow constraint on branch $branch_c is $rate_c."
    #     JuMP.@constraint(mod, mod[:pf0][findfirst(x -> x == branch_c, opf.branches)] == rate_c)
    # end
    solve_model!(mod)
    total_solve_time = solve_time(mod)
    termination_status(mod) != MOI.OPTIMAL && return mod, opf, pf, oplim, Pc, Pcc, Pccx, total_solve_time
    @debug "lower_bound = $(objective_value(mod))"

    # Set variables
    bd = benders(opf, mod)
    calc_θ!(pf, bd.Pᵢ)
    calc_Pline!(pf)
    overloads = zeros(length(opf.contingencies))

    contids = get_branch_bus_idx(opf.branches, opf.contingencies, bd.idx)
    flow = similar(pf.F)
    θ = similar(pf.θ)
    # bufs_F = [similar(pf.F) for _ in 1:Threads.nthreads()]
    # bufs_θ = [similar(pf.θ) for _ in 1:Threads.nthreads()]
    # num_islands = is_islanded.([pf], getindex.(contids, 2), getindex.(contids, 1))
    # if sum(num_islands) > 0.1 * length(opf.contingencies) 
    # LinearAlgebra.BLAS.set_num_threads(1)
    # mod_lock = Threads.ReentrantLock()
    for contid in contids
        # Threads.@threads for contid in contids
        # flow = bufs_F[Threads.threadid()]
        # θ = bufs_θ[Threads.threadid()]
        (c, cont) = contid
        # cont  = get_bus_idx(opf.contingencies[c], bd.idx)
        if !is_islanded(pf, cont, c)
            calculate_line_flows!(flow, pf, cont, c)
        else
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont, c, pf.slack)
            # Threads.lock(mod_lock) do
            basetype == SC::OPF && init_P!(Pc, opf, oplim, mod, bd.obj, bd.list, islands, island, c)
            type >= PCSC::OPF && init_P!(Pcc, opf, oplim, mod, bd.obj, bd.list, islands, island, c)
            type >= PCFSC::OPF && init_P!(Pccx, opf, oplim, mod, bd.obj, bd.list, islands, island, c)
            # end
            @debug "Island: Contingency line $(cont[1])-$(cont[2])-i_$(c)"
            if isempty(islands)
                fill!(flow, 0.0)
            else
                calculate_line_flows!(flow, pf, cont, c, nodes=islands[island], branches=island_b)
            end
        end

        # Calculate the power flow with the new outage and find if there are any overloads
        overload = filter_overload(flow, oplim.branch_rating * oplim.short_term_multi)

        if !isempty(overload)
            overloads[c] = maximum(x -> abs(x[2]), overload)
        end
    end

    total_solve_time = update_model!(mod, pf, bd, total_solve_time)
    termination_status(mod) != MOI.OPTIMAL && return mod, opf, pf, oplim, Pc, Pcc, Pccx, total_solve_time
    GC.safepoint()
    # end
    # LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())

    ΔPc = zeros(length(bd.Pg))
    ΔPcc = zeros(length(bd.Pg))
    ΔPccx = zeros(length(bd.Pg))
    pre = 0
    corr = 0

    ptdf = copy(pf.ϕ)
    X = copy(pf.X)
    B = copy(pf.B)
    θ = copy(pf.θ)
    flow = copy(pf.F)
    olc = Vector{Tuple{Int,Float64}}[]
    olcc = Vector{Tuple{Int,Float64}}[]
    olccx = Vector{Tuple{Int,Float64}}[]
    islands = Vector{Vector{Int}}[]
    island = 0
    island_b = Int[]

    permutation = sortperm(overloads, rev=true)
    cut_added = 1
    for iterations in 1:max_itr
        if cut_added == 0 # loops until no new cuts are added for the contingencies
            # @printf "\nEND: Total solve time %.4f.\n" total_solve_time
            return mod, opf, pf, oplim, Pc, Pcc, Pccx, total_solve_time
        end
        cut_added = 0
        @info "Iteration $iterations"

        for contid in contids[permutation]
            (c, cont) = contid
            empty!(olc)
            empty!(olcc)
            empty!(olccx)
            if is_islanded(pf, cont, c)
                islands, island, island_b = handle_islands(pf.B, pf.DA, cont, c, pf.slack)
            else
                empty!(islands)
                island = 0
                empty!(island_b)
            end

            olc = basetype == SC::OPF ? calculate_contingency_overload!(ΔPc, flow, (oplim.branch_rating * oplim.short_term_multi), θ, B, Pc, opf, mod, pf, bd, cont, c, islands, island, island_b) : Tuple{Int64,Float64}[]
            olcc = type >= PCSC::OPF ? calculate_contingency_overload!(ΔPcc, flow, (oplim.branch_rating * oplim.long_term_multi), θ, B, Pcc, opf, mod, pf, bd, cont, c, islands, island, island_b) : Tuple{Int64,Float64}[]
            olccx = type >= PCFSC::OPF ? calculate_contingency_overload!(ΔPccx, flow, (oplim.branch_rating * oplim.long_term_multi), θ, B, Pccx, opf, mod, pf, bd, cont, c, islands, island, island_b) : Tuple{Int64,Float64}[]
            if !isempty(olc) || !isempty(olcc) || !isempty(olccx) # ptdf calculation is more computational expensive than line flow
                if !is_islanded(pf, cont, c)
                    get_isf!(ptdf, X, pf.X, pf.B, pf.DA, cont, c)
                else
                    get_isf!(ptdf, pf.ϕ, islands[island], island_b)
                end


                if get(Pc, c, 0) == 0 || get(Pcc, c, 0) == 0 || get(Pccx, c, 0) == 0
                    # basetype == SC::OPF && get_ΔP!(ΔPc, mod, opf, bd.list, Pc, c)
                    type >= PCSC::OPF && get_ΔP!(ΔPcc, mod, opf, bd.list, Pcc, c)
                    type >= PCFSC::OPF && get_ΔP!(ΔPccx, mod, opf, bd.list, Pccx, c)
                end

                # # Calculate the power flow with the new outage and find if there are any overloads
                # overloads_c, overloads_cc, overloads_ccx = find_overloads(Val(type), flow, ptdf, bd.Pᵢ, oplim, ΔPc, ΔPcc, ΔPccx)
                # if !isempty(olc) || !isempty(overloads_c) 
                #     if isempty(olc)
                #         println([o[2] for o in overloads_c])
                #     elseif isempty(overloads_c) 
                #         println([o for o in olc])
                #     else
                #         println([(o1[2], o2[2], o1[2] - o2[2]) for (o1, o2) in zip(overloads_c, olc)])
                #     end
                #     # error(olc)
                # end
                # if !isempty(olcc) || !isempty(overloads_cc) 
                #     if isempty(olcc)
                #         println([o[2] for o in overloads_cc])
                #     elseif isempty(overloads_cc) 
                #         println([o for o in olcc])
                #     else
                #         println([(o1[2], o2[2], o1[2] - o2[2]) for (o1, o2) in zip(overloads_cc, olcc)])
                #     end
                #     # error(olcc)
                # end
            end

            # Cannot change the model before all data is exctracted!
            if !isempty(olc)
                cut_added, pre = add_cut(opf, mod, bd, ptdf, olc, cont, c, cut_added, lim, pre)
                # cut_added, pre = add_cut(Pc, opf, oplim, mod, bd, ΔPc, ptdf, olc, islands, island, cont, c, cut_added, lim, pre)
            end
            if !isempty(olcc)
                cut_added, corr = add_cut(Pcc, opf, oplim, mod, bd, ΔPcc, ptdf, olcc, islands, island, cont, c, cut_added, lim, corr)
            end
            if !isempty(olccx)
                cut_added = add_cut(Pccx, opf, oplim, mod, bd, ΔPccx, ptdf, olccx, islands, island, cont, c, cut_added, lim, corr)
            end
            if cut_added > 1
                total_solve_time = update_model!(mod, pf, bd, total_solve_time)
                cut_added = 1
            end
            termination_status(mod) != MOI.OPTIMAL && return mod, opf, pf, oplim, Pc, Pcc, Pccx, total_solve_time
        end

    end
    @warn "Reached $(length(opf.contingencies)) iterations without a stable solution."
    return mod, opf, pf, oplim, Pc, Pcc, Pccx, total_solve_time
end

""" Solve model and update the power flow object """
function update_model!(mod::Model, pf::DCPowerFlow, bd::Benders, total_solve_time::Real)
    # set_warm_start!(mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
    set_objective_function(mod, bd.obj)
    solve_model!(mod)
    total_solve_time += solve_time(mod)

    bd.Pᵢ = get_value(mod, :p0)
    @. bd.Pg = bd.Pᵢ - bd.Pd
    calc_θ!(pf, bd.Pᵢ)
    calc_Pline!(pf)
    return total_solve_time
end

function find_overloads(flow::Vector{<:Real}, ptdf::AbstractMatrix{<:Real}, Pᵢ::Vector{<:Real},
    ΔP::Vector{<:Real}, branch_rating::Vector{<:Real}
)
    LinearAlgebra.mul!(flow, ptdf, (Pᵢ .+ ΔP))
    return filter_overload(flow, branch_rating)
end

function calculate_contingency_line_flows!(ΔP::Vector{<:Real}, flow::Vector{<:Real}, θ::Vector{<:Real}, B::AbstractMatrix{<:Real}, 
    P::Dict{<:Integer, T}, opf::OPFsystem, mod::Model, pf::DCPowerFlow, bd::Benders, cont::Tuple{Real,Real}, c::Integer, islands::Vector, 
    island::Integer, island_b::Vector{<:Integer}
) where {T}
    if get(P, c, 0) == 0
        calculate_line_flows!(flow, pf, cont, c)
    else
        ΔP = get_ΔP!(ΔP, mod, opf, bd.list, P, c)
        if is_islanded(pf, cont, c)
            if iszero(ΔP)
                calculate_line_flows!(flow, pf, cont, c, nodes=islands[island], branches=island_b)
            else
                calculate_line_flows!(flow, pf, cont, c, Pᵢ=(bd.Pᵢ .+ ΔP), nodes=islands[island], branches=island_b)
            end
        else
            if iszero(ΔP)
                calculate_line_flows!(flow, pf, cont, c)
            else
                calculate_line_flows!(flow, pf, cont, c, Pᵢ=(bd.Pᵢ .+ ΔP))
            end
        end
    end
end
function calculate_contingency_overload!(ΔP::Vector{<:Real}, flow::Vector{<:Real}, branch_rating::Vector{<:Real},
    θ::Vector{<:Real}, B::AbstractMatrix{<:Real}, Pc::Dict, opf::OPFsystem, mod::Model, pf::DCPowerFlow, bd::Benders, 
    cont::Tuple{Real,Real}, c::Integer, islands::Vector, island::Integer, island_b::Vector{<:Integer}
)
    calculate_contingency_line_flows!(ΔP, flow, θ, B, Pc, opf, mod, pf, bd, cont, c, islands, island, island_b)
    return filter_overload(flow, branch_rating)
end

""" Return the short term power injection change at each node. """
function get_ΔP!(ΔP::Vector{T}, mod::Model, opf::OPFsystem, list::Vector{<:CTypes{Int}}, 
    Pc::Dict{<:Integer, Main.SCOPF.ExprC}, c::Integer
) where {T<:Real}
    fill!(ΔP, zero(T))
    x = get(Pc, c, 0)
    if x != 0
        pgc = get_value(mod, x.pgc)
        prc = get_value(mod, x.prc)
        lsc = get_value(mod, x.lsc)
        for (i, n) in enumerate(list)
            for g in n.ctrl_generation
                ΔP[i] -= pgc[g]
            end
            for g in n.renewables
                ΔP[i] -= prc[g]
            end
            for g in n.demands
                ΔP[i] += lsc[g]
            end
        end
    end
    return ΔP
end

""" Return the long term power injection change at each node. """
function get_ΔP!(ΔP::Vector{T}, mod::Model, opf::OPFsystem, list::Vector{<:CTypes{Int}}, 
    Pcc::Dict{<:Integer, Main.SCOPF.ExprCC}, c::Integer
) where {T<:Real}
    fill!(ΔP, zero(eltype(ΔP)))
    x = get(Pcc, c, 0)
    if x != 0
        pgu = get_value(mod, x.pgu)
        pgd = get_value(mod, x.pgd)
        pfdccc = get_value(mod, x.pfdccc)
        prcc = get_value(mod, x.prcc)
        lscc = get_value(mod, x.lscc)
        for (i, n) in enumerate(list)
            for g in n.ctrl_generation
                ΔP[i] += (pgu[g] - pgd[g])
            end
            for g in n.dc_branches
                ΔP[i] += beta(opf.nodes[i], opf.dc_branches[g]) * pfdccc[g]
            end
            for g in n.renewables
                ΔP[i] -= prcc[g]
            end
            for g in n.demands
                ΔP[i] += lscc[g]
            end
        end
    end
    return ΔP
end

""" Return the long term power injection change at each node. """
function get_ΔP!(ΔP::Vector{T}, mod::Model, opf::OPFsystem, list::Vector{<:CTypes{Int}}, 
    Pccx::Dict{<:Integer, Main.SCOPF.ExprCCX}, c::Integer
) where {T<:Real}
    fill!(ΔP, zero(eltype(ΔP)))
    x = get(Pccx, c, 0)
    if x != 0
        pgdx = get_value(mod, x.pgdx)
        lsccx = get_value(mod, x.lsccx)
        for (i, n) in enumerate(list)
            for g in n.ctrl_generation
                ΔP[i] -= pgdx[g]
            end
            for g in n.demands
                ΔP[i] += lsccx[g]
            end
        end
    end
    return ΔP
end

function get_P(P::Dict, opf::OPFsystem, oplim::Oplimits, mod::Model, obj::JuMP.AbstractJuMPScalar, list::Vector{<:CTypes{Int}},
    islands::Vector, island::Integer, c::Integer
)
    if get(P, c, 0) == 0  # If the contingency is not run before, a set of corrective variables is added
        init_P!(P, opf, oplim, mod, obj, list, islands, island, c)
    end
    return get(P, c, 0)
end

function add_overload_expr(expr::JuMP.AbstractJuMPScalar, mod::Model, ol::Real)
    if ol < 0
        JuMP.@constraint(mod, expr <= ol)
    else
        JuMP.@constraint(mod, expr >= ol)
    end
end

function add_cut(opf::OPFsystem, mod::Model, bd::Benders, ptdf::AbstractMatrix{<:Real}, overloads::Vector{<:Tuple{Integer,Real}},
    cont::Tuple{Integer,Integer}, c::Integer, cut_added::Integer, lim::Real, id::Integer
)
    for (i, ol) in overloads
        expr = JuMP.@expression(mod, sum((ptdf[i, inode] * (
                bd.Pg[inode] -
                sum((mod[:pg0][ctrl] for ctrl in sublist.ctrl_generation), init=0.0) -
                sum((beta(sublist.node, opf.dc_branches[d]) * mod[:pfdc0][d] for d in sublist.dc_branches), init=0.0) +
                sum((mod[:pr0][r] for r in sublist.renewables), init=0.0) -
                sum((mod[:ls0][d] for d in sublist.demands), init=0.0)
            ) for (inode, sublist) in enumerate(bd.list)), init=0.0
        ))

        id += 1
        @info @sprintf "Pre %d: Contingency line %d-%d-i_%d; overload on %s of %.4f" id cont[1] cont[2] c opf.branches[i].name ol
        @debug "Cut added: $(sprint_expr(expr,lim))\n"
        add_overload_expr(expr, mod, ol)
        cut_added = 2
    end
    return cut_added, id
end

function add_cut(Pc::Dict{<:Integer, Main.SCOPF.ExprC}, opf::OPFsystem, oplim::Oplimits, mod::Model, bd::Benders, ΔPc::Vector{<:Real},
    ptdf::AbstractMatrix{<:Real}, overloads::Vector{<:Tuple{Integer,Real}}, islands::Vector, island::Integer,
    cont::Tuple{Integer,Integer}, c::Integer, cut_added::Integer, lim::Real, id::Integer
)
    pc = get_P(Pc, opf, oplim, mod, bd.obj, bd.list, islands, island, c)
    for (i, ol) in overloads
        expr = JuMP.@expression(mod, sum((ptdf[i, inode] * (
                bd.Pg[inode] + ΔPc[inode] -
                sum((mod[:pg0][ctrl] - pc.pgc[ctrl] for ctrl in sublist.ctrl_generation), init=0.0) -
                sum((beta(sublist.node, opf.dc_branches[d]) * mod[:pfdc0][d] for d in sublist.dc_branches), init=0.0) +
                sum((pc.prc[r] for r in sublist.renewables), init=0.0) -
                sum((pc.lsc[d] for d in sublist.demands), init=0.0)
            ) for (inode, sublist) in enumerate(bd.list)), init=0.0
        ))

        id += 1
        @info @sprintf "Pre %d: Contingency line %d-%d-i_%d; overload on %s of %.4f" id cont[1] cont[2] c opf.branches[i].name ol
        @debug "Cut added: $(sprint_expr(expr,lim))\n"
        add_overload_expr(expr, mod, ol)
        cut_added = 2
    end
    return cut_added, id
end

function add_cut(Pcc::Dict{<:Integer, Main.SCOPF.ExprCC}, opf::OPFsystem, oplim::Oplimits, mod::Model, bd::Benders, ΔPcc::Vector{<:Real},
    ptdf::AbstractMatrix{<:Real}, overloads::Vector{<:Tuple{Integer,Real}}, islands::Vector, island::Integer,
    cont::Tuple{Integer,Integer}, c::Integer, cut_added::Integer, lim::Real, id::Integer
)
    pcc = get_P(Pcc, opf, oplim, mod, bd.obj, bd.list, islands, island, c)
    for (i, ol) in overloads
        # Finding and adding the Benders cut
        expr = JuMP.@expression(mod, sum((ptdf[i, inode] * (
                bd.Pg[inode] + ΔPcc[inode] -
                sum((mod[:pg0][ctrl] + pcc.pgu[ctrl] - pcc.pgd[ctrl] for ctrl in sublist.ctrl_generation), init=0.0) -
                sum((beta(sublist.node, opf.dc_branches[d]) * pcc.pfdccc[d] for d in sublist.dc_branches), init=0.0) +
                sum((pcc.prcc[r] for r in sublist.renewables), init=0.0) -
                sum((pcc.lscc[d] for d in sublist.demands), init=0.0)
            ) for (inode, sublist) in enumerate(bd.list)), init=0.0
        ))

        id += 1
        @info @sprintf "Corr %d: Contingency line %d-%d-i_%d; overload on %s of %.4f" id cont[1] cont[2] c opf.branches[i].name ol
        @debug "Cut added: $(sprint_expr(expr,lim))\n"
        add_overload_expr(expr, mod, ol)
        cut_added = 3
    end
    return cut_added, id
end

function add_cut(Pccx::Dict{<:Integer, Main.SCOPF.ExprCCX}, opf::OPFsystem, oplim::Oplimits, mod::Model, bd::Benders, ΔPccx::Vector{<:Real},
    ptdf::AbstractMatrix{<:Real}, overloads::Vector{<:Tuple{Integer,Real}}, islands::Vector, island::Integer,
    cont::Tuple{Integer,Integer}, c::Integer, cut_added::Integer, lim::Real, id::Integer
)
    if !isempty(overloads)
        pccx = get_P(Pccx, opf, oplim, mod, bd.obj, bd.list, islands, island, c)
        # sort!(overloads, rev = true, by = x -> abs(x[2]))
        for (i, ol) in overloads
            # Finding and adding the Benders cut
            expr = JuMP.@expression(mod, sum((ptdf[i, inode] * (
                    bd.Pg[inode] + ΔPccx[inode] -
                    sum((mod[:pg0][ctrl] - pccx.pgdx[ctrl] for ctrl in sublist.ctrl_generation), init=0.0) +
                    sum((pccx.prccx[r] for r in sublist.renewables), init=0.0) -
                    sum((pccx.lsccx[d] for d in sublist.demands), init=0.0)
                ) for (inode, sublist) in enumerate(bd.list)), init=0.0
            ))

            @info @sprintf "Corr.x %d: Contingency line %d-%d-i_%d; overload on %s of %.4f" id cont[1] cont[2] c opf.branches[i].name ol
            @debug "Cut added: $(sprint_expr(expr,lim))\n"
            add_overload_expr(expr, mod, ol)
            cut_added = 3
        end
    end
    return cut_added
end

" An AbstractJuMPScalar nicely formatted to a string "
sprint_expr(expr::AbstractJuMPScalar, lim=1e-6) =
    join(Printf.@sprintf("%s%5.2f %s ", (x[2] > 0 ? "+" : "-"), abs(x[2]), x[1])
         for x in expr.terms if abs(x[2]) > lim) *
    Printf.@sprintf("<= %s%.2f", (expr.constant > 0 ? "-" : " "), abs(expr.constant)
    )

function print_benders_results(opf::OPFsystem, mod::Model, Pc::Dict=Dict(), Pcc::Dict=Dict(), Pccx::Dict=Dict(), lim::Real=1e-6)
    function print_c(itr, symb::String, i_g::Int, lim::Real)
        for i in 1:length(opf.contingencies)
            c = get(itr, i, 0)
            if c != 0 && JuMP.value(getfield(c, Symbol(symb))[i_g]) > lim
                @printf("          c %12s: %s: %.3f\n", opf.contingencies[i].name, symb, JuMP.value(getfield(c, Symbol(symb))[i_g]))
            end
        end
    end
    for (i_g, g) in enumerate(opf.ctrl_generation)
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(mod[:pg0][i_g]), get_active_power_limits(g).max)
        print_c(Pc, "pgc", i_g, lim)
        print_c(Pcc, "pgu", i_g, lim)
        print_c(Pcc, "pgd", i_g, lim)
        print_c(Pccx, "pgdx", i_g, lim)
    end
    for (i_g, g) in enumerate(opf.dc_branches)
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(mod[:pfdc0][i_g]), get_active_power_limits(g).max)
        print_c(Pcc, "pfdccc", i_g, lim)
    end
    for (i_g, g) in enumerate(opf.renewables)
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(mod[:pr0][i_g]), get_active_power_limits(g).max)
        print_c(Pc, "prc", i_g, lim)
        print_c(Pcc, "prcc", i_g, lim)
        print_c(Pccx, "prccx", i_g, lim)
    end
    for (i_g, g) in enumerate(opf.demands)
        @printf("%12s: %5.3f (%.3f)\n", g.name, JuMP.value(mod[:ls0][i_g]), get_active_power(g))
        print_c(Pc, "lsc", i_g, lim)
        print_c(Pcc, "lscc", i_g, lim)
        print_c(Pccx, "lsccx", i_g, lim)
    end
end

function run_benders!(
    ::Val{SC::OPF}, 
    mod::Model, 
    opf::OPFsystem, 
    pf::DCPowerFlow, 
    oplim::Oplimits, 
    lim=1e-6,
    max_itr=length(opf.contingencies),
    branch_c=nothing, 
    rate_c=0.0, 
    debug=false
)
    solve_model!(mod)
    total_solve_time = solve_time(mod)
    termination_status(mod) != MOI.OPTIMAL && return mod, opf, pf, oplim, total_solve_time
    @debug "lower_bound = $(objective_value(mod))"

    # Set variables
    bd = benders(opf, mod)
    calc_θ!(pf, bd.Pᵢ)
    calc_Pline!(pf)
    overloads = zeros(length(opf.contingencies))

    contids = get_branch_bus_idx(opf.branches, opf.contingencies, bd.idx)
    flow = similar(pf.F)
    for contid in contids
        (c, cont) = contid
        if !is_islanded(pf, cont, c)
            calculate_line_flows!(flow, pf, cont, c)
        else
            islands, island, island_b = handle_islands(pf.B, pf.DA, cont, c, pf.slack)
            for in_vec in islands[1:end.!=island]
                for n in in_vec
                    add_to_expression!(bd.obj, opf.prob[c], sum((opf.voll[d] * get_active_power(opf.demands[d]) for d in bd.list[n].demands), init=0.0))
                end
            end
            continue
        end

        # Calculate the power flow with the new outage and find if there are any overloads
        overload = filter_overload(flow, oplim.branch_rating * oplim.short_term_multi)

        if !isempty(overload)
            overloads[c] = maximum(x -> abs(x[2]), overload)
        end
    end

    pre = 0
    ptdf = copy(pf.ϕ)
    X = copy(pf.X)
    olc = Vector{Tuple{Int,Float64}}[]

    permutation = sortperm(overloads, rev=true)
    cut_added = 1
    for iterations in 1:max_itr
        if cut_added == 0 # loops until no new cuts are added for the contingencies
            # @printf "\nEND: Total solve time %.4f.\n" total_solve_time
            return mod, opf, pf, oplim, total_solve_time
        end
        cut_added = 0
        @info "Iteration $iterations"

        for contid in contids[permutation]
            (c, cont) = contid
            empty!(olc)
            if is_islanded(pf, cont, c)
                continue
            end
            calculate_line_flows!(flow, pf, cont, c)
            olc = filter_overload(flow, (oplim.branch_rating * oplim.short_term_multi))
            if !isempty(olc) # ptdf calculation is more computational expensive than line flow
                get_isf!(ptdf, X, pf.X, pf.B, pf.DA, cont, c)
            end

            # Cannot change the model before all data is exctracted!
            if !isempty(olc)
                cut_added, pre = add_cut(opf, mod, bd, ptdf, olc, cont, c, cut_added, lim, pre)
            end
            if cut_added > 1
                total_solve_time = update_model!(mod, pf, bd, total_solve_time)
                cut_added = 1
            end
            termination_status(mod) != MOI.OPTIMAL && return mod, opf, pf, oplim, total_solve_time
        end

    end
    @warn "Reached $(length(opf.contingencies)) iterations without a stable solution."
    return mod, opf, pf, oplim, total_solve_time
end
