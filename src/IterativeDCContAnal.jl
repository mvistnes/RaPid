# Based on code by Sigurd Hofsmo Jakobsen, SINTEF Energy Research, 2022

mutable struct IterativeDCContAnal
    lodf::Array{<:Real, 3}
    contingencies::Vector{String}
    linerating::Vector{<:Real}
end

""" Note about LODF / PTDF:
In our current implementation we will store the PTDF for the base case (i.e. no outages) in LODF[:,:,1], while 
for all contingencies we will store the N first order contingencies at indices LODF[:,:,i] for i in 2:N. 
    Second/third order contingencies are not yet covered but will likely be indexed in a similar fashion. 

It is worth noting the terminology use of LODF / PTDF. The LODF for e.g.  branch 1 is the same as the PTDF when 
branch 1 is outaged. As such, it may not make sense that, in the code below, the LODF[:,:,1] calculated as the 
PTDF of the base case is not accurate use of this terminology, but rather a simplification for us to avoid 
having to create a separate PTDF variable for the base case (and thus increase the number of variables...).

"""
function iterative_cont_anal(sys::System, nodes::Vector{Bus}, branches::Vector{<:Branch}, contingencies::Vector{String}, idx::Dict{<:Any, <:Integer})
    lodf = zeros(length(branches), length(nodes), length(contingencies)+1)
    i_slack, slack = find_slack(nodes)
    lodf[:,:,1] = get_ptdf(branches, length(nodes), idx, i_slack)
    
    for (i,cont) in enumerate(contingencies)
        branches_cont = [b for b in branches if b.name != cont]
        islands = get_islands(sys, branches_cont)
        if length(islands) == 1
            drop_cont = [j for (j,b) in enumerate(branches) if b.name != cont]
            lodf[drop_cont, :, i+1] = get_ptdf(branches_cont, length(nodes), idx, i_slack)
                # The B-matrix can be altered instead of rebuilt from the ground up
        else
            # only the island with the slack bus is assumed stable
            for island in islands
                if slack.number ∈ get_number.(island[1]) # Nodes in island 1
                    drop_id = sort!(collect(values(get_nodes_idx(island[1]))))
                    drop_cont = [j for (j,b) in enumerate(branches) if b.name != cont && b.name ∈ get_name.(island[2])]
                    if length(drop_id) != 0 && length(drop_cont) != 0
                        lodf[drop_cont, drop_id, i+1] = get_ptdf(island[2], length(island[1]), idx, i_slack)
                    end
                    break
                end
            end
        end
    end
    return IterativeDCContAnal(lodf, contingencies, get_rate.(branches))
end

" Get the overload of all lines "
get_overload(contanal::IterativeDCContAnal, cont::Integer, Pᵢ::Vector{<:Real}) = 
    find_overload.(
            calculate_line_flows(get_cont_ptdf(contanal, cont), Pᵢ), 
            contanal.linerating
        )

        
get_overload(
    DA::AbstractMatrix{<:Real}, 
    B::AbstractMatrix{<:Real}, 
    X::AbstractMatrix{<:Real}, 
    δ::AbstractVector{<:Real}, 
    branch::Integer,
    cont::Tuple{Integer, Integer}, 
    slack::Integer
) = 
    find_overload.(calc_Pline(
                DA, 
                get_changed_angles(X, B, δ, cont[1], cont[2], slack),
                branch), 
            contanal.linerating
        )

" Get the PTDF corresponding to the contingency "
function get_cont_ptdf(contanal::IterativeDCContAnal, cont::Integer)
@assert 0 <= cont <= length(contanal.contingencies)  # "No entry in LODF matrix for the given contingency. "
# cont+1 as first element in LODF is PTDF for base case
return contanal.lodf[:,:, cont+1]  
end

function run_benders_cont_anal(system::System, voll, contingencies, prob, lim = 1e-6)

    opfm = scopf(SC, system, Gurobi.Optimizer, voll=voll)
    solve_model!(opfm.mod)
    lower_bound = objective_value(opfm.mod)

    # Set global variables
    nodes = get_sorted_nodes(opfm.sys)
    branches = get_sorted_branches(opfm.sys)
    idx = get_nodes_idx(nodes)
    list_ctrl = make_list(opfm, get_ctrl_generation)
    contanal = iterative_cont_anal(opfm.sys, nodes, branches, contingencies, idx)

    it = enumerate(contingencies)
    next = iterate(it)
    cut_added = false
    while next !== nothing || cut_added # loops until no new cuts are added for the contingencies
        if next === nothing
            next = iterate(it)
            cut_added = false
        end
        (c, cont), state = next
        overloads = get_overload(contanal, c, get_net_Pᵢ(opfm, nodes))
        overloads = [(i,ol) for (i,ol) in enumerate(overloads) if abs(ol) > lim]
        if length(overloads) > 0
            # sort!(overloads, rev = true, by = x -> abs(x[2]))
            (i,ol) = first(overloads)
            Pᵢ = get_Pᵢ(opfm, nodes)
            ptdf = get_cont_ptdf(contanal, c)
            
            expr = sum(
                    ptdf[i, in] * 
                    (Pᵢ[in] - sum((opfm.mod[:pg0][g.name] for g in list_ctrl[n.name]), init=0.0)) 
                    for (in, n) in enumerate(nodes) if abs(ptdf[i, in]) > lim
                )
            @info "Contingency on $(cont) resulted in overload on $(branches[i].name) of $(ol) \nCut added: $(sprint_expr(expr,lim))\n"
            set_warm_start!(opfm.mod, :pg0) # query of information then edit of model, else OptimizeNotCalled errors
            if ol < 0
                JuMP.@constraint(opfm.mod, expr <= ol)
            else
                JuMP.@constraint(opfm.mod, expr >= ol)
            end

            MOI.set(opfm.mod, MOI.Silent(), true) # supress output from the solver
            opfm.mod = solve_model!(opfm.mod)
            termination_status(opfm.mod) != MOI.OPTIMAL && return
            cut_added = true
        else
            next = iterate(it, state)
        end
    end
    return opfm, contanal
end
