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