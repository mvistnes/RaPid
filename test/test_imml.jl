using PowerSystems

# system = get_system("ELK14.json")
system = System("data\\ELK14\\A5.m")
opfm = scopf(SCOPF.SC, system, HiGHS.Optimizer)
solve_model!(opfm.mod)
nodes = get_sorted_nodes(opfm.sys)
branches = get_sorted_branches(opfm.sys)
idx = get_nodes_idx(nodes)
slack = find_slack(nodes)[1]

A = calc_A(branches, length(nodes), idx)
D = calc_D(branches)
B = calc_B(A, D)
X = calc_X(B, slack)
ϕ = calc_isf(A, D, X)
Pᵢ = get_net_Pᵢ(opfm, nodes, idx)
θ = run_pf(B, Pᵢ)

pf = SCOPF.DCPowerFlow(nodes, branches, idx, Pᵢ)
b1 = get_bus_idx(branches[1], idx)
@test calc_Pline(A, D, get_changed_angles(pf.X, pf.DA, pf.θ, b1[1], b1[2], slack)) ≈ 
    calculate_line_flows(calc_isf(D*A, get_changed_X(X, D*A, b1[1], b1[2])), Pᵢ) ≈
    calculate_line_flows(F, ptdf, X, B, angles, b1[1], b1[2], 1)

b3 = get_bus_idx(branches[3], idx)
@test get_changed_angles(pf.X, pf.DA, pf.θ, b2[1], b2[2], slack)
@test get_changed_X(X, D*A, b2[1], b2[2])
@test calculate_line_flows(F, ptdf, X, B, angles, b2[1], b2[2], 2)

@test θ ≈ pf.θ
@test D*A*θ ≈ calculate_line_flows(pf)

b5 = get_bus_idx(branches[5], idx)
@test_throws DivideError get_changed_anglespf.(X, pf.DA, pf.θ, b5[1], b5[2], slack)
@test_throws DivideError get_changed_X(X, D*A, b5[1], b5[2])
@test_throws DivideError calculate_line_flows(F, ptdf, X, B, angles, b5[1], b5[2], 2)

function test_imml_island_handling(pf::DCPowerFlow, bx::Vector{<:Tuple{Integer, Integer}})
    ptdf1 = copy(pf.ϕ)
    ptdf2 = copy(pf.ϕ)
    for (c,cont) in enumerate(bx)
        island, island_b = handle_islands(pf, cont, c)
        try
            flow = SCOPF.calculate_line_flows(pf, cont, c)
            if length(island) < 500
                println(c, cont)
            end
        catch DivideError
            if length(island) < 2
                @warn "Error and no islands!"
            end
        end
        try
            if isempty(island)
                ptdf1 = get_isf(pf.DA, pf.B, cont, c, pf.slack)
                ptdf2 = get_isf(pf, cont, c)
            end
            if any(abs.(ptdf1 .- ptdf2) .> 0.0001)
                @error "Diff $c $(cont[1])-$(cont[2])!"
            end
        catch
            @error "$c $(cont[1])-$(cont[2])!"
            break
        end
    end
end

test_imml_island_handling(pf, get_bus_idx.(branches, [idx]))
