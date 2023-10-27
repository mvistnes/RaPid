# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2023

"""Fast Decoupled Power Flow

it_max is the maximum number of iterations allowed

lim is the convergance limit

method: Standard Fast Decoupled Power Flow = "std",
Primal Fast Decoupled Power Flow = :pri and 
Dual Fast Decoupled Power Flow = "dual"

print_out = True -> state for each iteration"""
function fdlf(branches::Vector{<:Branch}, idx, slack, fixed, p0, q0, it_max::Int, lim::AbstractFloat, 
    method::Val; print_out = false
)
    deltaP = similar(p0)
    deltaQ = similar(q0)
    ybus = calc_Ybus(branches, idx, fixed)

    # building the two sub-matrices
    it = 1
    (dP, dQ) = init_d(branches, idx, fixed, ybus, slack, method)
    vang = zeros(length(p0))
    vmag = ones(length(p0))

    while true
        injected_power!(deltaP, deltaQ, p0, q0, ybus, vang, vmag)
        solve_d!(vang, vmag, deltaP, deltaQ, dP, dQ, method)
        if print_out
            println("\n-------------------------  Iteration", it, " -------------------------")
            println(vang)
            println(vmag)
            println(deltaP)
            println(deltaQ)
        end

        # check for convergance
        (all(isapprox.(deltaP, 0.0, atol=lim)) || all(isapprox.(deltaQ, 0.0, atol=lim)) || it >= it_max) && break
        it += 1
    end
    
    if it >= it_max
        println("\nNo convergance in ", it, " iterations\n")
    else
        println("\nConvergance in ", it, " iterations\n")
    end
    return it, deltaP, deltaQ
end

function init_d(branches, idx, fixed, ybus, slack, ::Val{:std})
    println("Standard Fast Decoupled Power Flow")
    dP = dQ = get_klu(calc_B_mark(branches, idx, fixed), slack)
    return dP, dQ
end

function init_d(branches, idx, fixed, ybus, slack, ::Val{:pri})
    println("Primal Fast Decoupled Power Flow")
    dP = get_klu(calc_B_doublemark(ybus), slack)
    dQ = get_klu(calc_B_mark(branches, idx, fixed), slack)
    return dP, dQ
end

function init_d(branches, idx, fixed, ybus, slack, ::Val{:dual})
    println("Dual Fast Decoupled Power Flow")
    dP = get_klu(calc_B_mark(branches, idx, fixed), slack)
    dQ = get_klu(calc_B_doublemark(ybus), slack)
    return dP, dQ
end

"""Solves the sub-matrices"""
function solve_d!(vang, vmag, deltaP, deltaQ, dP, dQ, ::Val{:std})
    ang = dP \ deltaP
    @. vang += ang
    mag = dQ \ deltaQ
    @. vmag += mag
end

solve_d!(vang, vmag, deltaP, deltaQ, dP, dQ, ::Val{:pri}) = 
    solve_d!(vang, vmag, deltaP, deltaQ, dP, dQ, Val(:std))

function solve_d!(vang, vmag, deltaP, deltaQ, dP, dQ, ::Val{:dual})
    ang = dQ \ deltaQ
    @. vmag += ang
    mag = dP \ deltaP
    @. vang += mag
end

"""Make a jacobi matrix out of the sub-matrices"""
function print_jacobi(dP, dQ)
    display(cat(dP, dQ, dims=(1,2)))
end  
