# CC BY 4.0 Matias Vistnes, Norwegian University of Science and Technology, 2022

module LinearProgramming

using LinearAlgebra
using Printf

LIM = 1*10^-10
IT_LIM = 20

"""
Optimize the objective function c with regards to the 
constraints of lhs. A and rhs. b using the Simplex 
primal method.

Input:
 - A: a matrix of coefficients for the constraints.
 - b: a vector of the value of the rhs. of the constraints.
 - c: a vector with the objective coefficients.
 - max: true if the objective should be maximized, false if
   it should be minimized.

Output:
 - x: a vector of the control variable values.
 - z: the optimal value of the objective function.
"""
function simplexmethod(A::Matrix{Float64}, b::Vector{Float64}, c::Vector{Float64}, max::Bool=true)
    T = maketableau(A, b, c, max)
    println("-------- Begin --------")
    iterations = 0
    while true
        # print the tableau
        display(T)
        println(" ")

        val, row, col = getpivot(T)
        # Condition for convergance
        val >= -LIM && break
        iterations >= IT_LIM && break

        iterations += 1
        println("Iteration ",iterations)
        T = pivotstep(T, row, col)
    end
    
    # Calculate the values of the control variables
    x = zeros(length(c))
    for i in 1:length(b)
        for j in 1:length(c)
            if !iszero(T[i+1,j+1])
                x[j] = T[i+1,end]/T[i+1,j+1]
                break
            end
        end
    end
    println("--------  End  --------")
    return x, ifelse(max, T[1,end], -T[1,end])
end

function maketableau(A::Matrix{Float64}, b::Vector{Float64}, c::Vector{Float64}, max::Bool)
    M, N = size(A)
    @assert rank(A) == M
    T = zeros(Float64, 1+N, 1+M+N)

    T[1,1] = 1
    T[1, 1:(M+1)] = ifelse(max, c, -c)
    T[1:(size(A,1)+1),1:(size(A,2)+1)] = A
    for i in 1:(N+1)
        T[i,i+M] = 1
    end
    T[1:end,end] = b
    return T
end

function getpivot(T::Matrix{Float64})
    # find min f value
    val, col = findmin(T[1,1:(end-1)])

    # Find smallest quotient
    row = 1
    x = T[1,end]/T[1,col]
    for i in 3:size(T,1)
        q = T[i,end]/T[i,col]
        if abs(q) < x
            x = q
            row = i
        end
    end
    return (val, row, col)
end

function pivotstep(T::Matrix{Float64}, row::Int64, col::Int64)
    # Row elimination
    for i in 1:size(T,1)
        if i != row
            ratio = T[i,col]/T[row,col]
            for j in 1:size(T,2)
                T[i,j] -= ratio * T[row,j]
            end
        end
    end
    return T
end

# A = [2 1;
#      2 3;
#      1 7]
# b = [12, 15, 21]
# c = [6, 14]
# x, z = @time simplexmethod(A,b,c,false)
# display(x)
# display(z)

# A = [2 8;
#      5 2]
# b = [60, 60]
# c = [-40, -88]
# println(simplexmethod(A,b,c,true))

# A = Float64[3 2 1;
#      2 5 3]
# b = Float64[10, 15]
# c = Float64[2, 3, 4]
# println(simplexmethod(A,b,c,false))

function interiorpointmethod(alpha::Float64, err::Float64, A::Matrix{Float64}, b::Vector{Float64}, b1::Vector{Float64}, c::Vector{Float64})
    if size(A,1) != length(b)+length(b1)
        println("'A' and 'b'+'b1' must have the same dimensions!")
        return
    end

    # Initiliazation
    num = length(b) + 2*length(b1) + length(c)
    A0 = zeros(length(b)+length(b1), num)
    A0[1:size(A,1),1:size(A,2)] = A
    for i in size(A,1):(length(b) + length(c))
        A0[i-size(A,1)+1,i] = 1
    end

    c1 = zeros(num)
    c1[1:length(c)] = c

    x1 = zeros(num)
    x1[1:length(b)] = repeat([0.001],length(b)) # inital value of variables
    for i in 1:length(c)
        x1[i+length(b)-1] = b[i] - sum(A[i,:])
    end
    for i in 1:length(b1)
        x1[i+2*(i-1)+length(b)+length(c)-1] = b1[i] - sum(A[i,:])
        x1[2*i+length(b)+length(c)-1] = -b1[i] + sum(A[i,:])
    end
    
    for i in 1:size(A0,1)
        for j in 1:size(A0,2)
            @printf("%7.1f ", A0[i,j])
        end
        println(' ')
    end
    println(c1, x1)

    # Solver
    println("-------- Begin --------")
    iterations = 0
    converge = false
    while !converge
        iterations += 1
        println("Iteration ", iterations)
        x0 = x1
        D = Diagonal(x1)
        A_hat = A0*D
        c_hat = D*c1
        P_hat = I - transpose(A_hat)*inv(A_hat*transpose(A_hat))*A_hat
        p0 = -P_hat*c_hat
        theta = -minimum(p0)
        x_hat1 = ones(num) + alpha/theta*p0
        x1 = D*x_hat1
        converge = true
        for i in 1:num
            if abs(x1[i] - x0[i]) > err
                converge = false
                break
            end
        end

        for i in x1
            @printf("%6.1f ", i)
        end
        println(transpose(c1)*x1)
    end


    println("--------  End  --------")
end


# # generators = ["G1", "G2", "G3", "LS"]
# A = Matrix{Float64}(undef, 10, 9)
# for i in 1:9 
#     A[10, i] = 1 
# end
# for i in 1:9 
#     A[i, i] = 1 
# end

# b = Vector{Float64}([40, 20, 40, 50, 50, 50, 60, 40, 50]) # capacity
# b1 = Vector{Float64}([140]) # load
# c = Vector{Float64}([ 1,  4,  6,  2,  5, 10,  3,  7,  9]) # cost

# const alpha = 0.9
# const err = 0.001

# # A = [2 1;
# #      2 3;
# #      1 7]
# # b = [12, 15, 21]
# # c = [-6, -14]
# @time interiorpointmethod(alpha, err, A, b, b1, c)


# Gaussian elimination to solve Ax = b:
# x = A \ b

end