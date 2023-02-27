# Very prelimerary code

module NonlinearProgramming

mutable struct System
    F # objective
    G # constraints
    L # lambda

end

function steepestdecent(alpha)
    U # control variables
    P # parameter values
    X # state variables
    W # dependent variables

    system = System(zeros(length(U)), zeros(length(P)), zeros(length(U)+length(P)))
    gradientL = ones(length(L))
    while max(gradientL) > 0
        gradientL = calcgradient(system, gradientL)
        system.U .-= alpha*gradientL
    end
end

function calcgradient(system, gradientL)
    for i in 1:length(system.L)
        gradientL[i] = 0
    end
    return gradientL
end

end