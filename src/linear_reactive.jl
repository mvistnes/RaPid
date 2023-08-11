# Equations from "Enchancement of Linear ATC Calculations by the Incorporation of Reactive Power Flows" by Grialva et al.

calc_Pjk_circ(Vj::Real, Gjk::Real) = Vj^2 * Gjk

calc_Qjk_circ(Vj::Real, Bjk::Real, Bjj::Real) = -1 *(Vj^2) * Bjk - (Vj^2) * Bjj

calc_Sjk_circ(Vj::Real, Vk::Real, Yjk::Real) = Vj * Vk * Yjk

calc_M2(Pjk_circ::Real, Qjk_circ::Real, Sjk_circ::Real) = Sjk_circ^2 - Pjk_circ^2 - Qjk_circ^2

calc_a(Pjk_circ::Real, Qjk_circ::Real) = Pjk_circ^2 + Qjk_circ^2

calc_b(Pjk_circ::Real, Qjk_circ::Real, Sjk_circ::Real, rate::Real) = 
    -Pjk_circ * (rate^2 - calc_M2(Pjk_circ, Qjk_circ, Sjk_circ))

calc_c(Pjk_circ::Real, Qjk_circ::Real, Sjk_circ::Real, rate::Real) = 
    0.25 * (rate^2 - calc_M2(Pjk_circ, Qjk_circ, Sjk_circ))^2 - Qjk_circ^2 * rate^2

calc_Pjk(a::Real, b::Real, c::Real, pm::Function) = pm(-b, sqrt(b^2 - 4 * a * c)) / (2 * a)

calc_Qjk(Pjk::Real, rate::Real) = sqrt(rate^2 - Pjk^2)

"""
Calculate the new rate of a branch accounting for reactive power flow in an 
active power rating assuming nearly constant voltage during the transfer.

Input:
    - phi: PTDF for the branch.
    - rate: rating for the branch.
    - v_fbus: voltage at the from end.
    - v_tbus: voltage at the to end.
    - y: admittance series component of the branch.
    - y_shunt: admittance shunt component of the branch.

Output:
    - Pjk: new active power rate
    - Qjk: the reactive power assumed
"""
function calc_linear_reactive_rate(phi::Real, rate::Real, v_fbus::Real, v_tbus::Real, 
        y::Complex{<:Real}, y_shunt::Real = 0.0)
    Pjk_circ = calc_Pjk_circ(v_fbus, real(y))
    Qjk_circ = calc_Qjk_circ(v_fbus, imag(y), y_shunt)
    Sjk_circ = calc_Sjk_circ(v_fbus, v_tbus, abs(y))
    a = calc_a(Pjk_circ, Qjk_circ)
    b = calc_b(Pjk_circ, Qjk_circ, Sjk_circ, rate)
    c = calc_c(Pjk_circ, Qjk_circ, Sjk_circ, rate)
    Pjk = signbit(phi) ? calc_Pjk(a, b, c, -) : calc_Pjk(a, b, c, +)
    Qjk = calc_Qjk(Pjk, rate)
    return abs(Pjk), abs(Qjk)
end