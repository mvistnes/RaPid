
function test_linear_reactive()
    x = [0.9, 0.37, 0.28]
    r = x .* 0.1
    z = r .+ x .* im
    y = 1 ./ z
    rate = [1.0, 1.3, 1.4]
    bx = [(1,2), (1,3), (2,3)]
    v = [1.0, 1.0, 1.0]
    Pl = [1.5, 0.6, 0.8]
    Pg = [1.5, 0.6, 0.8]
    pf = SCOPF.DCPowerFlow(bx, x, length(v), 1)
    to = calc_linear_reactive_rate.(pf.ϕ[1,:], rate, v[getindex.(bx, [1])], v[getindex.(bx, [2])], y)
    from = calc_linear_reactive_rate.(pf.ϕ[2,:], rate, v[getindex.(bx, [2])], v[getindex.(bx, [1])], y)
    return [maximum([t[1] f[1]]) for (t, f) in zip(to, from)]
end
