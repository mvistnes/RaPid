""" Draw a state from the probabilities p """
generate_state(p::Vector{<:Real}) = [p[i] < rand() for i in axes(p,1)]

generate_duration(λ::Real)::Int64 = round(-log(rand())*λ) 
generate_duration(λ::Real, state::AbstractVector{Bool}) = [x ? 0 : generate_duration(λ) for x in state]

""" Draw extra failures from p_co """
function generate_duration(λ::Real, state::AbstractVector{Bool}, p_co::Matrix{<:Real})
    duration = generate_duration(λ, state)
    for (i,s) in enumerate(state)
        if !s
            for j in axes(p_co,1)
                if p_co[i,j] > rand() && duration[j] < duration[i]
                    duration[j] = duration[i]
                end
            end
        end
    end
    return duration
end

function generate_states(p::Matrix{<:Real}, λ::Real)
    l = size(p,1)
    states = trues(l,size(p,2))
    for s in axes(states,1)
        state = generate_state(p[s,:])
        duration = generate_duration(λ, state)
        for (i,d) in enumerate(duration)
            if d > 0
                c_end = s+d > l ? l : s+d
                states[s:c_end,i] .= false
            end
        end
    end
    return states
end

struct Storm
    min::Float64
    max::Float64
    start::Int64
    duration::Int64
    ramp_up::Int64
    ramp_down::Int64
end
storm() = Storm(1.0, 1.0, 0, 0, 0, 0)

function generate_weather(size::Integer, s::Storm)
    @assert s.min < s.max
    @assert s.start - s.ramp_up > 0
    @assert s.start + s.duration + s.ramp_down <= size
    w = fill(s.min, size)
    for i in 1:(s.ramp_up)
        w[s.start-i] = s.max - i/s.ramp_up*(s.max-s.min)
    end
    w[s.start:s.start+s.duration] .= s.max
    for i in 1:(s.ramp_down)
        w[s.start+s.duration-1+i] = s.max - i/s.ramp_down*(s.max-s.min)
    end
    return w
end

function generate_p(prob::Vector{<:Real}, l::Integer)
    p = Matrix{Float64}(undef, l, length(prob))
    for i in axes(p,1)
        for j in axes(p,2)
            p[i,j] = prob[j] * (i < l/2 ? i : (l - i + 1))
        end
    end
    return p
end

function generate_p_from_weather(prob::Vector{<:Real}, size::Integer, weather::Vector{Storm})
    p = Matrix{Float64}(undef, size, length(prob))
    for i in axes(prob,1)
        if weather[i].duration > 0
            p[:,i] .= generate_weather(size, weather[i])
        else
            p[:,i] .= prob[i]
        end
    end
    return p
end

function generate_scenarios(p::Matrix{<:Real}, λ::Real, num::Int)
    res = [Dict() for _ in axes(p,1)]
    # counts = zeros(length(p)+1)
    for _ in 1:num
        states = generate_states(p, λ)
        # out = count(state)
        # counts[out+1] += 1
        for i in axes(states,1)
            if get(res[i], states[i,:], 0) == 0
                res[i][states[i,:]] = 1
            else
                res[i][states[i,:]] += 1
            end
        end
    end
    return res
end

function find_contingency_probabilities(states::Dict)
    sum_samples = 0
    prob = Vector{Float64}()
    cont = Vector{Pair{String, Vector{Int64}}}()
    for (state,num) in states
        push!(cont, "branch" => [i for (i,v) in enumerate(state) if !v])
        push!(prob, num)
        sum_samples += num
    end
    prob ./= sum_samples
    return cont, prob
end

function generate_contingency_probabilities(p::Matrix{<:Real}, λ::Real, num::Integer)
    states = generate_scenarios(p, λ, num)
    return find_contingency_probabilities.(states)
end
generate_contingency_probabilities(p::Vector{<:Real}, λ::Real, num::Integer, l::Integer) =
    generate_contingency_probabilities(generate_p(p, l), λ, num)
generate_contingency_probabilities(p::Vector{<:Real}, λ::Real, num::Integer, l::Integer, weather::Vector{Storm}) =
    generate_contingency_probabilities(generate_p_from_weather(p, l, weather), λ, num)