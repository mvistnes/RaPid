""" Draw a state from the probabilities p, then draw extra failures from p_co """
function generate_state(p::Vector{<:Real}, p_co::Matrix{<:Real}, λ::Real)
    @assert axes(p,1) == axes(p_co,1) == axes(p_co,2)
    state = generate_state(p)
    duration = [s ? generate_duration(λ) : 0.0 for s in state]
    state_co = trues(length(p))
    for i in axes(p,1)
        if !state[i]
            for j in axes(p_co,1)
                if p_co[i,j] > rand()
                    state_co[j] = false
                    if duration[j] < duration[i]
                        duration[j] = duration[i]
                    end
                end
            end
        end
    end
    return state .&& state_co, duration
end

generate_state(p::Vector{<:Real}) = [p[i] < rand() for i in axes(p,1)]

generate_duration(λ::Real)::Int64 = round(-log(rand())*λ) 
generate_duration(λ::Real, state::AbstractVector{Bool}) = [x ? 0 : generate_duration(λ) for x in state]

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

function generate_states(p::Matrix{<:Real}, p_co::Matrix{<:Real}, λ::Real)
    @assert size(p,2) == size(p_co,1) == size(p_co,2)
    l = size(p,1)
    states = trues(l,size(p,2))
    for s in axes(states,1)
        state = generate_state(p[s,:])
        duration = generate_duration(λ, state, p_co)
        for (i,d) in enumerate(duration)
            if d > 0
                c_end = s+d > l ? l : s+d
                states[s:c_end,i] .= false
            end
        end
    end
    return states
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

function generate_p_from_weather(prob::Vector{<:Real}, weather::Vector{<:Real}, sensitivity::Vector{<:Real})
    @assert length(prob) == length(sensitivity)
    p = Matrix{Float64}(undef, length(weather), length(prob))
    for i in axes(p,1)
        for j in axes(p,2)
            p[i,j] = prob[j] * (sensitivity[j] > 0.0 ? sensitivity[j] * weather[i] * (0.5 + rand()) : 1.0)
        end
    end
    return p
end

function generate_scenarios(p::Matrix{<:Real}, p_co::Matrix{<:Real}, λ::Real, num::Int)
    res = [Dict() for _ in axes(p,1)]
    # counts = zeros(length(p)+1)
    for _ in 1:num
        states = generate_states(p, p_co, λ)
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

generate_contingency_probabilities(p::Vector{<:Real}, p_co::Matrix{<:Real}, λ::Real, num::Integer, l::Integer) =
    generate_contingency_probabilities(generate_p(p, l), p_co, λ, num)
function generate_contingency_probabilities(p::Matrix{<:Real}, p_co::Matrix{<:Real}, λ::Real, num::Integer)
    states = generate_scenarios(p, p_co, λ, num)
    return find_contingency_probabilities.(states)
end