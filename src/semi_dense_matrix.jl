abstract type SemiDenseMatrix{TR} <: AbstractMatrix{TR} end

function Base.getindex(mx::SemiDenseMatrix, bus::Integer)
    if mx.filled[bus]
        return mx.data[bus]
    else
        mx.filled[bus] = true
        # println(sum(mx.filled))
        # if sum(mx.filled) > 5
        #     throw(@error "X")
        # end
        return _calc_vec(mx, bus)
    end
end

Base.getindex(mx::SemiDenseMatrix, i::Integer, ::Colon) = Base.getindex(mx, i)
Base.getindex(mx::SemiDenseMatrix, i::Integer, j::Integer) = Base.getindex(mx, i)[j]
Base.setindex!(mx::SemiDenseMatrix, val::Real, i::Integer, j::Integer) = mx.data[i][j] = val

# function *(A::SemiDenseMatrix{TR}, x::Vector) where {TR<:Real}
#     res = Vector{TR}(undef, size(A, 1))
#     @inbounds for i in axes(A,1)
#         res[i] = A[i]'x
#     end
#     return res
# end

mutable struct SemiDenseX{TI<:Integer, TR<:Real} <: SemiDenseMatrix{TR}
    data::Vector{Vector{TR}}
    filled::Vector{Bool}
    K::KLU.KLUFactorization{TR,TI}
    slack::TI
end

SemiDenseX(K::KLU.KLUFactorization{TR,TI}, slack::TI) where {TI<:Integer, TR<:Real} = 
    SemiDenseX{TI,TR}(fill(Vector{TR}(), size(K,1)), falses(size(K,1)), K, slack)

Base.size(mx::SemiDenseX) = size(mx.K)
Base.getindex(mx::SemiDenseX, ::Colon, j::Integer) = Base.getindex(mx, j)

function _calc_vec(mx::SemiDenseX{TI,TR}, bus::Integer) where {TI<:Integer, TR<:Real}
    x = zeros(TR, size(mx.K,1))
    x[bus] = one(TR)
    KLU.solve!(mx.K, x)
    x[mx.slack] = zero(TR)
    mx.data[bus] = x
    return x
end

mutable struct SemiDensePTDF{TI<:Integer, TR<:Real} <: SemiDenseMatrix{TR}
    data::Vector{Vector{TR}}
    filled::Vector{Bool}
    K::KLU.KLUFactorization{TR,TI}
    DA::SparseArrays.SparseMatrixCSC{TR,TI}
    slack::TI
    slack_array::Vector{TR}
end

function SemiDensePTDF(K::KLU.KLUFactorization{TR,TI}, DA::SparseArrays.SparseMatrixCSC{TR,TI}; slack::TI=1, 
    slack_array::AbstractVector{TR}=ones(TR,1), mgx::AbstractMatrix{TR}=ones(TR,1,1)
    ) where {TI<:Integer, TR<:Real}
    @assert sum(slack_array) ≈ one(TR)
    SemiDensePTDF{TI,TR}(fill(Vector{TR}(), size(DA,1)), falses(size(DA,1)), K, DA, slack, (slack_array' * mgx)')
end

Base.size(mx::SemiDensePTDF) = size(mx.DA)

function _calc_vec(mx::SemiDensePTDF{TI,TR}, branch::Integer) where {TI<:Integer, TR<:Real}
    vec = Vector(mx.DA[branch,:])
    KLU.solve!(mx.K, vec)
    set_tol_zero!(vec)
    vec[mx.slack] = zero(TR)
    if size(mx.slack_array) == size(vec)
        mx.data[branch] = vec .- (mx.slack_array'vec)'
    else
        mx.data[branch] = vec
    end
    return mx.data[branch]
end

function set_dist_slack!(mx::SemiDensePTDF{TI,TR}, mgx::AbstractMatrix{<:Real}, dist_slack::AbstractVector{<:Real}
    ) where {TI<:Integer, TR<:Real}
    @assert !iszero(sum(dist_slack))
    slack_array = dist_slack / sum(dist_slack)
    mx.slack_array = (slack_array' * mgx)'
    for i in axes(mx.data, 1)
        if !isempty(mx.data[i])
            mx.data[i] = mx.data[i] .- (mx.slack_array'mx.data[i])'
        end
    end
    return mx.slack_array
end