#=
Created on 07/12/2020 09:23:01
Last update: Monday 22 March 2021

@author: Michiel Stock
michielfmstock@gmail.com

Implementation of the smooth max operators and their gradients.
=#

# TODO: version that works with tuples

using LinearAlgebra: ⋅, norm
using StatsBase: mean
using ChainRulesCore
import ChainRulesCore: frule, rrule
using ChainRulesCore: NO_FIELDS
import Base: min, max, minimum, maximum

# Type hierarchy
# --------------

abstract type MaxOperator end

# ordinary max
struct Max <: MaxOperator end

# leaks a fraction p
struct LeakyMax{T<:AbstractFloat} <: MaxOperator
    p::Float64
    function LeakyMax(p::T=0.1) where {T<:AbstractFloat}
        @assert 0.0 < p < 1.0 "`p` has to be in ]0, 1["
        return new{T}(p)
    end
end

struct EntropyMax{T<:AbstractFloat} <: MaxOperator
    γ::T
    function EntropyMax(γ::T=1.0) where {T<:AbstractFloat}
        @assert γ > 0.0 "`γ` has to be greater than 0"
        return new{T}(γ)
    end
end

EntropyMax(γ::Integer) = EntropyMax(Float64(γ))

struct SquaredMax{T<:AbstractFloat} <: MaxOperator
    γ::T
    function SquaredMax(γ::T=1.0) where {T<:AbstractFloat}
        @assert γ > 0.0 "`γ` has to be greater than 0"
        return new{T}(γ)
    end
end

SquaredMax(γ::Integer) = SquaredMax(Float64(γ))

# Maximization functions
# ----------------------

smoothmax(::Max, x) = maximum(x)

function smoothmax(mo::LeakyMax{T}, x) where {T}
    p = mo.p
    return (one(T) - p) * maximum(x) + p * fin_mean(x)
end

function smoothmax(mo::EntropyMax, x)
    γ = mo.γ
    return logsumexp(x; γ)
end

function smoothmax(mo::SquaredMax, x)
    γ = mo.γ
    q = project_in_simplex(x ./ γ , one(γ))
    return q ⋅ x - (γ/2) * (q ⋅ q)
end

# Chain Rules
# -----------

function frule(::typeof(smoothmax), ::Max, x)
    i = argmax(x)
    m = x[i]
    q = zeros(eltype(x), length(x))
    q[i] = one(eltype(x))
    return m, q
end

function frule(::typeof(smoothmax), mo::LeakyMax, x)
    p = mo.p
    i = argmax(x)
    pcompl = (one(p) - p)
    m = pcompl * x[i] + p * fin_mean(x)
    n = count(el-> el > -Inf, x)
    q = Vector{typeof(p)}(undef, length(x))
    @. q = p / n * (x > -Inf)
    q[i] += pcompl
    return m, q
end

function frule(::typeof(smoothmax), mo::EntropyMax, x)
    γ = mo.γ
    m = logsumexp(x; γ)
    q = exp.((x ./ γ) .- m / γ)
    return m, q
end

function frule(::typeof(smoothmax), mo::SquaredMax, x)
    γ = mo.γ
    q = project_in_simplex(x ./ γ, one(γ))
    m = q ⋅ x - (γ/2) * (q ⋅ q)
    return m, q
end


function rrule(::typeof(smoothmax), mo::MaxOperator, x)
    m, q = max_argmax(mo, x)
    return m, ȳ -> (NO_FIELDS, Zero(), ȳ * q)
end

function rrule(::typeof(smoothmax), mo::EntropyMax, x)
    γ = mo.γ
    m = logsumexp(x; γ)
    q = exp.((x ./ γ) .- m / γ)
    return m, ȳ -> (NO_FIELDS, Zero(), ȳ * q)
end



"""
    smoothmax_argmax(mo::MaxOperator, x)

Returns the maximum and the argmaxim (i.e., the gradient of the max) of 
a vector `x` using a given `MaxOperator`. 
"""
smoothmax_argmax(mo::MaxOperator, x) = frule(smoothmax, mo, x)


# Minimum
# -------
smoothmin(::Max, x) = minimum(x)
smoothmin(mo::MaxOperator, x) = -smoothmax(mo, -x)

function frule(::typeof(smoothmin), mo::MaxOperator, x)
    m, q = frule(smoothmax, mo, -x)
    return -m, q
end

function rrule(::typeof(smoothmin), mo::MaxOperator, x)
    m, q = min_argmin(mo, x)
    return m, ȳ -> (NO_FIELDS, Zero(), ȳ * q)
end

smoothmin_argmin(mo::MaxOperator, x) = frule(smoothmin, mo, x)


# notational shortcuts

maxᵧ(args...) = smoothmax(args...) 
minᵧ(args...)  = smoothmin(args...) 
max_argmaxᵧ(args...)  = smoothmax_argmax(args...) 
min_argminᵧ(args...)  = smoothmin_argmin(args...) 