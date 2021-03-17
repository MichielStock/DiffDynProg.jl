#=
Created on 07/12/2020 09:23:01
Last update: Tuesday 2 March 2021

@author: Michiel Stock
michielfmstock@gmail.com

Implementation of the smooth max operators and their gradients.
=#

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

maximum(mo::Max, x::Vector{<:Number}) = maximum(x)

function maximum(mo::LeakyMax{T}, x::Vector{<:Number}) where {T}
    p = mo.p
    return (one(T) - p) * maximum(x) + p * fin_mean(x)
end

function maximum(mo::EntropyMax, x::Vector{<:Number})
    γ = mo.γ
    return logsumexp(x; γ)
end

function maximum(mo::SquaredMax, x::Vector{<:Number})
    γ = mo.γ
    q = project_in_simplex(x ./ γ , one(γ))
    return q ⋅ x - (γ/2) * norm(q)^2
end

# Chain Rules
# -----------

function frule(::typeof(maximum), ::Max, x::Vector{T}) where {T<:Number}
    i = argmax(x)
    m = x[i]
    q = zeros(T, length(x))
    q[i] = one(T)
    return m, q
end

function frule(::typeof(maximum), mo::LeakyMax, x::Vector{<:Number})
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

function frule(::typeof(maximum), mo::EntropyMax, x::Vector{<:Number})
    γ = mo.γ
    m = logsumexp(x; γ)
    q = exp.((x / γ) .- m / γ)
    return m, q
end

function frule(::typeof(maximum), mo::SquaredMax, x::Vector{<:Number})
    γ = mo.γ
    q = project_in_simplex(x ./ γ, one(γ))
    m = q ⋅ x - (γ/2) * (q ⋅ q)
    return m, q
end


function rrule(::typeof(maximum), mo::MaxOperator, x::Vector{T}) where {T<:Number}
    m, q = max_argmax(mo, x)
    return m, ȳ -> (NO_FIELDS, Zero(), ȳ * q)
end

function rrule(::typeof(maximum), mo::EntropyMax, x::Vector{T}) where {T<:Number}
    γ = mo.γ
    m = logsumexp(x; γ)
    q = exp.((x / γ) .- m / γ)
    return m, ȳ -> (NO_FIELDS, Zero(), ȳ * q)
end



"""
    max_argmax(mo::MaxOperator, x::Vector{<:Number})

Returns the maximum and the argmaxim (i.e., the gradient of the max) of 
a vector `x` using a given `MaxOperator`. 
"""
max_argmax(mo::MaxOperator, x::Vector{<:Number}) = frule(maximum, mo, x)


# Minimum
# -------

minimum(::Max, x::Vector{<:Number}) = minimum(x)
minimum(mo::MaxOperator, x::Vector{<:Number}) = -maximum(mo, -x)

function frule(::typeof(minimum), mo::MaxOperator, x::Vector{<:Number})
    m, q = frule(maximum, mo, -x)
    return -m, q
end

function rrule(::typeof(minimum), mo::MaxOperator, x::Vector{T}) where {T<:Number}
    m, q = min_argmin(mo, x)
    return m, ȳ -> (NO_FIELDS, Zero(), ȳ * q)
end

min_argmin(mo::MaxOperator, x::Vector{<:Number}) = frule(minimum, mo, x)
