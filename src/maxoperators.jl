
using LinearAlgebra: ⋅, norm
using StatsBase: mean

abstract type MaxOperator end

struct Max <: MaxOperator end

struct LeakyMax <: MaxOperator
    p::Float64
end

struct EntropyMax <: MaxOperator
    γ::Float64
end

EntropyMax() = EntropyMax(1.0)

struct SquaredMax <: MaxOperator
    γ::Float64
end

SquaredMax() = SquaredMax(1.0)

"""
    project_in_simplex(v::Vector, z::Number)

Projects a vector `v` in the simplex, optionally scale by `z`.
"""
function project_in_simplex(v::Vector, z::Number)
    n = length(v)
    μ = sort(v, rev=true)
    μcs = cumsum(μ)
    ρ = maximum((i for (i, μi) in enumerate(μ) if μi - (μcs[i] - z) / i > 0.0))
    θ = (μcs[ρ] - z) / ρ
    return max.(v .- θ, 0.0)
end

# Maximization functions

function max_argmax!(::Max, x::VecOrMat{T}) where {T}
    i = argmax(x)
    m = x[i]
    x .= zero(T)
    x[i] = one(T)
    return m, x
end

function max_argmax!(mo::LeakyMax, x::VecOrMat{T}) where {T}
    # if there are infinities, just use regular max
    !all(isfinite.(x)) && return max_argmax!(Max(), x)
    p = mo.p
    n = length(x)
    i = argmax(x)
    ρ = (1.0 - (n - 1) * p)
    m = p * sum(x) - (1.0 - (n - 1) * p) * x[i]
    x .= p
    x[i] = 1.0 - (n - 1) * p
    return m, x
end

function max_argmax!(mo::EntropyMax, x::VecOrMat{T}) where {T}
    # if there are infinities, just use regular max
    !all(isfinite.(x)) && return max_argmax!(Max(), x)
    γ = mo.γ
    xm = mean(x)
    # substract mean, scale and exponent
    @. x = exp((x - xm) / γ)
    se = sum(x)
    x ./= se
    return γ * log(se),  x
end

function max_argmax!(mo::SquaredMax, x::Vector)
    # if there are infinities, just use regular max
    !all(isfinite.(x)) && return max_argmax!(Max(), x)
    γ = mo.γ
    q = project_in_simplex(x ./ γ , 1.0)
    sqm = x ⋅ q - 0.5γ * norm(q)^2.0
    x .= q
    return sqm,  q
end

function max_argmax(mo::SquaredMax, x::Vector)
    γ = mo.γ
    q = project_in_simplex(x ./ γ , 1.0)
    sqm = x ⋅ q - 0.5γ * norm(q)^2.0
    return sqm,  q
end

max_argmax(mo::MaxOperator, x) = max_argmax!(mo, copy(x))

function min_argmin!(mo::MaxOperator, x)
    x .= -x
    m, am = max_argmax!(mo, x)
    return -m, am
end

min_argmin(mo::MaxOperator, x) = min_argmin!(mo, copy(x))
