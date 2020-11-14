
using LinearAlgebra: ⋅, norm

abstract type GenMax end

struct Max <: GenMax end

struct EntropyMax <: GenMax
    γ::Float64
end

EntropyMax() = EntropyMax(1.0)

struct SquaredMax <: GenMax
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

function max_argmax(::Max, x)
    i = argmax(x)
    am = zero(x)
    am[i] = one(eltype(x))
    return x[i],  am
end

function max_argmax(m::EntropyMax, x)
    γ = m.γ
    n = length(x)
    # substract mean, scale and exponent
    ex = (x .- sum(x) / n) ./ γ .|> exp
    se = sum(ex)
    return γ .* log(se),  ex ./ se
end

function max_argmax(m::SquaredMax, x)
    γ = m.γ
    q = project_in_simplex(x ./ γ , 1.0)
    sqm = x ⋅ q - 0.5γ * norm(q)^2.0
    return sqm,  q
end
