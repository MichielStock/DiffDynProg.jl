#=
Created on 07/12/2020 09:36:55
Last update: Monday 22 March 2021

@author: Michiel Stock
michielfmstock@gmail.com

General utilities used by the function of this package but
are themselves not the main product.
=#




# TODO: make this more efficient, special case for three dims, using a sorting network?
"""
    project_in_simplex(v::Vector, z::Number)

Projects a vector `v` in the simplex, optionally scale by `z`.
"""
function project_in_simplex(v::Vector{T}, z::Number) where {T<:Number}
    n = length(v)
    μ = sort(v, rev=true)
    μcs = cumsum(μ)
    ρ = maximum((i for (i, μi) in enumerate(μ) if μi - (μcs[i] - z) / i > 0.0))
    θ = (μcs[ρ] - z) / ρ
    return max.(v .- θ, zero(T))
end

project_in_simplex(v::Tuple, z) = project_in_simplex([v...], z)

function project_in_simplex((v₁, v₂, v₃)::NTuple{3,T}, z) where {T<:Number}
    μ₁, μ₂, μ₃ = v₁, v₂, v₃
    # three comparisions to sort
    if μ₁ < μ₂; μ₁, μ₂ = μ₂, μ₁; end
    if μ₂ < μ₃; μ₂, μ₃ = μ₃, μ₂; end
    if μ₁ < μ₂; μ₁, μ₂ = μ₂, μ₁; end
    # find theta
    θ = μ₁ - z
    μ₂ - 0.5 * (μ₁ + μ₂ - z) > 0.0 && (θ = 0.5 * (μ₁ + μ₂ - z))
    μ₃ - (1.0 / 3.0) * (μ₁ + μ₂ + μ₃ - z) > 0.0 && (θ = (1.0 / 3.0) * (μ₁ + μ₂ + μ₃ - z))
    # return vector
    return max.((v₁, v₂, v₃) .- θ, zero(T))
end

# dot product that only consider finite numbers
fin_dot(q, x) = sum(qᵢ * xᵢ for (qᵢ, xᵢ) in zip(q, x) if qᵢ > 0 && xᵢ > -Inf)
# finite mean that only considers finite numbers
fin_mean(x) = mean((xᵢ for xᵢ in x if xᵢ > -Inf))


gap_cost_matrix(n::Int, m::Int) = -(0:n-1) .- (0:m-1)'
gap_cost_matrix(cx::Vector, cy::Vector) = -cumsum(cx) .- cumsum(cy)'


function rrule(::typeof(gap_cost_matrix), cx::Vector, cy::Vector)
    n, m = length(cx), length(cy)
    csx = cumsum(cx)
    csy = cumsum(cy)
    return -csx .- csy',  ȳ -> (NO_FIELDS, -ȳ * (m:-1.0:1), -ȳ' * (n:-1.0:1))
end

# Generating Gumbel random values

randg() = - log(-log(rand()))

"""sample a vector or array of values from Gumbel(0, 1)"""
randg(n::Int...) = - log.(-log.(rand(n...)))

exprandg(n::Int) = 1.0 ./ -log.(rand(n))

function _logsumexp(x; γ=1)
    c = maximum(x)
    return c + γ * log(sum(exp.((x .- c) ./ γ)))
end

function logsumexp(X; dims::Union{Nothing,Int}=nothing, γ=1)
    dims isa Nothing && return _logsumexp(X; γ)
    c = maximum(X; dims)
    return c .+ γ * log.(sum(exp.((X .- c) ./ γ); dims))
end

function _softmax(x; γ=1)
    m = logsumexp(x; γ)
    return exp.((x .- m) ./ γ)
end

function softmax(X; dims::Union{Nothing,Int}=nothing, γ=1)
    dims isa Nothing && return _softmax(X; γ)
    m = logsumexp(X; dims, γ)
    return exp.((X .- m) ./ γ)
end

"""
    gumbel_softmax(lp::Vector; τ::Number=0.1)

Compute the Gumbel softmax approximation of sampling a one-hot-vector the logarithm
of an (unnormalized) probability vector. `τ` is the temperature parameter determining
the quality of the approximation.
"""
function gumbel_softmax(lp::Vector; τ::Number=0.1)
    z = lp .+ randg(length(lp))
    z = z .- logsumexp(z; γ=τ)
    return exp.(z ./ τ)
end

"""
    gumbel_softmax(lp::Array; τ::Number=0.1, dims=2)

Compute the Gumbel softmax approximation of sampling a one-hot-vector the logarithm
of an (unnormalized) probability matric (row-wise by default). `τ` is the temperature
parameter determining the quality of the approximation.
"""
function gumbel_softmax(lp::Array; τ::Number=0.1, dims=2)
	Z = lp .+ randg(size(lp)...)
	Z = Z .- logsumexp(Z; dims, γ=τ)
    return exp.(Z ./ τ)
end

