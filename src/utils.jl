#=
Created on 07/12/2020 09:36:55
Last update: 01/03/2021

@author: Michiel Stock
michielfmstock@gmail.com

General utilities used by the function of this package but
are themselves not the main product.
=#

# TODO: make this more efficient
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

function logsumexp(x; γ=1)
    c = maximum(x)
    return c + γ * log(sum(exp.((x .- c)/γ)))
end

# dot product that only consider finite numbers
fin_dot(q, x) = sum(qᵢ * xᵢ for (qᵢ, xᵢ) in zip(q, x) if qᵢ > 0 && xᵢ > -Inf)
# finite mean that only considers finite numbers
fin_mean(x) = mean((xᵢ for xᵢ in x if xᵢ > -Inf))


gap_cost_matrix(n::Int, m::Int) = -(0:n-1) .- (0:m-1)'
gap_cost_matrix(cx::Vector, cy::Vector) = -cumsum(cx) .- cumsum(cy)'

# TODO: check this gradient, I do not trust it 100%...
function rrule(::typeof(gap_cost_matrix), cx::Vector, cy::Vector)
    n, m = length(cx), length(cy)
    csx = cumsum(cx)
    csy = cumsum(cy)
    return -csx .- csy',  ȳ -> (NO_FIELDS, -ȳ * (m:-1.0:1), -ȳ' * (n:-1.0:1))
end

