#=
Created on 07/12/2020 09:36:55
Last update: -

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