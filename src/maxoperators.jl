
"""
    project_in_simplex(v::Vector, z::Number=1)

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