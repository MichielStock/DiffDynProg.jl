#=
Created on 06/03/2021 12:46:55
Last update: 17/03/2021

@author: Michiel Stock
michielfmstock@gmail.com

Data structure to store a dynamic programming matrix and its gradients.
=#

"""
	DP{T<:AbstractFloat}

Dynamic programming object. Contains all the matrices needed
for performing the dynamic programming and computing the gradient.
"""
struct DP{T<:AbstractFloat}
	D::Matrix{T}
	E::Matrix{T}
	Q::Array{T,3}
end

DP(T::Type{<:AbstractFloat}, n::Int, m::Int) = DP(Matrix{T}(undef, n+1, m+1),
											Matrix{T}(undef, n+2, m+2),
											Array{T}(undef, n+2, m+2, 3))
DP(n::Int, m::Int) = DP(Float64, n, m)
DP(θ::Matrix) = DP(eltype(θ), size(θ)...)

Base.size(dp::DP) = (size(dp.D, 1) - 1, size(dp.D, 2) - 1)
Base.size(dp::DP, dim::Integer) = Base.size(dp::DP)[dim]

Base.eltype(::DP{T}) where {T} = T

# getter functions

getE(dp::DP) = @view(dp.E[2:end-1, 2:end-1])
getD(dp::DP) = @view(dp.D[2:end, 2:end])
getQ(dp::DP) = @view(dp.Q[2:end-1, 2:end-1,:])

# subparts
getE(dp::DP, n, m) = @view(dp.E[2:n+1, 2:m+1])
getD(dp::DP, n, m) = @view(dp.D[2:n+1, 2:m+1])
getQ(dp::DP, n, m) = @view(dp.Q[2:n+1, 2:m+1,:])

# Recipes for plotting

# FIXME: set `yflip` => true
using RecipesBase

#=
Plots a heatmap of the gradient. Optionally give the size `(n, m)` if only a subset
is used. If you give the two sequences, these will be added to the labels.

The keyword argument `grad` can be specified depending on the type of gradient.

=#
@recipe function f(dp::DP, ns=nothing, mt=nothing; grad=nothing)
    seriestype := :heatmap
    yflip := true
    seriescolor --> :speed

    E = dp.E

    n = size(dp, 1)
    m = size(dp, 2)

    if ns isa AbstractString
        s = ns
        n = length(ns)
        yticks := (1:n, split(s, ""))
    elseif ns isa Integer
        n = ns
    end

    if mt isa AbstractString
        t = mt
        m = length(mt)
        xticks := (1:m, split(t, ""))
    elseif mt isa Integer
        m = mt
    end

    E = @view(E[2:n+1, 2:m+1])
    if grad == "pars" || grad isa typeof(∂NW_θ)
        E = getQ(dp)[1:n,1:m,2] .* E
    elseif grad == "gs" || grad isa typeof(∂NW_gs)
        E = getQ(dp)[1:n,1:m,1] .* E
    elseif grad == "gt" || grad isa typeof(∂NW_gt)
        E = getQ(dp)[1:n,1:m,3] .* E
    end
    E
end

@userplot GradientPlot

@recipe function f(h::GradientPlot)
    dp = h.args[1]

    
    a = 0.75  # scaling 
    c = 1e-2  # cutoff

    seriestype := :quiver
    seriescolor := "orange"
    #yflip := true

    Q = getQ(dp)
    n, m = size(dp)

    x = Int[]
    y = Int[]
    dirs = Tuple{Float64,Float64}[]
    for j in 1:m, i in 1:n
        q1 = (0.0, a * Q[i,j,1])
        if Q[i,j,1] > c 
            push!(x, j-1)
            push!(y, i-1)
            push!(dirs, q1)
        end
        q2 = (a * Q[i,j,2], a * Q[i,j,2])
        if Q[i,j,2] > c
            push!(x, j-1)
            push!(y, i-1)
            push!(dirs, q2)
        end
        q3 = (a * Q[i,j,3], 0.0)
        if Q[i,j,3] > c 
            push!(x, j-1)
            push!(y, i-1)
            push!(dirs, q3)
        end
    end

    quiver := dirs
    x, y
end