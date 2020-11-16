#=
Created on Friday 06 November 2020
Last update: Monday 16 November 2020

@author: Michiel Stock
michielfmstock@gmail.com

Implementation of dynamic time warping and its gradients
according to Mench and Blondel.
=#

# TODO: make part of a type tree for DP
struct DTW{T<:AbstractFloat}
	D::Matrix{T}
	E::Matrix{T}
	Q::Array{T,3}
end

DTW(T::Type{<:AbstractFloat}, n, m) = DTW(Matrix{T}(undef,n+1,m+1),Matrix{T}(undef,n+2,m+2),Array{T}(undef,n+2,m+2,3))
DTW(n, m) = DTW(Float64, n, m)
DTW(θ::Matrix{T} where {T<:AbstractFloat}) = DTW(T, size(θ)...)

function DPW!(D, θ)
    n, m = size(θ)
    @assert size(D) == (n+1, m+1) "The dimensions of the DP matrix and θ do not agree"
	D[:,1] .= Inf
	D[1,:] .= Inf
	D[1,1] = 0.0
	for i in 1:n, j in 1:m
	 	D[i+1,j+1] = min(D[i+1,j], D[i,j], D[i,j+1]) + θ[i,j]
	end
    return last(D)
end

DPW(θ) = DPW!(zeros(eltype(θ), size(θ,1)+1, size(θ,2)+1), θ)


function ∂DPW!(mo::MaxOperator, θ, D, E, Q)
    n, m = size(θ)
    fill!(Q, zero(eltype(Q)))
	Q[end, end, 2] = 1
	E[:,end] .= 0
	E[end,:] .= 0
	E[end,end] = 1
	D[:,1] .= Inf
	D[1,:] .= Inf
	D[1,1] = 0.0
    y = zeros(eltype(D), 3)
	for j in 1:m, i in 1:n
		y .= D[i+1,j], D[i,j], D[i,j+1]
		# caution, this overwrites y for performance purposes
		ymin, yargmin = min_argmin!(mo, y)
    	D[i+1,j+1] = ymin + θ[i,j]
       	Q[i+1,j+1,:] .= yargmin
	end
	for j in m:-1:1, i in n:-1:1
        E[i+1,j+1] = Q[i+1,j+2,1] * E[i+1,j+2] + Q[i+2,j+2,2] * E[i+2,j+2] + Q[i+2,j+1,3] * E[i+2,j+1]
	end
	return @view(D[2:end,2:end]), @view(E[2:end, 2:end])
end

∂DPW!(mo::MaxOperator, θ, dtw::DTW) = ∂DPW!(mo::MaxOperator, θ, dtw.D, dtw.E, dtw.Q) 

function ∂DPW(mo::MaxOperator, θ::Matrix{T} where {T})
    n, m = size(θ)
    Q = zeros(n+2, m+2, 3)
    E = zeros(n+2, m+2)
    D = zeros(n+1, m+1)
    return ∂DPW!(mo, θ, D, E, Q)
end

s = sin.(0:1:20pi)
t = cos.(1:1:10pi)

θ = (s .- t').^2


mo = EntropyMax(100.0)

dtw = DTW(θ)

D, E = ∂DPW!(mo::MaxOperator, θ, dtw::DTW)

using Plots

plot(heatmap(θ), heatmap(D), heatmap(E))




