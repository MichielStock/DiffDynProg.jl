#=
Created on Friday 06 November 2020
Last update: Tuesday 2 March 2021

@author: Michiel Stock
michielfmstock@gmail.com

Implementation of dynamic time warping and its gradients
according to Mench and Blondel.
=#

# TODO: make part of a type tree for DP

"""
	DTW{T<:AbstractFloat}

Dynamic time warping object. Contains all the matrices needed
for the dynamic programming and computing the gradient.
"""
struct DTW{T<:AbstractFloat}
	D::Matrix{T}
	E::Matrix{T}
	Q::Array{T,3}
end

getE(dtw::DTW) = @view(E[2:end-1, 2:end-1])
getD(dtw::DTW) = @view(D[2:end, 2:end])

DTW(T::Type{<:AbstractFloat}, n::Int, m::Int) = DTW(Matrix{T}(undef, n+1, m+1),
											Matrix{T}(undef, n+2, m+2),
											Array{T}(undef, n+2, m+2, 3))
DTW(n::Int, m::Int) = DTW(Float64, n, m)
DTW(θ::Matrix) = DTW(eltype(θ), size(θ)...)


function dynamic_time_warping(mo::MaxOperator, θ, D)
    n, m = size(θ)
    @assert size(D) == (n+1, m+1) "The dimensions of the DP matrix and θ do not agree"
	D[:,1] .= Inf
	D[1,:] .= Inf
	D[1,1] = 0.0
	y = zeros(eltype(D), 3)
	@inbounds for j in 1:m, i in 1:n
		y = zeros(eltype(D), 3)
	 	 D[i+1,j+1] = minimum(mo, y) + θ[i,j]
	end
    return last(D)
end

"""
	dynamic_time_warping(mo::MaxOperator, θ, dtw::DTW)

Performs dynamic time warping using a maximum operators `mo`, a weight matrix
`θ` and a dynamic time warping object `dtw`.
"""
dynamic_time_warping(mo::MaxOperator, θ, dtw::DTW) = dynamic_time_warping(mo, θ, dtw.D)

function ∂DPW(mo::MaxOperator, θ, D, E, Q)
    n, m = size(θ)
    fill!(Q, 0)
	Q[end, end, 2] = 1
	E[:,end] .= 0
	E[end,:] .= 0
	E[end,end] = 1
	D[:,1] .= Inf
	D[1,:] .= Inf
	D[1,1] = 0.0
    y = zeros(eltype(D), 3)
	@inbounds for j in 1:m, i in 1:n
		y .= D[i+1,j], D[i,j], D[i,j+1]
		# caution, this overwrites y for performance purposes
		ymin, yargmin = min_argmin(mo, y)
    	D[i+1,j+1] = ymin + θ[i,j]
       	Q[i+1,j+1,:] .= yargmin
	end
	@inbounds for j in m:-1:1, i in n:-1:1
        E[i+1,j+1] = Q[i+1,j+2,1] * E[i+1,j+2] + Q[i+2,j+2,2] * E[i+2,j+2] + Q[i+2,j+1,3] * E[i+2,j+1]
	end
	return @view(D[2:end,2:end]), @view(E[2:end-1, 2:end-1])
end

∂DPW(mo::MaxOperator, θ, dtw::DTW) = ∂DPW(mo::MaxOperator, θ, dtw.D, dtw.E, dtw.Q)

function rrule(::typeof(dynamic_time_warping), mo::MaxOperator, θ, dtw)
	D, E = ∂DPW(mo, θ, dtw)
	return last(D), ȳ -> (NO_FIELDS, Zero(), ȳ * E, Zero())
end

function ∂DPW(mo::MaxOperator, θ::Matrix{T} where {T})
    n, m = size(θ)
    Q = zeros(n+2, m+2, 3)
    E = zeros(n+2, m+2)
    D = zeros(n+1, m+1)
    return ∂DPW!(mo, θ, D, E, Q)
end


