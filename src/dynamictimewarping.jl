#=
Created on Friday 06 November 2020
Last update: Tuesday 09 March 2021

@author: Michiel Stock
michielfmstock@gmail.com

Implementation of dynamic time warping and its gradients
according to Mench and Blondel.
=#

function dynamic_time_warping(mo::MaxOperator, θ, D)
    n, m = size(θ)
    @assert size(D, 1) > n && size(D, 2) > m "The dimensions of the DP matrix `D`` and `θ` do not agree"
	D[:,1] .= maxintfloat(eltype(D))
	D[1,:] .= maxintfloat(eltype(D))
	D[1,1] = 0.0
	y = zeros(eltype(D), 3)
	@inbounds for j in 1:m, i in 1:n
		y .= D[i+1,j], D[i,j], D[i,j+1]
	 	 D[i+1,j+1] = minimum(mo, y) + θ[i,j]
	end
    return D[n+1, m+1]
end

"""
	dynamic_time_warping(mo::MaxOperator, θ, dp::DP)

Performs dynamic time warping using a maximum operators `mo`, a weight matrix
`θ` and a dynamic programming object `dp`.
"""
dynamic_time_warping(mo::MaxOperator, θ, dp::DP) = dynamic_time_warping(mo, θ, dp.D)

function ∂DTW(mo::MaxOperator, θ, D, E, Q)
    n, m = size(θ)
	@assert size(D, 1) > n && size(D, 2) > m "The dimensions of the DP matrix `D`` and `θ` do not agree"
    fill!(Q, 0)
	Q[n+2, m+2, 2] = 1
	E[:,m+2] .= 0
	E[n+2,:] .= 0
	E[n+2,m+2] = 1
	D[:,1] .= maxintfloat(eltype(D))
	D[1,:] .= maxintfloat(eltype(D))
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
	return @view(D[2:n+1,2:m+1]), @view(E[2:n+1,2:m+1])
end

∂DTW(mo::MaxOperator, θ, dp::DP) = ∂DTW(mo::MaxOperator, θ, dp.D, dp.E, dp.Q)


function ∂DTW(mo::MaxOperator, θ::Matrix{T} where {T})
    n, m = size(θ)
    Q = zeros(n+2, m+2, 3)
    E = zeros(n+2, m+2)
    D = zeros(n+1, m+1)
    return ∂DTW(mo, θ, D, E, Q)
end


function rrule(::typeof(dynamic_time_warping), mo::MaxOperator, θ, dp::DP)
	D, E = ∂DTW(mo, θ, dp)
	return last(D), ȳ -> (NO_FIELDS, Zero(), ȳ * E, Zero())
end

