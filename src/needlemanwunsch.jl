#=
Created on 05/03/2021 16:09:30
Last update: Monday 08 March 2021

@author: Michiel Stock
michielfmstock@gmail.com

Implementation of Needleman-Wunsch and its gradients
according to Mench and Blondel.
=#

# TODO: variable gap cost

function needleman_wunsch(mo::MaxOperator, θ::Matrix, g::Number, D::Matrix)
    n, m = size(θ)
    @assert size(D, 1) > n && size(D, 2) > m "The dimensions of the DP matrix `D`` and `θ` do not agree"
	D[:,1] .= range(0, step=-g, length=size(D,1))
	D[1,:] .= range(0, step=-g, length=size(D,2))
	@inbounds for j in 1:m, i in 1:n
	 	D[i+1,j+1] = maximum(mo, [D[i+1,j] - g, 
                                    D[i,j] + θ[i,j],
                                    D[i,j+1] - g])
	end
    return D[n+1, m+1]
end

function needleman_wunsch(mo::MaxOperator, θ::Matrix, g::Number)
    D = Matrix{eltype(θ)}(undef, size(θ,1)+1, size(θ,2)+1)
    return needleman_wunsch(mo, θ, g, D)
end

needleman_wunsch(mo::MaxOperator, θ::Matrix, g::Number, dp::DP) = needleman_wunsch(mo::MaxOperator, θ::Matrix, g::Number, dp.D) 

function ∂NW(mo::MaxOperator, θ::Matrix, g::Number, D::AbstractMatrix, E::AbstractMatrix, Q::AbstractArray)
    T = eltype(E)
    n, m = size(θ)
	@assert size(D, 1) > n && size(D, 2) > m "The dimensions of the DP matrix `D`` and `θ` do not agree"
    fill!(Q, zero(T))
	Q[n+2, m+2, 2] = 1
	E[:,m+2] .= 0
	E[n+2,:] .= 0
	E[n+2,m+2] = 1
	D[:,1] .= range(0, step=-g, length=size(D,1))
	D[1,:] .= range(0, step=-g, length=size(D,2))
	@inbounds for j in 1:m, i in 1:n
		v, q = max_argmax(mo, [D[i+1,j] - g, 
                                D[i,j] + θ[i,j],
                                D[i,j+1] - g])
    	D[i+1,j+1] = v
       	Q[i+1,j+1,:] .= q
	end
	@inbounds for j in m:-1:1, i in n:-1:1
        # change of D[n,m] by D[i,j]
        E[i+1,j+1] = Q[i+1,j+2,1] * E[i+1,j+2] + Q[i+2,j+2,2] * E[i+2,j+2] + Q[i+2,j+1,3] * E[i+2,j+1]
	end
    @inbounds for j in m:-1:1, i in n:-1:1
        # change of D[n,m] by θ[i,j]
        E[i+1,j+1] *=  Q[i,j,2]  # amount of step in this direction
	end
	return @view(D[2:n+1,2:m+1]), @view(E[2:n+1,2:m+1])
end

∂NW(mo::MaxOperator, θ::Matrix, g::Number, dp::DP) = ∂NW(mo::MaxOperator, θ::Matrix, g::Number, dp.D, dp.E, dp.Q)

∂NW(mo::MaxOperator, θ::Matrix, g::Number) = ∂NW(mo, θ, g, DP(θ))


# Gradients in ChainRulesCore

function rrule(::typeof(needleman_wunsch), mo::MaxOperator, θ::Matrix, g::Number, dp::DP)
	n, m = size(θ)
	D, E = ∂NW(mo, θ, -g, dp)
	return last(D), ȳ -> (NO_FIELDS, Zero(), ȳ * E, Zero(), Zero())
end

#=
This is a test
θ = [-3 -5 -3 0 -3;
    -1 7 -1 -3 -1;
    10 -1 10 -4 10.0]


n, m = size(θ)
g = 5
mo = EntropyMax()
Q = zeros(n+2, m+2, 3)
D = zeros(n+1, m+1)
E = zeros(n+2, m+2)
=#