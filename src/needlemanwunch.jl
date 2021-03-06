#=
Created on 05/03/2021 16:09:30
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Implementation of Needleman-Wunch and its gradients
according to Mench and Blondel.
=#

function needleman_wunch(mo::MaxOperator, θ::Matrix, g::Number, D)
    n, m = size(θ)
    @assert size(D, 1) > n && size(D, 2) > m "The dimensions of the DP matrix `D`` and `θ` do not agree"
	D[:,1] .= - g * 1:size(D, 1) 
	D[1,:] .= - g * 1:size(D, 2) 
	@inbounds for j in 1:m, i in 1:n
	 	D[i+1,j+1] = maximum(mo, [D[i+1,j] - g, 
                                    D[i,j] + θ[i,j],
                                    D[i,j+1] - g])
	end
    return D[n+1, m+1]
end

function ∂NW(mo::MaxOperator, θ::Matrix, g::Number, D, E, Q)
    T = eltype(E)
    n, m = size(θ)
	@assert size(D, 1) > n && size(D, 2) > m "The dimensions of the DP matrix `D`` and `θ` do not agree"
    fill!(Q, 0)
	Q[n+2, m+2, 2] = 1
	E[:,m+2] .= 0
	E[n+2,:] .= 0
	E[n+2,m+2] = 1
	D[:,1] .= - g * (1:size(D, 1)) 
	D[1,:] .= - g * (1:size(D, 2))
	D[1,1] = 0.0
	@inbounds for j in 1:m, i in 1:n
		v, q = max_argmax(mo, [D[i+1,j] - g, 
                                D[i,j] + θ[i,j],
                                D[i,j+1] - g])
    	D[i+1,j+1] = v
       	Q[i+1,j+1,:] .= q
	end
	@inbounds for j in m:-1:1, i in n:-1:1
        #TODO: check this
        # change of D[n,m] by D[i,j]
        E[i+1,j+1] = Q[i+1,j+2,1] * E[i+1,j+2] + Q[i+2,j+2,2] * E[i+2,j+2] + Q[i+2,j+1,3] * E[i+2,j+1]
	end
    @inbounds for j in m:-1:1, i in n:-1:1
        # change of D[n,m] by θ[i,j]
        E[i+1,j+1] *=  Q[i,j,2]  # amount of step in this direction
	end
	return @view(D[2:n+1,2:m+1]), @view(E[2:n+1,2:m+1])
end
