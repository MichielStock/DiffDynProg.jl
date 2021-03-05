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
	D[1,:] .= - g * 1:size(D, 1) 
	@inbounds for j in 1:m, i in 1:n
	 	D[i+1,j+1] = maximum(mo, [D[i,j] + θ[i,j],
                                    D[i+1,j] - g,
                                    D[i,j+1] - g])
	end
    return D[n+1, m+1]
end

