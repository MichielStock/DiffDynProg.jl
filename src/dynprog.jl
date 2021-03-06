#=
Created on 06/03/2021 12:46:55
Last update: -

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

# getter functions

getE(dp::DP) = @view(dp.E[2:end-1, 2:end-1])
getD(dp::DP) = @view(dp.D[2:end, 2:end])
getQ(dp::DP) = @view(dp.Q[2:end-1, 2:end-1,:])

# subparts
getE(dp::DP, n, m) = @view(dp.E[2:n+1, 2:m+1])
getD(dp::DP, n, m) = @view(dp.D[2:n+1, 2:m+1])
getQ(dp::DP, n, m) = @view(dp.Q[2:n+1, 2:m+1,:])