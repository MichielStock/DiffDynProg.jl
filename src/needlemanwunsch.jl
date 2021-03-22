#=
Created on 05/03/2021 16:09:30
Last update: wednesday 17 March 2021

@author: Michiel Stock
michielfmstock@gmail.com

Implementation of Needleman-Wunsch and its gradients
according to Mench and Blondel.
=#

# NW with fixed cost

function needleman_wunsch(mo::MaxOperator, θ::AbstractMatrix, g::Number, D::Matrix)
    n, m = size(θ)
    @assert size(D, 1) > n && size(D, 2) > m "The dimensions of the DP matrix `D`` and `θ` do not agree"
	D[:,1] .= range(0, step=-g, length=size(D,1))
	D[1,:] .= range(0, step=-g, length=size(D,2))
	@inbounds for j in 1:m, i in 1:n
	 	D[i+1,j+1] = maxᵧ(mo, (D[i+1,j] - g, 
                                    D[i,j] + θ[i,j],
                                    D[i,j+1] - g))
	end
    return D[n+1, m+1]
end

function needleman_wunsch(mo::MaxOperator, θ::AbstractMatrix, g::Number)
    D = Matrix{eltype(θ)}(undef, size(θ,1)+1, size(θ,2)+1)
    return needleman_wunsch(mo, θ, g, D)
end

needleman_wunsch(mo::MaxOperator, θ::Matrix, g, dp::DP) = needleman_wunsch(mo, θ, g, dp.D) 

# NW with variable cost

function needleman_wunsch(mo::MaxOperator, θ::Matrix, (gs, gt)::Tuple{<:AbstractVector,<:AbstractVector}, D::Matrix)
    n, m = size(θ)
    @assert size(D, 1) > n && size(D, 2) > m "The dimensions of the DP matrix `D`` and `θ` do not agree"
	D[1,1] = 0
	D[2:n+1,1] .= cumsum(gs[1:n])
	D[1,2:m+1] .= cumsum(gt[1:m])
	@inbounds for j in 1:m, i in 1:n
	 	D[i+1,j+1] = maxᵧ(mo, (D[i+1,j] - gs[i], 
                                    D[i,j] + θ[i,j],
                                    D[i,j+1] - gt[j]))
	end
    return D[n+1, m+1]
end

# NW with a substition matrix

needleman_wunsch(mo::MaxOperator, s::Vector, t::Vector, S::AbstractMatrix, args...) = needleman_wunsch(mo, S[s, t], args...)

# NW gradients

function ∂NW(mo::MaxOperator, θ::AbstractMatrix, g::Number, D::AbstractMatrix, E::AbstractMatrix, Q::AbstractArray)
    T = eltype(E)
    n, m = size(θ)
	@assert size(D, 1) > n && size(D, 2) > m "The dimensions of the DP matrix `D`` and `θ` do not agree"
    fill!(Q, zero(T))
	Q[n+2, m+2, 2] = one(T)
	E[:,m+2] .= zero(T)
	E[n+2,:] .= zero(T)
	E[n+2,m+2] = one(T)
	D[:,1] .= range(0, step=-g, length=size(D,1))
	D[1,:] .= range(0, step=-g, length=size(D,2))
	@inbounds for j in 1:m, i in 1:n
		v, q = max_argmaxᵧ(mo, (D[i+1,j] - g, 
                                D[i,j] + θ[i,j],
                                D[i,j+1] - g))
    	D[i+1,j+1] = v
       	Q[i+1,j+1,:] .= q
	end
	@inbounds for j in m:-1:1, i in n:-1:1
        # change of D[n,m] by D[i,j]
        E[i+1,j+1] = Q[i+1,j+2,1] * E[i+1,j+2] + Q[i+2,j+2,2] * E[i+2,j+2] + Q[i+2,j+1,3] * E[i+2,j+1]
	end
	return @view(D[2:n+1,2:m+1]), @view(E[2:n+1,2:m+1])
end

∂NW(mo::MaxOperator, θ::AbstractMatrix, g, dp::DP) = ∂NW(mo::MaxOperator, θ::Matrix, g, dp.D, dp.E, dp.Q)
∂NW(mo::MaxOperator, θ::AbstractMatrix, g) = ∂NW(mo, θ, g, DP(θ))

function ∂NW(mo::MaxOperator, θ::AbstractMatrix, (gs, gt)::Tuple{<:AbstractVector,<:AbstractVector}, D::AbstractMatrix, E::AbstractMatrix, Q::AbstractArray)
    T = eltype(E)
    n, m = size(θ)
	@assert size(D, 1) > n && size(D, 2) > m "The dimensions of the DP matrix `D`` and `θ` do not agree"
    fill!(Q, zero(T))
	Q[n+2, m+2, 2] = one(T)
	E[:,m+2] .= zero(T)
	E[n+2,:] .= zero(T)
	E[n+2,m+2] = one(T)
	D[1,1] = zero(T)
	D[2:n+1,1] .= -cumsum(gs[1:n])
	D[1,2:m+1] .= -cumsum(gt[1:m])
	@inbounds for j in 1:m, i in 1:n
		v, q = max_argmaxᵧ(mo, (D[i+1,j] - gs[i], 
                                D[i,j] + θ[i,j],
                                D[i,j+1] - gt[j]))
    	D[i+1,j+1] = v
       	Q[i+1,j+1,:] .= q
	end
	@inbounds for j in m:-1:1, i in n:-1:1
        # change of D[n,m] by D[i,j]
        E[i+1,j+1] = Q[i+1,j+2,1] * E[i+1,j+2] + Q[i+2,j+2,2] * E[i+2,j+2] + Q[i+2,j+1,3] * E[i+2,j+1]
	end
	return @view(D[2:n+1,2:m+1]), @view(E[2:n+1,2:m+1])
end

#∂NW(mo::MaxOperator, θ::AbstractMatrix, gs, gt, dp::DP) = ∂NW(mo::MaxOperator, θ::Matrix, gs, gt, dp.D, dp.E, dp.Q)
#∂NW(mo::MaxOperator, θ::AbstractMatrix, gs, gt) = ∂NW(mo::MaxOperator, θ::AbstractMatrix, gs, gt, DP(θ))

# specific gradients

# this is a dummy where one can specify the direction 

function ∂NW_(dir::Int, mo::MaxOperator, θ::AbstractMatrix, g, dp::DP)
	D, E = ∂NW(mo, θ, g, dp)
	dp.E .*= @view(dp.Q[:,:,dir])
	return D, E
end

∂NW_(dir::Int, mo::MaxOperator, θ::AbstractMatrix, g) = ∂NW_(dir, mo, θ, g, DP(θ))

∂NW_θ(args...) = ∂NW_(2, args...)
∂NW_gs(args...) = ∂NW_(1, args...)
∂NW_gt(args...) = ∂NW_(3, args...)

function ∂NW_all(mo::MaxOperator, θ::AbstractMatrix, g, dp::DP)
	n, m = size(θ)
	D, E = ∂NW(mo, θ, g, dp)
	dp.Q .*= dp.E
	return D, E, @view(dp.Q[2:n+1,2:m+1,2]), @view(dp.Q[2:n+1,2:m+1,1]), @view(dp.Q[2:n+1,2:m+1,3])
end

∂NW_all(mo::MaxOperator, θ::AbstractMatrix, g::Number) = ∂NW_all(mo, θ, g, DP(θ))

# Gradients in ChainRulesCore

function rrule(::typeof(needleman_wunsch), mo::MaxOperator, θ::AbstractMatrix, g::Number, dp::DP)
	n, m = size(θ)
	D, E = ∂NW(mo, θ, g, dp)
	# change of D[n,m] by θ[i,j]
	dp.E .*= @view(dp.Q[:,:,2])
	return last(D), ȳ -> (NO_FIELDS, Zero(), ȳ * E, Zero(), Zero())
end

function rrule(::typeof(needleman_wunsch), mo::MaxOperator, θ::AbstractMatrix, (gs, gt)::Tuple{<:AbstractVector,<:AbstractVector}, dp::DP)
	n, m = size(θ)
	D, E = ∂NW(mo, θ, (gs, gt), dp)
	# gradient wrt the two gap cost vectors
	dgs = similar(gs)
	for i in 1:n
		dgs[i] = @view(dp.Q[i+1,2:m+1,1]) ⋅ @view(dp.E[i+1,2:m+1])
	end
	dgt = similar(gt)
	for j in 1:m
		dgt[j] = @view(dp.Q[2:n+1,j+1,3]) ⋅ @view(dp.E[2:n+1,j+1])
	end
	# change of D[n,m] by θ[i,j]
	dp.E .*= @view(dp.Q[:,:,2])
	return last(D), ȳ -> (NO_FIELDS, Zero(), ȳ * E, (ȳ * dgs, ȳ * dgt), Zero())
end


#=
#This is a test
θ = 1.0 * [-3 -5 -3 0 -3;
    -1 7 -1 -3 -1;
    10 -1 10 -4 10]

s = "FGH"
t = "ABCDE"
n, m = size(θ)
g = 5
mo = EntropyMax()
Q = zeros(n+2, m+2, 3)
D = zeros(n+1, m+1)
E = zeros(n+2, m+2)
=#