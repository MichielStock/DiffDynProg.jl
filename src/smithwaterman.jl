#=
Created on 26/03/2021 14:08:31
Last update: Tuesday 30/03/2021

@author: Michiel Stock
michielfmstock@gmail.com

Implementation of Needleman-Wunsch and its gradients
following the principles of Mench and Blondel.
=#

# QUESTION: should I use a γ for the max?



function smith_waterman(mo::MaxOperator, θ::AbstractMatrix, (gs, gt)::TV, D::Matrix)
    n, m = size(θ)
    T = eltype(D)
    D[:,1] .= zero(T)
	D[1,:] .= zero(T)
    @assert size(D, 1) > n && size(D, 2) > m "The dimensions of the DP matrix `D`` and `θ` do not agree"
	@inbounds for j in 1:m, i in 1:n
	 	D[i+1,j+1] = maxᵧ(mo, (D[i+1,j] - gs[i], 
                                    D[i,j] + θ[i,j],
                                    D[i,j+1] - gt[j],
                                    zero(T)))
	end
    return logsumexp(@view(D[2:n+1,2:m+1]))
end

smith_waterman(mo::MaxOperator, θ::AbstractMatrix, g::TV, dp::DP) = smith_waterman(mo, θ, g, dp.D)
smith_waterman(mo::MaxOperator, θ::AbstractMatrix, g::TV) = smith_waterman(mo, θ, g, zeros(eltype(θ), size(θ,1)+1, size(θ,2)+1))
smith_waterman(mo::MaxOperator, θ::AbstractMatrix, g::Number, args...) = smith_waterman(mo, θ, g .* Ones.(size(θ)), args...)
smith_waterman(mo::MaxOperator, θ::AbstractMatrix, (gs, gt)::TN, args...) = smith_waterman(mo, θ, (gs*Ones(size(θ,1)), gt*Ones(size(θ,2))), args...)
# SW with a substition matrix
smith_waterman(mo::MaxOperator, s::Vector, t::Vector, S::AbstractMatrix, args...) = smith_waterman(mo, S[s, t], args...)



function ∂SW(mo::MaxOperator, θ::AbstractMatrix, (gs, gt)::TV, D::AbstractMatrix, E::AbstractMatrix, Q::AbstractArray)
    T = eltype(E)
    n, m = size(θ)
	@assert size(D, 1) > n && size(D, 2) > m "The dimensions of the DP matrix `D`` and `θ` do not agree"
    fill!(Q, zero(T))
	E[:,m+2] .= zero(T)
	E[n+2,:] .= zero(T)
    D[:,1] .= zero(T)
	D[1,:] .= zero(T)
	@inbounds for j in 1:m, i in 1:n
		v, q = max_argmaxᵧ(mo, (D[i+1,j] - gs[i], 
                                D[i,j] + θ[i,j],
                                D[i,j+1] - gt[j],
                                zero(T)))
    	D[i+1,j+1] = v
       	Q[i+1,j+1,:] .= q[1], q[2], q[3]  # last one is discarded
	end
    # find smooth max
    lse = logsumexp(@view(D[2:n+1,2:m+1]))
	@inbounds for j in m:-1:1, i in n:-1:1
        # change of D[n,m] by D[i,j]
        E[i+1,j+1] = exp.(D[i+1,j+1] - lse) +
                    Q[i+1,j+2,1] * E[i+1,j+2] +
                    Q[i+2,j+2,2] * E[i+2,j+2] +
                    Q[i+2,j+1,3] * E[i+2,j+1]
	end
	return lse, @view(E[2:n+1,2:m+1])
end

∂SW(mo::MaxOperator, θ::AbstractMatrix, g::TV, dp::DP) = ∂SW(mo, θ, g, dp.D, dp.E, dp.Q)
∂SW(mo::MaxOperator, θ::AbstractMatrix, g::TV) = ∂SW(mo, θ, g, DP(θ))
∂SW(mo::MaxOperator, θ::AbstractMatrix, g::Number, args...) = ∂SW(mo, θ, g .* Ones.(size(θ)), args...)

# specific gradients

# this is a dummy where one can specify the direction 

function ∂SW_(dir::Int, mo::MaxOperator, θ::AbstractMatrix, g, dp::DP)
	D, E = ∂SW(mo, θ, g, dp)
	dp.E .*= @view(dp.Q[:,:,dir])
	return D, E
end

∂SW_(dir::Int, mo::MaxOperator, θ::AbstractMatrix, g) = ∂SW_(dir, mo, θ, g, DP(θ))

∂SW_θ(args...) = ∂SW_(2, args...)
∂SW_gs(args...) = ∂SW_(1, args...)
∂SW_gt(args...) = ∂SW_(3, args...)

function ∂SW_all(mo::MaxOperator, θ::AbstractMatrix, g, dp::DP)
	n, m = size(θ)
	D, E = ∂SW(mo, θ, g, dp)
	dp.Q .*= dp.E
	return D, E, @view(dp.Q[2:n+1,2:m+1,2]), @view(dp.Q[2:n+1,2:m+1,1]), @view(dp.Q[2:n+1,2:m+1,3])
end

∂SW_all(mo::MaxOperator, θ::AbstractMatrix, g::Number) = ∂SW_all(mo, θ, g, DP(θ))

function rrule(::typeof(smith_waterman), mo::MaxOperator, θ::AbstractMatrix, g::Union{Number,TN}, dp::DP)
	n, m = size(θ)
	D, E = ∂SW(mo, θ, g, dp)
	# change of D[n,m] by θ[i,j]
	dp.E .*= @view(dp.Q[:,:,2])
	return last(D), ȳ -> (NO_FIELDS, Zero(), ȳ * E, Zero(), Zero())
end

function rrule(::typeof(smith_waterman), mo::MaxOperator, θ::AbstractMatrix, (gs, gt)::TV, dp::DP)
	n, m = size(θ)
	lse, E = ∂SW(mo, θ, (gs, gt), dp)
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
	return lse, ȳ -> (NO_FIELDS, Zero(), ȳ * E, (ȳ * dgs, ȳ * dgt), Zero())
end

s = [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0]
t = [0.1, 0.1, 1.2, 1.2, 0.8, 0, 0, -3, -2, -1]

θ = 1 .- abs.(s' .- t)
n, m = size(θ)
gs, gt = ones(n), ones(m)

mo = EntropyMax(0.1)
dp = DP(θ)