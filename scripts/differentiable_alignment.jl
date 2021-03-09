#=
Created on 05/03/2021 13:54:10
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Proof of concept of differentiable sequence alignment.

A sequence is just an array of integers. We simulate variations using a substitution matrix
and the probability of insertion/duplication.
=#

using DiffDynProg
using StatsBase
using Zygote

T = [0.9 0.05 0.02 0.03
    0.2 0.7 0.05 0.05;
    0.01 0.01 0.96 0.02;
    0.01 0.01 0.3 0.68]

l = 50
n = 25
k = size(T, 1)

progenitor_seq = rand(1:k, l)

function mutate(sequence, T; pdel=2e-2, pinsert=1e-2, pextend=0.6)
    seq_mutated = Int[]
    k = size(T, 1)
    l = length(sequence)
    for i in 1:l
        if rand() < pdel
            continue
        elseif rand() < pinsert
            while true
                push!(seq_mutated, rand(1:k))
                rand() > pextend && break
            end
        else
            c = sequence[i]
            push!(seq_mutated, sample(1:k, Weights(T[c,:])))
        end
    end
    return seq_mutated
end

const sequences = [mutate(progenitor_seq, T) for i in 1:n]

lmax = maximum(length.(sequences))
dp = DP(Float32, lmax, lmax)

S = zeros(Float32, k, k)
c = Float32(1)



function loss(S; λ=1e-1)
    cost = zero(eltype(S))
    n = length(S)
    lmax = maximum(length.(sequences))
    for (i, seq1) in enumerate(sequences)
        for (j, seq2) in enumerate(sequences[(i+1):end])
            θ = S[seq1,seq2]
            cost -= needleman_wunsch(EntropyMax(Float32(1)), θ, 1, dp)
        end
    end
    return cost / (n * (n-1) / 2) + λ * sum(abs2, S)
end

ΔS = zero(S)
momentum = 0.8
stepsize = 0.03

for t in 1:500
    println("step $t: loss=$(loss(S))")
    ΔS .= momentum .* ΔS + (1.0 - momentum) .* loss'(S)
    S .-= stepsize * (ΔS .+ ΔS')/2 
end