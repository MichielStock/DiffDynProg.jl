using DrWatson
@quickactivate "DiffDynProg"
DrWatson.greet()


import ChainRulesCore: frule, rrule, NO_FIELDS
using StatsBase

x = randn(5)

logsumexp(x; γ=1.0) = γ * (x ./ γ .|> exp |> sum |> log)

function softmax(x; γ=1.0)
    return (x .- mean(x)) ./ γ .|> exp |> x -> x ./ sum(x)
end

function rrule(::typeof(logsumexp), x; γ=1.0)
    println("rrule")
    return logsumexp(x; γ), ȳ -> (NO_FIELDS, ȳ * softmax(x; γ))
end

function frule(::typeof(logsumexp), x; γ=1.0)
    println("frule")
    return logsumexp(x; γ), softmax(x; γ)
end

using Zygote

logsumexp(x)

gradient(x-> logsumexp(x, γ=20), x)
