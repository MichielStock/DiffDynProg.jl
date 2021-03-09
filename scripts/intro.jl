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

# a more complex example

h(x, y; pow=2) = (2x + y^2)^pow

function rrule(::typeof(h), x, y; pow=2)
    println("rrule")
    return h(x,y;pow), (ȳ) -> (NO_FIELDS, ȳ*2pow * (2x+y^2)^(pow-1), y*2y*pow * (2x+y^2)^(pow-1))
end

gradient(h, 1.0, 2.0)

gradient((x,y)->h(x,y, pow=3), 1.0, 2.0)

g(x, y, a=2) = (2x + y^2)^a

function rrule(::typeof(g), x, y, a=3)
    println("rrule")
    return g(x,y,a), (ȳ) -> (NO_FIELDS, ȳ*a * (2x+y^2)^(a-1), y*2y*a * (2x+y^2)^(a-1), Zero())
end