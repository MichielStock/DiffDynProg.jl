#=
Created on 26/02/2021 15:21:24
Last update: 09 March 2021

@author: Michiel Stock
michielfmstock@gmail.com

Example of differentiable dynamic time warping. Uses a simple alignment
in time.
=#

using DrWatson
@quickactivate "DiffDynProg"

using DiffDynProg
using Plots, StatsBase
using Zygote

dx = 0.4
σ = 5e-2

s = 1.5cos.(0:dx:4π)
t = sin.(-1:dx:(pi+1))

s .+= σ * randn(length(s))
t .+= σ * randn(length(t))

n, m = length(s), length(t)

ts = scatter(s, label="s", xlabel="time")
scatter!(t, label="t")

θ = (s .- t').^2

hθ = heatmap(θ, title="weight matrix", xlabel="i", ylabel="j")

# build a DTW
dtw = DP(θ)

# align
c = dynamic_time_warping(Max(), θ, dtw)
# gradient, here the gradient is just the mapping
c, E_max = ∂DTW(Max(), θ, dtw)
# get filled DP matrix
D = dtw.D[2:end, 2:end]
hD = heatmap(D, title="D (regular max)")
hmax = heatmap(E_max, title="Max")

# entropy max
c, E_em = ∂DTW(EntropyMax(0.1), θ, dtw)
hem = heatmap(E_em, title="Entropy max")

# squared max
c, E_sm = ∂DTW(SquaredMax(.1), θ, dtw)
hsm = heatmap(E_sm, title="Squared max")

for (c, E) in zip(["red"], [E_em])
    for i in 1:n
        for j in 1:m
            E[i,j] > 1e-3 && plot!(ts, [i,j], [s[i], t[j]],
                            color=c, alpha=0.4E[i,j], label="")
        end
    end
end
ts

# combine plots
plot(ts, hθ, hD, hmax, hem, hsm, size=(1200, 600))

savefig(plotsdir("dtw_test.png"))


# computing a loss function

L(θ) = dynamic_time_warping(EntropyMax(0.1), θ, dtw) + 0.01sum(abs2, θ)

L(θ)  # function
L'(θ)  # gradient

# why not optimizing the sequences?

L(s, t) = dynamic_time_warping(EntropyMax(0.1), (s .- t').^2, dtw) + 0.01sum(abs2, θ) +
                     0.001sum(abs2, s.-mean(s)) +  # reg on s and t
                     0.001sum(abs2, t.-mean(t))

∇s, ∇t = gradient(L, s, t)

scatter(s, label="s", title="time series + gradient")
quiver!(1:n,s, quiver=(zeros(n), ∇s))
scatter!(t, label="t")
quiver!(t, quiver=(zeros(m), ∇t))
savefig(plotsdir("dtw_gradient.png"))

# with a gapcost
Lfixedgap(θ, cost) = dynamic_time_warping(EntropyMax(1.0), θ + cost * gap_cost_matrix(n, m), dtw)

∇θ, dc = gradient(Lfixedgap, θ, 0.4)
heatmap(∇θ)

# with a variabel gap cost
Lvargap(θ, cs, ct) = dynamic_time_warping(EntropyMax(.1), θ + gap_cost_matrix(cs, ct), dtw)
cs, ct = 0.01ones(n), 0.01*ones(m)
∇θ, dcs, dct = gradient(Lvargap, θ, cs, ct)