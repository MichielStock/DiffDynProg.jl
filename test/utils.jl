

@testset "Utilities" begin

    @testset "Gumbel" begin
        using DiffDynProg: randg

        @test randg() isa Number
        @test randg(2, 3) isa Matrix{Float64}
    end

    @testset "project in simplex" begin
        using DiffDynProg: project_in_simplex

        q = project_in_simplex([1.0, 2.0, 3.0], 1.0)
        @test sum(q) ≈ 1 && all(0 .≤ q .≤ 1)
        @test q ≈ project_in_simplex([1.0, 2.0, 3.0] .+1, 1.0)
    end

    @testset "softmax etc" begin
        using DiffDynProg: softmax, logsumexp, gumbel_softmax

        z = randn(8)
        Z = randn(10, 20)
        γ = 0.1
        τ = 0.1

        # vector
        @test  maximum(z) ≤ logsumexp(z; γ) ≤ maximum(z) + γ * log(length(z))
        q = softmax(z; γ)
        @test sum(q) ≈ 1 && all(0 .≤ q .≤ 1)
        y = gumbel_softmax(z; τ)
        @test sum(y) ≈ 1 && all(0 .≤ y .≤ 1)

        # matrix
        @test all(maximum(Z, dims=2) .≤ logsumexp(Z; γ, dims=2) .≤ maximum(Z, dims=2) .+ γ * log(size(Z, 1)))
        Q = softmax(Z, dims=2)
        @test all(sum(Q, dims=2) .≈ 1) && all(0 .≤ Q .≤ 1)
        Y = gumbel_softmax(Z, dims=2; τ)
        @test all(sum(Y, dims=2) .≈ 1) && all(0 .≤ Y .≤ 1)
    end

    @testset "gap costs" begin
        C = gap_cost_matrix(3, 4)
        @test all(diff(C, dims=1) .== -1)
        @test all(diff(C, dims=2) .== -1)
        @test C[1,1] == C[2,2] + 2

        C = gap_cost_matrix(ones(5), [1, 2, 3])
        @test all(diff(C, dims=1) .== -1)
        @test all(diff(C, dims=2)[1,:] .== -[2, 3])
    end
end