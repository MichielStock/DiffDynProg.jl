

@testset "Dynamic Programming" begin
    n, m = 10, 8
    s, t = randn(n), randn(m)
    θ = s .- t' .|> abs2

    dp = DP(θ)

    @testset "DP data structure" begin
        
        @test size(dp) == (n, m)
        @test eltype(dp) == Float64

        @test size(getD(dp)) == (n, m)
        @test size(getE(dp)) == (n, m)
        @test size(getQ(dp)) == (n, m, 3)

        @test size(getD(dp, 2, 3)) == (2, 3)
        @test size(getE(dp, 2, 3)) == (2, 3)
        @test size(getQ(dp, 2, 3)) == (2, 3, 3)   


        @test eltype(DP(Float16, n, m)) == Float16
    end

    @testset "dynamic time warping" begin
        s1 = dynamic_time_warping(Max(), θ, dp)
        s2 = dynamic_time_warping(EntropyMax(1.0), θ, dp)
        
        D, E = ∂DTW(Max(), θ, dp)

        @test s1 ≈ last(D)
        @test size(E) == (n, m)

        D, E = ∂DTW(EntropyMax(1.0), θ)
        @test s2 ≈ last(D)
        @test size(E) == (n, m)
    end

    @testset "Needleman-Wunsch" begin
        s1 = needleman_wunsch(Max(), θ, 1, dp)
        s2 = needleman_wunsch(EntropyMax(1.0), θ, 1, dp)
        s3 = needleman_wunsch(EntropyMax(1.0), θ, (ones(n), ones(m)), dp)
        
        D, E = ∂NW(Max(), θ, 1, dp)

        @test s1 ≈ last(D)
        @test size(E) == (n, m)

        D, E = ∂NW(EntropyMax(1.0), θ, 1)
        @test s2 ≈ last(D)
        @test size(E) == (n, m)
    end

    @testset "Needleman-Wunsch subtitution matrix" begin
        s = rand(1:10, n)
        t = rand(1:10, m)

        S = randn(10, 10)
        S .+= S'

        θ = S[s, t]
        @test needleman_wunsch(EntropyMax(), s, t, S, 1) ≈ needleman_wunsch(EntropyMax(), θ, 1) 
        @test needleman_wunsch(SquaredMax(), s, t, S, 1, dp) ≈ needleman_wunsch(SquaredMax(), θ, 1, dp) 
    end

    @testset "Smith-Waterman" begin
        using DiffDynProg: logsumexp
        s1 = smith_waterman(Max(), θ, 1, dp)
        s2 = smith_waterman(EntropyMax(1.0), θ, 1, dp)
        s3 = smith_waterman(EntropyMax(1.0), θ, (ones(n), ones(m)), dp)
        
        d, E = ∂SW(Max(), θ, 1, dp)

        @test s1 ≈ d
        @test size(E) == (n, m)

        d, E = ∂SW(EntropyMax(1.0), θ, 1)
        @test s2 ≈ d
        @test size(E) == (n, m)
    end

    @testset "Smith-Waterman subtitution matrix" begin
        s = rand(1:10, n)
        t = rand(1:10, m)

        S = randn(10, 10)
        S .+= S'

        θ = S[s, t]
        @test smith_waterman(EntropyMax(), s, t, S, 1) ≈ smith_waterman(EntropyMax(), θ, 1) 
        @test smith_waterman(SquaredMax(), s, t, S, 1, dp) ≈ smith_waterman(SquaredMax(), θ, 1, dp) 
    end

end