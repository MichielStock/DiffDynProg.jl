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