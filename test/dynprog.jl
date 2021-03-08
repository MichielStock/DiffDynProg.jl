n, m = 10, 8
s, t = randn(n), randn(m)
θ = s .- t' .|> abs2

@testset "DP data structure" begin
    s, t = randn(10), randn(8)
    dp = DP(θ)
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