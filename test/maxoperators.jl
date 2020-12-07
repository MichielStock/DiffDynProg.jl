@testset "maxoperators" begin


    x = [1.0, 2.0, 5.0]

    for mo in [Max(), LeakyMax(), EntropyMax(), SquaredMax()]
        @test maximum(mo, x) isa Number
        @test maximum(mo, x .+ 1) ≈ maximum(mo, x) + 1
        @test isfinite(maximum(mo, [-Inf, 2, 3, 9]))
    end
end