@testset "maxoperators" begin

    @testset "maximum" begin 
        x = [1.0, 2.0, 5.0]
        y = [-1.2, 3.3, 5, -Inf, 0.1]


        for mo in [Max(), LeakyMax(), EntropyMax(), SquaredMax()]
            @test maximum(mo, x) isa Number
            @test maximum(mo, x .+ 1) ≈ maximum(mo, x) + 1

            @test maximum(mo, y) isa Number
            @test maximum(mo, y.+ 1) ≈ maximum(mo, y) + 1
            @test isfinite(maximum(mo, y))

            # max and argmax
            m, q = max_argmax(mo, x)
            @test m isa Number
            @test m ≈ maximum(mo, x)
            @test sum(q) ≈ 1.0 && all(q .≥ 0.0)
        end
    end

    @testset "minimum" begin 
    x = -[1.0, 2.0, 5.0]
    y = -[-1.2, 3.3, 5, -Inf, 0.1]


    for mo in [Max(), LeakyMax(), EntropyMax(), SquaredMax()]
        @test minimum(mo, x) isa Number
        @test minimum(mo, x .+ 1) ≈ minimum(mo, x) + 1

        @test minimum(mo, y) isa Number
        @test minimum(mo, y.+ 1) ≈ minimum(mo, y) + 1
        @test isfinite(minimum(mo, y))

        # max and argmax
        m, q = min_argmin(mo, x)
        @test m isa Number
        @test m ≈ minimum(mo, x)
        @test sum(q) ≈ 1.0 && all(q .≥ 0.0)
    end
end


end