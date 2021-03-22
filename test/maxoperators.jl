@testset "maxoperators" begin

    @testset "maximum" begin 
        x = [1.0, 2.0, 5.0]
        y = [-1.2, 3.3, 5, -maxintfloat(Float64), 0.1]


        for mo in [Max(), LeakyMax(), EntropyMax(), SquaredMax()]
            @test maxᵧ(mo, x) isa Number
            @test maxᵧ(mo, x .+ 1) ≈ maxᵧ(mo, x) + 1

            @test maxᵧ(mo, y) isa Number
            @test maxᵧ(mo, y.+ 1) ≈ maxᵧ(mo, y) + 1
            @test isfinite(maxᵧ(mo, y))

            # max and argmax
            m, q = max_argmaxᵧ(mo, x)
            @test m isa Number
            @test m ≈ maxᵧ(mo, x)
            @test sum(q) ≈ 1.0 && all(q .≥ 0.0)
        end
    end

    @testset "maximum (tuple0" begin 
        x = (1.0, 2.0, 5.0)
        y = (-1.2, 3.3, 5, -maxintfloat(Float64), 0.1)


        for mo in [Max(), LeakyMax(), EntropyMax(), SquaredMax()]
            @test maxᵧ(mo, x) isa Number
            @test maxᵧ(mo, x .+ 1) ≈ maxᵧ(mo, x) + 1

            @test maxᵧ(mo, y) isa Number
            @test maxᵧ(mo, y.+ 1) ≈ maxᵧ(mo, y) + 1
            @test isfinite(maxᵧ(mo, y))

            # max and argmax
            m, q = max_argmaxᵧ(mo, x)
            @test m isa Number
            @test m ≈ maxᵧ(mo, x)
            @test sum(q) ≈ 1.0 && all(q .≥ 0.0)
        end
    end


    @testset "minimum" begin 
    x = -[1.0, 2.0, 5.0]
    y = -[-1.2, 3.3, 5, maxintfloat(Float64), 0.1]


    for mo in [Max(), LeakyMax(), EntropyMax(), SquaredMax()]
        @test minᵧ(mo, x) isa Number
        @test minᵧ(mo, x .+ 1) ≈ minᵧ(mo, x) + 1

        @test minᵧ(mo, y) isa Number
        @test minᵧ(mo, y.+ 1) ≈ minᵧ(mo, y) + 1
        @test isfinite(minᵧ(mo, y))

        # max and argmax
        m, q = min_argminᵧ(mo, x)
        @test m isa Number
        @test m ≈ minᵧ(mo, x)
        @test sum(q) ≈ 1.0 && all(q .≥ 0.0)
    end
end


end