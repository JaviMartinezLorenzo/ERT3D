using Test
using ERT3D

@testset "ERT3D" begin

    @testset "Grid" begin

        grid = Grid(32)

        @test grid.N == 32
        @test grid.Δx ≈ 2π/32
        @test length(grid.x) == 32
        @test length(grid.y) == 32
        @test length(grid.z) == 32
        @test length(grid.k) == 32

    end

    @testset "Taylor-Green IC" begin

        grid = Grid(32)
        params = Parameters(1.4)

        state = taylor_green_ic(grid, params, 0.07)

        primitive = primitive_variables(state, params.gamma)

        @test maximum(abs.(primitive.rho .- 1.0)) < 1e-14
        @test maximum(abs.(primitive.w)) < 1e-14

    end

    @testset "Primitive/Conserved conversion" begin

        grid = Grid(16)
        params = Parameters(1.4)

        state = taylor_green_ic(grid, params, 0.07)
        primitive = primitive_variables(state, params.gamma)
        recovered = conserved_variables(primitive, params.gamma)

        @test recovered.rho ≈ state.rho
        @test recovered.rhou ≈ state.rhou
        @test recovered.rhov ≈ state.rhov
        @test recovered.rhow ≈ state.rhow
        @test recovered.rhoE ≈ state.rhoE
    end


end
