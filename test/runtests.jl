using Test
using ERT3D

@testset "ERT3D" begin

    @testset "Grid" begin
        # TODO: Δx = 2π/N for a few N, coordinate array lengths, etc.
    end

    @testset "ConservativeFE flux" begin
        # TODO: verify eq. (13) reduces to expected flux for a simple
        #       hand-computable case (e.g. constant field -> zero flux divergence)
    end

    @testset "ExplicitRK3 vs ImplicitMidpoint structure" begin
        # TODO: verify ImplicitMidpoint is self-adjoint to machine precision
        #       on a trivial linear test problem (step forward then backward,
        #       check ~1e-14 recovery) BEFORE trusting it on the full nonlinear system
    end

    @testset "Energy conservation (correctness gate, small case)" begin
        # TODO: small-N, short-time synthetic turbulence run,
        #       assert kinetic energy drift below a tolerance
    end

end
