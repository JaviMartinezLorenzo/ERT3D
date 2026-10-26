using ERT3D
using Test

function run_precision_test(::Type{T}) where {T<:AbstractFloat}
    @testset "ERT3D precision: $T" begin
        N = 32

        # --------------------------------------------------------
        # Setup
        # --------------------------------------------------------

        grid = Grid(N; T=T)

        params = Parameters(
            gamma = 1.4,
            Mt0 = 0.1,
            T = T,
        )

        derivative = Central4(grid)
        scheme = KennedyGruber()
        integrator = ExplicitRK3()

        state = State(grid)
        initialize!(state, TaylorGreen(), grid, params)

        workspace = RK3Workspace(grid)

        # --------------------------------------------------------
        # Verify precision before simulation
        # --------------------------------------------------------

        @test eltype(state.rho) == T
        @test eltype(state.rhou) == T
        @test eltype(state.rhov) == T
        @test eltype(state.rhow) == T
        @test eltype(state.rhoE) == T

        @test eltype(derivative.xwork) == T
        @test eltype(derivative.ywork) == T
        @test eltype(derivative.zwork) == T

        @test eltype(params.gamma) == T
        @test eltype(workspace.stage1.rho) == T
        @test eltype(workspace.stage2.rho) == T
        @test eltype(workspace.rhs.rho) == T

        # --------------------------------------------------------
        # Simulation
        # --------------------------------------------------------

        sim = Simulation(
            state = state,
            grid = grid,
            params = params,
            derivative = derivative,
            scheme = scheme,
            integrator = integrator,
            workspace = workspace,
        )

        diagnostics = Diagnostics()

        hooks = [
            OutputHook(
                "diagnostics",
                1,
                (sim, step) -> record!(diagnostics, sim),
            ),
        ]

        dt = 1e-3
        t_end = 0.01

        run!(
            sim,
            t_end;
            dt = dt,
            hooks = hooks,
            verbose = true,
        )

        # --------------------------------------------------------
        # Verify precision and numerical validity after simulation
        # --------------------------------------------------------

        @test sim.t ≈ t_end
        @test eltype(sim.state.rho) == T

        for field in (
            sim.state.rho,
            sim.state.rhou,
            sim.state.rhov,
            sim.state.rhow,
            sim.state.rhoE,
        )
            @test all(isfinite, field)
        end

        @test all(sim.state.rho .> zero(T))

        println("Precision test passed for $T")
        println("Final time: ", sim.t)
        println("State type: ", eltype(sim.state.rho))
    end
end

run_precision_test(Float64)
run_precision_test(Float32)