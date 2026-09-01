using ERT3D
using Printf


# ================================================================
# 1. Constant-state preservation
# ================================================================

function test_constant_state()

    println()
    println("1. Constant-state preservation")
    println("================================")

    for N in (16, 32, 64)

        grid = Grid(N)
        params = ERT3D.Parameters(1.4, 0.07)

        state = State(grid)
        out = State(grid)
        primitive = PrimitiveState(grid)

        # Uniform physical state
        @. primitive.rho = 1.0
        @. primitive.u   = 1.0
        @. primitive.v   = 0.5
        @. primitive.w   = -0.25
        @. primitive.p   = 1.0 / params.gamma

        state = conserved_variables(primitive, params)

        D = Central4(grid)

        spatial_operator!(
            out,
            Pirozzoli(),
            state,
            D,
            grid,
            params,
        )

        max_error = maximum([
            maximum(abs.(out.rho)),
            maximum(abs.(out.rhou)),
            maximum(abs.(out.rhov)),
            maximum(abs.(out.rhow)),
            maximum(abs.(out.rhoE)),
        ])

        @printf(
            "N = %3d    max spatial residual = %.6e\n",
            N,
            max_error,
        )
    end
end


# ================================================================
# 2. Global conservation
# ================================================================

function test_global_conservation()

    println()
    println("2. Global conservation")
    println("=======================")

    for N in (16, 32, 64)

        grid = Grid(N)
        params = Parameters(1.4, 0.07)

        # Smooth periodic TGV state
        state = taylor_green_ic(grid, params)
        out = State(grid)

        D = Central4(grid)

        spatial_operator!(
            out,
            Pirozzoli(),
            state,
            D,
            grid,
            params,
        )

        sums = (
            sum(out.rho),
            sum(out.rhou),
            sum(out.rhov),
            sum(out.rhow),
            sum(out.rhoE),
        )

        max_sum = maximum(abs.(collect(sums)))

        @printf(
            "N = %3d    max |sum(spatial residual)| = %.6e\n",
            N,
            max_sum,
        )
    end
end


# ================================================================
# 3. Direct vs Pirozzoli
# ================================================================

function test_direct_vs_pirozzoli()

    println()
    println("3. Direct vs Pirozzoli")
    println("=======================")

    for N in (16, 32, 64)

        grid = Grid(N)
        params = Parameters(1.4, 0.07)

        state = taylor_green_ic(grid, params)

        direct = State(grid)
        pirozzoli = State(grid)

        D = Central4(grid)

        spatial_operator!(
            direct,
            Direct(),
            state,
            D,
            grid,
            params,
        )

        spatial_operator!(
            pirozzoli,
            Pirozzoli(),
            state,
            D,
            grid,
            params,
        )

        difference = maximum([
            maximum(abs.(direct.rho   .- pirozzoli.rho)),
            maximum(abs.(direct.rhou  .- pirozzoli.rhou)),
            maximum(abs.(direct.rhov  .- pirozzoli.rhov)),
            maximum(abs.(direct.rhow  .- pirozzoli.rhow)),
            maximum(abs.(direct.rhoE  .- pirozzoli.rhoE)),
        ])

        @printf(
            "N = %3d    ||Direct - Pirozzoli||∞ = %.6e\n",
            N,
            difference,
        )
    end
end


# ================================================================
# Run tests
# ================================================================

test_constant_state()
test_global_conservation()
test_direct_vs_pirozzoli()