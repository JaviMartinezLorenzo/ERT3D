using ERT3D
using Printf


# ================================================================
# 1. Constant-state preservation
# ================================================================

function test_constant_state(scheme; name)

    println()
    println("1. Constant-state preservation — $name")
    println("========================================")

    for N in (16, 32, 64)

        grid = Grid(N)
        params = Parameters(1.4, 0.07)

        primitive = PrimitiveState(grid)

        # Uniform physical state
        @. primitive.rho = 1.0
        @. primitive.u   = 1.0
        @. primitive.v   = 0.5
        @. primitive.w   = -0.25
        @. primitive.p   = 1.0 / params.gamma

        state = conserved_variables(primitive, params)
        out = State(grid)

        D = Central4(grid)

        spatial_operator!(
            out,
            scheme,
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

function test_global_conservation(scheme; name)

    println()
    println("2. Global conservation — $name")
    println("==============================")

    for N in (16, 32, 64)

        grid = Grid(N)
        params = Parameters(1.4, 0.07)

        # Smooth periodic TGV state
        state = taylor_green_ic(grid, params)
        out = State(grid)

        D = Central4(grid)

        spatial_operator!(
            out,
            scheme,
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
# 3. Direct vs scheme
# ================================================================

function test_direct_vs_scheme(scheme; name)

    println()
    println("3. Direct vs $name")
    println("==================")

    for N in (16, 32, 64)

        grid = Grid(N)
        params = Parameters(1.4, 0.07)

        state = taylor_green_ic(grid, params)

        direct = State(grid)
        split = State(grid)

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
            split,
            scheme,
            state,
            D,
            grid,
            params,
        )

        difference = maximum([
            maximum(abs.(direct.rho  .- split.rho)),
            maximum(abs.(direct.rhou .- split.rhou)),
            maximum(abs.(direct.rhov .- split.rhov)),
            maximum(abs.(direct.rhow .- split.rhow)),
            maximum(abs.(direct.rhoE .- split.rhoE)),
        ])

        @printf(
            "N = %3d    ||Direct - %s||∞ = %.6e\n",
            N,
            name,
            difference,
        )
    end
end


# ================================================================
# Run tests
# ================================================================

test_constant_state(
    KennedyGruber();
    name = "KennedyGruber",
)

test_global_conservation(
    KennedyGruber();
    name = "KennedyGruber",
)

test_direct_vs_scheme(
    KennedyGruber();
    name = "KennedyGruber",
)


test_constant_state(
    Feiereisen();
    name = "Feiereisen",
)

test_global_conservation(
    Feiereisen();
    name = "Feiereisen",
)

test_direct_vs_scheme(
    Feiereisen();
    name = "Feiereisen",
)