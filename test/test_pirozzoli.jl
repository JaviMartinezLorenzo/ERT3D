using ERT3D
using CairoMakie
using Printf


# ================================================================
# Test fields
# ================================================================

function test_fields(grid::Grid)

    X = reshape(grid.x, :, 1, 1)
    Y = reshape(grid.y, 1, :, 1)
    Z = reshape(grid.z, 1, 1, :)

    rho =
        1.0 .+
        0.1 .* sin.(X .+ 0.0 .* Y .+ 0.0 .* Z)

    u =
        1.0 .+
        0.2 .* cos.(X .+ 0.0 .* Y .+ 0.0 .* Z)

    v =
        1.0 .+
        0.2 .* sin.(0.0 .* X .+ Y .+ 0.0 .* Z)

    w =
        1.0 .+
        0.2 .* cos.(0.0 .* X .+ 0.0 .* Y .+ Z)

    phi = 1.0 .+ 0.1 .* sin.(X .+ Y .+ Z)

    return rho, u, v, w, phi, X, Y, Z
end


# ================================================================
# Exact continuous Pirozzoli operator
#
# For smooth continuous fields, the Pirozzoli split form is
# algebraically equivalent to the derivative of
#
#     rho*u_j*phi
#
# summed over j = x,y,z.
#
# The numerical Pirozzoli operator is NOT replaced by this form;
# this is only the analytical reference used for the test.
# ================================================================

function exact_pirozzoli(
    rho,
    u,
    v,
    w,
    phi,
    X,
    Y,
    Z,
)

    S = X .+ Y .+ Z

    # ------------------------------------------------------------
    # x-direction
    # ------------------------------------------------------------

    drho_dx = 0.1 .* cos.(X)
    du_dx   = -0.2 .* sin.(X)
    dphi_dx = 0.1 .* cos.(S)

    dx =
        drho_dx .* u .* phi .+
        rho .* du_dx .* phi .+
        rho .* u .* dphi_dx


    # ------------------------------------------------------------
    # y-direction
    # ------------------------------------------------------------

    dv_dy   = 0.2 .* cos.(Y)
    dphi_dy = 0.1 .* cos.(S)

    dy =
        rho .* dv_dy .* phi .+
        rho .* v .* dphi_dy


    # ------------------------------------------------------------
    # z-direction
    # ------------------------------------------------------------

    dw_dz   = -0.2 .* sin.(Z)
    dphi_dz = 0.1 .* cos.(S)

    dz =
        rho .* dw_dz .* phi .+
        rho .* w .* dphi_dz


    return dx .+ dy .+ dz
end


# ================================================================
# Test continuity specialization
#
# Compare the specialized phi = 1 implementation against the
# general Pirozzoli implementation with phi represented by a
# three-dimensional array of ones.
# ================================================================

function continuity_specialization_error(
    grid::Grid,
    D::DerivativeOperator,
)

    rho, u, v, w, _, _, _, _ =
        test_fields(grid)

    phi = ones(size(rho))

    specialized = zeros(size(rho))
    general     = zeros(size(rho))

    ERT3D.pirozzoli!(
        specialized,
        rho,
        u,
        v,
        w,
        1.0,
        D,
        grid,
    )

    ERT3D.pirozzoli!(
        general,
        rho,
        u,
        v,
        w,
        phi,
        D,
        grid,
    )

    return maximum(abs.(specialized .- general))
end


# ================================================================
# Main convergence study
# ================================================================

function convergence_study()

    resolutions = [16, 32, 64, 128]

    errors = Float64[]
    continuity_errors = Float64[]

    println()
    println("Pirozzoli split-form verification")
    println("=================================")
    println()

    for N in resolutions

        # --------------------------------------------------------
        # Grid and derivative operator
        # --------------------------------------------------------

        grid = Grid(N)
        D = Central8(grid)

        # --------------------------------------------------------
        # Test fields
        # --------------------------------------------------------

        rho, u, v, w, phi, X, Y, Z =
            test_fields(grid)

        # --------------------------------------------------------
        # Numerical Pirozzoli operator
        # --------------------------------------------------------

        numerical = zeros(size(rho))

        ERT3D.pirozzoli!(
            numerical,
            rho,
            u,
            v,
            w,
            phi,
            D,
            grid,
        )

        # --------------------------------------------------------
        # Exact analytical result
        # --------------------------------------------------------

        exact = exact_pirozzoli(
            rho,
            u,
            v,
            w,
            phi,
            X,
            Y,
            Z,
        )

        error = maximum(abs.(numerical .- exact))

        push!(errors, error)

        # --------------------------------------------------------
        # Continuity specialization
        # --------------------------------------------------------

        continuity_error =
            continuity_specialization_error(grid, D)

        push!(continuity_errors, continuity_error)

        @printf(
            "N = %3d    general error = %.6e    continuity difference = %.6e\n",
            N,
            error,
            continuity_error,
        )
    end


    # ============================================================
    # Observed convergence order
    # ============================================================

    println()
    println("Observed convergence order")
    println("==========================")

    for i in 1:length(resolutions)-1

        p =
            log(errors[i] / errors[i+1]) /
            log(2.0)

        @printf(
            "N = %3d → %3d:    p = %.6f\n",
            resolutions[i],
            resolutions[i+1],
            p,
        )
    end


    # ============================================================
    # Convergence plot
    # ============================================================

    fig = Figure(size = (900, 650))

    ax = Axis(
        fig[1, 1],
        xlabel = "Grid resolution N",
        ylabel = "L∞ error",
        xscale = log2,
        yscale = log10,
        title = "Pirozzoli split-form convergence — Central4",
    )

    lines!(
        ax,
        resolutions,
        errors,
        linewidth = 2,
        label = "Pirozzoli",
    )

    scatter!(
        ax,
        resolutions,
        errors,
    )

    axislegend(
        ax,
        position = :lb,
    )

    save(
        "test_pirozzoli_convergence.png",
        fig,
        px_per_unit = 2,
    )

    println()
    println("Figure saved to:")
    println("test_pirozzoli_convergence.png")
    println()

    return resolutions, errors, continuity_errors
end


# ================================================================
# Run
# ================================================================

convergence_study()