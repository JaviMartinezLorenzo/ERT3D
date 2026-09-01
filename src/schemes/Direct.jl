"""
    spatial_operator!(out, ::Direct, state, D, grid, params)

Compute the spatial contribution of the compressible Euler equations
using their direct conservative formulation.

The semi-discrete equations are written as

    ∂Q/∂t = -spatial_operator!(...)

where `out` contains the spatial terms for

    rho
    rhou
    rhov
    rhow
    rhoE

The derivative operator `D` is independent of the spatial formulation.
It may be a `Central4`, `Central6`, `Central8`, or `Spectral` operator.
"""
function spatial_operator!(
    out::State,
    ::Direct,
    state::State,
    D::DerivativeOperator,
    grid::Grid,
    params::Parameters,
)

    # ------------------------------------------------------------
    # Primitive variables
    # ------------------------------------------------------------

    primitive = primitive_variables(state, params)

    rho = primitive.rho
    u = primitive.u
    v = primitive.v
    w = primitive.w
    p = primitive.p

    rhoE = state.rhoE


    # ------------------------------------------------------------
    # Workspace
    # ------------------------------------------------------------

    tmp = similar(rho)


    # ------------------------------------------------------------
    # Continuity
    #
    # Rrho =
    #     ∂(rhou)/∂x
    #   + ∂(rhov)/∂y
    #   + ∂(rhow)/∂z
    # ------------------------------------------------------------

    derivative_x!(
        out.rho,
        rho .* u,
        D,
        grid,
    )

    derivative_y!(
        tmp,
        rho .* v,
        D,
        grid,
    )

    out.rho .+= tmp

    derivative_z!(
        tmp,
        rho .* w,
        D,
        grid,
    )

    out.rho .+= tmp


    # ------------------------------------------------------------
    # x-momentum
    #
    # Rrhou =
    #     ∂(rhou² + p)/∂x
    #   + ∂(rhouv)/∂y
    #   + ∂(rhouw)/∂z
    # ------------------------------------------------------------

    derivative_x!(
        out.rhou,
        rho .* u .* u .+ p,
        D,
        grid,
    )

    derivative_y!(
        tmp,
        rho .* u .* v,
        D,
        grid,
    )

    out.rhou .+= tmp

    derivative_z!(
        tmp,
        rho .* u .* w,
        D,
        grid,
    )

    out.rhou .+= tmp


    # ------------------------------------------------------------
    # y-momentum
    #
    # Rrhov =
    #     ∂(rhouv)/∂x
    #   + ∂(rhov² + p)/∂y
    #   + ∂(rhovw)/∂z
    # ------------------------------------------------------------

    derivative_x!(
        out.rhov,
        rho .* u .* v,
        D,
        grid,
    )

    derivative_y!(
        tmp,
        rho .* v .* v .+ p,
        D,
        grid,
    )

    out.rhov .+= tmp

    derivative_z!(
        tmp,
        rho .* v .* w,
        D,
        grid,
    )

    out.rhov .+= tmp


    # ------------------------------------------------------------
    # z-momentum
    #
    # Rrhow =
    #     ∂(rhouw)/∂x
    #   + ∂(rhovw)/∂y
    #   + ∂(rhow² + p)/∂z
    # ------------------------------------------------------------

    derivative_x!(
        out.rhow,
        rho .* u .* w,
        D,
        grid,
    )

    derivative_y!(
        tmp,
        rho .* v .* w,
        D,
        grid,
    )

    out.rhow .+= tmp

    derivative_z!(
        tmp,
        rho .* w .* w .+ p,
        D,
        grid,
    )

    out.rhow .+= tmp


    # ------------------------------------------------------------
    # Energy
    #
    # RE =
    #     ∂((rhoE + p)u)/∂x
    #   + ∂((rhoE + p)v)/∂y
    #   + ∂((rhoE + p)w)/∂z
    # ------------------------------------------------------------

    derivative_x!(
        out.rhoE,
        (rhoE .+ p) .* u,
        D,
        grid,
    )

    derivative_y!(
        tmp,
        (rhoE .+ p) .* v,
        D,
        grid,
    )

    out.rhoE .+= tmp

    derivative_z!(
        tmp,
        (rhoE .+ p) .* w,
        D,
        grid,
    )

    out.rhoE .+= tmp


    return out
end