"""
    spatial_operator!(
        out,
        ::Direct,
        state,
        D,
        grid,
        params,
        workspace,
    )

Compute the spatial contribution of the compressible Euler equations
using the direct conservative formulation.

The semi-discrete equations are written as

    ∂Q/∂t = -spatial_operator!(...)

where `out` contains the spatial contributions for the five conserved
variables:

    rho
    rhou
    rhov
    rhow
    rhoE

The convective fluxes are discretized directly as

    ∂(ρu_j φ)/∂x_j

using the supplied `DerivativeOperator`.

The spatial derivative operator `D` is independent of the spatial
formulation and may be a `Central4`, `Central6`, `Central8`, or
`Spectral` operator.

The `workspace` provides reusable storage for primitive variables and
temporary arrays, avoiding repeated allocation during spatial
operator evaluations.
"""
function spatial_operator!(
    out::State,
    ::Direct,
    state::State,
    D::DerivativeOperator,
    grid::Grid,
    params::Parameters,
    workspace::SpatialWorkspace,
)

    # ------------------------------------------------------------
    # Primitive variables
    # ------------------------------------------------------------

    primitive_variables!(
        workspace.primitive,
        state,
        params,
    )

    primitive = workspace.primitive

    rho = primitive.rho
    u   = primitive.u
    v   = primitive.v
    w   = primitive.w
    p   = primitive.p

    rhoE = state.rhoE


    # ------------------------------------------------------------
    # Workspace
    # ------------------------------------------------------------

    tmp  = workspace.tmp
    flux = workspace.flux


    # ============================================================
    # Continuity
    # ============================================================

    @. flux = rho * u
    derivative_x!(out.rho, flux, D, grid)

    @. flux = rho * v
    derivative_y!(tmp, flux, D, grid)
    @. out.rho += tmp

    @. flux = rho * w
    derivative_z!(tmp, flux, D, grid)
    @. out.rho += tmp


    # ============================================================
    # x-momentum
    # ============================================================

    @. flux = rho * u * u + p
    derivative_x!(out.rhou, flux, D, grid)

    @. flux = rho * u * v
    derivative_y!(tmp, flux, D, grid)
    @. out.rhou += tmp

    @. flux = rho * u * w
    derivative_z!(tmp, flux, D, grid)
    @. out.rhou += tmp


    # ============================================================
    # y-momentum
    # ============================================================

    @. flux = rho * u * v
    derivative_x!(out.rhov, flux, D, grid)

    @. flux = rho * v * v + p
    derivative_y!(tmp, flux, D, grid)
    @. out.rhov += tmp

    @. flux = rho * v * w
    derivative_z!(tmp, flux, D, grid)
    @. out.rhov += tmp


    # ============================================================
    # z-momentum
    # ============================================================

    @. flux = rho * u * w
    derivative_x!(out.rhow, flux, D, grid)

    @. flux = rho * v * w
    derivative_y!(tmp, flux, D, grid)
    @. out.rhow += tmp

    @. flux = rho * w * w + p
    derivative_z!(tmp, flux, D, grid)
    @. out.rhow += tmp


    # ============================================================
    # Energy
    # ============================================================

    @. flux = (rhoE + p) * u
    derivative_x!(out.rhoE, flux, D, grid)

    @. flux = (rhoE + p) * v
    derivative_y!(tmp, flux, D, grid)
    @. out.rhoE += tmp

    @. flux = (rhoE + p) * w
    derivative_z!(tmp, flux, D, grid)
    @. out.rhoE += tmp

    @. out.rho  = -out.rho
    @. out.rhou = -out.rhou
    @. out.rhov = -out.rhov
    @. out.rhow = -out.rhow
    @. out.rhoE = -out.rhoE
    
    return out
end