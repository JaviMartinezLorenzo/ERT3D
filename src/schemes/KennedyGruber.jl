"""
    spatial_operator!(
        out,
        ::KennedyGruber,
        state,
        D,
        grid,
        params,
        workspace,
    )

Compute the spatial contribution of the compressible Euler equations
using the Direct Kennedy–Gruber split-form spatial discretization.

The semi-discrete equations are written as

    ∂Q/∂t = -spatial_operator!(...)

where `out` contains the spatial contributions for the five conserved
variables:

    rho
    rhou
    rhov
    rhow
    rhoE

The convective terms are formulated using the generalized split form

    ∂(ρuⱼφ)/∂xⱼ =
        1/4 ∂(ρuⱼφ)/∂xⱼ
      + 1/4 (uⱼ ∂(ρφ)/∂xⱼ
           + ρ ∂(uⱼφ)/∂xⱼ
           + φ ∂(ρuⱼ)/∂xⱼ)
      + 1/4 (ρuⱼ ∂φ/∂xⱼ
           + ρφ ∂uⱼ/∂xⱼ
           + uⱼφ ∂ρ/∂xⱼ).

Here `j` denotes the spatial directions `x`, `y`, and `z`, and the
summation over `j` is evaluated explicitly in the implementation.

The scalar `φ` is chosen according to the conserved equation:

- `φ = 1`        for continuity
- `φ = u`        for x-momentum
- `φ = v`        for y-momentum
- `φ = w`        for z-momentum
- `φ = H`        for energy

where `H` is the total specific enthalpy,

    H = E + p/ρ
      = (ρE + p)/ρ.

For the momentum equations, the pressure contribution is added
separately as

    ∂(p δᵢⱼ)/∂xⱼ = ∂p/∂xᵢ.

The `workspace` provides reusable storage for primitive variables
and temporary arrays.
"""
function KennedyGruber!(
    out::Array{Float64,3},
    rho::Array{Float64,3},
    u::Array{Float64,3},
    v::Array{Float64,3},
    w::Array{Float64,3},
    phi::Array{Float64,3},
    D::DerivativeOperator,
    grid::Grid,
    workspace::SpatialWorkspace,
)

    tmp  = workspace.tmp
    flux = workspace.flux

    # ============================================================
    # x-direction
    # ============================================================

    @. flux = rho * u * phi
    derivative_x!(tmp, flux, D, grid)
    @. out = 0.25 * tmp

    @. flux = rho * phi
    derivative_x!(tmp, flux, D, grid)
    @. out += 0.25 * u * tmp

    @. flux = u * phi
    derivative_x!(tmp, flux, D, grid)
    @. out += 0.25 * rho * tmp

    @. flux = rho * u
    derivative_x!(tmp, flux, D, grid)
    @. out += 0.25 * phi * tmp

    derivative_x!(tmp, phi, D, grid)
    @. out += 0.25 * rho * u * tmp

    derivative_x!(tmp, u, D, grid)
    @. out += 0.25 * rho * phi * tmp

    derivative_x!(tmp, rho, D, grid)
    @. out += 0.25 * u * phi * tmp


    # ============================================================
    # y-direction
    # ============================================================

    @. flux = rho * v * phi
    derivative_y!(tmp, flux, D, grid)
    @. out += 0.25 * tmp

    @. flux = rho * phi
    derivative_y!(tmp, flux, D, grid)
    @. out += 0.25 * v * tmp

    @. flux = v * phi
    derivative_y!(tmp, flux, D, grid)
    @. out += 0.25 * rho * tmp

    @. flux = rho * v
    derivative_y!(tmp, flux, D, grid)
    @. out += 0.25 * phi * tmp

    derivative_y!(tmp, phi, D, grid)
    @. out += 0.25 * rho * v * tmp

    derivative_y!(tmp, v, D, grid)
    @. out += 0.25 * rho * phi * tmp

    derivative_y!(tmp, rho, D, grid)
    @. out += 0.25 * v * phi * tmp


    # ============================================================
    # z-direction
    # ============================================================

    @. flux = rho * w * phi
    derivative_z!(tmp, flux, D, grid)
    @. out += 0.25 * tmp

    @. flux = rho * phi
    derivative_z!(tmp, flux, D, grid)
    @. out += 0.25 * w * tmp

    @. flux = w * phi
    derivative_z!(tmp, flux, D, grid)
    @. out += 0.25 * rho * tmp

    @. flux = rho * w
    derivative_z!(tmp, flux, D, grid)
    @. out += 0.25 * phi * tmp

    derivative_z!(tmp, phi, D, grid)
    @. out += 0.25 * rho * w * tmp

    derivative_z!(tmp, w, D, grid)
    @. out += 0.25 * rho * phi * tmp

    derivative_z!(tmp, rho, D, grid)
    @. out += 0.25 * w * phi * tmp

    return out
end


function KennedyGruber!(
    out::Array{Float64,3},
    rho::Array{Float64,3},
    u::Array{Float64,3},
    v::Array{Float64,3},
    w::Array{Float64,3},
    phi::Float64,
    D::DerivativeOperator,
    grid::Grid,
    workspace::SpatialWorkspace,
)

    @assert phi == 1.0

    tmp  = workspace.tmp
    flux = workspace.flux

    # ------------------------------------------------------------
    # x-direction
    # ------------------------------------------------------------

    @. flux = rho * u
    derivative_x!(tmp, flux, D, grid)
    @. out = 0.5 * tmp

    derivative_x!(tmp, rho, D, grid)
    @. out += 0.5 * u * tmp

    derivative_x!(tmp, u, D, grid)
    @. out += 0.5 * rho * tmp


    # ------------------------------------------------------------
    # y-direction
    # ------------------------------------------------------------

    @. flux = rho * v
    derivative_y!(tmp, flux, D, grid)
    @. out += 0.5 * tmp

    derivative_y!(tmp, rho, D, grid)
    @. out += 0.5 * v * tmp

    derivative_y!(tmp, v, D, grid)
    @. out += 0.5 * rho * tmp


    # ------------------------------------------------------------
    # z-direction
    # ------------------------------------------------------------

    @. flux = rho * w
    derivative_z!(tmp, flux, D, grid)
    @. out += 0.5 * tmp

    derivative_z!(tmp, rho, D, grid)
    @. out += 0.5 * w * tmp

    derivative_z!(tmp, w, D, grid)
    @. out += 0.5 * rho * tmp

    return out
end


function spatial_operator!(
    out::State,
    ::KennedyGruber,
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

    tmp = workspace.tmp
    H   = workspace.H

    @. H = (rhoE + p) / rho


    # ------------------------------------------------------------
    # Continuity φ = 1
    # ------------------------------------------------------------

    KennedyGruber!(
        out.rho,
        rho,
        u,
        v,
        w,
        1.0,
        D,
        grid,
        workspace,
    )


    # ------------------------------------------------------------
    # x-momentum φ = u
    # ------------------------------------------------------------

    KennedyGruber!(
        out.rhou,
        rho,
        u,
        v,
        w,
        u,
        D,
        grid,
        workspace,
    )

    derivative_x!(
        tmp,
        p,
        D,
        grid,
    )

    @. out.rhou += tmp


    # ------------------------------------------------------------
    # y-momentum φ = v
    # ------------------------------------------------------------

    KennedyGruber!(
        out.rhov,
        rho,
        u,
        v,
        w,
        v,
        D,
        grid,
        workspace,
    )

    derivative_y!(
        tmp,
        p,
        D,
        grid,
    )

    @. out.rhov += tmp


    # ------------------------------------------------------------
    # z-momentum φ = w
    # ------------------------------------------------------------

    KennedyGruber!(
        out.rhow,
        rho,
        u,
        v,
        w,
        w,
        D,
        grid,
        workspace,
    )

    derivative_z!(
        tmp,
        p,
        D,
        grid,
    )

    @. out.rhow += tmp


    # ------------------------------------------------------------
    # Energy φ = H
    # ------------------------------------------------------------

    KennedyGruber!(
        out.rhoE,
        rho,
        u,
        v,
        w,
        H,
        D,
        grid,
        workspace,
    )

   
    @. out.rho  = -out.rho
    @. out.rhou = -out.rhou
    @. out.rhov = -out.rhov
    @. out.rhow = -out.rhow
    @. out.rhoE = -out.rhoE

    return out
end