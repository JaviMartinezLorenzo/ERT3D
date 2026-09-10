"""
    spatial_operator!(
        out,
        ::Feiereisen,
        state,
        D,
        grid,
        params,
        workspace,
    )

Compute the spatial contribution of the compressible Euler equations
using the Feiereisen et al. split formulation.

The generic convective term is written as

    ∂(ρuⱼφ)/∂xⱼ =
        1/2 ∂(ρuⱼφ)/∂xⱼ
      + 1/2 φ ∂(ρuⱼ)/∂xⱼ
      + 1/2 ρuⱼ ∂φ/∂xⱼ.

The scalar `φ` is chosen as:

- `φ = 1`        for continuity
- `φ = u, v, w`  for the three momentum equations
- `φ = H`        for the energy equation

where `H` is the total specific enthalpy,

    H = E + p/ρ.

The spatial derivative operator is supplied independently through
a `DerivativeOperator`.

The `workspace` provides reusable storage for primitive variables
and temporary arrays, avoiding repeated allocation during spatial
operator evaluations.
"""
function feiereisen!(
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

    # ------------------------------------------------------------
    # x-direction
    # ------------------------------------------------------------

    @. flux = rho * u * phi
    derivative_x!(tmp, flux, D, grid)
    @. out = 0.5 * tmp

    @. flux = rho * u
    derivative_x!(tmp, flux, D, grid)
    @. out += 0.5 * phi * tmp

    @. flux = phi
    derivative_x!(tmp, flux, D, grid)
    @. out += 0.5 * rho * u * tmp


    # ------------------------------------------------------------
    # y-direction
    # ------------------------------------------------------------

    @. flux = rho * v * phi
    derivative_y!(tmp, flux, D, grid)
    @. out += 0.5 * tmp

    @. flux = rho * v
    derivative_y!(tmp, flux, D, grid)
    @. out += 0.5 * phi * tmp

    @. flux = phi
    derivative_y!(tmp, flux, D, grid)
    @. out += 0.5 * rho * v * tmp


    # ------------------------------------------------------------
    # z-direction
    # ------------------------------------------------------------

    @. flux = rho * w * phi
    derivative_z!(tmp, flux, D, grid)
    @. out += 0.5 * tmp

    @. flux = rho * w
    derivative_z!(tmp, flux, D, grid)
    @. out += 0.5 * phi * tmp

    @. flux = phi
    derivative_z!(tmp, flux, D, grid)
    @. out += 0.5 * rho * w * tmp

    return out
end


function feiereisen!(
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

    @. flux = rho * u * phi
    derivative_x!(tmp, flux, D, grid)
    @. out = 0.5 * tmp

    @. flux = rho * u
    derivative_x!(tmp, flux, D, grid)
    @. out += 0.5 * phi * tmp


    # ------------------------------------------------------------
    # y-direction
    # ------------------------------------------------------------

    @. flux = rho * v * phi
    derivative_y!(tmp, flux, D, grid)
    @. out += 0.5 * tmp

    @. flux = rho * v
    derivative_y!(tmp, flux, D, grid)
    @. out += 0.5 * phi * tmp


    # ------------------------------------------------------------
    # z-direction
    # ------------------------------------------------------------

    @. flux = rho * w * phi
    derivative_z!(tmp, flux, D, grid)
    @. out += 0.5 * tmp

    @. flux = rho * w
    derivative_z!(tmp, flux, D, grid)
    @. out += 0.5 * phi * tmp

    return out
end


function spatial_operator!(
    out::State,
    ::Feiereisen,
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

    feiereisen!(
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

    feiereisen!(
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

    feiereisen!(
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

    feiereisen!(
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

    feiereisen!(
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