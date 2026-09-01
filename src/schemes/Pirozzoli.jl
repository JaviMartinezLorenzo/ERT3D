"""
    spatial_operator!(out, ::Pirozzoli, state, D, grid, params)

Compute the spatial contribution of the compressible Euler equations
using the Pirozzoli split-form spatial discretization.

The semi-discrete equations are written as

    ∂Q/∂t = -spatial_operator!(...)

where `out` contains the spatial contributions for the five conserved
variables:

    rho
    rhou
    rhov
    rhow
    rhoE

The convective terms are formulated using the generic split form

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
"""
function pirozzoli!(
    out::Array{Float64,3},
    rho::Array{Float64,3},
    u::Array{Float64,3},
    v::Array{Float64,3},
    w::Array{Float64,3},
    phi::Array{Float64,3},
    D::DerivativeOperator,
    grid::Grid,
)

    # Pre-allocate working arrays to eliminate broadcast allocations
    tmp = similar(rho)
    flux = similar(rho)

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

function pirozzoli!(
    out::Array{Float64,3},
    rho::Array{Float64,3},
    u::Array{Float64,3},
    v::Array{Float64,3},
    w::Array{Float64,3},
    phi::Float64,
    D::DerivativeOperator,
    grid::Grid,
)

    @assert phi == 1.0

    # Pre-allocate working arrays for the continuity equation
    tmp = similar(rho)
    flux = similar(rho)

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
    ::Pirozzoli,
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
    u   = primitive.u
    v   = primitive.v
    w   = primitive.w
    p   = primitive.p

    rhoE = state.rhoE
    
    # Pre-allocate arrays used directly in spatial_operator!
    tmp = similar(rho)
    H   = similar(rho)

    # In-place broadcast for enthalpy calculation
    @. H = (rhoE + p) / rho


    # ------------------------------------------------------------
    # Continuity φ = 1
    # ------------------------------------------------------------

    pirozzoli!(
        out.rho,
        rho,
        u,
        v,
        w,
        1.0,
        D,
        grid,
    )


    # ------------------------------------------------------------
    # x-momentum φ = u
    # ------------------------------------------------------------

    pirozzoli!(
        out.rhou,
        rho,
        u,
        v,
        w,
        u,
        D,
        grid,
    )

    # Pressure contribution
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

    pirozzoli!(
        out.rhov,
        rho,
        u,
        v,
        w,
        v,
        D,
        grid,
    )

    # Pressure contribution
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

    pirozzoli!(
        out.rhow,
        rho,
        u,
        v,
        w,
        w,
        D,
        grid,
    )

    # Pressure contribution
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

    pirozzoli!(
        out.rhoE,
        rho,
        u,
        v,
        w,
        H,
        D,
        grid,
    )

    return out
end