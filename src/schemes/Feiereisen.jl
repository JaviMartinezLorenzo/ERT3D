"""
    spatial_operator!(out, ::Feiereisen, state, D, grid, params)

Compute the spatial contribution of the compressible Euler equations
using the Feiereisen et al. split formulation.

The generic convective term is written as

    ∂(ρuⱼφ)/∂xⱼ =
        1/2 ∂(ρuⱼφ)/∂xⱼ
      + 1/2 φ ∂(ρuⱼ)/∂xⱼ
      + 1/2 ρuⱼ ∂φ/∂xⱼ.

The scalar φ is chosen as:

- φ = 1        for continuity
- φ = u, v, w  for the three momentum equations
- φ = H        for the energy equation

where H is the total specific enthalpy,

    H = E + p/ρ.

The spatial derivative operator is supplied independently through
a `DerivativeOperator`.
"""
function feiereisen!(
    out::Array{Float64,3},
    rho::Array{Float64,3},
    u::Array{Float64,3},
    v::Array{Float64,3},
    w::Array{Float64,3},
    phi::Union{Float64,Array{Float64,3}},
    D::DerivativeOperator,
    grid::Grid,
)

    tmp = similar(rho)

    # ------------------------------------------------------------
    # x-direction
    # ------------------------------------------------------------

    derivative_x!(tmp, rho .* u .* phi, D, grid)
    out .= 0.5 .* tmp

    derivative_x!(tmp, rho .* u, D, grid)
    out .+= 0.5 .* phi .* tmp

    derivative_x!(tmp, phi, D, grid)
    out .+= 0.5 .* rho .* u .* tmp


    # ------------------------------------------------------------
    # y-direction
    # ------------------------------------------------------------

    derivative_y!(tmp, rho .* v .* phi, D, grid)
    out .+= 0.5 .* tmp

    derivative_y!(tmp, rho .* v, D, grid)
    out .+= 0.5 .* phi .* tmp

    derivative_y!(tmp, phi, D, grid)
    out .+= 0.5 .* rho .* v .* tmp


    # ------------------------------------------------------------
    # z-direction
    # ------------------------------------------------------------

    derivative_z!(tmp, rho .* w .* phi, D, grid)
    out .+= 0.5 .* tmp

    derivative_z!(tmp, rho .* w, D, grid)
    out .+= 0.5 .* phi .* tmp

    derivative_z!(tmp, phi, D, grid)
    out .+= 0.5 .* rho .* w .* tmp

    return out
end

function spatial_operator!(
    out::State,
    ::Feiereisen,
    state::State,
    D::DerivativeOperator,
    grid::Grid,
    params::Parameters,
)
    # implementation
end