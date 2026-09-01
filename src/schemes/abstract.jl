"""
    FluxScheme

Abstract type for spatial formulations of the compressible Euler
equations.

Concrete subtypes define how the convective terms are formulated.
The spatial derivative operator is supplied separately through a
`DerivativeOperator`, allowing the flux formulation and derivative
discretization to be varied independently.

Implemented formulations include:

- `Direct`      — direct discretization of the original Euler equations.
- `Pirozzoli`   — split formulation following Pirozzoli.
- `Feiereisen`  — formulation following Feiereisen et al.

A concrete subtype must implement:

    spatial_operator!(out, scheme, state, D, grid, params)

where `out` contains the spatial contribution to the five conserved
equations.
"""
abstract type FluxScheme end


"""
    spatial_operator!

Interface for the spatial discretization of the compressible Euler
equations.

The concrete method is selected by the type of `scheme`.
"""
function spatial_operator! end


struct Direct <: FluxScheme end

struct Pirozzoli <: FluxScheme end

struct Feiereisen <: FluxScheme end