"""
    FluxScheme

Abstract type for spatial discretization of the convective (flux) terms.
Every concrete subtype must implement:

    compute_flux(scheme::FluxScheme, state::State, grid::Grid) -> dstate/dt contribution

Concrete subtypes should each carry their own stencil order `L` as a field,
so order-of-accuracy sweeps (L=1,2,3,...) don't require new types.

Adding a new formulation later (e.g. C-KG-SF, eq. 16) means:
  1. define a new `struct ConservativeKG <: FluxScheme` with its fields
  2. implement `compute_flux` for it
  3. nothing else in the codebase needs to change
"""
abstract type FluxScheme end

# TODO: function compute_flux(scheme::FluxScheme, state, grid) end
#       (defined generically here only as a documented interface;
#        each scheme file provides its own method)
