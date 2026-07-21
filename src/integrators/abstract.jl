"""
    TimeIntegrator

Abstract type for temporal advancement schemes. Every concrete subtype
must implement:

    step!(state::State, dt::Float64, scheme::FluxScheme, integrator::TimeIntegrator, grid::Grid)

which advances `state` in place by one step of size `dt`, calling
`compute_flux(scheme, state, grid)` internally — the integrator must not
depend on which FluxScheme is passed in; that's what keeps the two axes
decoupled.

Key structural distinction for this project:
  - ExplicitRK3 / ExplicitRK4  -> NOT self-adjoint; Φ_{-Δt} ∘ Φ_{Δt} ≠ id exactly
  - ImplicitMidpoint            -> exactly self-adjoint (up to round-off);
                                   this symmetry is the actual mechanism
                                   being tested against reversibility
"""
abstract type TimeIntegrator end

# TODO: function step!(state, dt, scheme, integrator::TimeIntegrator, grid) end
#       (interface documented here; each integrator file provides its own method)
