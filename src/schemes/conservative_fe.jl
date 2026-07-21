"""
    ConservativeFE(L)

Locally conservative flux-difference form of the Feiereisen split (FE-SF),
Pirozzoli eq. (13):
    f̂_{j+1/2} = 2 * Σ_{l=1}^{L} a_l * Σ_{m=0}^{l-1} (ρg·u)_{j-m,l}

using the two-point averaging operator (eq. 10):
    (f̃,g)_{j,l} = (1/4)(f_j + f_{j+l})(g_j + g_{j+l})

This is the primary scheme for the main project — implement this first
(before ConservativeBL / ConservativeKG), get correctness verified via:
  1. synthetic-turbulence correctness gate (Pirozzoli §3.1, M_t0=0.07)
  2. TGV validation against Brachet et al. (1983)
before building the reversibility comparison on top of it.

φ mapping per Pirozzoli's generic-scalar notation:
  - continuity: φ=1, flux uses (ρ, u)
  - momentum:   φ=u_i
  - energy:     φ=H  (total enthalpy, NOT E — see eq. 13 vs energy flux term)
"""
struct ConservativeFE <: FluxScheme
    L::Int
end

# TODO: compute_flux(scheme::ConservativeFE, state, grid)
#       implement eq. (10), (13) — precompute/cache the two-point averages
#       per Pirozzoli's efficiency note (Appendix A) before assembling flux differences
