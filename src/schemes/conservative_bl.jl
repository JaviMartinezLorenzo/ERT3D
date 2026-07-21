"""
    ConservativeBL(L)

Locally conservative flux-difference form of the Blaisdell split (BL-SF),
Pirozzoli eq. (14). Same averaging-operator machinery as ConservativeFE,
different (f,g) pairing in the flux formula.

Implement after ConservativeFE is verified — mostly a variant of the same
machinery, good second data point for the "does split-form choice affect
reversibility" comparison (secondary to the main conservative-vs-
non-conservative and RK3-vs-implicit-midpoint axes).
"""
struct ConservativeBL <: FluxScheme
    L::Int
end

# TODO: compute_flux(scheme::ConservativeBL, state, grid)
