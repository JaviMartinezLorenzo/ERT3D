"""
    NonConservative(L)

Direct (non-locally-conservative) application of central differences to the
split convective form, Pirozzoli eq. (5):
    D_S(fg)_j = (1/2) D(fg)_j + (1/2) f_j D(g)_j + (1/2) g_j D(f)_j

Corresponds to D-FE-SF / D-BL-SF depending on which (f,g) pairing is chosen
upstream — decide whether to make FE/BL a type parameter or a field.

L = stencil half-width (L=1 -> 2nd order, L=2 -> 4th order, L=3 -> 6th order).

This scheme is expected to be less robust (Pirozzoli notes it can diverge
under strong density variation) — useful as a baseline/negative control,
not for the main reversibility comparison.
"""
struct NonConservative <: FluxScheme
    L::Int
end

# TODO: compute_flux(scheme::NonConservative, state, grid)
#       implement eq. (5)/(6) directly (not in flux-difference form)
