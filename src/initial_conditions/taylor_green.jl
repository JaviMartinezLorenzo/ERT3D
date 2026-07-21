"""
    taylor_green_ic(grid; ρ0=1.0, p0=..., U0=1.0)

Closed-form 3D Taylor-Green Vortex velocity field:
    u =  U0 sin(x) cos(y) cos(z)
    v = -U0 cos(x) sin(y) cos(z)
    w = 0
initialized with uniform density and a pressure field satisfying the
low-Mach compressible closure (verify exact standard TGV pressure formula
before implementing — several conventions exist in the literature).

Primary IC for the main project (reversibility comparison). Keep initial
Mach number low (see project scoping notes — avoid shocklet formation,
which would introduce genuine thermodynamic irreversibility unrelated to
the numerical comparison being made).

Validate resulting energy-decay/enstrophy curves against Brachet et al.
(1983) before trusting this for anything downstream.
"""
function taylor_green_ic end
# TODO: implement, returning a State

# TODO: helper to rescale amplitude U0 to hit a target M_t0 = u_rms / c
