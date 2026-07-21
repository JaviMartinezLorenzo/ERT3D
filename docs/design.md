# ERT3D — design notes

## Motivation

Time-reversibility of the compressible Euler equations is a structural
property of the *time integrator*, not the spatial discretization —
established for FD2/FD4/pseudo-spectral schemes on TGV. This project
extends that question to Pirozzoli's family of locally-conservative
split-flux discretizations (which nobody has paired with a genuinely
symmetric time integrator), asking specifically:

- Does the "lower order reverses better" trend still hold within a
  properly conservative flux family?
- Does an exactly self-adjoint integrator (implicit midpoint) change the
  picture relative to the asymmetric integrators (RK3, RK4) used in the
  existing literature?

Scoped deliberately small: 2-3 values of stencil order L, two integrators,
TGV at low initial Mach number (avoid shocklets — a real thermodynamic
irreversibility that would contaminate the numerical comparison).
Explicitly NOT attempting a full mechanistic explanation or an exhaustive
order sweep — see project scoping notes (this is a portfolio/verification
project, with a possible small-journal submission as a stretch outcome,
not a guaranteed publication).

## Architecture

Two orthogonal axes, decoupled via multiple dispatch:
  - `FluxScheme`     — spatial discretization (src/schemes/)
  - `TimeIntegrator` — temporal advancement (src/integrators/)

Neither depends on the other's implementation; `step!` dispatches on both,
calling `compute_flux` internally. Adding a new scheme or integrator later
means adding one file, not touching existing code.

## Correctness gates (in order)

1. **Synthetic turbulence, M_t0=0.07** vs. Pirozzoli Fig. 2 — verifies the
   ConservativeFE flux implementation itself, independent of TGV/reversibility.
2. **TGV forward-only** vs. Brachet et al. (1983) — verifies the flow/IC
   before reversibility is layered on top.
3. **ImplicitMidpoint self-adjointness** on a trivial linear problem —
   verify the integrator's structural symmetry holds to ~machine precision
   before trusting it on the full nonlinear system.

## Known scope cuts (see project planning discussion)

- No hand-rolled JFNK/GMRES — NLsolve.jl handles the implicit nonlinear solve.
- No KG-SF (density-weighted) scheme initially — only relevant at higher
  Mach/strong density variation, which is out of scope (low-Mach only,
  to avoid shocklets).
- No full order sweep (L up to 8-10) — 2-3 values of L only.
- No mechanistic explanation of the order-vs-reversibility trend — report
  the observed trend honestly as an open question if it doesn't get fully
  explained within the timeframe.
