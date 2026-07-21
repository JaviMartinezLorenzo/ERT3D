"""
    synthetic_turbulence_ic(grid; k0=6, Mt0=0.07)

Synthetic isotropic solenoidal turbulence field matching Pirozzoli §3.1
(originally Honein & Moin), used ONLY as a correctness gate — NOT part of
the main reversibility comparison (see project scoping discussion: this
IC is harder to build than TGV and not what the reversibility literature
uses, so keep it scoped to "does my flux implementation reproduce the
known energy-conservation result", nothing more).

Spectrum shape:
    E(k) = A (k/k0)^4 exp(-2 (k/k0)^2),  peaking at k0

Steps:
  1. generate random field in spectral space shaped to E(k)
  2. Helmholtz-project to enforce ∇·u = 0
  3. inverse FFT to physical space
  4. rescale amplitude so u_rms/c = Mt0
  5. zero initial density/temperature fluctuations (per Pirozzoli)

Target check: run with ConservativeFE, compare resulting K(t) and ρ' rms
curves qualitatively against Pirozzoli Fig. 2 (M_t0=0.07 case only —
skip the M_t0=0.3/KG-SF case, out of scope per project decisions).
"""
function synthetic_turbulence_ic end
# TODO: implement using FFTW, returning a State
