"""
Week 1 correctness gate.

Reproduce Pirozzoli §3.1 (M_t0 = 0.07 case only): synthetic turbulence IC,
ConservativeFE scheme, L=1 and L=3, verify:
  - kinetic energy stays ~constant (within a few % over a couple eddy
    turnover times)
  - density rms levels off near ρ'/ρ0/Mt0² ≈ 0.35

This is a throwaway-quality script, not part of the reusable Experiment
sweep machinery — its only job is "is my flux operator correct".
"""

# TODO: using ERT3D
# TODO: build grid, synthetic_turbulence_ic, run for L=1 and L=3, plot K(t) and ρ'_rms
