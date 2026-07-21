"""
The core result of the project.

Sweep over:
  schemes:     [ConservativeFE(1), ConservativeFE(3)]   (start small; extend if time allows)
  integrators: [ExplicitRK3(), ImplicitMidpoint(1e-10)]
  (optionally NonConservative(L) as a negative control)

For each combination, build an Experiment, run_experiment(), collect:
  - L2 reconstruction error at t_reversal
  - K(t) drift curves (forward and reversed branches)
  - E(k) spectra at t_forward and t_reversal

Output: comparison plots (L2 error vs. L, per integrator) — this is the
key figure for the write-up / portfolio page / possible paper.
"""

# TODO: using ERT3D
# TODO: build the Experiment list, loop, collect results, save via JLD2, plot
