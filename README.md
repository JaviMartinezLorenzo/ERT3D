# ERT3D — Euler Reversibility Testbench, 3D

[TODO: hero figure here — e.g. Q-criterion before/after reversal, or the
key L2-error-vs-order comparison plot, once available]

A benchmark for studying how spatial flux formulation (conservative vs.
non-conservative split forms, at variable order) and temporal integration
(symmetric/implicit vs. asymmetric/explicit) each affect discrete
time-reversibility of the compressible Euler equations, using the
inviscid Taylor-Green Vortex as the primary test case.

Built to be extensible: adding a new flux scheme or time integrator means
adding one file (see `docs/design.md`), not modifying existing code.

## Status

Early development — see `docs/design.md` for scope and design rationale.

## Quick start

```julia
# TODO once implemented:
# using ERT3D
# exp = Experiment(scheme=ConservativeFE(3), integrator=ImplicitMidpoint(1e-10),
#                  N=64, dt=..., t_forward=10.0, t_reversal=20.0, label="cfe3_impmid")
# result = run_experiment(exp)
```

## Repository structure

```
src/
  ERT3D.jl                    — module entry point, wiring the two axes together
  state.jl, grid.jl           — core data containers
  schemes/                    — FluxScheme implementations (Axis 1)
  integrators/                — TimeIntegrator implementations (Axis 2)
  initial_conditions/         — TGV (main) and synthetic turbulence (correctness gate)
  experiment.jl, metrics.jl, io.jl
scripts/                      — orchestration scripts (correctness gates, main sweep)
test/                         — automated tests, incl. structural checks
docs/                         — design rationale
```

## References

- Pirozzoli, S. "Generalized conservative approximations of split
  convective derivative operators." *J. Comput. Phys.* 229 (2010): 7180-7190.
- Brachet, M.E. et al. "Small-scale structure of the Taylor-Green vortex."
  *J. Fluid Mech.* 130 (1983): 411-452.
- [TODO: add the 2008 reversibility-benchmark reference once confirmed]

## License

[TODO: MIT recommended]
