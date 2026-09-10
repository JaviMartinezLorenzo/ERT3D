````markdown
# ERT3D — Euler Reversibility Testbench, 3D

[![Julia](https://img.shields.io/badge/Julia-1.12-%239558B2.svg)](https://julialang.org/)

A research code for studying **discrete time-reversibility in the compressible Euler equations**.

ERT3D investigates how spatial discretization and temporal integration influence the preservation of reversibility in numerical simulations. The primary benchmark is the **three-dimensional inviscid Taylor–Green vortex**.

The code is designed around two independent numerical components:

- **Spatial discretization:** conservative and split-form representations of the convective terms, combined with configurable derivative operators.
- **Time integration:** explicit and symmetric/implicit time-integration schemes.

The main objective is to quantify how these choices affect conservation, accuracy, stability, and ultimately the error accumulated when a simulation is integrated forward and then backward in time.

> **Status:** Early development. The core spatial discretization, state representation, derivative operators, and initial time-integration components are currently being developed and validated.

---

## Overview

For the compressible Euler equations,

\[
\frac{\partial Q}{\partial t} + \nabla \cdot F(Q) = 0,
\]

a continuous solution is formally time-reversible under

\[
t \rightarrow -t,
\qquad
\mathbf{u} \rightarrow -\mathbf{u}.
\]

A discrete numerical method does not necessarily preserve this property.

ERT3D provides a controlled framework for studying the resulting **reversibility error** and its dependence on:

- spatial derivative order,
- flux formulation,
- split-form representation,
- temporal integration scheme,
- spatial resolution,
- timestep size.

The Taylor–Green vortex provides the primary test problem because it generates increasingly complex three-dimensional flow structures while remaining well suited to controlled numerical experiments.

---

## Numerical methods

### Spatial discretization

The current spatial discretization framework supports:

- Central finite-difference derivative operators
- Fourth-, sixth-, and eighth-order central differences
- Direct conservative formulation
- Feiereisen split form
- Kennedy–Gruber split form

The formulation is designed so that additional spatial schemes can be introduced without modifying the existing solver infrastructure.

### Time integration

The current time-integration framework includes:

- Shu–Osher third-order explicit Runge–Kutta
- Classical fourth-order explicit Runge–Kutta
- Implicit midpoint
- Two-stage fourth-order Gauss–Legendre Runge–Kutta

Explicit methods are being implemented first, followed by symmetric implicit methods relevant to reversibility studies.

---

## Design

The code separates the main numerical concerns into independent abstractions:

```text
                    ERT3D
                      │
          ┌───────────┴───────────┐
          │                       │
   Spatial discretization     Time integration
          │                       │
    DerivativeOperator       TimeIntegrator
          │                       │
      FluxScheme                step!
          │                       │
          └───────────┬───────────┘
                      │
                    State
                      │
              Euler semi-discrete
                  equations
````

Adding a new derivative operator, flux formulation, or time integrator should require implementing the corresponding interface rather than modifying the rest of the solver.

The detailed design rationale is documented in [`docs/design.md`](docs/design.md).

---

## Project structure

```text
src/
├── ERT3D.jl                     # Module entry point
├── parameters.jl                # Physical and numerical parameters
├── grid.jl                      # Computational grid
├── state.jl                     # Conserved and primitive states
├── physics.jl                   # Euler-variable conversions and physics
├── workspace.jl                 # Reusable solver workspaces
├── derivatives/
│   ├── abstract.jl              # Derivative operator interface
│   └── central.jl               # Central finite differences
├── schemes/
│   ├── abstract.jl              # Flux scheme interface
│   ├── Direct.jl                # Direct conservative formulation
│   ├── Feiereisen.jl            # Feiereisen split form
│   └── KennedyGruber.jl         # Kennedy–Gruber split form
├── integrators/
│   ├── abstract.jl              # Time-integrator interface
│   ├── explicit_rk3.jl          # Shu–Osher RK3
│   ├── explicit_rk4.jl          # Classical RK4
│   └── implicit_midpoint.jl     # Implicit midpoint
├── initial_conditions/
│   ├── taylor_green.jl          # Taylor–Green vortex
│   └── synthetic_turbulence.jl  # Synthetic turbulence test case
├── experiment.jl                # Experiment orchestration
├── metrics.jl                   # Error and diagnostic metrics
└── io.jl                        # Output and data handling

test/
├── derivatives/
├── schemes/
├── integrators/
└── ...

docs/
└── design.md                    # Design and implementation rationale

scripts/                         # Reproducible simulations and parameter sweeps
```

---

## Development

The project is currently under active development. Numerical components are being validated independently before being combined into complete simulation workflows.

Current validation includes:

* manufactured derivative convergence tests,
* discrete skew-symmetry of centered derivatives,
* constant-state preservation,
* global conservation,
* spatial convergence of split-form operators,
* reusable workspace and state-operation tests,
* Runge–Kutta stage and timestep validation.

The next stage is to expose these components through a high-level simulation interface and run complete Taylor–Green vortex simulations.

---

## Planned workflow

The intended user-facing interface is:

```julia
using ERT3D

sim = Simulation(
    grid,
    params;
    derivative = Central8(grid),
    scheme = KennedyGruber(),
    integrator = ExplicitRK3(),
)

initialize!(sim, TaylorGreen())

run!(sim, Tfinal, dt)
```

This interface is **planned and will evolve during development**.

---

## Reversibility benchmark

The central experiment will follow the sequence:

```text
Initial condition
       │
       ▼
 Forward integration
       │
       ▼
   State at T
       │
       ▼
 Backward integration
       │
       ▼
   State at 0
       │
       ▼
 Reversibility error
```

The accumulated error is measured relative to the original initial condition.

This makes it possible to compare spatial and temporal discretizations on the same physical problem while isolating the contribution of each numerical choice.

---

## References

Pirozzoli, S.
“Generalized conservative approximations of split convective derivative operators.”
*Journal of Computational Physics*, 229 (2010), 7180–7190.

Brachet, M. E., et al.
“Small-scale structure of the Taylor–Green vortex.”
*Journal of Fluid Mechanics*, 130 (1983), 411–452.

Duponcheel, M., Orlandi, P., Winckelmans, G.
“Time-reversibility of the Euler equations as a benchmark for energy conserving schemes.”
*Journal of Computational Physics*, 227 (2008), 8736–8752.

```
```
