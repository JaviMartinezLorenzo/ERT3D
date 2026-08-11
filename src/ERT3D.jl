"""
    ERT3D

Euler Reversibility Testbench, 3D.

Benchmark for studying how spatial flux formulation (conservative vs.
non-conservative split forms, at variable order) and temporal integration
(symmetric/implicit vs. asymmetric/explicit) each affect discrete
time-reversibility of the compressible Euler equations, using the
Taylor-Green Vortex as the primary test case.

Two orthogonal axes, kept fully decoupled via multiple dispatch:
  - FluxScheme      (src/schemes/*.jl)   — spatial discretization of convective terms
  - TimeIntegrator  (src/integrators/*.jl) — temporal advancement

See docs/design.md for the architecture rationale.
"""
module ERT3D

using FFTW
using NLsolve
using StaticArrays
using JLD2
using WriteVTK
using ProgressMeter

# ---- Core abstractions ------------------------------------------------
include("parameters.jl")
include("grid.jl")             # Grid: 2π-periodic Cartesian grid, Δx, wavenumbers
include("state.jl")            # State: the conserved-variable container (ρ, ρu, ρE)
include("physics.jl")          # Primitive variables calculations


# ---- Axis 1: spatial flux formulations --------------------------------
include("schemes/abstract.jl")         # FluxScheme abstract type + compute_flux interface
include("schemes/nonconservative.jl")  # D-FE-SF / D-BL-SF, eq. (5)
include("schemes/conservative_fe.jl")  # C-FE-SF, eq. (13)
include("schemes/conservative_bl.jl")  # C-BL-SF, eq. (14)
# include("schemes/conservative_kg.jl") # C-KG-SF, eq. (16) — add later if needed

# ---- Axis 2: time integrators ------------------------------------------
include("integrators/abstract.jl")       # TimeIntegrator abstract type + step! interface
include("integrators/explicit_rk3.jl")   # Shu-Osher TVD RK3 (matches Pirozzoli's paper)
include("integrators/explicit_rk4.jl")   # standard RK4 (asymmetric baseline)
include("integrators/implicit_midpoint.jl") # symmetric/time-reversible integrator

# ---- Initial conditions -------------------------------------------------
include("initial_conditions/taylor_green.jl")     # closed-form TGV velocity field
include("initial_conditions/synthetic_turbulence.jl") # Pirozzoli §3.1 correctness-gate IC

# ---- Experiment orchestration --------------------------------------------
include("experiment.jl")   # Experiment struct + run_experiment(): forward → reverse → metrics
include("metrics.jl")      # L2 reconstruction error, K(t) drift, E(k) spectra
include("io.jl")           # JLD2 checkpointing, VTK export for ParaView/Q-criterion

export Experiment, run_experiment
export FluxScheme, NonConservative, ConservativeFE, ConservativeBL
export TimeIntegrator, ExplicitRK3, ExplicitRK4, ImplicitMidpoint
export taylor_green_ic, synthetic_turbulence_ic

end # module ERT3D
