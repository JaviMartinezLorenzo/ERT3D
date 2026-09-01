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

# ---- Spatial derivative operators --------------------------------------
include("derivatives/abstract.jl")
include("derivatives/central.jl")
# include("derivatives/spectral.jl")  # add when implemented

# ---- Axis 1: spatial flux formulations --------------------------------
include("schemes/abstract.jl")         # FluxScheme abstract type + compute_flux interface
include("schemes/direct.jl") 
include("schemes/Pirozzoli.jl")
include("schemes/Feiereisen.jl") 

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

# ---- Core types ----------------------------------------------------

export Parameters
export Grid
export State
export PrimitiveState

# ---- Spatial derivatives -------------------------------------------

export DerivativeOperator
export CentralDifference
export Central4
export Central6
export Central8
export Spectral

export derivative_x!
export derivative_y!
export derivative_z!

# ---- Spatial formulations ------------------------------------------

export FluxScheme
export Direct
export Pirozzoli
export Feiereisen

export spatial_operator!

# ---- Time integration ----------------------------------------------

export TimeIntegrator
export ExplicitRK3
export ExplicitRK4
export ImplicitMidpoint

# ---- Initial conditions --------------------------------------------

export taylor_green_ic
export synthetic_turbulence_ic

# ---- Diagnostics / experiments ------------------------------------

export primitive_variables
export conserved_variables
export rms_velocity
export run_experiment

end # module ERT3D
