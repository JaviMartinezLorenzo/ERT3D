"""
    Experiment

A single configuration of the testbench: one flux scheme + one time
integrator + resolution + timing parameters. Designed so a full sweep
(scheme × integrator × L) is a loop over `Experiment` values rather than
copy-pasted scripts.

Fields:
    scheme::FluxScheme
    integrator::TimeIntegrator
    N::Int              # grid resolution (points per direction)
    dt::Float64
    t_forward::Float64  # e.g. 10.0 — run forward to here
    t_reversal::Float64 # e.g. 20.0 — after u→-u, run forward again to here
    label::String        # human-readable tag for filenames/plots
"""
Base.@kwdef struct Experiment
    scheme::FluxScheme
    integrator::TimeIntegrator
    N::Int
    dt::Float64
    t_forward::Float64
    t_reversal::Float64
    label::String
end

"""
    run_experiment(exp::Experiment) -> ExperimentResult

Orchestration:
  1. build Grid(exp.N)
  2. initialize State via taylor_green_ic (main project) — swap IC function
     here for correctness-gate runs
  3. step! forward from t=0 to exp.t_forward, recording K(t) each step
  4. apply the reversal operator u → -u (see state.jl: reverse_velocity!)
  5. step! forward from "0" to (exp.t_reversal - exp.t_forward), recording K(t)
  6. compute metrics.jl diagnostics (L2 reconstruction error, spectra, etc.)
  7. return a result struct + optionally checkpoint via io.jl

Keep this function IO/plotting-free — return data, let scripts/ handle
figures, so this stays reusable for both the correctness gate and the
main sweep.
"""
function run_experiment end
# TODO: implement per steps above
