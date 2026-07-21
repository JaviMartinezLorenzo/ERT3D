"""
    ExplicitRK3()

Third-order TVD Runge-Kutta (Shu-Osher, 1988) — the integrator used in
Pirozzoli (2010) for the Euler-turbulence test case. Included primarily
to reproduce his setup exactly for the correctness gate, and to serve as
one of the two asymmetric/non-reversible baselines in the main comparison
(alongside ExplicitRK4).

NOT self-adjoint — expect approximate-but-not-exact reversibility, bounded
by O(dt^3) local truncation error per step.
"""
struct ExplicitRK3 <: TimeIntegrator end

# TODO: step!(state, dt, scheme, integrator::ExplicitRK3, grid)
#       standard 3-stage Shu-Osher SSP-RK3 update
