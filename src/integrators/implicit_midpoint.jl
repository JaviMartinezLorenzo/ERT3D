"""
    ImplicitMidpoint(tol)

Implicit midpoint rule:
    state_{n+1} = state_n + dt * F((state_n + state_{n+1})/2)

Exactly self-adjoint (time-symmetric) up to floating-point round-off —
this is the structural property the whole reversibility comparison hinges
on. `tol` sets the nonlinear-solve tolerance (passed to NLsolve).

Implementation notes:
  - Use NLsolve.jl for the per-step nonlinear solve; do NOT hand-roll
    Newton-Krylov/GMRES (see project scoping discussion — this was
    explicitly de-risked by using the library instead).
  - Verify empirically in week 1 that NLsolve's default method converges
    reasonably for this system size before committing further design
    around it; if not, revisit tolerance/method choice early.
"""
struct ImplicitMidpoint <: TimeIntegrator
    tol::Float64
end

# TODO: step!(state, dt, scheme, integrator::ImplicitMidpoint, grid)
#       build residual function R(state_{n+1}) = state_{n+1} - state_n - dt*F(midpoint)
#       solve via NLsolve.nlsolve
