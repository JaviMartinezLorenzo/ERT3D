"""
    ExplicitRK4()

Standard 4-stage explicit Runge-Kutta. Second asymmetric baseline for
the reversibility comparison (in addition to ExplicitRK3) — track its
O(dt^4) reversibility error and high-wavenumber energy damping separately
from RK3's, since they're expected to diverge differently even though
neither is self-adjoint.
"""
struct ExplicitRK4 <: TimeIntegrator end

# TODO: step!(state, dt, scheme, integrator::ExplicitRK4, grid)
