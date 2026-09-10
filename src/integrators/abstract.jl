"""
    TimeIntegrator

Abstract supertype for temporal integration operators.
Concrete subtypes define the numerical method used to advance the
semi-discrete system in time.

"""
abstract type TimeIntegrator end

"""
    ExplicitRK3

Third-order Shu-Osher explicit Runge-Kutta time integration scheme.

Uses three explicit stages to advance the semi-discrete solution by
one time step.

"""
struct ExplicitRK3 <: TimeIntegrator end

"""
    ExplicitRK4

Fourth-order classical explicit Runge-Kutta time integration scheme.

Uses four explicit stages to advance the semi-discrete solution by
one time step.

"""
struct ExplicitRK4 <: TimeIntegrator end

"""
    ImplicitMidpoint

Second-order implicit midpoint time integration scheme.

A one-stage implicit Runge-Kutta method with a self-adjoint time
discretization.

"""
struct ImplicitMidpoint <: TimeIntegrator end

"""
    GaussLegendre4

Fourth-order two-stage Gauss-Legendre implicit Runge-Kutta time
integration scheme.

A self-adjoint implicit method based on the two-stage Gauss-Legendre
collocation scheme.

"""
struct GaussLegendre4 <: TimeIntegrator end