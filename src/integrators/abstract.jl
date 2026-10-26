"""
    TimeIntegrator

Abstract supertype for temporal integration operators.
Concrete subtypes define the numerical method used to advance the
semi-discrete system in time.
"""
abstract type TimeIntegrator end

"""
    ExplicitRK3(; formulation=:shu_osher, compensated=false)

Third-order explicit Runge-Kutta time integrator.

`formulation` may be `:shu_osher` (the original SSP Shu-Osher
implementation) or `:increment` (accumulate the RK increment and add it
to the solution once at the end of the step).

`compensated=true` enables Kahan compensated accumulation of the final
increment. It is available only with `formulation=:increment`.

The default constructor preserves the original Shu-Osher formulation.
"""
struct ExplicitRK3 <: TimeIntegrator
    formulation::Symbol
    compensated::Bool

    function ExplicitRK3(;
        formulation::Symbol = :shu_osher,
        compensated::Bool = false,
    )
        formulation in (:shu_osher, :increment) ||
            throw(ArgumentError("formulation must be :shu_osher or :increment"))
        compensated && formulation != :increment &&
            throw(ArgumentError("compensated=true requires formulation=:increment"))
        return new(formulation, compensated)
    end
end

"""
    ExplicitRK4

Fourth-order classical explicit Runge-Kutta time integration scheme.
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
"""
struct GaussLegendre4 <: TimeIntegrator end
