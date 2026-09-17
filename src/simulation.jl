"""

Top-level orchestration object for compressible Euler simulations.

Holds the solution state, mesh geometry, physical parameters, spatial discretizations,
time integrators, and preallocated scratch buffers. `Simulation` coordinates the
execution of time steps without implementing the numerical flux or integration
algorithms directly.

"""

Base.@kwdef mutable struct Simulation
    state::State
    grid::Grid
    params::Parameters

    derivative::DerivativeOperator
    scheme::FluxScheme
    integrator::TimeIntegrator

    workspace::TimeIntegratorWorkspace

    t::Float64 = 0.0
end

"""
    OutputHook

A scheduled action fired during run!, based on step count not simulation time. 
This is deliberate: letting a hook clip dt to land on a specific time would make the step-size sequence differ
between the forward and post-reversal legs of a reversibility run, introducing an error source that would 
contaminate the actual metrics the experiment are trying to isolate.
"""
struct OutputHook
    name::String
    every_n_steps::Int
    action::Function   # action(sim::Simulation, step::Int) -> nothing
end

function fire_due_hooks!(hooks::Vector{OutputHook}, sim::Simulation, step::Int)
    for h in hooks
        @assert h.every_n_steps > 0 "every_n_steps must be positive"
        
        if step % h.every_n_steps == 0
            h.action(sim, step)
        end
    end
end



function run!(
    sim::Simulation,
    t_end::Float64;
    dt::Float64,
    hooks::Vector{OutputHook} = OutputHook[],
    verbose::Bool = false,
    cfl_limit::Float64 = 1.5,
)
    @assert t_end >= sim.t
    @assert dt > 0.0

    n_steps = round(Int, (t_end - sim.t) / dt)

    @assert isapprox(
        sim.t + n_steps * dt,
        t_end;
        atol = 1e-10,
    ) """
        t_end - sim.t must be (very nearly) an exact multiple of dt —
        choose dt so this divides evenly.
    """

    fire_due_hooks!(hooks, sim, 0)

    for step in 1:n_steps

        cfl = current_cfl(sim, dt)

        if cfl > cfl_limit
            error("""
                CFL limit exceeded at step $step (t=$(sim.t)):
                CFL=$(round(cfl, digits=3)) > $cfl_limit.
                Reduce dt or check for a developing instability.
                """)
        end

        step!(
            sim.state,
            sim.integrator,
            dt,
            sim.scheme,
            sim.derivative,
            sim.grid,
            sim.params,
            sim.workspace,
        )

        sim.t += dt

        fire_due_hooks!(hooks, sim, step)

        verbose && step % 50 == 0 && println(
            "step $step  t=$(round(sim.t, digits=3))  " *
            "CFL=$(round(cfl, digits=3))"
        )
    end

    return sim
end