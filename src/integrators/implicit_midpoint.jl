"""
    step!(state, ::ImplicitMidpoint, dt, scheme, D, grid, params, workspace)

Advance the compressible Euler solution by one time step `dt` using
the implicit midpoint method.

The semi-discrete system is

    dQ/dt = R(Q),

where `R(Q)` is evaluated by `spatial_operator!`.

The implicit update is

    Qⁿ⁺¹ = Qⁿ + dt R((Qⁿ + Qⁿ⁺¹) / 2)

The nonlinear system is solved using fixed-point iteration.
All intermediate states and right-hand-side storage are supplied
by the pre-allocated `ImplicitMidpointWorkspace`.
"""
function step!(
    state::State,
    ::ImplicitMidpoint,
    dt::Float64,
    scheme::FluxScheme,
    D::DerivativeOperator,
    grid::Grid,
    params::Parameters,
    workspace::ImplicitMidpointWorkspace,
)
    Qk = workspace.Qk
    Qnew = workspace.Qnew
    Qmid = workspace.Qmid
    rhs = workspace.R
    spatial = workspace.spatial

    copy_state!(Qk, state)

    T = eltype(state.rho)

    qscale = max(
        maximum(abs, state.rho),
        maximum(abs, state.rhou),
        maximum(abs, state.rhov),
        maximum(abs, state.rhow),
        maximum(abs, state.rhoE),
        one(T),
    )

    rtol = 50eps(T)
    atol = zero(T)
    tolerance = atol + rtol * qscale
    delta = T(Inf)

    for k in 1:workspace.maxiter
        # Qmid = (Qn + Qk) / 2
        linear_combination!(Qmid, T(0.5), state, T(0.5), Qk)

        # R = R(Qmid)
        spatial_operator!(
            rhs, scheme, Qmid, D, grid, params, spatial,
        )

        # Qnew = Qn + dt * R(Qmid)
        copy_state!(Qnew, state)
        axpy!(Qnew, dt, rhs)

        if !all(isfinite, Qnew.rho)  ||
           !all(isfinite, Qnew.rhou) ||
           !all(isfinite, Qnew.rhov) ||
           !all(isfinite, Qnew.rhow) ||
           !all(isfinite, Qnew.rhoE)
            error(
                "Implicit midpoint iteration produced non-finite " *
                "state values at iteration $k."
            )
        end

        # Allocation-free update norm:
        # delta = ||Qnew - Qk||∞
        delta = max(
            max_abs_diff(Qnew.rho,  Qk.rho),
            max_abs_diff(Qnew.rhou, Qk.rhou),
            max_abs_diff(Qnew.rhov, Qk.rhov),
            max_abs_diff(Qnew.rhow, Qk.rhow),
            max_abs_diff(Qnew.rhoE, Qk.rhoE),
        )

        # Swap buffers rather than copying the new iterate.
        Qk, Qnew = Qnew, Qk

        if delta <= tolerance
            copy_state!(state, Qk)
            return state
        end
    end

    error(
        "Implicit midpoint fixed-point iteration failed to converge " *
        "after $(workspace.maxiter) iterations (update norm = $delta, " *
        "tolerance = $tolerance).",
    )
end
