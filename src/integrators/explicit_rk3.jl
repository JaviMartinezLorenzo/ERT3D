"""
    step!(state, ::ExplicitRK3, dt, scheme, D, grid, params, workspace)

Advance the compressible Euler solution by one time step `dt` using
the third-order TVD Runge-Kutta scheme of Shu and Osher.

The semi-discrete system is

    dQ/dt = R(Q),

where `R(Q)` is evaluated by `spatial_operator!`.

The three stages are

    Q¹ = Qⁿ + dt R(Qⁿ)

    Q² = 3/4 Qⁿ
         + 1/4 [Q¹ + dt R(Q¹)]

    Qⁿ⁺¹ = 1/3 Qⁿ
           + 2/3 [Q² + dt R(Q²)]

All intermediate states and right-hand-side storage are supplied by
the pre-allocated `RK3Workspace`.
"""
function step!(
    state::State,
    ::ExplicitRK3,
    dt::Float64,
    scheme::FluxScheme,
    D::DerivativeOperator,
    grid::Grid,
    params::Parameters,
    workspace::RK3Workspace,
)

    stage1 = workspace.stage1
    stage2 = workspace.stage2
    rhs    = workspace.rhs
    spatial = workspace.spatial

    # ------------------------------------------------------------
    # Stage 1
    #
    # Q1 = Qn + dt * R(Qn)
    # ------------------------------------------------------------

    spatial_operator!(
        rhs,
        scheme,
        state,
        D,
        grid,
        params,
        spatial,
    )

    copy_state!(stage1, state)
    axpy!(stage1, dt, rhs)


    # ------------------------------------------------------------
    # Stage 2
    #
    # Q2 = 3/4 Qn + 1/4 [Q1 + dt R(Q1)]
    # ------------------------------------------------------------

    spatial_operator!(
        rhs,
        scheme,
        stage1,
        D,
        grid,
        params,
        spatial,
    )

    axpy!(stage1, dt, rhs)

    linear_combination!(
        stage2,
        0.75,
        state,
        0.25,
        stage1,
    )


    # ------------------------------------------------------------
    # Stage 3
    #
    # Qn+1 = 1/3 Qn + 2/3 [Q2 + dt R(Q2)]
    # ------------------------------------------------------------

    spatial_operator!(
        rhs,
        scheme,
        stage2,
        D,
        grid,
        params,
        spatial,
    )

    axpy!(stage2, dt, rhs)

    linear_combination!(
        state,
        1.0 / 3.0,
        state,
        2.0 / 3.0,
        stage2,
    )

    return state
end