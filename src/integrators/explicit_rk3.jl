"""
    step!(state, ::ExplicitRK3, dt, scheme, D, grid, params, workspace)

Advance the compressible Euler solution by one time step using third-order
explicit Runge-Kutta.

The default `:shu_osher` formulation preserves the original SSP update.
The optional `:increment` formulation computes

    ΔQ = dt/6 * (R(Qⁿ) + R(Q¹) + 4R(Q²))
    Qⁿ⁺¹ = Qⁿ + ΔQ

and can optionally use Kahan compensation for the final state update.
Compensation only reduces rounding in that final addition; it does not
compensate the spatial operator or RK stage calculations.
"""
function step!(
    state::State,
    integrator::ExplicitRK3,
    dt::Float64,
    scheme::FluxScheme,
    D::DerivativeOperator,
    grid::Grid,
    params::Parameters,
    workspace::RK3Workspace,
)
    if integrator.formulation === :shu_osher
        return _step_rk3_shu_osher!(
            state, dt, scheme, D, grid, params, workspace,
        )
    else
        return _step_rk3_increment!(
            state, dt, integrator.compensated,
            scheme, D, grid, params, workspace,
        )
    end
end

function _step_rk3_shu_osher!(
    state::State,
    dt::Float64,
    scheme::FluxScheme,
    D::DerivativeOperator,
    grid::Grid,
    params::Parameters,
    workspace::RK3Workspace,
)
    stage1 = workspace.stage1
    stage2 = workspace.stage2
    rhs = workspace.rhs
    spatial = workspace.spatial
    T = eltype(state.rho)
    h = T(dt)

    spatial_operator!(rhs, scheme, state, D, grid, params, spatial)
    copy_state!(stage1, state)
    axpy!(stage1, h, rhs)

    spatial_operator!(rhs, scheme, stage1, D, grid, params, spatial)
    axpy!(stage1, h, rhs)
    linear_combination!(stage2, T(0.75), state, T(0.25), stage1)

    spatial_operator!(rhs, scheme, stage2, D, grid, params, spatial)
    axpy!(stage2, h, rhs)
    linear_combination!(state, T(1 / 3), state, T(2 / 3), stage2)

    return state
end

function _step_rk3_increment!(
    state::State,
    dt::Float64,
    compensated::Bool,
    scheme::FluxScheme,
    D::DerivativeOperator,
    grid::Grid,
    params::Parameters,
    workspace::RK3Workspace,
)
    stage1 = workspace.stage1
    stage2 = workspace.stage2
    rhs = workspace.rhs
    rhs1 = workspace.rhs1
    rhs2 = workspace.rhs2
    increment = workspace.increment
    spatial = workspace.spatial
    T = eltype(state.rho)
    h = T(dt)

    # R1 = R(Qn); Q1 = Qn + dt*R1
    spatial_operator!(rhs1, scheme, state, D, grid, params, spatial)
    copy_state!(stage1, state)
    axpy!(stage1, h, rhs1)

    # R2 = R(Q1); Q2 = Qn + dt/4*(R1 + R2)
    spatial_operator!(rhs2, scheme, stage1, D, grid, params, spatial)
    copy_state!(stage2, state)
    axpy!(stage2, h / T(4), rhs1)
    axpy!(stage2, h / T(4), rhs2)

    # R3 = R(Q2)
    spatial_operator!(rhs, scheme, stage2, D, grid, params, spatial)

    # ΔQ = dt/6*(R1 + R2 + 4R3)
    copy_state!(increment, rhs1)
    axpy!(increment, one(T), rhs2)
    axpy!(increment, T(4), rhs)
    # Scale the weighted residual into the final increment without needing
    # an additional state-scaling helper.
    _scale_state!(increment, h / T(6))

    if compensated
        _kahan_add_state!(state, increment, workspace.compensation)
    else
        axpy!(state, one(T), increment)
    end

    return state
end

"""
Apply a Kahan-compensated addition `state += increment` field by field.

The compensation state is persistent workspace storage. It should be
zeroed whenever the solution is externally reset or replaced (for example,
when initializing a new run or manually reversing the state).
"""
function _kahan_add_state!(
    state::State{T},
    increment::State{T},
    compensation::State{T},
) where {T<:AbstractFloat}
    for field in (:rho, :rhou, :rhov, :rhow, :rhoE)
        q = getfield(state, field)
        dq = getfield(increment, field)
        c = getfield(compensation, field)
        @inbounds for i in eachindex(q, dq, c)
            y = dq[i] - c[i]
            updated = q[i] + y
            c[i] = (updated - q[i]) - y
            q[i] = updated
        end
    end
    return state
end

function _scale_state!(state::State{T}, α::T) where {T<:AbstractFloat}
    for field in (:rho, :rhou, :rhov, :rhow, :rhoE)
        q = getfield(state, field)
        @inbounds for i in eachindex(q)
            q[i] *= α
        end
    end
    return state
end
