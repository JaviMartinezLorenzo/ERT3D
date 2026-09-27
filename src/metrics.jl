using Statistics


# ================================================================
# Root-mean-square velocity
# ================================================================

"""
    rms_velocity(state)

Compute the root-mean-square fluctuating velocity,

    u_rms = sqrt(<u'² + v'² + w'²>),

where `u' = u - <u>` and `<·>` denotes the spatial average.
"""
function rms_velocity(state::State)

    u = state.rhou ./ state.rho
    v = state.rhov ./ state.rho
    w = state.rhow ./ state.rho

    u .-= mean(u)
    v .-= mean(v)
    w .-= mean(w)

    return sqrt(
        (
            sum(abs2, u) +
            sum(abs2, v) +
            sum(abs2, w)
        ) / length(u)
    )
end


# ================================================================
# Velocity magnitude
# ================================================================

"""
    velocity_magnitude(state)

Compute the velocity magnitude,

    |u| = sqrt(u² + v² + w²),

at every grid point.
"""
function velocity_magnitude(state::State)

    u = state.rhou ./ state.rho
    v = state.rhov ./ state.rho
    w = state.rhow ./ state.rho

    return sqrt.(u.^2 .+ v.^2 .+ w.^2)
end


# ================================================================
# Kinetic energy
# ================================================================

"""
    kinetic_energy(state, grid)

Compute the total kinetic energy,

    K = ∫ 1/2 ρ |u|² dV,

over the periodic computational domain.
"""
function kinetic_energy(
    state::State,
    grid::Grid,
)

    u = state.rhou ./ state.rho
    v = state.rhov ./ state.rho
    w = state.rhow ./ state.rho

    kinetic = 0.5 .* state.rho .* (
        u.^2 .+ v.^2 .+ w.^2
    )

    N = length(grid.x)
    dV = (2π / N)^3

    return sum(kinetic) * dV
end


# ================================================================
# Density RMS
# ================================================================

"""
    density_rms(state)

Compute the root-mean-square density fluctuation,

    ρ'_rms = sqrt(<(ρ - <ρ>)²>).
"""
function density_rms(state::State)

    rho = state.rho
    rho_fluctuating = rho .- mean(rho)

    return sqrt(
        sum(abs2, rho_fluctuating) / length(rho)
    )
end


# ================================================================
# L2 reconstruction error
# ================================================================

"""
    l2_reconstruction_error(state0, state_reversed, grid)

Compute the normalized L2 reconstruction error between the initial
and time-reversed conserved states.
"""
function l2_reconstruction_error(
    state0::State,
    state_reversed::State,
)

    error² =
        sum(abs2, state_reversed.rho  .- state0.rho)  +
        sum(abs2, state_reversed.rhou .- state0.rhou) +
        sum(abs2, state_reversed.rhov .- state0.rhov) +
        sum(abs2, state_reversed.rhow .- state0.rhow) +
        sum(abs2, state_reversed.rhoE .- state0.rhoE)

    norm² =
        sum(abs2, state0.rho)  +
        sum(abs2, state0.rhou) +
        sum(abs2, state0.rhov) +
        sum(abs2, state0.rhow) +
        sum(abs2, state0.rhoE)

    return sqrt(error² / norm²)
end

function total_energy(
    state::State,
    grid::Grid,
)
    N = length(grid.x)
    dV = (2π / N)^3

    return sum(state.rhoE) * dV
end

"""
    current_cfl(sim::Simulation, dt::Float64) -> Float64

Compute the maximum local CFL number for the current state and timestep `dt`.
"""
function current_cfl(sim::Simulation, dt::Float64)

    state = sim.state
    gamma = sim.params.gamma
    inv_dx = 1.0 / sim.grid.dx

    max_cfl = 0.0

    @inbounds for k in axes(state.rho, 3)
        for j in axes(state.rho, 2)
            for i in axes(state.rho, 1)

                rho  = state.rho[i, j, k]
                rhou = state.rhou[i, j, k]
                rhov = state.rhov[i, j, k]
                rhow = state.rhow[i, j, k]
                rhoE = state.rhoE[i, j, k]

                u = rhou / rho
                v = rhov / rho
                w = rhow / rho

                e = rhoE / rho -
                    0.5 * (u^2 + v^2 + w^2)

                c = sqrt(gamma * (gamma - 1.0) * e)

                local_cfl = dt * inv_dx * (
                    abs(u) + abs(v) + abs(w) + 3.0 * c
                )

                max_cfl = max(max_cfl, local_cfl)
            end
        end
    end

    return max_cfl
end

"""
    velocity_gradient_tensor(state, params, D, grid)

Compute the full velocity gradient tensor ∂u_i/∂x_j at every grid point.

Returned as a 3×3 tuple of tuples: grad[i][j] is the array holding
∂u_i/∂x_j everywhere. Built from the existing derivative operator —
used by q_criterion and by the Brachet et al. early-time gradient check.
"""
function velocity_gradient_tensor(state::State, params::Parameters,
                                    D::DerivativeOperator, grid::Grid)
    primitive = primitive_variables(state, params)
    u, v, w = primitive.u, primitive.v, primitive.w

    dudx = similar(u); dudy = similar(u); dudz = similar(u)
    dvdx = similar(u); dvdy = similar(u); dvdz = similar(u)
    dwdx = similar(u); dwdy = similar(u); dwdz = similar(u)

    derivative_x!(dudx, u, D, grid)
    derivative_y!(dudy, u, D, grid)
    derivative_z!(dudz, u, D, grid)

    derivative_x!(dvdx, v, D, grid)
    derivative_y!(dvdy, v, D, grid)
    derivative_z!(dvdz, v, D, grid)

    derivative_x!(dwdx, w, D, grid)
    derivative_y!(dwdy, w, D, grid)
    derivative_z!(dwdz, w, D, grid)

    return ((dudx, dudy, dudz),
            (dvdx, dvdy, dvdz),
            (dwdx, dwdy, dwdz))
end


"""
    q_criterion(state, params, derivative, grid)

Compute the Q-criterion field, Q = 0.5*(||Ω||² - ||S||²), from the
velocity gradient tensor. Positive Q marks vortex cores — the standard
diagnostic for isosurface visualization of vortical structures.

Reuses velocity_gradient_tensor (same machinery as the Brachet
early-time gradient check).
"""
function q_criterion(state::State, params::Parameters, D::DerivativeOperator, grid::Grid)
    grad = velocity_gradient_tensor(state, params, D, grid)
    # grad[i][j] = ∂u_i/∂x_j

    Q = similar(grad[1][1])
    @inbounds for idx in eachindex(Q)
        s2 = 0.0  # ||S||²
        w2 = 0.0  # ||Ω||²
        for i in 1:3, j in 1:3
            gij = grad[i][j][idx]
            gji = grad[j][i][idx]
            s = 0.5 * (gij + gji)
            w = 0.5 * (gij - gji)
            s2 += s^2
            w2 += w^2
        end
        Q[idx] = 0.5 * (w2 - s2)
    end
    return Q
end

# ================================================================
# TODO
# ================================================================

# function energy_spectrum(state, grid) end