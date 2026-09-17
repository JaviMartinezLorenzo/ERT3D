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
    grid::Grid,
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

# ================================================================
# TODO
# ================================================================

# function energy_spectrum(state, grid) end