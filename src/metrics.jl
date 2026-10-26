
using Statistics

# ================================================================
# Root-mean-square velocity
# ================================================================

"""
    rms_velocity(state)

Compute the root-mean-square fluctuating velocity,

    u_rms = sqrt(<u'² + v'² + w'²>),

where `u' = u - <u>` and `<·>` denotes the spatial average.

Returns a Float64 scalar.
"""
function rms_velocity(state::State{T}) where {T<:AbstractFloat}
    rho = state.rho
    n = length(rho)

    # Compute velocity means and fluctuations in Float64.
    u_mean = sum(Float64(state.rhou[i] / rho[i]) for i in eachindex(rho)) / n
    v_mean = sum(Float64(state.rhov[i] / rho[i]) for i in eachindex(rho)) / n
    w_mean = sum(Float64(state.rhow[i] / rho[i]) for i in eachindex(rho)) / n

    sum_sq = 0.0
    @inbounds for i in eachindex(rho)
        u = Float64(state.rhou[i] / rho[i]) - u_mean
        v = Float64(state.rhov[i] / rho[i]) - v_mean
        w = Float64(state.rhow[i] / rho[i]) - w_mean

        sum_sq += u^2 + v^2 + w^2
    end

    return sqrt(sum_sq / n)
end


# ================================================================
# Velocity magnitude
# ================================================================

"""
    velocity_magnitude(state)

Compute the velocity magnitude,

    |u| = sqrt(u² + v² + w²),

at every grid point. Returns an array with the same element type
as the state.
"""
function velocity_magnitude(state::State{T}) where {T<:AbstractFloat}
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

Returns a Float64 scalar.
"""
function kinetic_energy(
    state::State{T},
    grid::Grid{T},
) where {T<:AbstractFloat}

    N = length(grid.x)
    dV = (2π / N)^3
    total = 0.0

    @inbounds for i in eachindex(state.rho)
        rho = Float64(state.rho[i])
        u = Float64(state.rhou[i]) / rho
        v = Float64(state.rhov[i]) / rho
        w = Float64(state.rhow[i]) / rho

        total += 0.5 * rho * (u^2 + v^2 + w^2)
    end

    return total * dV
end


# ================================================================
# Density RMS
# ================================================================

"""
    density_rms(state)

Compute the root-mean-square density fluctuation,

    ρ'_rms = sqrt(<(ρ - <ρ>)²>).

Returns a Float64 scalar.
"""
function density_rms(state::State{T}) where {T<:AbstractFloat}
    rho = state.rho
    n = length(rho)

    rho_mean = sum(Float64(x) for x in rho) / n

    sum_sq = 0.0
    @inbounds for i in eachindex(rho)
        fluctuation = Float64(rho[i]) - rho_mean
        sum_sq += fluctuation^2
    end

    return sqrt(sum_sq / n)
end


# ================================================================
# Normalized L2 reconstruction / solution error
# ================================================================

"""
    l2_reconstruction_error(state0, state_reversed)

Compute the normalized L2 error between two conserved states.

Supports comparisons between different floating-point precisions.
All differences, norms and accumulations are evaluated in Float64.

For a cross-precision comparison, ensure both states represent the
same physical time and use compatible grids.
"""
function l2_reconstruction_error(
    state0::State{T1},
    state_reversed::State{T2},
) where {T1<:AbstractFloat,T2<:AbstractFloat}

    fields0 = (
        state0.rho, state0.rhou, state0.rhov,
        state0.rhow, state0.rhoE,
    )
    fields1 = (
        state_reversed.rho, state_reversed.rhou, state_reversed.rhov,
        state_reversed.rhow, state_reversed.rhoE,
    )

    error² = 0.0
    norm² = 0.0

    for (a, b) in zip(fields0, fields1)
        length(a) == length(b) ||
            throw(DimensionMismatch("States must have matching dimensions"))

        @inbounds for i in eachindex(a, b)
            ai = Float64(a[i])
            diff = Float64(b[i]) - ai

            error² += diff^2
            norm² += ai^2
        end
    end

    norm² > 0.0 ||
        throw(ArgumentError("Cannot normalize the error: initial state has zero norm"))

    return sqrt(error² / norm²)
end


# ================================================================
# Total energy
# ================================================================

"""
    total_energy(state, grid)

Compute the total energy,

    E = ∫ ρE dV,

over the periodic computational domain.

Returns a Float64 scalar.
"""
function total_energy(
    state::State{T},
    grid::Grid{T},
) where {T<:AbstractFloat}

    N = length(grid.x)
    dV = (2π / N)^3

    total = sum(Float64(x) for x in state.rhoE)

    return total * dV
end


# ================================================================
# CFL number
# ================================================================

"""
    current_cfl(sim::Simulation, dt::Float64) -> Float64

Compute the maximum local CFL number for the current state and
timestep. The calculation and result use Float64.
"""
function current_cfl(sim::Simulation, dt::Float64)
    state = sim.state
    gamma = Float64(sim.params.gamma)
    inv_dx = 1.0 / Float64(sim.grid.dx)

    max_cfl = 0.0

    @inbounds for i in eachindex(state.rho)
        rho = Float64(state.rho[i])
        rhou = Float64(state.rhou[i])
        rhov = Float64(state.rhov[i])
        rhow = Float64(state.rhow[i])
        rhoE = Float64(state.rhoE[i])

        u = rhou / rho
        v = rhov / rho
        w = rhow / rho

        e = rhoE / rho - 0.5 * (u^2 + v^2 + w^2)
        c² = gamma * (gamma - 1.0) * e

        c² >= 0.0 ||
            throw(DomainError(c², "Negative sound-speed squared; check state positivity"))

        c = sqrt(c²)

        local_cfl = dt * inv_dx * (
            abs(u) + abs(v) + abs(w) + 3.0 * c
        )

        max_cfl = max(max_cfl, local_cfl)
    end

    return max_cfl
end


# ================================================================
# Velocity gradient tensor
# ================================================================

"""
    velocity_gradient_tensor(state, params, D, grid)

Compute the full velocity gradient tensor ∂uᵢ/∂xⱼ at every grid point.

Returns a 3×3 tuple of tuples, where `grad[i][j]` stores ∂uᵢ/∂xⱼ.
The derivative fields retain the state precision.
"""
function velocity_gradient_tensor(
    state::State{T},
    params::Parameters{T},
    D::DerivativeOperator,
    grid::Grid{T},
) where {T<:AbstractFloat}

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

    return (
        (dudx, dudy, dudz),
        (dvdx, dvdy, dvdz),
        (dwdx, dwdy, dwdz),
    )
end


# ================================================================
# Q-criterion
# ================================================================

"""
    q_criterion(state, params, derivative, grid)

Compute the Q-criterion field,

    Q = 0.5 * (||Ω||² - ||S||²),

from the velocity gradient tensor. Positive Q identifies regions
where rotation dominates strain.

The returned field retains the state precision.
"""
function q_criterion(
    state::State{T},
    params::Parameters{T},
    D::DerivativeOperator,
    grid::Grid{T},
) where {T<:AbstractFloat}

    grad = velocity_gradient_tensor(state, params, D, grid)
    Q = similar(grad[1][1])

    @inbounds for idx in eachindex(Q)
        s2 = 0.0
        w2 = 0.0

        for i in 1:3, j in 1:3
            gij = Float64(grad[i][j][idx])
            gji = Float64(grad[j][i][idx])

            s = 0.5 * (gij + gji)
            w = 0.5 * (gij - gji)

            s2 += s^2
            w2 += w^2
        end

        Q[idx] = T(0.5 * (w2 - s2))
    end

    return Q
end

# ================================================================
# TODO
# ================================================================

# function energy_spectrum(state, grid) end

# ================================================================
# Field-wise normalized state errors
# ================================================================

"""
    state_errors(reference, state)

Compute normalized errors for density, momentum and total energy
density, using the spatial fluctuations of the reference fields
as normalization scales.

Returns a named tuple containing:
    all      : legacy normalized L2 error
    rho      : normalized density fluctuation error
    momentum : normalized momentum fluctuation error
    rhoE     : normalized total energy density fluctuation error
    worst    : maximum finite field-group error

Each field group is normalized by its own reference fluctuation RMS.
Momentum combines rhou, rhov and rhow.

If a reference group has exactly zero fluctuations, its error is NaN.
The `worst` metric ignores such groups.

All diagnostic calculations are accumulated in Float64.
"""
function state_errors(
    reference::State{T1},
    state::State{T2},
) where {T1<:AbstractFloat,T2<:AbstractFloat}

    rho_error = _group_fluctuation_error(
        (reference.rho,),
        (state.rho,),
    )

    momentum_error = _group_fluctuation_error(
        (reference.rhou, reference.rhov, reference.rhow),
        (state.rhou, state.rhov, state.rhow),
    )

    rhoE_error = _group_fluctuation_error(
        (reference.rhoE,),
        (state.rhoE,),
    )

    valid_errors = filter(
        isfinite,
        (rho_error, momentum_error, rhoE_error),
    )

    worst_error = isempty(valid_errors) ? NaN : maximum(valid_errors)

    return (
        all      = l2_reconstruction_error(reference, state),
        rho      = rho_error,
        momentum = momentum_error,
        rhoE     = rhoE_error,
        worst    = worst_error,
    )
end


"""
    state_error(reference, state; measure=:worst)

Return one field-wise normalized state error.

Available measures:
    :all      - legacy combined conservative-state error
    :rho      - density fluctuation error
    :momentum - momentum fluctuation error
    :rhoE     - total energy density fluctuation error
    :worst    - maximum finite field-group error
"""
function state_error(
    reference::State,
    state::State;
    measure::Symbol = :worst,
)
    errors = state_errors(reference, state)

    hasproperty(errors, measure) ||
        throw(ArgumentError("Unknown error measure: $measure"))

    return getproperty(errors, measure)
end


# ================================================================
# Internal helpers
# ================================================================

"""
Compute the RMS error divided by the reference fluctuation RMS
for a group of fields.

All fields are treated as components of one group. Each component
is centered by its own spatial mean before calculating the
reference fluctuation RMS.
"""
function _group_fluctuation_error(reference_fields, state_fields)

    length(reference_fields) == length(state_fields) ||
        throw(DimensionMismatch("Reference and state field groups differ"))

    error_squared = 0.0
    fluctuation_squared = 0.0
    count = 0

    for (ref, val) in zip(reference_fields, state_fields)
        length(ref) == length(val) ||
            throw(DimensionMismatch("Reference and state fields must match"))

        n = length(ref)
        n > 0 || throw(ArgumentError("State fields cannot be empty"))

        # Calculate the reference mean in Float64.
        reference_mean = sum(Float64(x) for x in ref) / n

        @inbounds for i in eachindex(ref, val)
            r = Float64(ref[i])
            q = Float64(val[i])

            difference = q - r
            fluctuation = r - reference_mean

            error_squared += difference^2
            fluctuation_squared += fluctuation^2
        end

        count += n
    end

    error_rms = sqrt(error_squared / count)
    fluctuation_rms = sqrt(fluctuation_squared / count)

    # A spatially uniform reference has no fluctuation scale.
    fluctuation_rms == 0.0 && return NaN

    return error_rms / fluctuation_rms
end


"""
    state_error_scales(reference)

Report the RMS magnitude and fluctuation RMS of each conservative
field group in the reference state.

For momentum, the RMS values combine rhou, rhov and rhow.
All returned values are Float64.
"""
function state_error_scales(reference::State)

    return (
        rho = _group_rms_scales((reference.rho,)),
        momentum = _group_rms_scales((
            reference.rhou,
            reference.rhov,
            reference.rhow,
        )),
        rhoE = _group_rms_scales((reference.rhoE,)),
    )
end


function _group_rms_scales(fields)

    total_squared = 0.0
    fluctuation_squared = 0.0
    count = 0

    for field in fields
        n = length(field)
        n > 0 || throw(ArgumentError("State fields cannot be empty"))

        mean_value = sum(Float64(x) for x in field) / n

        @inbounds for x in field
            value = Float64(x)
            total_squared += value^2
            fluctuation_squared += (value - mean_value)^2
        end

        count += n
    end

    return (
        rms = sqrt(total_squared / count),
        fluctuation_rms = sqrt(fluctuation_squared / count),
    )
end

function max_abs_diff(a, b)
    @assert axes(a) == axes(b)

    delta = zero(eltype(a))

    @inbounds @simd for i in eachindex(a, b)
        delta = max(delta, abs(a[i] - b[i]))
    end

    return delta
end