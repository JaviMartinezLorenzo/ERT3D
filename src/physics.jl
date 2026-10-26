"""
    primitive_variables(state::State{T}, params::Parameters{T})

Convert a conserved state to primitive variables using the calorically
perfect ideal-gas equation of state:

    p = (γ - 1)ρe

Returns a `PrimitiveState{T}`.
"""
function primitive_variables(
    state::State{T},
    params::Parameters{T},
) where {T<:AbstractFloat}

    gamma = params.gamma
    rho = state.rho

    u = state.rhou ./ rho
    v = state.rhov ./ rho
    w = state.rhow ./ rho

    kinetic = T(0.5) .* (u.^2 .+ v.^2 .+ w.^2)
    e = state.rhoE ./ rho .- kinetic
    p = (gamma - one(T)) .* rho .* e

    return PrimitiveState(rho, u, v, w, p)
end


"""
    conserved_variables(primitive::PrimitiveState{T}, params::Parameters{T})

Convert a primitive state to conserved variables.

Returns a `State{T}`.
"""
function conserved_variables(
    primitive::PrimitiveState{T},
    params::Parameters{T},
) where {T<:AbstractFloat}

    gamma = params.gamma
    rho = primitive.rho
    u   = primitive.u
    v   = primitive.v
    w   = primitive.w
    p   = primitive.p

    rhou = rho .* u
    rhov = rho .* v
    rhow = rho .* w

    kinetic = T(0.5) .* (u.^2 .+ v.^2 .+ w.^2)
    rhoE = p ./ (gamma - one(T)) .+ rho .* kinetic

    return State(rho, rhou, rhov, rhow, rhoE)
end


"""
    reverse_velocity!(state::State)

Apply the velocity-reversal operator:

    u → -u

Density and total energy remain unchanged.
"""
function reverse_velocity!(state::State{T}) where {T<:AbstractFloat}
    state.rhou .*= -one(T)
    state.rhov .*= -one(T)
    state.rhow .*= -one(T)
    return state
end


function primitive_variables!(
    out::PrimitiveState{T},
    state::State{T},
    params::Parameters{T},
) where {T<:AbstractFloat}

    gamma = params.gamma

    rho  = state.rho
    rhou = state.rhou
    rhov = state.rhov
    rhow = state.rhow
    rhoE = state.rhoE

    @. out.rho = rho
    @. out.u   = rhou / rho
    @. out.v   = rhov / rho
    @. out.w   = rhow / rho

    @. out.p = (gamma - one(T)) * (
        rhoE -
        T(0.5) * (
            rhou * out.u +
            rhov * out.v +
            rhow * out.w
        )
    )

    return out
end