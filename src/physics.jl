"""
    primitive_variables(state::State, gamma)

Convert a conserved state to primitive variables using the calorically
perfect ideal-gas equation of state

    p = (γ - 1)ρe

Returns a `PrimitiveState`.
"""
function primitive_variables(state::State, params::Parameters)

    gamma = params.gamma
    rho = state.rho

    u = state.rhou ./ rho
    v = state.rhov ./ rho
    w = state.rhow ./ rho

    kinetic = 0.5 .* (u.^2 .+ v.^2 .+ w.^2)

    e = state.rhoE ./ rho .- kinetic

    p = (gamma - 1) .* rho .* e

    return PrimitiveState(rho, u, v, w, p)
end

"""
    conserved_variables(primitive::PrimitiveState, gamma)

Convert a primitive state to conserved variables.

Returns a `State`.
"""
function conserved_variables(primitive::PrimitiveState, params::Parameters)

    gamma = params.gamma
    rho = primitive.rho
    u   = primitive.u
    v   = primitive.v
    w   = primitive.w
    p   = primitive.p

    rhou = rho .* u
    rhov = rho .* v
    rhow = rho .* w

    kinetic = 0.5 .* (u.^2 .+ v.^2 .+ w.^2)

    rhoE = p ./ (gamma - 1) .+ rho .* kinetic

    return State(rho, rhou, rhov, rhow, rhoE)
end

"""
    reverse_velocity!(state)

Apply the velocity-reversal operator

    u → -u

used by the Euler Reversibility Test.

Density and total energy remain unchanged.
"""
function reverse_velocity!(state::State)
    state.rhou .*= -1
    state.rhov .*= -1
    state.rhow .*= -1
    return state
end