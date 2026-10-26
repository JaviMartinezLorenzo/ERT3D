"""
    initialize!(state, ::TaylorGreen, grid, params)

Initialize `state` with the compressible Taylor–Green vortex.

The velocity amplitude is `V0 = 2*Mt0`, giving
`sqrt(<u² + v² + w²>) = Mt0` for the present TGV convention.

The reference pressure is `p0 = 1/gamma`, consistent with the
project nondimensionalization.
"""
function initialize!(
    state::State{T},
    ::TaylorGreen,
    grid::Grid{T},
    params::Parameters{T},
) where {T<:AbstractFloat}

    Mt0 = params.Mt0
    gamma = params.gamma
    V0 = T(2) * Mt0

    X = reshape(grid.x, :, 1, 1)
    Y = reshape(grid.y, 1, :, 1)
    Z = reshape(grid.z, 1, 1, :)

    rho  = state.rho
    rhou = state.rhou
    rhov = state.rhov
    rhow = state.rhow
    rhoE = state.rhoE

    @. rho = one(T)

    @. rhou = rho * V0 * sin(X) * cos(Y) * cos(Z)
    @. rhov = -rho * V0 * cos(X) * sin(Y) * cos(Z)
    @. rhow = zero(T)

    @. rhoE = (
        one(T) / gamma +
        (Mt0^2 / T(4)) *
        (cos(T(2) * X) + cos(T(2) * Y)) *
        (cos(T(2) * Z) + T(2))
    ) / (gamma - one(T)) +
    T(0.5) * (
        rhou^2 +
        rhov^2 +
        rhow^2
    ) / rho

    return state
end