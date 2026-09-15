"""
    initialize!(state, ::TaylorGreen, grid, params)

Initialize `state` with the compressible Taylor–Green vortex.

The velocity amplitude is `V0 = 2*Mt0`, giving
`sqrt(<u² + v² + w²>) = Mt0` for the present TGV convention.

The reference pressure is `p0 = 1/gamma`, consistent with the
project nondimensionalization.
"""
function initialize!(
    state::State,
    ::TaylorGreen,
    grid::Grid,
    params::Parameters,
)
    Mt0 = params.Mt0
    gamma = params.gamma
    V0 = 2.0 * Mt0

    X = reshape(grid.x, :, 1, 1)
    Y = reshape(grid.y, 1, :, 1)
    Z = reshape(grid.z, 1, 1, :)

    rho = state.rho
    rhou = state.rhou
    rhov = state.rhov
    rhow = state.rhow
    rhoE = state.rhoE

    @. rho = 1.0

    @. rhou = rho * V0 * sin(X) * cos(Y) * cos(Z)
    @. rhov = -rho * V0 * cos(X) * sin(Y) * cos(Z)
    @. rhow = 0.0

    @. state.rhoE =
        (
            1.0 / gamma +
            (Mt0^2 / 4.0) *
            (cos(2.0 * X) + cos(2.0 * Y)) *
            (cos(2.0 * Z) + 2.0)
        ) / (gamma - 1.0) +
        0.5 * (
            rhou^2 +
            rhov^2 +
            rhow^2
        ) / rho

    return state
end