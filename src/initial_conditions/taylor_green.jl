"""
    taylor_green_ic(grid, params)

Generate the compressible Taylor-Green vortex initial condition as a State.

Velocity amplitude V0 = 2*Mt0 is chosen so that rms_velocity(state) returns
exactly Mt0 (see rms_velocity docstring for the convention:
sqrt(<u²+v²+w²>), matching Honein & Moin 2004's Mt0 = sqrt(3)*u_rms,component
convention — verified algebraically: <u²+v²+w²> = V0²/4 for this field).

Pressure closure matches the standard compressible TGV benchmark, with
reference pressure p0 = 1/gamma consistent with the project's nondimensionalization
(rho0 = T0 = 1 at reference state — see equations note, Section 4).
"""
function taylor_green_ic(
    grid::Grid,
    params::Parameters
)
    Mt0 = params.Mt0
    gamma = params.gamma
    V0 = 2.0 * Mt0
    

    # Reshape the 1D coordinate vectors so Julia broadcasts them to a full
    # (N × N × N) Cartesian grid without explicitly storing X, Y and Z.
    X = reshape(grid.x, :, 1, 1)
    Y = reshape(grid.y, 1, :, 1)
    Z = reshape(grid.z, 1, 1, :)

    u =  V0 .* sin.(X) .* cos.(Y) .* cos.(Z)
    v = -V0 .* cos.(X) .* sin.(Y) .* cos.(Z)
    w = zeros(size(u))

    rho = ones(size(u))

    p = 1.0 / gamma .+
        (Mt0^2 / 4.0) .* (cos.(2 .* X) .+ cos.(2 .* Y)) .* (cos.(2 .* Z) .+ 2.0)

    primitive = PrimitiveState(rho, u, v, w, p)

    return conserved_variables(primitive, params)

end