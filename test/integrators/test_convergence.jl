using Test

println()
println("Explicit RK3")
println("============")

N = 16

grid = Grid(N)
params = Parameters(1.4, 0.5)
D = Central4(grid)

scheme = KennedyGruber()
integrator = ExplicitRK3()

state = State(grid)
workspace = ERT3D.RK3Workspace(grid)

# -------------------------------------------------------------------------
# Test 1: constant state must remain constant
# -------------------------------------------------------------------------

state.rho  .= 1.0
state.rhou .= 0.5
state.rhov .= 0.25
state.rhow .= -0.125
state.rhoE .= 2.5

state_initial = state

dt = 1.0e-3

step!(
    state,
    integrator,
    dt,
    scheme,
    D,
    grid,
    params,
    workspace,
)

err = maximum([
    maximum(abs.(state.rho  .- state_initial.rho)),
    maximum(abs.(state.rhou .- state_initial.rhou)),
    maximum(abs.(state.rhov .- state_initial.rhov)),
    maximum(abs.(state.rhow .- state_initial.rhow)),
    maximum(abs.(state.rhoE .- state_initial.rhoE)),
])

println("1. Constant-state preservation")
println("   max error = ", err)

@test err < 1.0e-14


# -------------------------------------------------------------------------
# Test 2: one-step consistency
#
# For a smooth state, compare one RK3 step against
# Q + dt * RHS(Q). The difference should be O(dt^2).
# -------------------------------------------------------------------------

state = State(grid)

x = grid.x
y = grid.y
z = grid.z

for k in 1:N, j in 1:N, i in 1:N
    xi = x[i]
    yj = y[j]
    zk = z[k]

    rho = 1.0 + 0.05 * sin(xi) * cos(yj)
    u   = 0.2 + 0.03 * cos(zk)
    v   = 0.1 + 0.02 * sin(xi)
    w   = 0.05 + 0.02 * cos(yj)
    p   = 1.0 + 0.05 * cos(xi + yj + zk)

    state.rho[i, j, k]  = rho
    state.rhou[i, j, k] = rho * u
    state.rhov[i, j, k] = rho * v
    state.rhow[i, j, k] = rho * w
    state.rhoE[i, j, k] =
        p / (params.gamma - 1.0) +
        0.5 * rho * (u^2 + v^2 + w^2)
end

state_initial = state
rhs = State(grid)

spatial = workspace.spatial

spatial_operator!(
    rhs,
    scheme,
    state_initial,
    D,
    grid,
    params,
    spatial,
)

dt = 1.0e-6

expected = state_initial
axpy!(expected, dt, rhs)

step!(
    state,
    integrator,
    dt,
    scheme,
    D,
    grid,
    params,
    workspace,
)

err = maximum([
    maximum(abs.(state.rho  .- expected.rho)),
    maximum(abs.(state.rhou .- expected.rhou)),
    maximum(abs.(state.rhov .- expected.rhov)),
    maximum(abs.(state.rhow .- expected.rhow)),
    maximum(abs.(state.rhoE .- expected.rhoE)),
])

println("2. One-step consistency")
println("   max error = ", err)

@test err < 1.0e-10

println()
println("All RK3 tests passed.")