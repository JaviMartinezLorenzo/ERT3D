using ERT3D

N =32
T = Float32

grid = Grid(N; T=T)

params = Parameters(
    gamma = T(1.4),
    Mt0 = T(0.1),
    T = T,
)

derivative = Central4(grid)
scheme = KennedyGruber()
integrator = ExplicitRK3()

state = State(grid)

initialize!(
    state,
    TaylorGreen(),
    grid,
    params,
)

workspace = RK3Workspace(grid)

sim = Simulation(
    state = state,
    grid = grid,
    params = params,
    derivative = derivative,
    scheme = scheme,
    integrator = integrator,
    workspace = workspace,
)

diagnostics = Diagnostics()

hooks = [
    OutputHook(
        "diagnostics",
        1,
        (sim, step) -> record!(diagnostics, sim),
    ),
]

dt = 0.125/2
t_end = dt * 55.0

run!(
    sim,
    t_end;
    dt = dt,
    hooks = hooks,
    verbose = true,
)
#=
save_diagnostics(
    diagnostics,
    "tgv_diagnostics.jld2",
)
=#