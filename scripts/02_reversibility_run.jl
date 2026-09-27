using ERT3D

N = 64
dt = 2.5e-2
t_rev = 8.0

grid = Grid(N)

params = Parameters(
    gamma = 1.4,
    Mt0 = 0.07,
)

derivative = Central6(grid)
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

s0 = State(grid)
copy_state!(s0, state)

diagnostics = Diagnostics()
coll = VTKCollection("data/processed/reversibility")

hooks = [
    OutputHook(
        "diagnostics",
        1,
        (sim, step) -> record!(diagnostics, sim),
    ),
    OutputHook(
        "vtk",
        2,
        (sim, step) -> add_snapshot!(coll, sim.state, sim.grid, sim.params, sim.derivative, sim.t),
    ),
]


# ------------------------------------------------------------
# Forward leg: 0 -> t_rev
# ------------------------------------------------------------

println("Forward leg: t = 0 -> $t_rev")

run!(
    sim,
    t_rev;
    dt = dt,
    hooks = hooks,
    verbose = true,
)

save_checkpoint(sim, "data/raw/reversibility/ckpt_reversal.jld2")

# ------------------------------------------------------------
# Reversal: flip velocity, run forward again for the same duration
# ------------------------------------------------------------

reverse_velocity!(sim.state)

println("Reversed leg: t = $t_rev -> $(2*t_rev)")

run!(
    sim,
    2 * t_rev;
    dt = dt,
    hooks = hooks,
    verbose = true,
)

close_collection!(coll)

# ------------------------------------------------------------
# Reconstruction error and diagnostics
# ------------------------------------------------------------

error = l2_reconstruction_error(s0, sim.state)
println()
println("L2 reconstruction error at t = $(2*t_rev): $error")

save_diagnostics(
    diagnostics,
    "data/raw/reversibility/diagnostics.jld2",
)

println("Run complete. Final t = $(sim.t), $(length(diagnostics.t)) diagnostic records saved.")