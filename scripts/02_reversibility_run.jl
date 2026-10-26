using ERT3D

# ============================================================
# Simulation parameters
# ============================================================

N = 32
dt_target =  0.06
t_rev = 8

dt = compatible_dt(t_rev, dt_target)
nsteps = round(Int, t_rev / dt)

grid = Grid(N)

params = Parameters(
    gamma = 1.4,
    Mt0 = 0.07,
)

derivative = Central6(grid)
scheme = KennedyGruber()
integrator = ImplicitMidpoint()

state = State(grid)

initialize!(
    state,
    TaylorGreen(),
    grid,
    params,
)
#=
coll = VTKCollection("data/raw/initial_condition")

add_snapshot!(
    coll,
    sim.state,
    sim.grid,
    sim.params,
    sim.derivative,
    sim.t,
)

close_collection!(coll)
=#

workspace = ImplicitMidpointWorkspace(grid)

sim = Simulation(
    state = state,
    grid = grid,
    params = params,
    derivative = derivative,
    scheme = scheme,
    integrator = integrator,
    workspace = workspace,
)

# Preserve the initial state for the final reconstruction test
s0 = State(grid)
copy_state!(s0, state)

# ============================================================
# Output configuration
# ============================================================

diagnostics = Diagnostics()

# Optional VTK output for visualization:
#
# coll = VTKCollection("data/processed/reversibility")
#
# hooks = [
#     OutputHook(
#         "diagnostics",
#         1,
#         (sim, step) -> record!(diagnostics, sim),
#     ),
#     OutputHook(
#         "vtk",
#         2,
#         (sim, step) -> add_snapshot!(
#             coll,
#             sim.state,
#             sim.grid,
#             sim.params,
#             sim.derivative,
#             sim.t,
#         ),
#     ),
# ]

hooks = [
    OutputHook(
        "diagnostics",
        1,
        (sim, step) -> record!(diagnostics, sim),
    ),
]

# ============================================================
# Simulation information
# ============================================================

println()
println("============================================================")
println("ERT3D — TIME REVERSIBILITY TEST")
println("============================================================")
println("Initial condition : Taylor–Green vortex")
println("Grid              : $(N)^3")
println("Mach number       : $(params.Mt0)")
println("Spatial scheme    : $(nameof(typeof(scheme)))")
println("Derivative        : $(nameof(typeof(derivative)))")
println("Time integrator   : $(nameof(typeof(integrator)))")
println("dt                : $(dt)")
println("Reversal time     : $(t_rev)")
println("Steps per leg     : $(nsteps)")
println("Total steps       : $(2 * nsteps)")
println("============================================================")
println()

# ============================================================
# Forward leg
# ============================================================

println("Forward evolution: 0 → $t_rev")

run!(
    sim,
    t_rev;
    dt = dt,
    hooks = hooks,
    verbose = true,
)

# Optional checkpoint at the reversal point:
#
# save_checkpoint(
#     sim,
#     "data/raw/reversibility/ckpt_reversal.jld2",
# )

# ============================================================
# Velocity reversal + reversed leg
# ============================================================

println()
println("Velocity reversal at t = $t_rev")
println("Reversed evolution: $t_rev → $(2 * t_rev)")

reverse_velocity!(sim.state)

run!(
    sim,
    2 * t_rev;
    dt = dt,
    hooks = hooks,
    verbose = true,
)

# Optional VTK output:
#
# close_collection!(coll)

# ============================================================
# Reconstruction error
# ============================================================

# Restore the original velocity direction before comparing
# the final state with the initial state.
reverse_velocity!(sim.state)

error = l2_reconstruction_error(s0, sim.state)

save_diagnostics(
    diagnostics,
    "data/raw/IM/diagnostics.jld2",
)

println()
println("============================================================")
println("REVERSIBILITY RESULT")
println("============================================================")
println("Final time        : $(sim.t)")
println("L2 reconstruction : $(error)")
println("Diagnostic records: $(length(diagnostics.t))")
println("============================================================")
println()