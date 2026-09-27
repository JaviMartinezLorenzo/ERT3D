using ERT3D

# ============================================================
# Pirozzoli / Honein-Moin isotropic turbulence test
# ============================================================

N = 32

Mt0 = 0.07
k0  = 6.0

# Eddy turnover time from Honein & Moin
tau = 2sqrt(3) / (k0 * Mt0)

# CFL ~ 1.
dt = 5.5e-2

# Reproduce approximately 60 turnover times
t_final = 60.0 * tau
nsteps = round(Int, t_final/dt)
t_final = nsteps * dt



println("============================================================")
println("Pirozzoli isotropic turbulence test")
println("============================================================")
println("N       = $N")
println("Mt0     = $Mt0")
println("k0      = $k0")
println("tau     = $tau")
println("t_final = $t_final")
println("dt      = $dt")
println("t_final/tau = $(t_final/tau)")
println()

# ------------------------------------------------------------
# Grid and parameters
# ------------------------------------------------------------

grid = Grid(N)

params = Parameters(
    gamma = 1.4,
    Mt0   = Mt0,
    k0    = k0,
)

# Pirozzoli's C-KG-SF:

derivative = Central6(grid)
scheme = Feiereisen()
integrator = ExplicitRK3()

# ------------------------------------------------------------
# Initial condition
# ------------------------------------------------------------

state = State(grid)

initialize!(
    state,
    SyntheticTurbulence(),
    grid,
    params,
)

# ------------------------------------------------------------
# Simulation
# ------------------------------------------------------------

workspace = RK3Workspace(grid)

sim = Simulation(
    state      = state,
    grid       = grid,
    params     = params,
    derivative = derivative,
    scheme     = scheme,
    integrator = integrator,
    workspace  = workspace,
)

# ------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------

diagnostics = Diagnostics()


hooks = [

    # Record every timestep
    OutputHook(
        "diagnostics",
        1,
        (sim, step) -> record!(diagnostics, sim),
    ),

]

# ------------------------------------------------------------
# Run
# ------------------------------------------------------------

println("Forward run: t = 0 -> $t_final")
println("Target: approximately 60 eddy turnover times")
println()

run!(
    sim,
    t_final;
    dt = dt,
    hooks = hooks,
    verbose = true,
)

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------
mkpath("data/raw/C6_F")

save_diagnostics(
    diagnostics,
    "data/raw/C6_F/diagnostics.jld2",
)

println()
println("============================================================")
println("Run complete")
println("============================================================")
println("Final t       = $(sim.t)")
println("Final t/tau   = $(sim.t / tau)")
println("Records       = $(length(diagnostics.t))")