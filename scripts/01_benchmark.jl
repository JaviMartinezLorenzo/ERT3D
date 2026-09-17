using ERT3D
using BenchmarkTools
using Profile

N = 32

grid = Grid(N)

params = Parameters(
    gamma = 1.4,
    Mt0 = 0.1,
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

dt = 1e-3

@btime begin
    step!(
        $state,
        $integrator,
        $dt,
        $scheme,
        $derivative,
        $grid,
        $params,
        $workspace,
    )
    nothing

end

# ------------------------------------------------------------
# Profiling: where is the time actually going inside step!
# ------------------------------------------------------------

# Warm up once outside the profiler, so JIT compilation time
# (which is real but one-time, not representative of steady-state
# per-step cost) doesn't pollute the profile.
step!(state, integrator, dt, scheme, derivative, grid, params, workspace)

Profile.clear()
@profile for _ in 1:200
    step!(state, integrator, dt, scheme, derivative, grid, params, workspace)
end

println()
println("Profile (sorted by time, top entries):")
println("========================================")
Profile.print(mincount=20, format=:flat, sortedby=:count)