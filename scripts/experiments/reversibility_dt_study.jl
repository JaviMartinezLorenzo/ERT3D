using ERT3D
using JLD2

# Time-step convergence study: fixed reversal time, varying dt.
# Run from the repository root with: julia --project=. scripts/04_reversibility_dt_study.jl

const N = 32
const T_REV = 8.0
const DT_TARGETS = [0.06, 0.05, 0.025, 0.01, 0.005, 0.001]
const INTEGRATORS = ["ImplicitMidpoint"]
const OUTPUT_DIR = "data/raw/reversibility_studies_dt"
const OUTPUT_CSV = joinpath(OUTPUT_DIR, "dt_study.csv")
const OUTPUT_JLD2 = joinpath(OUTPUT_DIR, "dt_study.jld2")

function make_simulation(integrator_name::String, N::Int)
    grid = Grid(N)
    params = Parameters(gamma=1.4, Mt0=0.07)
    derivative = Central6(grid)
    scheme = KennedyGruber()
    state = State(grid)
    initialize!(state, TaylorGreen(), grid, params)

    if integrator_name == "RK3"
        integrator = ExplicitRK3()
        workspace = RK3Workspace(grid)
    elseif integrator_name == "ImplicitMidpoint"
        integrator = ImplicitMidpoint()
        workspace = ImplicitMidpointWorkspace(grid)
    else
        error("Unknown integrator: $integrator_name")
    end

    sim = Simulation(
        state=state,
        grid=grid,
        params=params,
        derivative=derivative,
        scheme=scheme,
        integrator=integrator,
        workspace=workspace,
    )
    state0 = State(grid)
    copy_state!(state0, state)
    return sim, state0
end

function run_reversibility_case(integrator_name::String, N::Int,
                                t_rev::Float64, dt_target::Float64)
    # Ensure both legs land exactly on their requested end times.
    dt = compatible_dt(t_rev, dt_target)
    nsteps = round(Int, t_rev / dt)
    sim, state0 = make_simulation(integrator_name, N)

    elapsed = @elapsed begin
        run!(sim, t_rev; dt=dt, verbose=false)
        reverse_velocity!(sim.state)
        run!(sim, 2t_rev; dt=dt, verbose=false)
        # Restore the original velocity orientation before comparing states.
        reverse_velocity!(sim.state)
    end

    error = l2_reconstruction_error(state0, sim.state)
    return (integrator=integrator_name, N=N, t_rev=t_rev,
            dt_target=dt_target, dt=dt, nsteps_per_leg=nsteps,
            total_steps=2nsteps, reconstruction_error=error,
            runtime_s=elapsed, final_time=sim.t)
end

function csv_row(r)
    return join((r.integrator, r.N, r.t_rev, r.dt_target, r.dt,
                 r.nsteps_per_leg, r.total_steps, r.reconstruction_error,
                 r.runtime_s, r.final_time), ",")
end

mkpath(OUTPUT_DIR)
open(OUTPUT_CSV, "w") do io
    println(io, "integrator,N,t_rev,dt_target,dt_actual,nsteps_per_leg,total_steps,reconstruction_error,runtime_s,final_time")
end

results = NamedTuple[]
for dt_target in DT_TARGETS
    for integrator_name in INTEGRATORS
        println("Running dt study: integrator=$integrator_name, N=$N, T_rev=$T_REV, dt_target=$dt_target")
        result = run_reversibility_case(integrator_name, N, T_REV, dt_target)
        push!(results, result)
        open(OUTPUT_CSV, "a") do io
            println(io, csv_row(result))
        end
        println("  dt=$(result.dt), steps/leg=$(result.nsteps_per_leg), error=$(result.reconstruction_error), runtime=$(result.runtime_s) s")
    end
end

# Save structured native Julia data as well as a CSV for plotting.
jldsave(OUTPUT_JLD2;
    integrator=[r.integrator for r in results],
    N=[r.N for r in results],
    t_rev=[r.t_rev for r in results],
    dt_target=[r.dt_target for r in results],
    dt_actual=[r.dt for r in results],
    nsteps_per_leg=[r.nsteps_per_leg for r in results],
    total_steps=[r.total_steps for r in results],
    reconstruction_error=[r.reconstruction_error for r in results],
    runtime_s=[r.runtime_s for r in results],
    final_time=[r.final_time for r in results],
)

println("\nSaved results to:\n  $OUTPUT_CSV\n  $OUTPUT_JLD2")
