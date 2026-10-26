using ERT3D
using JLD2

# Reversal-time study: fixed timestep, varying reversal time.
# Run from the repository root with: julia --project=. scripts/05_reversibility_trev_study.jl

const N = 32
const DT_FIXED = 0.025
# Every T_REV value must be an integer multiple of DT_FIXED because run!
# requires exact step alignment. These values include the T_rev=8 baseline.
const T_REVS = [1.0, 2.0, 4.0, 8.0, 16.0, 20.0, 32.0, 40.0, 60.0, 70.0, 80.0]
const INTEGRATORS = ["ImplicitMidpoint"]
const OUTPUT_DIR = "data/raw/reversibility_studies"
const OUTPUT_CSV = joinpath(OUTPUT_DIR, "trev_study.csv")
const OUTPUT_JLD2 = joinpath(OUTPUT_DIR, "trev_study.jld2")

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
                                t_rev::Float64, dt::Float64)
    nsteps_float = t_rev / dt
    nsteps = round(Int, nsteps_float)
    isapprox(nsteps * dt, t_rev; atol=1e-10, rtol=0.0) ||
        error("T_REV=$t_rev must be an integer multiple of DT_FIXED=$dt")

    sim, state0 = make_simulation(integrator_name, N)
    elapsed = @elapsed begin
        run!(sim, t_rev; dt=dt, verbose=false)
        reverse_velocity!(sim.state)
        run!(sim, 2t_rev; dt=dt, verbose=false)
        reverse_velocity!(sim.state)
    end

    error = l2_reconstruction_error(state0, sim.state)
    return (integrator=integrator_name, N=N, t_rev=t_rev,
            dt=dt, nsteps_per_leg=nsteps, total_steps=2nsteps,
            reconstruction_error=error, runtime_s=elapsed,
            final_time=sim.t)
end

function csv_row(r)
    return join((r.integrator, r.N, r.t_rev, r.dt,
                 r.nsteps_per_leg, r.total_steps, r.reconstruction_error,
                 r.runtime_s, r.final_time), ",")
end

mkpath(OUTPUT_DIR)
open(OUTPUT_CSV, "w") do io
    println(io, "integrator,N,t_rev,dt_actual,nsteps_per_leg,total_steps,reconstruction_error,runtime_s,final_time")
end

results = NamedTuple[]
for t_rev in T_REVS
    for integrator_name in INTEGRATORS
        println("Running T_rev study: integrator=$integrator_name, N=$N, T_rev=$t_rev, dt=$DT_FIXED")
        result = run_reversibility_case(integrator_name, N, t_rev, DT_FIXED)
        push!(results, result)
        open(OUTPUT_CSV, "a") do io
            println(io, csv_row(result))
        end
        println("  steps/leg=$(result.nsteps_per_leg), error=$(result.reconstruction_error), runtime=$(result.runtime_s) s")
    end
end

jldsave(OUTPUT_JLD2;
    integrator=[r.integrator for r in results],
    N=[r.N for r in results],
    t_rev=[r.t_rev for r in results],
    dt_actual=[r.dt for r in results],
    nsteps_per_leg=[r.nsteps_per_leg for r in results],
    total_steps=[r.total_steps for r in results],
    reconstruction_error=[r.reconstruction_error for r in results],
    runtime_s=[r.runtime_s for r in results],
    final_time=[r.final_time for r in results],
)

println("\nSaved results to:\n  $OUTPUT_CSV\n  $OUTPUT_JLD2")
