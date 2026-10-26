using ERT3D
using Statistics
using Printf
using Random

# ================================================================
# ERT3D — Integrator performance benchmark
#
# Compares explicit RK3 (Shu-Osher, increment, increment + Kahan) with
# implicit midpoint. The integrators are timed through `step!` directly,
# so the cost of `run!` (CFL evaluation, hooks) is not attributed to
# either of them.
#
# What is reported per (precision, dt):
#   - time per step (median and minimum) and the relative spread
#   - allocations per step
#   - ns per cell and step
#   - residual-evaluation equivalents per step: time per step divided by
#     the time of one `spatial_operator!` call. RK3 should be close to 3;
#     for the midpoint it is the effective iteration count, which grows
#     with the CFL number.
#   - the difference of the final state from RK3 Shu-Osher (a check that
#     all variants solve the same problem)
#
# Assumed from the existing code base (adapt if they differ):
#   ERT3D.step!(state, integrator, dt, scheme, derivative, grid, params, workspace)
#   ERT3D.spatial_operator!(rhs, scheme, Q, derivative, grid, params, spatial_ws)
#   Simulation fields: state, grid, params, derivative, scheme, integrator, workspace
#   workspace.spatial holds the spatial-operator workspace
#   l2_reconstruction_error(state_a, state_b)
# ================================================================

# ----------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------

const N = 32
const T_END = 0.2                    # simulated time per run
const DT_VALUES = [1e-3, 5e-3, 2.5e-2]
const PRECISIONS = (Float64, Float32)
const NTRIALS = 20
const NWARMUP_STEPS = 3

const GAMMA = 1.4
const MT0 = 0.07

const CSV_PATH = "data/processed/benchmark_integrators.csv"

const VARIANTS = [
    (name = "RK3 Shu-Osher",         kind = :rk3,
     make = () -> ExplicitRK3()),
    (name = "RK3 increment",         kind = :rk3,
     make = () -> ExplicitRK3(formulation = :increment)),
    (name = "RK3 increment + Kahan", kind = :rk3,
     make = () -> ExplicitRK3(formulation = :increment, compensated = true)),
    (name = "Implicit midpoint",     kind = :midpoint,
     make = () -> ImplicitMidpoint()),
]

const BASELINE = VARIANTS[1].name

# ----------------------------------------------------------------
# Setup
# ----------------------------------------------------------------

function build_simulation(variant, T)
    grid = Grid(N; T = T)
    params = Parameters(gamma = T(GAMMA), Mt0 = T(MT0), T = T)

    state = State(grid)
    initialize!(state, TaylorGreen(), grid, params)

    workspace = variant.kind === :midpoint ?
        ImplicitMidpointWorkspace(grid) :
        RK3Workspace(grid)

    return Simulation(
        state = state,
        grid = grid,
        params = params,
        derivative = Central6(grid),
        scheme = KennedyGruber(),
        integrator = variant.make(),
        workspace = workspace,
    )
end

# Advance `nsteps` steps of size `dt` through step! only.
function advance!(sim, dt, nsteps)
    for _ in 1:nsteps
        ERT3D.step!(
            sim.state, sim.integrator, dt, sim.scheme,
            sim.derivative, sim.grid, sim.params, sim.workspace,
        )
    end
    return sim
end

function all_finite(state)
    return all(isfinite, state.rho) && all(isfinite, state.rhou) &&
           all(isfinite, state.rhov) && all(isfinite, state.rhow) &&
           all(isfinite, state.rhoE)
end

# ----------------------------------------------------------------
# Measurements
# ----------------------------------------------------------------

# One timed run. The simulation is built before the clock starts.
function measure_run(variant, T, dt, nsteps)
    sim = build_simulation(variant, T)
    GC.gc()

    stats = @timed advance!(sim, dt, nsteps)

    all_finite(sim.state) || error("non-finite state after $nsteps steps")

    return (
        time = stats.time,
        bytes = stats.bytes,
        gctime = stats.gctime,
        state = sim.state,
    )
end

# Median time of one spatial-operator call (seconds), or NaN if the
# assumed interface does not match.
function rhs_cost(T; nrep = 30)
    try
        sim = build_simulation(VARIANTS[1], T)
        spatial = sim.workspace.spatial
        rhs = deepcopy(sim.state)

        call() = ERT3D.spatial_operator!(
            rhs, sim.scheme, sim.state, sim.derivative,
            sim.grid, sim.params, spatial,
        )

        call(); call()                       # compile and warm caches
        GC.gc()
        return median([(@elapsed call()) for _ in 1:nrep])
    catch err
        @warn "Could not time the residual evaluation; RHS-equivalents will be NaN" exception = err
        return NaN
    end
end

# ----------------------------------------------------------------
# Output helpers
# ----------------------------------------------------------------

function print_header()
    println("=" ^ 92)
    println("ERT3D INTEGRATOR PERFORMANCE BENCHMARK")
    println("=" ^ 92)

    version = try
        string(pkgversion(ERT3D))
    catch
        "unknown"
    end

    total_steps = length(PRECISIONS) * length(VARIANTS) * NTRIALS *
                  sum(round(Int, T_END / dt) for dt in DT_VALUES)

    @printf("Julia:            %s\n", VERSION)
    @printf("Threads:          %d Julia threads, %d CPU threads\n",
            Threads.nthreads(), Sys.CPU_THREADS)
    @printf("CPU:              %s\n", Sys.cpu_info()[1].model)
    @printf("ERT3D:            %s\n", version)
    @printf("Grid:             N = %d (%d cells), Central6 + Kennedy-Gruber\n", N, N^3)
    @printf("Simulated time:   %.3f per run\n", T_END)
    @printf("Timesteps dt:     %s\n", string(DT_VALUES))
    @printf("Precisions:       %s\n", join(string.(PRECISIONS), ", "))
    @printf("Trials:           %d per configuration (order shuffled)\n", NTRIALS)
    @printf("Timed step! calls: about %d in total\n", total_steps)
    println()
end

function print_table(rows, base_state_name)
    @printf("  %-24s %9s %9s %8s %10s %12s %10s %9s  %s\n",
            "Variant", "ms/step", "min", "spread", "KiB/step",
            "ns/cell-step", "RHS evals", "vs RK3", "diff vs RK3")
    println("  ", "-" ^ 104)

    for r in rows
        if r.failed !== nothing
            println("  ", rpad(r.name, 24), " FAILED: ", r.failed)
            continue
        end
        diffstr = r.name == base_state_name ? "—" : @sprintf("%.2e", r.diff)
        @printf("  %-24s %9.3f %9.3f %7.1f%% %10.2f %12.1f %10.2f %8.2f×  %s\n",
                r.name, r.med_ms, r.min_ms, r.spread, r.kib,
                r.ns_cell, r.evals, r.ratio, diffstr)
    end
end

# ----------------------------------------------------------------
# Main
# ----------------------------------------------------------------

function main()
    print_header()
    rng = MersenneTwister(1)

    mkpath(dirname(CSV_PATH))
    open(CSV_PATH, "w") do io
        println(io, "precision,variant,dt,trial,nsteps,time_s,bytes,gctime_s")

        for T in PRECISIONS
            println("=" ^ 92)
            println("Precision: $T")
            println("=" ^ 92)

            t_rhs = rhs_cost(T)
            @printf("One residual evaluation: %.3f ms\n", 1e3 * t_rhs)

            # Compile every code path before any timed run
            for variant in VARIANTS
                try
                    advance!(build_simulation(variant, T), DT_VALUES[1], NWARMUP_STEPS)
                catch err
                    @warn "Warm-up failed for $(variant.name)" exception = err
                end
            end

            for dt in DT_VALUES
                nsteps = round(Int, T_END / dt)
                println("\n$T, dt = $dt ($nsteps steps)")

                times = Dict(v.name => Float64[] for v in VARIANTS)
                bytes = Dict(v.name => Int[] for v in VARIANTS)
                gcs = Dict(v.name => Float64[] for v in VARIANTS)
                final_states = Dict{String,Any}()
                failures = Dict{String,String}()

                print("  trials: ")
                for trial in 1:NTRIALS
                    for variant in shuffle(rng, VARIANTS)
                        haskey(failures, variant.name) && continue
                        try
                            r = measure_run(variant, T, dt, nsteps)
                            push!(times[variant.name], r.time)
                            push!(bytes[variant.name], r.bytes)
                            push!(gcs[variant.name], r.gctime)
                            final_states[variant.name] = r.state
                            println(io, join((T, variant.name, dt, trial, nsteps,
                                              r.time, r.bytes, r.gctime), ","))
                        catch err
                            failures[variant.name] = sprint(showerror, err)
                        end
                    end
                    print(".")
                end
                println()

                # Summary rows
                base_time = isempty(times[BASELINE]) ? NaN : median(times[BASELINE])
                rows = []
                for variant in VARIANTS
                    name = variant.name
                    if haskey(failures, name) || isempty(times[name])
                        push!(rows, (name = name,
                                     failed = get(failures, name, "no successful run")))
                        continue
                    end

                    t = times[name]
                    med = median(t)
                    diff = (name == BASELINE || !haskey(final_states, BASELINE)) ? NaN :
                           l2_reconstruction_error(final_states[BASELINE], final_states[name])

                    push!(rows, (
                        name = name,
                        failed = nothing,
                        med_ms = 1e3 * med / nsteps,
                        min_ms = 1e3 * minimum(t) / nsteps,
                        spread = 100 * (maximum(t) - minimum(t)) / med,
                        kib = median(bytes[name]) / nsteps / 1024,
                        ns_cell = 1e9 * med / (nsteps * N^3),
                        evals = (med / nsteps) / t_rhs,
                        ratio = med / base_time,
                        diff = diff,
                    ))
                end
                print_table(rows, BASELINE)
            end
        end
    end

    println("\nRaw per-trial data: $CSV_PATH")
    println("=" ^ 92)
end

main()