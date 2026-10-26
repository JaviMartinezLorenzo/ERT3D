using ERT3D
using Statistics
using Test

function make_sim(::Type{T}, N) where {T<:AbstractFloat}
    grid = Grid(N; T=T)

    params = Parameters(
        gamma = 1.4,
        Mt0 = 0.1,
        T = T,
    )

    derivative = Central4(grid)
    scheme = KennedyGruber()
    integrator = ExplicitRK3(
        formulation=:increment,
        compensated=true,
    )

    state = State(grid)
    initialize!(state, TaylorGreen(), grid, params)

    workspace = RK3Workspace(grid)

    return Simulation(
        state=state,
        grid=grid,
        params=params,
        derivative=derivative,
        scheme=scheme,
        integrator=integrator,
        workspace=workspace,
    )
end

function benchmark_precision(
    ::Type{T};
    N=64,
    dt=0.03125,
    t_end=1.0,
    trials=5,
) where {T<:AbstractFloat}

    @assert dt isa Float64

    nsteps = round(Int, t_end / dt)
    @assert isapprox(nsteps * dt, t_end; atol=1e-12)

    # Warm-up: compile this precision before measuring.
    warmup = make_sim(T, N)
    run!(warmup, t_end; dt=dt, verbose=false)

    @test eltype(warmup.state.rho) == T
    @test isapprox(warmup.t, t_end; atol=1e-12)
    @test all(isfinite, warmup.state.rho)

    # Repeated trials. Setup is excluded from timing.
    times = Float64[]

    for _ in 1:trials
        sim = make_sim(T, N)

        @test eltype(sim.state.rho) == T

        elapsed = @elapsed run!(
            sim,
            t_end;
            dt=dt,
            verbose=false,
        )

        @test isapprox(sim.t, t_end; atol=1e-12)
        @test all(isfinite, sim.state.rho)
        @test all(sim.state.rho .> zero(T))

        push!(times, elapsed)
    end

    median_time = median(times)

    return (
        precision = T,
        median_seconds = median_time,
        seconds_per_step = median_time / nsteps,
        steps = nsteps,
        trial_times = times,
    )
end

N = 64
dt = 0.03125
t_end = 1.0
trials = 5

println("Float32 vs Float64 RK3 benchmark")
println("================================")
println("N = $N")
println("dt = $dt ($(typeof(dt)))")
println("t_end = $t_end")
println("Expected steps = $(round(Int, t_end / dt))")
println("Trials = $trials")
println("Integrator = Increment + Kahan")
println()

dp = benchmark_precision(
    Float64; N=N, dt=dt, t_end=t_end, trials=trials,
)

sp = benchmark_precision(
    Float32; N=N, dt=dt, t_end=t_end, trials=trials,
)

speedup = dp.median_seconds / sp.median_seconds
reduction = 100 * (1 - sp.median_seconds / dp.median_seconds)

println("\nResults")
println("----------------------------------------")
println("Float64 median:      $(round(dp.median_seconds, digits=4)) s")
println("Float32 median:      $(round(sp.median_seconds, digits=4)) s")
println("Float64 per step:    $(round(dp.seconds_per_step * 1e3, digits=4)) ms")
println("Float32 per step:    $(round(sp.seconds_per_step * 1e3, digits=4)) ms")
println("Speedup:             $(round(speedup, digits=3))x")
println("Runtime reduction:   $(round(reduction, digits=2))%")
println("\nFloat64 trial times: ", dp.trial_times)
println("Float32 trial times: ", sp.trial_times)
