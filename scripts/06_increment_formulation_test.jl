using ERT3D
using Statistics

function check_finite(state)
    @assert all(isfinite, state.rho)
    @assert all(isfinite, state.rhou)
    @assert all(isfinite, state.rhov)
    @assert all(isfinite, state.rhow)
    @assert all(isfinite, state.rhoE)
end

function run_case(integrator, T, N, t_end, dt_values)
    params = Parameters(gamma=T(1.4), Mt0=T(0.07), T=T)
    grid = Grid(N; T=T)
    derivative = Central6(grid)
    scheme = KennedyGruber()

    final_states = []

    for dt in dt_values
        println("  Running dt = $dt")

        state = State(grid)
        initialize!(state, TaylorGreen(), grid, params)

        workspace = RK3Workspace(grid)

        sim = Simulation(
            state=state,
            grid=grid,
            params=params,
            derivative=derivative,
            scheme=scheme,
            integrator=integrator,
            workspace=workspace,
        )

        run!(sim, t_end; dt=dt, verbose=false)

        @assert isapprox(sim.t, t_end; atol=1e-12) (
            "Simulation did not reach t_end: " *
            "t = $(sim.t), expected $t_end"
        )

        check_finite(sim.state)
        push!(final_states, deepcopy(sim.state))
    end

    # Pairwise differences between successive timestep refinements.
    # Each result compares dt_values[i] against dt_values[i+1].
    pairwise_metrics = [
        state_errors(
            final_states[i + 1],
            final_states[i],
        )
        for i in 1:(length(final_states) - 1)
    ]

    # Preserve the legacy whole-state L2 errors for comparability.
    errors = [metrics.all for metrics in pairwise_metrics]

    # Three timesteps produce two pairwise errors and one observed order.
    ratios = errors[1:end-1] ./ errors[2:end]
    orders = log2.(ratios)

    # Calculate observed orders for each field-group diagnostic.
    diagnostic_names = (:rho, :momentum, :rhoE, :worst)

    diagnostic_orders = Dict{Symbol, Vector{Float64}}()

    for metric in diagnostic_names
        metric_errors = [
            getproperty(metrics, metric)
            for metrics in pairwise_metrics
        ]

        diagnostic_orders[metric] = log2.(
            metric_errors[1:end-1] ./ metric_errors[2:end]
        )
    end

    return (;
        final_states,
        pairwise_metrics,
        errors,
        ratios,
        orders,
        diagnostic_orders,
    )
end

function print_metrics(metrics)
    println("    Legacy whole-state L2 : $(metrics.all)")
    println("    Density fluctuation  : $(metrics.rho)")
    println("    Momentum fluctuations: $(metrics.momentum)")
    println("    Energy fluctuation   : $(metrics.rhoE)")
    println("    Worst group error    : $(metrics.worst)")
end

function main()
    # --------------------------------------------------
    # Configuration: three timesteps only
    # --------------------------------------------------
    N = 32
    t_end = 1.0

    dt_values = [
        0.0625,
        0.03125,
        0.015625,
    ]

    precisions = (Float32, Float64)

    configurations = [
        ("Shu-Osher", ExplicitRK3()),
        ("Increment", ExplicitRK3(formulation=:increment)),
        (
            "Increment + Kahan",
            ExplicitRK3(
                formulation=:increment,
                compensated=true,
            ),
        ),
    ]

    println("ERT3D: RK3 increment and compensation experiment")
    println("=" ^ 60)
    println("N = $N")
    println("Spatial derivative = Central6")
    println("Scheme = Kennedy-Gruber")
    println("t_end = $t_end")
    println("Timesteps = $dt_values")
    println()

    for T in precisions
        println()
        println("=" ^ 60)
        println("Precision: $T")
        println("=" ^ 60)

        results = Dict{String, Any}()

        for (name, integrator) in configurations
            println()
            println("Integrator: $name")
            println("-" ^ 40)

            result = run_case(
                integrator,
                T,
                N,
                t_end,
                dt_values,
            )

            results[name] = result

            println("\nPairwise temporal errors:")

            for i in eachindex(result.pairwise_metrics)
                println(
                    "\n  dt = $(dt_values[i]) vs $(dt_values[i+1])"
                )
                print_metrics(result.pairwise_metrics[i])
            end

            println("\nObserved convergence orders:")

            println("  Legacy whole-state L2: $(only(result.orders))")

            for metric in (:rho, :momentum, :rhoE, :worst)
                println(
                    "  $(metric): " *
                    "$(only(result.diagnostic_orders[metric]))"
                )
            end
        end

        # --------------------------------------------------
        # Compare against Shu-Osher at identical timesteps.
        # --------------------------------------------------
        baseline = results["Shu-Osher"].final_states

        println()
        println("Difference from Shu-Osher baseline")
        println("-" ^ 40)

        for name in ("Increment", "Increment + Kahan")
            states = results[name].final_states

            println("\n$name")

            for i in eachindex(dt_values)
                metrics = state_errors(
                    baseline[i],
                    states[i],
                )

                println("\n  dt = $(dt_values[i])")
                print_metrics(metrics)
            end
        end
    end

    println("\nExperiment completed.")
end

main()