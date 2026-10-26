using ERT3D

function main()
    # --------------------------------------------------
    # Configuration
    # --------------------------------------------------
    N = 16
    t_end = 0.2
    dt_values =  [1e-2, 5e-3, 2.5e-3]
    T = Float64

    println("Temporal convergence")
    println("--------------------------------------")
    println("Precision = $T")
    println("N = $N")
    println("t_end = $t_end")
    println("Timesteps = $dt_values")
    println()

    # --------------------------------------------------
    # Common numerical setup
    # --------------------------------------------------
    params = Parameters(gamma=T(1.4), Mt0=T(0.1), T=T)
    grid = Grid(N; T=T)
    derivative = Central4(grid)
    scheme = KennedyGruber()
    integrator = ExplicitRK3()

    final_states = State[]

    # --------------------------------------------------
    # Run simulations with successively halved dt
    # --------------------------------------------------
    for dt in dt_values
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

        run!(sim, t_end; dt=dt)

        @assert isapprox(sim.t, t_end; atol=1e-12) (
            "Simulation did not reach t_end: " *
            "t = $(sim.t), expected $t_end"
        )

        @assert all(isfinite, sim.state.rho)
        @assert all(isfinite, sim.state.rhou)
        @assert all(isfinite, sim.state.rhov)
        @assert all(isfinite, sim.state.rhow)
        @assert all(isfinite, sim.state.rhoE)

        push!(final_states, deepcopy(sim.state))

        println("dt = $dt: completed successfully")
    end

    # --------------------------------------------------
    # Estimate temporal convergence order
    # --------------------------------------------------
    errors = [
        l2_reconstruction_error(final_states[2], final_states[1]),
        l2_reconstruction_error(final_states[3], final_states[2]),
    ]

    observed_order = log2(errors[1] / errors[2])

    println()
    println("Results")
    println("--------------------------------------")
    println("Error(dt, dt/2)   = $(errors[1])")
    println("Error(dt/2, dt/4) = $(errors[2])")
    println("Error ratio       = $(errors[1] / errors[2])")
    println("Observed order    = $observed_order")
end

main()
