using ERT3D
using CairoMakie


function convergence_study()

    # ------------------------------------------------------------
    # Configuration
    # ------------------------------------------------------------

    Ns = [16, 32, 64, 128]

    operators = [
        ("Central4", Central4, 4),
        ("Central6", Central6, 6),
        ("Central8", Central8, 8),
    ]

    errors = Dict{String,Vector{Float64}}()


    # ------------------------------------------------------------
    # Compute errors
    # ------------------------------------------------------------

    for (name, Operator, theoretical_order) in operators

        E = Float64[]

        for N in Ns

            grid = Grid(N)
            x = grid.x

            # Test function:  f(x,y,z) = sin(x)
            # Exact derivative: df/dx = cos(x)
            f = reshape(sin.(x), N, 1, 1) .* ones(N, N, N)
            exact = reshape(cos.(x), N, 1, 1) .* ones(N, N, N)

            D = Operator(grid)
            out = similar(f)

            derivative_x!(out, f, D, grid)

            error = maximum(abs.(out .- exact))
            push!(E, error)
        end

        errors[name] = E
    end


    # ------------------------------------------------------------
    # Print observed orders
    # ------------------------------------------------------------

    println()
    println("Observed convergence orders")
    println("============================")

    for (name, Operator, theoretical_order) in operators
        E = errors[name]
        println()
        println("$name (theoretical order = $theoretical_order)")
        for n in 1:length(Ns)-1
            p = log(E[n] / E[n+1]) / log(2)
            println("N = $(Ns[n]) → $(Ns[n+1]):  p = $(round(p, digits=6))")
        end
    end


    # ------------------------------------------------------------
    # Grid spacing
    # ------------------------------------------------------------

    dx = [Grid(N).Δx for N in Ns]


    # ------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------

    fig = Figure(size = (1000, 700), fontsize = 18)

    ax = Axis(
        fig[1, 1],
        xscale = log10,
        yscale = log10,
        xlabel = L"\Delta x",
        ylabel = L"L_\infty\ \mathrm{error}",
    )

    # ------------------------------------------------------------
    # Numerical convergence curves — capture handles for the legend
    # ------------------------------------------------------------

    curve_entries = Vector{Any}()
    curve_labels  = String[]

    for (name, _, _) in operators
        l = lines!(ax, dx, errors[name], linewidth = 2)
        s = scatter!(ax, dx, errors[name], markersize = 10)
        push!(curve_entries, [l, s])   # group line+marker as one legend icon
        push!(curve_labels, name)
    end

    # ------------------------------------------------------------
    # Theoretical convergence slopes
    # ------------------------------------------------------------

    ref_entries = Vector{Any}()
    ref_labels  = String[]

    for order in (4, 6, 8)
        E = errors["Central$order"]
        C = E[1] / dx[1]^order          # anchor to the first numerical point
        reference = C .* dx .^ order

        r = lines!(ax, dx, reference, linestyle = :dash, linewidth = 1.5)
        push!(ref_entries, r)
        push!(ref_labels, "O(Δx^$order)")
    end

    # ------------------------------------------------------------
    # Legend — now matches entries to labels correctly
    # ------------------------------------------------------------

    Legend(
        fig[1, 2],
        vcat(curve_entries, ref_entries),
        vcat(curve_labels, ref_labels),
        "Scheme / reference",
    )

    # ------------------------------------------------------------
    # Save
    # ------------------------------------------------------------

    mkpath("figures")
    save("figures/derivative_convergence.pdf", fig)
    save("figures/derivative_convergence.png", fig, px_per_unit = 3)

    println()
    println("Figure saved:")
    println("  figures/derivative_convergence.pdf")
    println("  figures/derivative_convergence.png")

    return fig
end

convergence_study()