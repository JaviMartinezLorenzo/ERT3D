using ERT3D
using CairoMakie


function convergence_study()

    # ============================================================
    # Configuration
    # ============================================================

    Ns = [16, 32, 64, 128]

    operators = [
        ("Central 4", Central4, 4),
        ("Central 6", Central6, 6),
        ("Central 8", Central8, 8),
    ]

    errors = Dict{String, Vector{Float64}}()

    # ============================================================
    # Manufactured-solution test
    # ============================================================

    for (name, Operator, theoretical_order) in operators

        E = Float64[]

        for N in Ns

            grid = Grid(N)
            x = grid.x

            # Test function:
            #     f(x,y,z) = sin(x)
            #
            # Exact derivative:
            #     df/dx = cos(x)

            f = reshape(sin.(x), N, 1, 1) .* ones(N, N, N)
            exact = reshape(cos.(x), N, 1, 1) .* ones(N, N, N)

            D = Operator(grid)
            numerical = similar(f)

            derivative_x!(numerical, f, D, grid)

            err = maximum(abs.(numerical .- exact))
            push!(E, err)
        end

        errors[name] = E
    end

    # ============================================================
    # Grid spacing
    # ============================================================

    dx_values = [Grid(N).dx for N in Ns]

    # Plot from coarse grid -> fine grid:
    # large Δx on the left, small Δx on the right.
    order = sortperm(dx_values, rev = true)

    dx = dx_values[order]

    for name in keys(errors)
        errors[name] = errors[name][order]
    end

    # ============================================================
    # Observed convergence orders
    # ============================================================

    println()
    println("Observed convergence orders")
    println("============================")

    for (name, _, theoretical_order) in operators

        E = errors[name]

        observed = [
            log(E[i] / E[i + 1]) / log(2)
            for i in 1:length(E)-1
        ]

        println()
        println("$name (theoretical order = $theoretical_order)")

        for i in eachindex(observed)
            println(
                "Δx: ",
                round(dx[i], sigdigits = 5),
                " → ",
                round(dx[i + 1], sigdigits = 5),
                "    p = ",
                round(observed[i], digits = 3),
            )
        end
    end

    # ============================================================
    # Figure
    # ============================================================

    fig = Figure(
        size = (900, 650),
        fontsize = 20,
    )

    ax = Axis(
        fig[1, 1],

        xscale = log10,
        yscale = log10,

        xlabel = L"\Delta x",
        ylabel = L"L_\infty\ \mathrm{error}",

        # FIX 1: Match the tick labels to the correctly sorted dx values.
        xticks = (dx, string.(Ns[order])),

        xlabelsize = 21,
        ylabelsize = 21,

        xticklabelsize = 17,
        yticklabelsize = 17,

        xminorticksvisible = false,
        yminorticksvisible = false,

        xgridvisible = true,
        ygridvisible = true,

        xminorgridvisible = false,
        yminorgridvisible = false,

        spinewidth = 1.2,
        xgridwidth = 0.8,
        ygridwidth = 0.8,
    )

    # ============================================================
    # Numerical results
    # ============================================================

    for (name, _, _) in operators

        lines!(
            ax,
            dx,
            errors[name],
            linewidth = 2.8,
            label = name,
        )

        scatter!(
            ax,
            dx,
            errors[name],
            markersize = 9,
        )
    end

    # ============================================================
    # Reference convergence slopes
    # ============================================================

    xref = [minimum(dx) * 0.52, minimum(dx) * 0.82]

    reference_y = Dict(
        4 => 2e-7,
        6 => 2e-10,
        8 => 2e-13,
    )

    for order in (4, 6, 8)

        y0 = reference_y[order]

        C = y0 / xref[1]^order
        yref = C .* xref .^ order

        lines!(
            ax,
            xref,
            yref,
            linestyle = :dash,
            linewidth = 1.6,
        )

        # FIX 2: Use %$ for variable interpolation inside L"..." macro.
        text!(
            ax,
            xref[2],
            yref[2],
            text = L"\mathcal{O}(\Delta x^{%$order})",
            fontsize = 14,
            align = (:left, :center),
            offset = (8, 0),
        )
    end

    # ============================================================
    # Axis limits
    # ============================================================

    xlims!(
        ax,
        minimum(dx) * 0.4,
        maximum(dx) * 1.15,
    )

    # ============================================================
    # Legend
    # ============================================================

    axislegend(
        ax,
        position = :rb,
        framevisible = true,
        labelsize = 16,
        patchsize = (30, 2),
        padding = (10, 10, 10, 10),
        rowgap = 4,
    )

    # ============================================================
    # Save
    # ============================================================

    mkpath("figures")

    save(
        "figures/derivative_convergence.png",
        fig,
        px_per_unit = 3,
    )

    println()
    println("Figure saved:")
    println("  figures/derivative_convergence.png")

    return fig
end


convergence_study()