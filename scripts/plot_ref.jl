using JLD2
using CairoMakie
using XLSX

# ================================================================
# Pirozzoli isotropic-turbulence diagnostics, several orders together,
# with digitized reference data from Pirozzoli (2010) overlaid.
#
#   top:    T / T0
#   bottom: (rho' / rho0) / Mt0^2
# ================================================================

# ----------------------------------------------------------------
# User settings
# ----------------------------------------------------------------

output_file = "data/processed/pirozzoli_orders_comparison.png"
reference_file = "data/reference/pirozzoli_density_reference_processed.xlsx"

Mt0 = 0.07
tau = 8.247860988423223
rho0 = 1.0

# Number of markers drawn per simulation curve. The line always uses every record.
n_markers = 55

# Edit the file names to match your runs. Markers follow the paper:
# L = 2 triangle up, L = 3 triangle down, L = 4 diamond.
# sheet is the name of the matching tab in reference_file.
runs = [
    (label = "Central 4th order (L = 2)",
     file = "data/raw/pirozzoli_C4/diagnostics.jld2",
     marker = :utriangle,
     sheet = "L2"),
    (label = "Central 6th order (L = 3)",
     file = "data/raw/pirozzoli_C6/diagnostics.jld2",
     marker = :dtriangle,
     sheet = "L3"),
    (label = "Central 8th order (L = 4)",
     file = "data/raw/pirozzoli_C8/diagnostics.jld2",
     marker = :diamond,
     sheet = "L4"),
]

@assert tau > 0.0 "tau must be positive."
@assert Mt0 > 0.0 "Mt0 must be positive."

# ----------------------------------------------------------------
# Figure
# ----------------------------------------------------------------

fig = Figure(
    size = (800, 650),
    fontsize = 20,
)

ax1 = Axis(
    fig[1, 1],
    ylabel = L"T/T_0",
    limits = (0, 60, 0.90, 1.10),
    yticks = 0.90:0.05:1.10,
    yticklabelspace = 40.0,
    xtickalign = 1,
    ytickalign = 1,
    xgridvisible = false,
    ygridvisible = false,
)

ax2 = Axis(
    fig[2, 1],
    xlabel = L"t/\tau",
    ylabel = L"\rho'/\rho_0\,/\,M_{\mathrm{t0}}^2",
    limits = (0, 60, 0, 2.0),
    yticks = 0:0.5:2.0,
    yticklabelspace = 40.0,
    xtickalign = 1,
    ytickalign = 1,
    xgridvisible = false,
    ygridvisible = false,
)

hidexdecorations!(ax1, grid = false)

# ----------------------------------------------------------------
# Plot one simulation run (function, so loop variables stay local)
# ----------------------------------------------------------------

function plot_run!(ax1, ax2, run)
    d = load(run.file)

    t = d["t"]
    T = d["kinetic_energy"]
    rho_rms = d["density_rms"]

    @assert !isempty(t) "No diagnostic records found in $(run.file)."
    @assert length(t) == length(T) == length(rho_rms) "Inconsistent lengths in $(run.file)."
    @assert T[1] > 0.0 "Initial kinetic energy must be positive in $(run.file)."

    t_norm = t ./ tau
    T_norm = T ./ T[1]
    rho_norm = (rho_rms ./ rho0) ./ Mt0^2

    stride = max(1, cld(length(t_norm), n_markers))
    idx = unique(vcat(1:stride:length(t_norm), length(t_norm)))

    println("$(run.label): records = $(length(t)), t_final/tau = $(t_norm[end]), ",
            "T/T0 final = $(T_norm[end]), rho final = $(rho_norm[end])")

    lines!(ax1, t_norm, T_norm, color = :black, linewidth = 1.2)
    scatter!(ax1, t_norm[idx], T_norm[idx],
        marker = run.marker, color = :white,
        strokecolor = :black, strokewidth = 1.2, markersize = 9)

    lines!(ax2, t_norm, rho_norm, color = :black, linewidth = 1.2)
    scatter!(ax2, t_norm[idx], rho_norm[idx],
        marker = run.marker, color = :white,
        strokecolor = :black, strokewidth = 1.2, markersize = 9,
        label = run.label)
end

# ----------------------------------------------------------------
# Plot one digitized reference curve from the paper
# ----------------------------------------------------------------

function plot_reference!(ax2, reference_file, run)
    table = XLSX.readtable(reference_file, run.sheet)
    t_ref = Float64.(table.data[findfirst(==(:t_tau), table.column_labels)])
    rho_ref = Float64.(table.data[findfirst(==(:rho_norm), table.column_labels)])

    @assert !isempty(t_ref) "No rows found in sheet $(run.sheet) of $(reference_file)."

    scatter!(ax2, t_ref, rho_ref,
        marker = run.marker, color = :black, markersize = 7,
        label = "Pirozzoli (2010), $(run.label)")
end

# ----------------------------------------------------------------
# Load and plot everything
# ----------------------------------------------------------------

println("============================================================")
println("Pirozzoli comparison of central schemes")
println("============================================================")

for run in runs
    plot_run!(ax1, ax2, run)
    plot_reference!(ax2, reference_file, run)
end

println("============================================================")

# The bottom right corner of the lower panel is empty, so the legend fits there.
axislegend(ax2, position = :rb, framevisible = false, labelsize = 14)

# Panel labels inside the axes.
text!(ax1, 0.015, 0.95, text = "(a)", space = :relative, align = (:left, :top), fontsize = 20)
text!(ax2, 0.015, 0.95, text = "(b)", space = :relative, align = (:left, :top), fontsize = 20)

save(output_file, fig, px_per_unit = 2)

display(fig)

println("Saved plot to: $output_file")
