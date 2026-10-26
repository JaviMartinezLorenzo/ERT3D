using JLD2
using CairoMakie
using XLSX

# ================================================================
# Pirozzoli (2010) isotropic-turbulence benchmark: ERT3D vs reference
#
#   (a) kinetic energy           T / T0
#   (b) density fluctuations     (rho' / rho0) / Mt0^2
#
# Lines + open markers: ERT3D runs.  Filled markers: digitized
# reference data from Pirozzoli (2010).
# ================================================================

# ----------------------------------------------------------------
# User settings
# ----------------------------------------------------------------

output_png     = "data/processed/pirozzoli_orders_comparison.png"
output_pdf     = "data/processed/pirozzoli_orders_comparison.pdf"   # vector version for reports/papers
reference_file = "data/reference/pirozzoli_density_reference_processed.xlsx"

Mt0  = 0.07
tau  = 8.247860988423223
rho0 = 1.0

t_max = 60

# Footer caption. Edit to match the runs that produced the data.
config_note = "Synthetic solenoidal initialization  ·  k₀ = 6  ·  Mt₀ = 0.07  ·  N = 32³"

# Number of markers drawn per simulation curve. The line always uses every record.
n_markers = 55

# Markers follow the paper: L = 2 triangle up, L = 3 triangle down, L = 4 diamond.
# `sheet` is the tab of `reference_file` that holds the matching digitized curve.
# `color` is black by default (paper style); use distinct colours, e.g.
# Makie.wong_colors()[1:3], if the curves are hard to tell apart.
runs = [
    (label = "ERT3D - 4th order (L = 2)",
     file = "data/raw/pirozzoli_C4/diagnostics.jld2",
     marker = :utriangle, color = :black, sheet = "L2"),
    (label = "ERT3D - 6th order (L = 3)",
     file = "data/raw/pirozzoli_C6/diagnostics.jld2",
     marker = :dtriangle, color = :black, sheet = "L3"),
    (label = "ERT3D - 8th order (L = 4)",
     file = "data/raw/pirozzoli_C8/diagnostics.jld2",
     marker = :diamond, color = :black, sheet = "L4"),
]

@assert tau > 0.0 "tau must be positive."
@assert Mt0 > 0.0 "Mt0 must be positive."

# ----------------------------------------------------------------
# Figure and axes
# ----------------------------------------------------------------

set_theme!(theme_latexfonts())

fig = Figure(size = (1150, 720), fontsize = 22)

# Settings shared by both panels
axis_style = (
    xticks = 0:10:t_max,
    xminorticksvisible = true,
    yminorticksvisible = true,
    xminorticks = IntervalsBetween(5),
    yminorticks = IntervalsBetween(5),
    xtickalign = 1,
    ytickalign = 1,
    xminortickalign = 1,
    yminortickalign = 1,
    xgridvisible = false,
    ygridvisible = false,
    yticklabelspace = 45.0,            # same width in both panels -> aligned y labels
    titlealign = :left,
    titlesize = 22,
    titlefont = :regular,
    spinewidth = 1.2,
)

ax1 = Axis(
    fig[1, 1];
    axis_style...,
    title = "(a)  Kinetic energy",
    ylabel = L"T/T_0",
    limits = (0, t_max, 0.90, 1.10),
    yticks = (0.90:0.05:1.10, ["0.90", "0.95", "1.00", "1.05", "1.10"]),
)

ax2 = Axis(
    fig[2, 1];
    axis_style...,
    title = "(b)  RMS density fluctuation",
    xlabel = L"t/\tau",
    ylabel = L"(\rho'/\rho_0)\,/\,M_{\mathrm{t0}}^{2}",
    limits = (0, t_max, 0, 2.0),
    yticks = 0:0.5:2.0,
)

linkxaxes!(ax1, ax2)
hidexdecorations!(ax1; ticks = false, minorticks = false, grid = false)

# ----------------------------------------------------------------
# Plot one simulation run
# ----------------------------------------------------------------

function plot_run!(ax1, ax2, run)
    d = load(run.file)

    t = d["t"]
    T = d["kinetic_energy"]
    rho_rms = d["density_rms"]

    @assert !isempty(t) "No diagnostic records found in $(run.file)."
    @assert length(t) == length(T) == length(rho_rms) "Inconsistent lengths in $(run.file)."
    @assert T[1] > 0.0 "Initial kinetic energy must be positive in $(run.file)."

    t_norm   = t ./ tau
    T_norm   = T ./ T[1]
    rho_norm = (rho_rms ./ rho0) ./ Mt0^2

    stride = max(1, cld(length(t_norm), n_markers))
    idx = unique(vcat(1:stride:length(t_norm), length(t_norm)))

    println("$(run.label): records = $(length(t)), t_final/tau = $(t_norm[end]), ",
            "T/T0 final = $(T_norm[end]), rho final = $(rho_norm[end])")

    for (ax, y) in ((ax1, T_norm), (ax2, rho_norm))
        lines!(ax, t_norm, y; color = run.color, linewidth = 1.4)
        scatter!(ax, t_norm[idx], y[idx];
            marker = run.marker, color = :white,
            strokecolor = run.color, strokewidth = 1.4, markersize = 10)
    end
end

# ----------------------------------------------------------------
# Plot one digitized reference curve from the paper
# ----------------------------------------------------------------

function plot_reference!(ax2, reference_file, run)
    table = XLSX.readtable(reference_file, run.sheet)
    t_ref   = Float64.(table.data[findfirst(==(:t_tau), table.column_labels)])
    rho_ref = Float64.(table.data[findfirst(==(:rho_norm), table.column_labels)])

    @assert !isempty(t_ref) "No rows found in sheet $(run.sheet) of $(reference_file)."

    scatter!(ax2, t_ref, rho_ref;
        marker = run.marker, color = run.color, markersize = 8)
end

# ----------------------------------------------------------------
# Load and plot everything (references last, so they stay on top)
# ----------------------------------------------------------------

println("============================================================")
println("Pirozzoli comparison of central schemes")
println("============================================================")

foreach(run -> plot_run!(ax1, ax2, run), runs)
foreach(run -> plot_reference!(ax2, reference_file, run), runs)

println("============================================================")

# ----------------------------------------------------------------
# Legend (outside the axes, grouped: marker shape = scheme, style = data source)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# Legend (Semantic grouping: Shape = Order, Style = Formulation)
# ----------------------------------------------------------------

# Group 1: Spatial Order (mapped to marker shape, neutral fill)
order_entries = [
    MarkerElement(marker = :utriangle, color = :transparent, strokecolor = :black, strokewidth = 1.4, markersize = 14),
    MarkerElement(marker = :dtriangle, color = :transparent, strokecolor = :black, strokewidth = 1.4, markersize = 14),
    MarkerElement(marker = :diamond,   color = :transparent, strokecolor = :black, strokewidth = 1.4, markersize = 14)
]
order_labels = ["4th order (L = 2)", "6th order (L = 3)", "8th order (L = 4)"]

# Group 2: Numerical Formulation (mapped to fill style and lines)
method_entries = [
    # ERT3D: Line + open (white) marker
    [LineElement(color = :black, linewidth = 1.4),
     MarkerElement(marker = :circle, color = :white, strokecolor = :black, 
                   strokewidth = 1.4, markersize = 11)],
    # Pirozzoli: Solid black marker
    MarkerElement(marker = :circle, color = :black, markersize = 10)
]
method_labels = [
    "ERT3D (D-KG-SF)",
    "Pirozzoli 2010 (C-KG-SF)"
]

Legend(
    fig[1:2, 2],
    [order_entries, method_entries],
    [order_labels, method_labels],
    ["Spatial Order", "Numerical Formulation"];
    framevisible = false,
    labelsize = 16,
    titlesize = 18,
    gridshalign = :left,
    titlehalign = :left
)

# ----------------------------------------------------------------
# Titles, footer, layout
# ----------------------------------------------------------------

Label(fig[0, 1:2], "Isotropic Turbulence: Reference Comparison\nDirect vs. Conservative Kennedy-Gruber Splitting"; 
    fontsize = 26, font = :bold, justification = :center)
Label(fig[3, 1:2], config_note; fontsize = 18, color = (:black, 0.6))

rowgap!(fig.layout, 1, 12)
colgap!(fig.layout, 30)

# ----------------------------------------------------------------
# Save
# ----------------------------------------------------------------

mkpath(dirname(output_png))
save(output_png, fig; px_per_unit = 2)
save(output_pdf, fig)

display(fig)

println("Saved plot to: $output_png and $output_pdf")
