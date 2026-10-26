using JLD2
using CairoMakie

# ================================================================
# Pirozzoli isotropic-turbulence diagnostics
#
# Reproduces the two quantities shown in Pirozzoli (2010), Fig. 2:
#
#   top:    T / T0
#   bottom: (rho' / rho0) / Mt0^2
#
# versus normalized time t / tau.
#
# The JLD2 file is produced by ERT3D's save_diagnostics().
# ================================================================

# ----------------------------------------------------------------
# User settings
# ----------------------------------------------------------------

diagnostics_file = "data/raw/C6_22222/diagnostics.jld2"
output_file      = "data/processed/pirozzoli_st_diagnostics.png"

Mt0 = 0.07
tau = 8.247860988423223

# ----------------------------------------------------------------
# Load diagnostics
# ----------------------------------------------------------------

d = load(diagnostics_file)

t = d["t"]
T = d["kinetic_energy"]
rho_rms = d["density_rms"]

@assert !isempty(t) "No diagnostic records found."
@assert length(t) == length(T) == length(rho_rms) "Diagnostic arrays have inconsistent lengths."
@assert T[1] > 0.0 "Initial kinetic energy must be positive."
@assert tau > 0.0 "tau must be positive."
@assert Mt0 > 0.0 "Mt0 must be positive."

# ----------------------------------------------------------------
# Pirozzoli normalization
# ----------------------------------------------------------------
#
# Pirozzoli defines total kinetic energy as
#
#     T = ∫ rho q² dV
#
# whereas ERT3D stores
#
#     kinetic_energy = ∫ 1/2 rho q² dV.
#
# The factor 1/2 cancels in T/T0, so the stored diagnostic can be
# normalized directly by its initial value.
#
# rho_rms is already the RMS density fluctuation:
#
#     rho' = sqrt(<rho²> - <rho>²)
#
# and rho0 = 1 for our initialization.
#
T_norm = T ./ T[1]
rho0 = 1.0
rho_rms_norm = (rho_rms ./ rho0) ./ Mt0^2

t_norm = t ./ tau

# ----------------------------------------------------------------
# Basic consistency information
# ----------------------------------------------------------------

println("============================================================")
println("Pirozzoli diagnostic plot")
println("============================================================")
println("File       = $diagnostics_file")
println("Mt0        = $Mt0")
println("tau        = $tau")
println("Records    = $(length(t))")
println("t_final    = $(t[end])")
println("t_final/tau = $(t[end] / tau)")
println("T/T0 final = $(T_norm[end])")
println("rho'/rho0/Mt0² final = $(rho_rms_norm[end])")
println("============================================================")

# ----------------------------------------------------------------
# Figure
# ----------------------------------------------------------------
#
# Pirozzoli Fig. 2 uses two vertically stacked panels with the same
# horizontal coordinate t/tau and x-range approximately 0--60.
# ----------------------------------------------------------------

fig = Figure(
    size = (900, 700),
    fontsize = 18,
)

ax1 = Axis(
    fig[1, 1],
    ylabel = L"T/T_0",
    xlabelvisible = false,
    limits = (0, 60, 0.90, 1.10),
    yticks = 0.90:0.05:1.10,
)

ax2 = Axis(
    fig[2, 1],
    xlabel = L"t/\tau",
    ylabel = L"\rho'/\rho_0\,/\,M_{t0}^2",
    limits = (0, 60, 0, 2.0),
    yticks = 0:0.5:2.0,
)

n_markers = 55
stride = max(1, cld(length(t_norm), n_markers))
idx = unique(vcat(1:stride:length(t_norm), length(t_norm)))

# Same fixed space for tick labels so both panels line up.
ax1.yticklabelspace = 40.0
ax2.yticklabelspace = 40.0

lines!(ax1, t_norm, T_norm, color = :black, linewidth = 1.2)
scatter!(ax1, t_norm[idx], T_norm[idx],
    color = :white, strokecolor = :black, strokewidth = 1.2, markersize = 8)

lines!(ax2, t_norm, rho_rms_norm, color = :black, linewidth = 1.2)
scatter!(ax2, t_norm[idx], rho_rms_norm[idx],
    color = :white, strokecolor = :black, strokewidth = 1.2, markersize = 8)

# Panel labels placed inside the axes, so they cannot hit the tick numbers.
text!(ax1, 0.015, 0.95, text = "(a)", space = :relative, align = (:left, :top), fontsize = 18)
text!(ax2, 0.015, 0.95, text = "(b)", space = :relative, align = (:left, :top), fontsize = 18)


# Keep the visual style compact and paper-like.
hidexdecorations!(ax1, grid = false)
ax1.xgridvisible = false
ax2.xgridvisible = false
ax1.ygridvisible = false
ax2.ygridvisible = false


save(output_file, fig, px_per_unit = 2)

display(fig)

println("Saved plot to: $output_file")
