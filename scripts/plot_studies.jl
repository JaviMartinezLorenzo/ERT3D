using CairoMakie

# ================================================================
# ERT3D — Reversibility studies (Δt and T_rev dependence)
# ================================================================

output_png = "data/processed/reversibility_studies.png"
output_pdf = "data/processed/reversibility_studies.pdf"   # vector version for reports/papers

# Footer caption. Edit to match the runs that produced the data.
config_note = "N = 32  ·  Central6 + Kennedy–Gruber  ·  Shu–Osher RK3  ·  Mt₀ = 0.07"

# Baseline configuration shared by both studies
dt_base, Trev_base = 0.025, 8

# ----------------------------------------------------------------
# Data
# ----------------------------------------------------------------

# Time-step study (T_rev fixed at Trev_base)
dt = [0.001, 0.005, 0.01, 0.025, 0.05, 0.06015037594]
dt_error = [
    3.4412395415279355e-12,
    4.10929539854871e-10,
    3.285690542534291e-9,
    5.110110925019787e-8,
    3.997287947629732e-7,
    6.849337599870788e-7,
]

# Reversal-time study (Δt fixed at dt_base)
Trev = [1, 2, 4, 8, 16, 20, 32, 40, 70, 80]
Trev_error = [
    6.3807813764489065e-9,
    1.2795498188008415e-8,
    2.557088295633877e-8,
    5.110110925019787e-8,
    1.0196922270200766e-7,
    1.273255066287274e-7,
    2.031080063254848e-7,
    2.530851885718827e-7,
    7.825932115912968e-7,
    2.7258678167069525e-6,
]

# Sanity check: both studies must share the same baseline point
err_base_dt = dt_error[findfirst(==(dt_base), dt)]
err_base_T  = Trev_error[findfirst(==(Trev_base), Trev)]
@assert isapprox(err_base_dt, err_base_T; rtol = 1e-12) "baseline point differs between studies"
err_base = err_base_dt

# ----------------------------------------------------------------
# Fitted slopes (least squares in log-log space)
# ----------------------------------------------------------------

function loglog_fit(x, y)
    A = hcat(log10.(x), ones(length(x)))
    p, c = A \ log10.(y)          # y ≈ 10^c * x^p
    return p, c
end

Trev_fit_max = 32                  # upper end of the linear-growth regime
lin = Trev .<= Trev_fit_max

p_dt, c_dt = loglog_fit(dt, dt_error)
p_T,  c_T  = loglog_fit(Trev[lin], Trev_error[lin])

println("Fitted slope in Δt:                 ", round(p_dt, digits = 3))
println("Fitted slope in T_rev (≤ $Trev_fit_max):   ", round(p_T, digits = 3))

# ----------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------

# Data curve: black line with open circles
function curve!(ax, x, y)
    lines!(ax, x, y; color = :black, linewidth = 2.0)
    scatter!(ax, x, y;
        color = :white, strokecolor = :black, strokewidth = 1.8, markersize = 11)
end

# Reference slope drawn parallel to the fit but shifted down by `shift`,
# so it never sits on top of the data.
function slope_ref!(ax, x1, x2, c, p, shift, label)
    y1 = 10^(c + p * log10(x1)) / shift
    y2 = y1 * (x2 / x1)^p
    lines!(ax, [x1, x2], [y1, y2]; color = :gray, linewidth = 2.0, linestyle = :dash)
    text!(ax, x2, y2; text = label, color = :gray, fontsize = 22,
        align = (:left, :top), offset = (8, -4))
end

# Highlight the shared baseline point
function baseline!(ax, x, y)
    scatter!(ax, [x], [y]; color = :black, markersize = 16)
    text!(ax, x, y; text = "baseline", fontsize = 18,
        align = (:right, :bottom), offset = (-10, 8))
end

# ----------------------------------------------------------------
# Figure
# ----------------------------------------------------------------

set_theme!(theme_latexfonts())

fig = Figure(size = (1400, 640), fontsize = 22)

ax_dt = Axis(
    fig[1, 1],
    xscale = log10,
    yscale = log10,
    title = "Time-step dependence",
    subtitle = rich("T", subscript("rev"), " = $Trev_base"),
    titlesize = 28,
    subtitlesize = 20,
    xlabel = "Δt",
    ylabel = "L₂ reconstruction error",
    # fewer ticks: the previous set (0.025 / 0.05 / 0.06) overlapped
    xticks = ([0.001, 0.0025, 0.005, 0.01, 0.025, 0.05],
              ["0.001", "0.0025", "0.005", "0.01", "0.025", "0.05"]),
    limits = (8e-4, 8e-2, 5e-13, 5e-6),
    yminorticksvisible = true,
    yminorticks = IntervalsBetween(9),
    xgridvisible = false,
    ygridvisible = false,
    xtickalign = 1,
    ytickalign = 1,
    yminortickalign = 1,
    spinewidth = 1.2,
)

ax_trev = Axis(
    fig[1, 2],
    xscale = log10,
    yscale = log10,
    title = "Reversal-time dependence",
    subtitle = "Δt = $dt_base",
    titlesize = 28,
    subtitlesize = 20,
    xlabel = rich("T", subscript("rev")),
    ylabel = "L₂ reconstruction error",
    # powers of two: no overlap between 32/40 or 70/80
    xticks = ([1, 2, 4, 8, 16, 32, 64], ["1", "2", "4", "8", "16", "32", "64"]),
    limits = (0.8, 100, 1e-9, 1e-5),
    yminorticksvisible = true,
    yminorticks = IntervalsBetween(9),
    xgridvisible = false,
    ygridvisible = false,
    xtickalign = 1,
    ytickalign = 1,
    yminortickalign = 1,
    spinewidth = 1.2,
)

# ----------------------------------------------------------------
# Panel (a): Δt study
# ----------------------------------------------------------------

curve!(ax_dt, dt, dt_error)
slope_ref!(ax_dt, 0.0015, 0.008, c_dt, 3, 5, rich("∝ Δt", superscript("3")))
baseline!(ax_dt, dt_base, err_base)

text!(ax_dt, 0.97, 0.05;
    text = "fitted slope ≈ $(round(p_dt, digits = 2))",
    space = :relative, align = (:right, :bottom), fontsize = 20)

# ----------------------------------------------------------------
# Panel (b): T_rev study
# ----------------------------------------------------------------

curve!(ax_trev, Trev, Trev_error)
slope_ref!(ax_trev, 2.0, 10.0, c_T, 1, 4, rich("∝ T", subscript("rev")))
baseline!(ax_trev, Trev_base, err_base)

text!(ax_trev, 0.97, 0.05;
    text = rich("fitted slope (T", subscript("rev"), " ≤ $Trev_fit_max) ≈ $(round(p_T, digits = 2))"),
    space = :relative, align = (:right, :bottom), fontsize = 20)

# ----------------------------------------------------------------
# Panel labels and footer
# ----------------------------------------------------------------

for (ax, lab) in ((ax_dt, "(a)"), (ax_trev, "(b)"))
    text!(ax, 0.03, 0.97; text = lab, space = :relative,
        align = (:left, :top), fontsize = 28, font = :bold)
end

Label(fig[2, 1:2], config_note; fontsize = 18, color = (:black, 0.6))

colgap!(fig.layout, 60)
rowgap!(fig.layout, 10)

# ----------------------------------------------------------------
# Save
# ----------------------------------------------------------------

mkpath(dirname(output_png))
save(output_png, fig; px_per_unit = 2)
save(output_pdf, fig)

display(fig)
println("Saved plot to: $output_png and $output_pdf")
