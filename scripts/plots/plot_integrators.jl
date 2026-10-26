using CairoMakie

# ================================================================
# ERT3D — Reversibility studies: explicit RK3 vs implicit midpoint
# ================================================================

output_png = "data/processed/reversibility_studies.png"

# Footer caption
config_note = "N = 32  ·  Central6 + Kennedy–Gruber  ·  Mt₀ = 0.07  ·  Float64  " 

# Reference configuration for plot subtitles
dt_base, Trev_base = 0.025, 8

# Integrator colours. For a monochrome version use c_mid = :gray45.
const c_rk3 = :black
const c_mid = "#0072B2"

# ----------------------------------------------------------------
# Data
# ----------------------------------------------------------------

# Time-step study (T_rev fixed at Trev_base)
dt = [0.001, 0.005, 0.01, 0.025, 0.05, 0.06015037594]

dt_error_rk3 = [
    3.4412395415279355e-12,
    4.10929539854871e-10,
    3.285690542534291e-9,
    5.110110925019787e-8,
    3.997287947629732e-7,
    6.849337599870788e-7,
]

dt_error_midpoint = [
    2.276917799423883e-15,
    1.7186879275056663e-15,
    2.1350472903123794e-15,
    1.7452732948179986e-15,
    6.029927477543876e-14,
    8.882270882761238e-15,
]

# Reversal-time study (Δt fixed at dt_base)
Trev = [1, 2, 4, 8, 16, 20, 32, 40, 60, 70, 80]

Trev_error_rk3 = [
    6.3807813764489065e-9,
    1.2795498188008415e-8,
    2.557088295633877e-8,
    5.110110925019787e-8,
    1.0196922270200766e-7,
    1.273255066287274e-7,
    2.031080063254848e-7,
    2.530851885718827e-7,
    4.514385827391257e-7,
    7.825932115912968e-7,
    2.7258678167069525e-6,
]

Trev_error_midpoint = [
    4.202842869649727e-16,
    5.501121214432302e-16,
    1.1148169089247136e-15,
    1.7452732948179986e-15,
    1.7454278151549984e-15,
    3.567345743120666e-15,
    8.448212644779348e-15,
    8.477825263267002e-15,
    8.750052553745903e-15,
    9.6448879503603e-15,
    1.7853486440338448e-14,
]

@assert length(dt) == length(dt_error_rk3) == length(dt_error_midpoint)
@assert length(Trev) == length(Trev_error_rk3) == length(Trev_error_midpoint)

# Sanity check: both studies must share the same baseline point
i_dt = findfirst(==(dt_base), dt)
i_T  = findfirst(==(Trev_base), Trev)
@assert isapprox(dt_error_rk3[i_dt], Trev_error_rk3[i_T]; rtol = 1e-12) "RK3 baseline differs between studies"
if !isapprox(dt_error_midpoint[i_dt], Trev_error_midpoint[i_T]; rtol = 1e-6)
    @warn "Midpoint baseline differs between studies" dt_error_midpoint[i_dt] Trev_error_midpoint[i_T]
end

# ----------------------------------------------------------------
# Fitted slopes: explicit RK3 only (least squares in log-log space)
# ----------------------------------------------------------------

function loglog_fit(x, y)
    A = hcat(log10.(x), ones(length(x)))
    p, c = A \ log10.(y)          # y ≈ 10^c * x^p
    return p, c
end

Trev_fit_max = 32                  # upper end of the linear-growth regime
lin = Trev .<= Trev_fit_max

p_dt, c_dt = loglog_fit(dt, dt_error_rk3)
p_T,  c_T  = loglog_fit(Trev[lin], Trev_error_rk3[lin])

println("RK3 fitted slope in Δt:               ", round(p_dt, digits = 3))
println("RK3 fitted slope in T_rev (≤ $Trev_fit_max): ", round(p_T, digits = 3))

# ----------------------------------------------------------------
# Plot helpers
# ----------------------------------------------------------------

# Data curve: line with open markers; `label` feeds the shared legend
function curve!(ax, x, y; color, marker, label)
    lines!(ax, x, y; color, linewidth = 2.0, label)
    scatter!(ax, x, y; color = :white, strokecolor = color, strokewidth = 1.8,
        marker, markersize = 11, label)
end

# Reference slope drawn parallel to the fit but shifted down by `shift`,
# so it never sits on top of the data.
function slope_ref!(ax, x1, x2, c, p, shift, label)
    y1 = 10^(c + p * log10(x1)) / shift
    y2 = y1 * (x2 / x1)^p
    lines!(ax, [x1, x2], [y1, y2]; color = :gray, linewidth = 2.0, linestyle = :dash)
    text!(ax, x2, y2; text = label, color = :gray, fontsize = 20,
        align = (:left, :top), offset = (8, -4))
end

# Vertical connector between the two curves at the baseline, with the ratio
function gap_annotation!(ax, x, y_hi, y_lo)
    lines!(ax, [x, x], [y_lo, y_hi]; color = (:gray, 0.9), linewidth = 1.6, linestyle = :dot)
    r = y_hi / y_lo
    e = floor(Int, log10(r))
    m = round(r / 10.0^e; digits = 1)
    if m >= 10
        m /= 10
        e += 1
    end
    text!(ax, x, sqrt(y_hi * y_lo);
        text = rich("≈ ", string(m), "×10", superscript(string(e)), " lower"),
        color = :gray, fontsize = 20, align = (:left, :center), offset = (10, 0))
end

# Machine epsilon reference
function eps_line!(ax, xhi)
    hlines!(ax, [eps(Float64)]; color = :gray, linewidth = 1.6, linestyle = :dot)
    text!(ax, xhi, eps(Float64); text = rich("ε", subscript("64")),
        color = :gray, fontsize = 20, align = (:right, :bottom), offset = (-8, 4))
end

# ----------------------------------------------------------------
# Figure
# ----------------------------------------------------------------

set_theme!(theme_latexfonts())

all_err = vcat(dt_error_rk3, Trev_error_rk3, dt_error_midpoint, Trev_error_midpoint,
               [eps(Float64)])
ylo = 10.0^(floor(log10(minimum(all_err))) - 0.3)
yhi = 10.0^(ceil(log10(maximum(all_err))) + 0.3)

xlim_dt   = (8e-4, 8e-2)
xlim_trev = (0.8, 100.0)

fig = Figure(size = (1400, 720), fontsize = 22)

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
    # the last point (Δt = 0.0602) is plotted but left unlabeled: its label
    # used to overprint the 0.05 tick
    xticks = ([0.001, 0.0025, 0.005, 0.01, 0.025, 0.05],
              ["0.001", "0.0025", "0.005", "0.01", "0.025", "0.05"]),
    limits = (xlim_dt..., ylo, yhi),
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
    xticks = ([1, 2, 4, 8, 16, 32, 64], ["1", "2", "4", "8", "16", "32", "64"]),
    limits = (xlim_trev..., ylo, yhi),
    yminorticksvisible = true,
    yminorticks = IntervalsBetween(9),
    xgridvisible = false,
    ygridvisible = false,
    xtickalign = 1,
    ytickalign = 1,
    yminortickalign = 1,
    spinewidth = 1.2,
)

# Same vertical range in both panels, so the gap is directly comparable;
# the right panel keeps its ticks but drops the repeated labels.
linkyaxes!(ax_dt, ax_trev)
hideydecorations!(ax_trev; ticks = false, minorticks = false, grid = false)

for (ax, x, rk3, mid, xbase, xhi) in (
        (ax_dt,   dt,   dt_error_rk3,   dt_error_midpoint,   dt_base,   xlim_dt[2]),
        (ax_trev, Trev, Trev_error_rk3, Trev_error_midpoint, Trev_base, xlim_trev[2]),
    )
    eps_line!(ax, xhi)
    curve!(ax, x, rk3; color = c_rk3, marker = :circle,  label = "Explicit RK3")
    curve!(ax, x, mid; color = c_mid, marker = :diamond, label = "Implicit midpoint")

end

# Show the gap annotation only on the reversal-time plot.
i = findfirst(==(Trev_base), Trev)
gap_annotation!(ax_trev, Trev_base, Trev_error_rk3[i], Trev_error_midpoint[i])

# Reference slopes for explicit RK3 only
slope_ref!(ax_dt, 0.0015, 0.008, c_dt, 3, 5,
    rich("∝ Δt", superscript("3"), "  (fit ≈ $(round(p_dt, digits = 2)))"))
slope_ref!(ax_trev, 2.0, 10.0, c_T, 1, 4,
    rich("∝ T", subscript("rev"), "  (fit ≈ $(round(p_T, digits = 1)), T", subscript("rev"),
         " ≤ $Trev_fit_max)"))

# One shared legend below the panels
Legend(fig[2, 1:2], ax_dt; orientation = :horizontal, framevisible = false,
    labelsize = 22, merge = true, unique = true, tellheight = true)

Label(fig[3, 1:2], config_note; fontsize = 18, color = (:black, 0.6))

colgap!(fig.layout, 30)
rowgap!(fig.layout, 10)

# ----------------------------------------------------------------
# Save
# ----------------------------------------------------------------

mkpath(dirname(output_png))
save(output_png, fig; px_per_unit = 2)

display(fig)
println("Saved plot to: $output_png")