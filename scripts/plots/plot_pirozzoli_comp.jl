using JLD2
using CairoMakie
using XLSX

# ================================================================
# C6 comparison: ERT3D (explicit RK), implicit midpoint, Pirozzoli
#
#   top:    T / T0
#   bottom: (rho' / rho0) / Mt0^2
# ================================================================

output_file = "data/processed/pirozzoli_C6_IM_comp.png"

Mt0  = 0.07
k0   = 6
N    = 32
tau  = 8.247860988423223
rho0 = 1.0

n_markers = 55

const c_mid = "#0072B2"

# Pirozzoli reference curve (L = 3), digitized from his Fig. 2.
# One sheet per order (L2, L3, L4), columns t_tau and rho_norm.
pirozzoli_file  = "data/reference/pirozzoli_density_reference_processed.xlsx"
pirozzoli_sheet = "L3"

runs = [
    (label = "ERT3D C6-KG-SF",
     file = "data/raw/pirozzoli_C6/diagnostics.jld2",
     color = :black),
    (label = "ERT3D C6-KG-SF-IM",
     file = "data/raw/C6_KG_IM/diagnostics.jld2",
     color = c_mid),
]

marker = :dtriangle   # L = 3 in the paper's convention

@assert tau > 0.0 "tau must be positive."
@assert Mt0 > 0.0 "Mt0 must be positive."

# ----------------------------------------------------------------
# Figure
# ----------------------------------------------------------------

fig = Figure(size = (800, 720), fontsize = 20)

ax1 = Axis(fig[1, 1],
    ylabel = L"T/T_0",
    limits = (0, 60, 0.90, 1.10),
    yticks = 0.90:0.05:1.10,
    yticklabelspace = 40.0,
    xtickalign = 1, ytickalign = 1,
    xgridvisible = false, ygridvisible = false)

ax2 = Axis(fig[2, 1],
    xlabel = L"t/\tau",
    ylabel = L"\rho'/\rho_0\,/\,M_{\mathrm{t0}}^2",
    limits = (0, 60, 0, 2.0),
    yticks = 0:0.5:2.0,
    yticklabelspace = 40.0,
    xtickalign = 1, ytickalign = 1,
    xgridvisible = false, ygridvisible = false)

hidexdecorations!(ax1, grid = false)

Label(fig[0, 1], "Isotropic Turbulence: Reference Comparison",
    fontsize = 22, font = :bold, tellwidth = false)

Label(fig[3, 1],
    rich("Synthetic solenoidal initialization  ·  k", subscript("0"), " = $k0  ·  M",
         subscript("t0"), " = $Mt0  ·  N = $N", superscript("3")),
    fontsize = 14, color = :gray40, tellwidth = false)

# ----------------------------------------------------------------
# Pirozzoli reference data
# ----------------------------------------------------------------

# Only the density is digitized, so the reference appears in panel (b) only.
# The first row is (0, 0), an axis-origin artifact of the digitization, so skip it.
sheet = XLSX.readxlsx(pirozzoli_file)[pirozzoli_sheet]
data = sheet[:]                      # row 1 is the header, row 2 is the (0, 0) artifact
ref_t   = Float64.(data[3:end, 1])
ref_rho = Float64.(data[3:end, 2])

# ----------------------------------------------------------------
# ERT3D runs
# ----------------------------------------------------------------

for run in runs
    d = load(run.file)

    t       = d["t"]
    T       = d["kinetic_energy"]
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

    lines!(ax1, t_norm, T_norm, color = run.color, linewidth = 1.2)
    scatter!(ax1, t_norm[idx], T_norm[idx],
        marker = marker, color = :white,
        strokecolor = run.color, strokewidth = 1.2, markersize = 9)

    lines!(ax2, t_norm, rho_norm, color = run.color, linewidth = 1.2)
    scatter!(ax2, t_norm[idx], rho_norm[idx],
        marker = marker, color = :white,
        strokecolor = run.color, strokewidth = 1.2, markersize = 9,
        label = run.label)
end

# Reference drawn last so the ERT3D markers and lines do not cover it.
# Filled circles, as in the paper's legend, slightly smaller than the open markers.
scatter!(ax2, ref_t, ref_rho,
    marker = :circle, color = :black, markersize = 6,
    label = "Pirozzoli 2010 (C6-KG-SF)")

axislegend(ax2, position = :rb, framevisible = false)

text!(ax1, 0.015, 0.95, text = "(a)", space = :relative, align = (:left, :top), fontsize = 20)
text!(ax2, 0.015, 0.95, text = "(b)", space = :relative, align = (:left, :top), fontsize = 20)

save(output_file, fig, px_per_unit = 2)
display(fig)

println("Saved plot to: $output_file")