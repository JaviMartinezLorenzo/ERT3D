using CairoMakie

# ================================================================
# ERT3D — Temporal convergence order (reads the CSV written by
# scripts/07_convergence_order.jl: precision,variant,dt,error)
# ================================================================

csv_path   = "data/processed/convergence_order.csv"
output_png = "data/processed/convergence_order.png"

# Footer caption. Edit to match the runs that produced the data.
config_note = "N = 32  ·  Central6 + Kennedy–Gruber  ·  Mt₀ = 0.07  ·  t_end = 1  ·  reference: Float64, Δt = 2⁻¹¹"

# Appearance per (precision, variant). Unknown series fall back to gray.
styles = Dict(
    ("Float64", "Increment")         => (label = "Float64",                 color = "#000000", marker = :circle),
    ("Float32", "Shu-Osher")         => (label = "Float32 · Shu–Osher",     color = "#c0392b", marker = :rect),
    ("Float32", "Increment")         => (label = "Float32 · increment",     color = "#2471a3", marker = :utriangle),
    ("Float32", "Increment + Kahan") => (label = "Float32 · increment + Kahan", color = "#1e8449", marker = :diamond),
)
fallback_style = (label = "", color = "#7f7f7f", marker = :circle)

# ----------------------------------------------------------------
# Read CSV (no extra dependencies)
# ----------------------------------------------------------------

function read_series(path)
    isfile(path) || error("CSV not found: $path — run scripts/07_convergence_order.jl first")
    data  = Dict{Tuple{String,String},Tuple{Vector{Float64},Vector{Float64}}}()
    order = Tuple{String,String}[]
    for line in Iterators.drop(eachline(path), 1)
        isempty(strip(line)) && continue
        prec, variant, dt, err = split(strip(line), ',')
        key = (String(prec), String(variant))
        if !haskey(data, key)
            data[key] = (Float64[], Float64[])
            push!(order, key)
        end
        push!(data[key][1], parse(Float64, dt))
        push!(data[key][2], parse(Float64, err))
    end
    # sort every series by dt (ascending)
    for key in order
        x, y = data[key]
        p = sortperm(x)
        data[key] = (x[p], y[p])
    end
    return data, order
end

data, order = read_series(csv_path)

# ----------------------------------------------------------------
# Fits and helpers
# ----------------------------------------------------------------

# Least-squares slope in log-log space: y ≈ 10^c * x^p
function loglog_fit(x, y)
    A = hcat(log10.(x), ones(length(x)))
    p, c = A \ log10.(y)
    return p, c
end

# Intercept of a reference line with fixed slope p through the series
anchor(x, y, p) = sum(log10.(y) .- p .* log10.(x)) / length(x)

function slope_ref!(ax, x1, x2, c, p, mult, label; valign = :top)
    y1 = mult * 10^(c + p * log10(x1))
    y2 = mult * 10^(c + p * log10(x2))
    lines!(ax, [x1, x2], [y1, y2]; color = :gray, linewidth = 2.0, linestyle = :dash)
    text!(ax, x2, y2; text = label, color = :gray, fontsize = 22,
        align = (:left, valign), offset = (8, valign == :top ? -4 : 4))
end

const SUP = Dict('0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
                 '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹', '-' => '⁻')
pow2_label(x) = "2" * join(SUP[c] for c in string(round(Int, log2(x))))

println("Fitted slopes (all points of each series):")
fits = Dict{Tuple{String,String},Float64}()
for key in order
    x, y = data[key]
    p, _ = loglog_fit(x, y)
    fits[key] = p
    println("  ", rpad(join(key, " · "), 34), round(p, digits = 3))
end

# ----------------------------------------------------------------
# Figure
# ----------------------------------------------------------------

set_theme!(theme_latexfonts())

all_dt  = sort(unique(vcat((data[k][1] for k in order)...)))
all_err = vcat((data[k][2] for k in order)...)
ylo = 10.0^(floor(log10(minimum(all_err))) - 0.5)
yhi = 10.0^(ceil(log10(maximum(all_err)))  + 0.5)

fig = Figure(size = (1300, 720), fontsize = 22)

ax = Axis(
    fig[1, 1],
    xscale = log10,
    yscale = log10,
    title = "Temporal convergence order",
    subtitle = rich("relative L", subscript("2"), " error against a Float64 reference"),
    titlesize = 28,
    subtitlesize = 20,
    xlabel = "Δt",
    ylabel = "Relative L₂ error",
    xticks = (all_dt, pow2_label.(all_dt)),
    limits = (minimum(all_dt) / 1.5, maximum(all_dt) * 1.5, ylo, yhi),
    yminorticksvisible = true,
    yminorticks = IntervalsBetween(9),
    xgridvisible = false,
    ygridvisible = false,
    xtickalign = 1,
    ytickalign = 1,
    yminortickalign = 1,
    spinewidth = 1.2,
)

# Curves: line + open markers, one legend entry per series
for key in order
    st = get(styles, key, (label = join(key, " · "), color = fallback_style.color,
                           marker = fallback_style.marker))
    x, y = data[key]

    lab = st.label
    # report the fitted slope only where it is meaningful
    if key == ("Float64", "Increment") || key == ("Float32", "Shu-Osher")
        lab *= "  (slope ≈ $(round(fits[key], digits = 1)))"
    end

    lines!(ax, x, y; color = st.color, linewidth = 2.0, label = lab)
    scatter!(ax, x, y; color = :white, strokecolor = st.color, strokewidth = 1.8,
        marker = st.marker, markersize = 12, label = lab)
end

# Reference slopes: Δt³ below the Float64 curve, Δt⁻¹ above the Shu–Osher curve
x1, x2 = all_dt[min(2, end)], all_dt[min(4, end)]

if haskey(data, ("Float64", "Increment"))
    x, y = data[("Float64", "Increment")]
    slope_ref!(ax, x1, x2, anchor(x, y, 3), 3, 1 / 8, rich("∝ Δt", superscript("3")))
end
if haskey(data, ("Float32", "Shu-Osher"))
    x, y = data[("Float32", "Shu-Osher")]
    slope_ref!(ax, x1, x2, anchor(x, y, -1), -1, 4, rich("∝ Δt", superscript("−1")); valign = :bottom)
end

Legend(fig[1, 2], ax; framevisible = false, labelsize = 20, merge = true, unique = true)

Label(fig[2, 1:2], config_note; fontsize = 18, color = (:black, 0.6))

colgap!(fig.layout, 30)
rowgap!(fig.layout, 10)

# ----------------------------------------------------------------
# Save
# ----------------------------------------------------------------

mkpath(dirname(output_png))
save(output_png, fig; px_per_unit = 2)

display(fig)
println("Saved plot to: $output_png")