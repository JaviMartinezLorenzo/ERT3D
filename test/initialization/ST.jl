using ERT3D
using FFTW
using Statistics
using Printf

# ================================================================
# Configuration
# ================================================================

N = 64
gamma = 1.4
Mt0 = 0.07
k0 = 6.0
seed = 12345

grid = Grid(N)

params = Parameters(
    gamma = gamma,
    Mt0   = Mt0,
    k0    = k0,
)

state = State(grid)

ic = SyntheticTurbulence()

initialize!(
    state,
    ic,
    grid,
    params,
)

# ================================================================
# 1. Basic physical checks
# ================================================================

println()
println("============================================================")
println("Synthetic turbulence correctness gate")
println("============================================================")
println("N     = ", N)
println("Mt0   = ", Mt0)
println("k0    = ", k0)
println("seed  = ", seed)
println()

# Mean velocity
u = state.rhou ./ state.rho
v = state.rhov ./ state.rho
w = state.rhow ./ state.rho

mean_u = mean(u)
mean_v = mean(v)
mean_w = mean(w)

println("Mean velocity:")
@printf("  <u> = %.6e\n", mean_u)
@printf("  <v> = %.6e\n", mean_v)
@printf("  <w> = %.6e\n", mean_w)

# RMS velocity
urms = rms_velocity(state)

println()
@printf("Velocity RMS = %.8f\n", urms)
@printf("Target Mt0   = %.8f\n", Mt0)

# ================================================================
# 2. Divergence-free check in Fourier space
# ================================================================

uh = fft(u)
vh = fft(v)
wh = fft(w)

KX, KY, KZ, Kmag = ERT3D.wavenumber_grids(grid)

div_hat =
    KX .* uh .+
    KY .* vh .+
    KZ .* wh

# Relative divergence error.
div_norm = sqrt(sum(abs2, div_hat))
vel_norm = sqrt(
    sum(abs2, uh) +
    sum(abs2, vh) +
    sum(abs2, wh)
)

relative_divergence = div_norm / vel_norm

println()
@printf("Relative spectral divergence = %.6e\n", relative_divergence)

# ================================================================
# 3. Empirical energy spectrum
# ================================================================

"""
Compute the kinetic-energy spectrum of the physical velocity field.

The returned E[k+1] is the kinetic energy per unit volume contained
in the integer shell k.

Because FFTW's forward transform is unnormalised,

    u_hat = FFT(u),

Parseval gives

    <u²> = sum(|u_hat|²) / N^6.

Therefore

    E(k) = 1/(2 N^6)
           sum_shell (|u_hat|² + |v_hat|² + |w_hat|²).
"""
function empirical_spectrum(state::State, grid::Grid)

    u = state.rhou ./ state.rho
    v = state.rhov ./ state.rho
    w = state.rhow ./ state.rho

    uh = fft(u)
    vh = fft(v)
    wh = fft(w)

    _, _, _, Kmag = ERT3D.wavenumber_grids(grid)

    shells = ERT3D.shell_indices(Kmag)

    kmax = maximum(shells)

    E = zeros(Float64, kmax + 1)

    normalization = 2.0 * grid.N^6

    @inbounds for idx in eachindex(shells)

        k = shells[idx]

        E[k + 1] += (
            abs2(uh[idx]) +
            abs2(vh[idx]) +
            abs2(wh[idx])
        ) / normalization

    end

    return 0:kmax, E
end

k_values, E = empirical_spectrum(state, grid)

# ================================================================
# 4. Compare against target spectrum
# ================================================================

target = zeros(Float64, length(k_values))

for (i, k) in enumerate(k_values)

    target[i] = ERT3D.target_spectrum(
        Float64(k),
        k0,
    )

end

# Ignore k = 0 because both target and turbulence field should
# contain zero mean velocity.
valid = k_values .> 0

# Normalize both spectra by their total energy.
E_normalized =
    E ./ sum(E[valid])

target_normalized =
    target ./ sum(target[valid])

# Relative L2 error of the SPECTRUM SHAPE.
spectrum_error = sqrt(
    sum(
        (E_normalized[valid] .-
         target_normalized[valid]).^2
    )
    /
    sum(target_normalized[valid].^2)
)

println()
println("Spectrum comparison:")
@printf(
    "  Relative spectral L2 error = %.6e\n",
    spectrum_error
)

# ================================================================
# 5. Peak location
# ================================================================

# Don't consider k=0.
peak_index = argmax(E[valid])

valid_k = k_values[valid]

measured_peak = valid_k[peak_index]

println()
println("Spectrum peak:")
println("  Target k0    = ", k0)
println("  Measured peak = ", measured_peak)

# ================================================================
# 6. Energy / RMS consistency
# ================================================================

# The spectrum is energy per unit volume.
#
# Therefore:
#
#     total kinetic energy density
#         = sum(E)
#         = 1/2 <u²+v²+w²>
#
# and therefore
#
#     urms = sqrt(2 sum(E))

spectral_urms = sqrt(
    2.0 * sum(E)
)

println()
println("Energy consistency:")
@printf(
    "  RMS from state    = %.8f\n",
    urms
)

@printf(
    "  RMS from spectrum = %.8f\n",
    spectral_urms
)

@printf(
    "  Difference         = %.6e\n",
    abs(urms - spectral_urms)
)

# ================================================================
# 7. Density and pressure checks
# ================================================================

primitive = primitive_variables(
    state,
    params,
)

rho_min = minimum(primitive.rho)
rho_max = maximum(primitive.rho)

p_min = minimum(primitive.p)
p_max = maximum(primitive.p)

target_pressure = 1.0 / gamma

println()
println("Thermodynamic state:")
@printf("  rho range = [%.6e, %.6e]\n", rho_min, rho_max)
@printf("  p range   = [%.6e, %.6e]\n", p_min, p_max)
@printf("  p target  = %.6e\n", target_pressure)

# ================================================================
# 8. Automated correctness gate
# ================================================================

println()
println("============================================================")
println("CORRECTNESS GATE")
println("============================================================")

# Tolerances
rms_tol = 1e-12
mean_tol = 1e-12
div_tol = 1e-10
thermo_tol = 1e-12

# Spectrum tolerance should NOT be extremely tight.
# The spectrum is intentionally generated on a discrete FFT grid
# and shell-binned, so a tolerance of a few percent is reasonable.
spectrum_tol = 5e-2

rms_ok =
    abs(urms - Mt0) < rms_tol

mean_ok =
    abs(mean_u) < mean_tol &&
    abs(mean_v) < mean_tol &&
    abs(mean_w) < mean_tol

div_ok =
    relative_divergence < div_tol

thermo_ok =
    abs(rho_min - 1.0) < thermo_tol &&
    abs(rho_max - 1.0) < thermo_tol &&
    abs(p_min - target_pressure) < thermo_tol &&
    abs(p_max - target_pressure) < thermo_tol

spectrum_ok =
    spectrum_error < spectrum_tol

println("  RMS normalization : ", rms_ok ? "PASS" : "FAIL")
println("  Zero mean velocity : ", mean_ok ? "PASS" : "FAIL")
println("  Divergence-free    : ", div_ok ? "PASS" : "FAIL")
println("  Thermodynamic IC   : ", thermo_ok ? "PASS" : "FAIL")
println("  Target spectrum    : ", spectrum_ok ? "PASS" : "FAIL")

all_ok =
    rms_ok &&
    mean_ok &&
    div_ok &&
    thermo_ok &&
    spectrum_ok

println()

if all_ok
    println("✓ SYNTHETIC TURBULENCE CORRECTNESS GATE PASSED")
else
    error("✗ SYNTHETIC TURBULENCE CORRECTNESS GATE FAILED")
end