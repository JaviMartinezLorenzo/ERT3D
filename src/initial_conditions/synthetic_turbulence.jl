using FFTW
using Random
using Statistics


# ================================================================
# Wavenumber grids
# ================================================================

"""
    wavenumber_grids(grid)

Construct the three-dimensional Fourier wavenumber grids and their
magnitude for the periodic computational domain.

The wavenumbers are taken directly from `grid.k`, which must contain
the FFT-compatible integer wavenumbers for the 2π-periodic domain.
"""
function wavenumber_grids(grid::Grid)

    KX = reshape(grid.k, :, 1, 1)
    KY = reshape(grid.k, 1, :, 1)
    KZ = reshape(grid.k, 1, 1, :)

    Kmag = sqrt.(
        KX.^2 .+
        KY.^2 .+
        KZ.^2
    )

    return KX, KY, KZ, Kmag
end


# ================================================================
# Target spectrum
# ================================================================

"""
    target_spectrum(k, k0)

Target modal spectrum shape,

    E(k) ∝ (k/k0)^4 exp(-2(k/k0)^2).

The overall amplitude is imposed separately by matching the desired
rms velocity `Mt0`.
"""
@inline function target_spectrum(k::Float64, k0::Float64)

    k == 0.0 && return 0.0

    q = k / k0

    return q^4 * exp(-2.0 * q^2)
end


"""
    target_mode_amplitude(k, k0)

Return the Fourier-mode amplitude corresponding to the target
spectrum.

For non-zero wavenumbers,

    |û| ∝ sqrt(E(k)) / k.

The absolute normalization is fixed later by the requested rms
velocity.
"""
@inline function target_mode_amplitude(
    k::Float64,
    k0::Float64,
)

    k == 0.0 && return 0.0

    return sqrt(target_spectrum(k, k0)) / k
end


# ================================================================
# Divergence-free projection
# ================================================================

"""
    project_divergence_free!(ux, uy, uz, KX, KY, KZ, Kmag)

Project the Fourier-space velocity field onto the divergence-free
subspace,

    û ← û - k (k·û)/|k|².

The zero mode is explicitly set to zero.
"""
function project_divergence_free!(
    ux,
    uy,
    uz,
    KX,
    KY,
    KZ,
    Kmag,
)

    @inbounds for I in eachindex(ux)

        k2 = Kmag[I]^2

        if k2 == 0.0

            ux[I] = 0.0
            uy[I] = 0.0
            uz[I] = 0.0

        else

            kdotu =
                (
                    ux[I] * KX[I] +
                    uy[I] * KY[I] +
                    uz[I] * KZ[I]
                ) / k2

            ux[I] -= kdotu * KX[I]
            uy[I] -= kdotu * KY[I]
            uz[I] -= kdotu * KZ[I]

        end
    end

    return nothing
end


# ================================================================
# Spectrum shaping
# ================================================================

"""
    shape_to_spectrum!(ux, uy, uz, Kmag, k0)

Rescale each Fourier velocity mode so that its amplitude follows the
prescribed target spectrum.

The direction and phase of each non-zero mode are retained.

The zero mode is set to zero.
"""
function shape_to_spectrum!(
    ux,
    uy,
    uz,
    Kmag,
    k0::Float64,
)

    @inbounds for I in eachindex(ux)

        k = Kmag[I]

        if k == 0.0

            ux[I] = 0.0
            uy[I] = 0.0
            uz[I] = 0.0

            continue
        end

        magnitude = sqrt(
            abs2(ux[I]) +
            abs2(uy[I]) +
            abs2(uz[I])
        )

        if magnitude == 0.0

            ux[I] = 0.0
            uy[I] = 0.0
            uz[I] = 0.0

        else

            desired = target_mode_amplitude(k, k0)
            scale = desired / magnitude

            ux[I] *= scale
            uy[I] *= scale
            uz[I] *= scale

        end
    end

    return nothing
end


# ================================================================
# Initial condition
# ================================================================

"""
    initialize!(state, ::SyntheticTurbulence, grid, params)

Initialize the state with synthetic isotropic turbulence.

The construction is:

1. Generate reproducible Gaussian random velocity fields.
2. Transform them to Fourier space.
3. Project onto the divergence-free subspace.
4. Impose the prescribed modal spectrum.
5. Project again to remove numerical compressible components.
6. Transform back to physical space.
7. Remove numerical mean velocity.
8. Rescale the velocity so that `rms_velocity(state) = params.Mt0`.
9. Set uniform density and thermodynamic pressure.

This initial condition is intended as a correctness gate for the
turbulence implementation, not as the primary reversibility case.
"""
function initialize!(
    state::State,
    ic::SyntheticTurbulence,
    grid::Grid,
    params::Parameters,
)
    Random.seed!(ic.seed)
    N = grid.N
    gamma = params.gamma
    k0 = params.k0

    # ------------------------------------------------------------
    # Reproducible random field
    # ------------------------------------------------------------

    rng = MersenneTwister(12345)

    ux = fft(randn(rng, N, N, N))
    uy = fft(randn(rng, N, N, N))
    uz = fft(randn(rng, N, N, N))

    # ------------------------------------------------------------
    # Fourier wavenumbers
    # ------------------------------------------------------------

    KX, KY, KZ, Kmag = wavenumber_grids(grid)

    # ------------------------------------------------------------
    # Solenoidal projection
    # ------------------------------------------------------------

    project_divergence_free!(
        ux,
        uy,
        uz,
        KX,
        KY,
        KZ,
        Kmag,
    )

    # ------------------------------------------------------------
    # Impose target spectrum
    # ------------------------------------------------------------

    shape_to_spectrum!(
        ux,
        uy,
        uz,
        Kmag,
        k0,
    )

    # Projection once more removes any numerical longitudinal
    # component introduced by the spectrum operation.
    project_divergence_free!(
        ux,
        uy,
        uz,
        KX,
        KY,
        KZ,
        Kmag,
    )

    # ------------------------------------------------------------
    # Back to physical space
    # ------------------------------------------------------------

    u = real.(ifft(ux))
    v = real.(ifft(uy))
    w = real.(ifft(uz))

    # ------------------------------------------------------------
    # Remove numerical mean
    # ------------------------------------------------------------

    u .-= mean(u)
    v .-= mean(v)
    w .-= mean(w)

    # ------------------------------------------------------------
    # Normalize velocity amplitude
    # ------------------------------------------------------------

    current_rms = sqrt(
        (
            sum(abs2, u) +
            sum(abs2, v) +
            sum(abs2, w)
        ) / length(u)
    )

    @assert current_rms > 0.0

    amplitude = params.Mt0 / current_rms

    u .*= amplitude
    v .*= amplitude
    w .*= amplitude

    # ------------------------------------------------------------
    # Conservative variables
    # ------------------------------------------------------------

    rho = state.rho
    rhou = state.rhou
    rhov = state.rhov
    rhow = state.rhow
    rhoE = state.rhoE

    @. rho = 1.0

    @. rhou = rho * u
    @. rhov = rho * v
    @. rhow = rho * w

    # Uniform reference pressure:
    #
    #     p0 = 1 / gamma
    #
    # therefore
    #
    #     rho*e = p0/(gamma-1).

    @. rhoE =
        (1.0 / gamma) / (gamma - 1.0) +
        0.5 * (
            rhou^2 +
            rhov^2 +
            rhow^2
        ) / rho

    return state
end