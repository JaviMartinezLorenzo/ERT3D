"""
Synthetic isotropic turbulence initial condition.

Generates a random, divergence-free velocity field with target spectrum

    E(k) ∝ (k/k0)^4 exp[-2(k/k0)^2]

as used in the Pirozzoli correctness test.

The field is normalized so that

    u_rms = Mt0

with rho = 1 and uniform initial pressure.
"""

using FFTW
using Random
using Statistics


# ================================================================
# Wavenumber grids
# ================================================================

"""
    wavenumber_grids(grid)

Construct the 3D Cartesian wavenumber grids and their magnitude.
"""
function wavenumber_grids(grid::Grid)
    N = grid.N
    KX = reshape(grid.k, :, 1, 1) .* ones(1, N, N)
    KY = reshape(grid.k, 1, :, 1) .* ones(N, 1, N)
    KZ = reshape(grid.k, 1, 1, :) .* ones(N, N, 1)
    Kmag = sqrt.(KX.^2 .+ KY.^2 .+ KZ.^2)
    return KX, KY, KZ, Kmag
end


# ================================================================
# Target spectrum
# ================================================================

"""
    target_spectrum(k, k0)

Unnormalised target shell spectrum

    E(k) ∝ (k/k0)^4 exp[-2(k/k0)^2].
"""
@inline function target_spectrum(k, k0)

    k == 0 && return 0.0

    q = k / k0

    return q^4 * exp(-2.0 * q^2)
end


# ================================================================
# Shell indexing
# ================================================================

"""
    shell_indices(Kmag)

Assign every Fourier mode to an integer wavenumber shell.
"""
function shell_indices(Kmag)

    return round.(Int, Kmag)

end


"""
    shell_counts(shells)

Count the number of discrete Fourier modes in every shell.
"""
function shell_counts(shells)

    kmax = maximum(shells)

    counts = zeros(Int, kmax + 1)

    @inbounds for s in shells
        counts[s + 1] += 1
    end

    return counts
end


# ================================================================
# Divergence-free projection
# ================================================================

"""
    project_divergence_free!(ux, uy, uz, KX, KY, KZ, Kmag)

Project the Fourier-space velocity field onto the solenoidal
(divergence-free) subspace:

    u_hat <- u_hat - k (k · u_hat) / |k|²

The zero-wavenumber mode is explicitly removed.
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

    @inbounds for idx in eachindex(Kmag)

        k2 = Kmag[idx]^2

        if k2 == 0.0

            ux[idx] = 0.0
            uy[idx] = 0.0
            uz[idx] = 0.0

        else

            kdotu =
                ux[idx] * KX[idx] +
                uy[idx] * KY[idx] +
                uz[idx] * KZ[idx]

            factor = kdotu / k2

            ux[idx] -= factor * KX[idx]
            uy[idx] -= factor * KY[idx]
            uz[idx] -= factor * KZ[idx]

        end
    end

    return nothing
end


# ================================================================
# Spectrum shaping
# ================================================================

"""
    shape_to_spectrum!(ux, uy, uz, Kmag, k0)

Rescale the Fourier-space velocity field so that the total energy
contained in each discrete shell follows the prescribed target
spectrum.
"""
function shape_to_spectrum!(
    ux,
    uy,
    uz,
    Kmag,
    k0,
)

    shells = shell_indices(Kmag)
    counts = shell_counts(shells)

    kmax = maximum(shells)

    for k in 1:kmax

        n = counts[k + 1]

        # No modes in this shell.
        n == 0 && continue

        target = target_spectrum(k, k0)

        target == 0.0 && continue

        current_energy = 0.0

        @inbounds for idx in eachindex(shells)

            if shells[idx] == k

                current_energy +=
                    abs2(ux[idx]) +
                    abs2(uy[idx]) +
                    abs2(uz[idx])

            end
        end

        current_energy == 0.0 && continue

        scale = sqrt(target / current_energy)

        @inbounds for idx in eachindex(shells)

            if shells[idx] == k

                ux[idx] *= scale
                uy[idx] *= scale
                uz[idx] *= scale

            end
        end
    end

    # Remove mean mode explicitly.
    ux[1, 1, 1] = 0.0
    uy[1, 1, 1] = 0.0
    uz[1, 1, 1] = 0.0

    return nothing
end


# ================================================================
# Synthetic turbulence initialization
# ================================================================

"""
    initialize!(state, ::SyntheticTurbulence, grid, params)

Initialize a synthetic isotropic solenoidal turbulent velocity field.

The procedure is:

1. Generate deterministic Gaussian random velocity.
2. Transform to Fourier space.
3. Project onto the divergence-free subspace.
4. Impose the target shell spectrum.
5. Transform back to physical space.
6. Normalize the velocity so that `rms_velocity(state) == Mt0`.
7. Set uniform density and pressure.
"""
function initialize!(
    state::State,
    ic::SyntheticTurbulence,
    grid::Grid,
    params::Parameters,
)

    N = grid.N

    gamma = params.gamma
    Mt0 = params.Mt0
    k0 = params.k0

    # ------------------------------------------------------------
    # 1. Deterministic random velocity
    # ------------------------------------------------------------

    rng = MersenneTwister(ic.seed)

    u0 = randn(rng, N, N, N)
    v0 = randn(rng, N, N, N)
    w0 = randn(rng, N, N, N)

    # ------------------------------------------------------------
    # 2. Fourier transform
    # ------------------------------------------------------------

    uh = fft(u0)
    vh = fft(v0)
    wh = fft(w0)

    # ------------------------------------------------------------
    # 3. Wavenumber grids
    # ------------------------------------------------------------

    KX, KY, KZ, Kmag =
        wavenumber_grids(grid)

    # ------------------------------------------------------------
    # 4. Divergence-free projection
    # ------------------------------------------------------------

    project_divergence_free!(
        uh,
        vh,
        wh,
        KX,
        KY,
        KZ,
        Kmag,
    )

    # ------------------------------------------------------------
    # 5. Impose target spectrum
    # ------------------------------------------------------------

    shape_to_spectrum!(
        uh,
        vh,
        wh,
        Kmag,
        k0,
    )

    # ------------------------------------------------------------
    # 6. Transform back to physical space
    # ------------------------------------------------------------

    u = real.(ifft(uh))
    v = real.(ifft(vh))
    w = real.(ifft(wh))

    # ------------------------------------------------------------
    # 7. Set uniform density
    # ------------------------------------------------------------

    @. state.rho = 1.0

    # ------------------------------------------------------------
    # 8. Store velocity as momentum
    # ------------------------------------------------------------

    @. state.rhou = u
    @. state.rhov = v
    @. state.rhow = w

    # ------------------------------------------------------------
    # 9. Normalize to desired turbulent Mach number
    # ------------------------------------------------------------

    current_rms = rms_velocity(state)

    @assert current_rms > 0.0

    amplitude = Mt0 / current_rms

    @. state.rhou *= amplitude
    @. state.rhov *= amplitude
    @. state.rhow *= amplitude

    # ------------------------------------------------------------
    # 10. Uniform thermodynamic state
    #
    # p0 = 1/gamma
    # rho0 = 1
    #
    # rhoE = p/(gamma-1) + 1/2 rho |u|²
    # ------------------------------------------------------------

    @. state.rhoE =
        (1.0 / gamma) / (gamma - 1.0) +
        0.5 * (
            state.rhou^2 +
            state.rhov^2 +
            state.rhow^2
        ) / state.rho

    return state
end