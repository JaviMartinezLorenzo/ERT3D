"""
Diagnostics used across correctness gates and the main reversibility study:

- kinetic_energy(state, grid)         -> K(t), global kinetic energy
                                          (correctness gate: should stay
                                          ~constant for conservative schemes,
                                          per Pirozzoli Fig. 2/3)
- l2_reconstruction_error(state0, state_reversed, grid)
                                       -> ||U_2T - U_0|| for the reversibility metric
- energy_spectrum(state, grid)        -> E(k) via 3D FFT (FFTW), for
                                          high-wavenumber dissipation diagnosis
- density_rms(state)                  -> ρ' rms, second correctness-gate metric
                                          (compare against Pirozzoli's
                                          ρ'/ρ0/Mt0² ≈ 0.35 plateau)
"""

"""
    rms_velocity(state::State)

Compute the instantaneous root-mean-square velocity

    u_rms = sqrt(<u² + v² + w²>)

where <> denotes the spatial average over the computational domain.
"""
function rms_velocity(state::State)

    u = state.rhou ./ state.rho
    v = state.rhov ./ state.rho
    w = state.rhow ./ state.rho

    return sqrt(
        (
            sum(abs2, u) +
            sum(abs2, v) +
            sum(abs2, w)
        ) / length(u)
    )

end

"""
    velocity_magnitude(state::State)

Compute the velocity magnitude

    |u| = sqrt(u² + v² + w²)

at every grid point.
"""
function velocity_magnitude(state::State)

    u = state.rhou ./ state.rho
    v = state.rhov ./ state.rho
    w = state.rhow ./ state.rho

    return sqrt.(u.^2 .+ v.^2 .+ w.^2)

end

# TODO: function l2_reconstruction_error(state0, state_reversed, grid) end
# TODO: function energy_spectrum(state, grid) end
# TODO: function density_rms(state) end
