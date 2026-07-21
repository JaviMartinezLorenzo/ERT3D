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

# TODO: function kinetic_energy(state, grid) end
# TODO: function l2_reconstruction_error(state0, state_reversed, grid) end
# TODO: function energy_spectrum(state, grid) end
# TODO: function density_rms(state) end
