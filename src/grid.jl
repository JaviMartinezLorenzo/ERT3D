"""
    Grid

2π-periodic Cartesian grid, matching the isotropic-turbulence/TGV convention
(domain = (2π)^3, so integer wavenumbers fall out naturally for FFT-based
diagnostics and comparison against Brachet et al. / Pirozzoli).

Fields to include:
- N::Int              # points per direction (e.g. 32, 64)
- Δx::Float64         # = 2π / N
- x, y, z             # coordinate vectors, 0:Δx:2π-Δx
- k                   # wavenumber vector for FFT-based spectra (if using FFTW)

Note: this is metadata only — no solution data lives here (that's State).
"""
struct Grid
    # TODO
end

# TODO: Grid(N::Int) constructor computing Δx, coordinate arrays, wavenumbers
