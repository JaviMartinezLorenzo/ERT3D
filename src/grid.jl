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
    N::Int
    Δx::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    k::Vector{Float64}
end

function Grid(N::Int)
    Δx = 2π / N
    x = collect(0:Δx:(2π - Δx))
    y = copy(x)
    z = copy(x)
    k = fftfreq(N, 1/Δx) .* 2π
    return Grid(N, Δx, x, y, z, k)
end

"""
    meshgrid(grid::Grid)

Construct the Cartesian coordinate arrays associated with `grid`.

Returns

    (X, Y, Z)

where each array has dimensions (N,N,N).
"""
function meshgrid(grid::Grid)

    X = [x for x in grid.x, y in grid.y, z in grid.z]
    Y = [y for x in grid.x, y in grid.y, z in grid.z]
    Z = [z for x in grid.x, y in grid.y, z in grid.z]

    return X, Y, Z

end