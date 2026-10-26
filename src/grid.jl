
"""
    Grid(N::Int; T::Type{<:AbstractFloat}=Float64)

Construct a 2π-periodic Cartesian grid with N points per direction.

The domain is (2π)^3, matching the isotropic-turbulence/TGV convention.
The floating-point precision T defaults to Float64 and can be set to
Float32 for single-precision simulations.

Fields:
- N::Int       # points per direction
- dx::T        # grid spacing, 2π/N
- x, y, z      # coordinate vectors
- k            # wavenumber vector for FFT-based diagnostics
"""
struct Grid{T<:AbstractFloat}
    N::Int
    dx::T
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    k::Vector{T}
end

function Grid(N::Int; T::Type{<:AbstractFloat}=Float64)
    N > 0 || throw(ArgumentError("N must be positive"))

    L = T(2π)
    dx = L / T(N)

    x = collect(T(0):dx:(L - dx))
    y = copy(x)
    z = copy(x)

    k = T.(fftfreq(N, one(T) / dx)) .* L

    return Grid{T}(N, dx, x, y, z, k)
end

"""
    meshgrid(grid::Grid)

Construct the Cartesian coordinate arrays associated with `grid`.

Returns `(X, Y, Z)`, each with dimensions `(N, N, N)`.
"""
function meshgrid(grid::Grid{T}) where {T<:AbstractFloat}
    X = [x for x in grid.x, y in grid.y, z in grid.z]
    Y = [y for x in grid.x, y in grid.y, z in grid.z]
    Z = [z for x in grid.x, y in grid.y, z in grid.z]

    return X, Y, Z
end
