
"""
    Central4(grid::Grid{T}) where {T<:AbstractFloat}

Construct a fourth-order central finite-difference operator.

Uses a five-point stencil and two periodic ghost cells on each side
of the differentiated direction. The precision is inferred from the grid.
"""
function Central4(grid::Grid{T}) where {T<:AbstractFloat}
    N = grid.N
    L = 2

    coeffs = (T(2) / T(3), -T(1) / T(12))

    xwork = zeros(T, N + 2L, N, N)
    ywork = zeros(T, N, N + 2L, N)
    zwork = zeros(T, N, N, N + 2L)

    return Central4{T}(coeffs, xwork, ywork, zwork)
end


"""
    Central6(grid::Grid{T}) where {T<:AbstractFloat}

Construct a sixth-order central finite-difference operator.

Uses a seven-point stencil and three periodic ghost cells on each side
of the differentiated direction. The precision is inferred from the grid.
"""
function Central6(grid::Grid{T}) where {T<:AbstractFloat}
    N = grid.N
    L = 3

    coeffs = (T(3) / T(4), -T(3) / T(20), T(1) / T(60))

    xwork = zeros(T, N + 2L, N, N)
    ywork = zeros(T, N, N + 2L, N)
    zwork = zeros(T, N, N, N + 2L)

    return Central6{T}(coeffs, xwork, ywork, zwork)
end


"""
    Central8(grid::Grid{T}) where {T<:AbstractFloat}

Construct an eighth-order central finite-difference operator.

Uses a nine-point stencil and four periodic ghost cells on each side
of the differentiated direction. The precision is inferred from the grid.
"""
function Central8(grid::Grid{T}) where {T<:AbstractFloat}
    N = grid.N
    L = 4

    coeffs = (
        T(4) / T(5),
        -T(1) / T(5),
        T(4) / T(105),
        -T(1) / T(280),
    )

    xwork = zeros(T, N + 2L, N, N)
    ywork = zeros(T, N, N + 2L, N)
    zwork = zeros(T, N, N, N + 2L)

    return Central8{T}(coeffs, xwork, ywork, zwork)
end









"""
    derivative_x!(out, f, D::CentralDifference, grid::Grid)

Compute the central finite-difference approximation of ∂f/∂x.

The result is written in-place to `out`. Periodic ghost cells stored
in `D.xwork` are updated before evaluating the stencil.
"""
function derivative_x!(
    out::Array{T,3},
    f::Array{T,3},
    D::CentralDifference,
    grid::Grid
) where {T<:AbstractFloat}

    L = length(D.coeffs)
    N = grid.N

    D.xwork[L+1:L+N, :, :] .= f
    D.xwork[1:L, :, :] .= @view f[end-L+1:end, :, :]
    D.xwork[L+N+1:end, :, :] .= @view f[1:L, :, :]

    inv_dx = one(T) / T(grid.dx)

    @inbounds for k in 1:N, j in 1:N, i in 1:N
        idx = i + L
        acc = zero(T)

        for l in 1:L
            acc += D.coeffs[l] *
                (D.xwork[idx+l, j, k] -
                 D.xwork[idx-l, j, k])
        end

        out[i, j, k] = acc * inv_dx
    end

    return out
end


"""
    derivative_y!(out, f, D::CentralDifference, grid::Grid)

Compute the central finite-difference approximation of ∂f/∂y.

The result is written in-place to `out`.
"""
function derivative_y!(
    out::Array{T,3},
    f::Array{T,3},
    D::CentralDifference,
    grid::Grid
) where {T<:AbstractFloat}

    L = length(D.coeffs)
    N = grid.N

    D.ywork[:, L+1:L+N, :] .= f
    D.ywork[:, 1:L, :] .= @view f[:, end-L+1:end, :]
    D.ywork[:, L+N+1:end, :] .= @view f[:, 1:L, :]

    inv_dy = one(T) / T(grid.dx)

    @inbounds for k in 1:N, j in 1:N, i in 1:N
        jdx = j + L
        acc = zero(T)

        for l in 1:L
            acc += D.coeffs[l] *
                (D.ywork[i, jdx+l, k] -
                 D.ywork[i, jdx-l, k])
        end

        out[i, j, k] = acc * inv_dy
    end

    return out
end


"""
    derivative_z!(out, f, D::CentralDifference, grid::Grid)

Compute the central finite-difference approximation of ∂f/∂z.

The result is written in-place to `out`.
"""
function derivative_z!(
    out::Array{T,3},
    f::Array{T,3},
    D::CentralDifference,
    grid::Grid
) where {T<:AbstractFloat}

    L = length(D.coeffs)
    N = grid.N

    D.zwork[:, :, L+1:L+N] .= f
    D.zwork[:, :, 1:L] .= @view f[:, :, end-L+1:end]
    D.zwork[:, :, L+N+1:end] .= @view f[:, :, 1:L]

    inv_dz = one(T) / T(grid.dx)

    @inbounds for k in 1:N, j in 1:N, i in 1:N
        kdx = k + L
        acc = zero(T)

        for l in 1:L
            acc += D.coeffs[l] *
                (D.zwork[i, j, kdx+l] -
                 D.zwork[i, j, kdx-l])
        end

        out[i, j, k] = acc * inv_dz
    end

    return out
end
