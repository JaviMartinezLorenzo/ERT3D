"""
    Central4(grid::Grid)

Construct a fourth-order central finite-difference operator.

Uses a five-point stencil and two periodic ghost cells on each side
of the differentiated direction.
"""
function Central4(grid::Grid)
    N = grid.N
    L = 2

    coeffs = (2/3, -1/12)

    xwork = zeros(N + 2L, N, N)
    ywork = zeros(N, N + 2L, N)
    zwork = zeros(N, N, N + 2L)

    return Central4(coeffs, xwork, ywork, zwork)
end


"""
    Central6(grid::Grid)

Construct a sixth-order central finite-difference operator.

Uses a seven-point stencil and three periodic ghost cells on each side
of the differentiated direction.
"""
function Central6(grid::Grid)
    N = grid.N
    L = 3

    coeffs = (3/4, -3/20, 1/60)

    xwork = zeros(N + 2L, N, N)
    ywork = zeros(N, N + 2L, N)
    zwork = zeros(N, N, N + 2L)

    return Central6(coeffs, xwork, ywork, zwork)
end


"""
    Central8(grid::Grid)

Construct an eighth-order central finite-difference operator.

Uses a nine-point stencil and four periodic ghost cells on each side
of the differentiated direction.
"""
function Central8(grid::Grid)
    N = grid.N
    L = 4

    coeffs = (4/5, -1/5, 4/105, -1/280)

    xwork = zeros(N + 2L, N, N)
    ywork = zeros(N, N + 2L, N)
    zwork = zeros(N, N, N + 2L)

    return Central8(coeffs, xwork, ywork, zwork)
end


"""
    derivative_x!(out, f, D::CentralDifference, grid::Grid)

Compute the central finite-difference approximation of `∂f/∂x`.

The result is written in-place to `out`. Periodic ghost cells stored
in `D.xwork` are updated before evaluating the stencil.
"""
function derivative_x!(
    out::Array{Float64,3},
    f::Array{Float64,3},
    D::CentralDifference,
    grid::Grid
)
    L = length(D.coeffs)
    N = grid.N

    D.xwork[L+1:L+N, :, :] .= f
    D.xwork[1:L, :, :] .= f[end-L+1:end, :, :]
    D.xwork[L+N+1:end, :, :] .= f[1:L, :, :]

    inv_dx = 1 / grid.dx

    @inbounds for k in 1:N
        for j in 1:N
            for i in 1:N
                idx = i + L
                acc = 0.0

                for l in 1:L
                    acc += D.coeffs[l] *
                        (D.xwork[idx+l, j, k] -
                         D.xwork[idx-l, j, k])
                end

                out[i, j, k] = acc * inv_dx
            end
        end
    end

    return out
end


"""
    derivative_y!(out, f, D::CentralDifference, grid::Grid)

Compute the central finite-difference approximation of `∂f/∂y`.

The result is written in-place to `out`.
"""
function derivative_y!(
    out::Array{Float64,3},
    f::Array{Float64,3},
    D::CentralDifference,
    grid::Grid
)
    L = length(D.coeffs)
    N = grid.N

    D.ywork[:, L+1:L+N, :] .= f
    D.ywork[:, 1:L, :] .= f[:, end-L+1:end, :]
    D.ywork[:, L+N+1:end, :] .= f[:, 1:L, :]

    inv_dy = 1 / grid.dx

    @inbounds for k in 1:N
        for j in 1:N
            for i in 1:N
                jdx = j + L
                acc = 0.0

                for l in 1:L
                    acc += D.coeffs[l] *
                        (D.ywork[i, jdx+l, k] -
                         D.ywork[i, jdx-l, k])
                end

                out[i, j, k] = acc * inv_dy
            end
        end
    end

    return out
end


"""
    derivative_z!(out, f, D::CentralDifference, grid::Grid)

Compute the central finite-difference approximation of `∂f/∂z`.

The result is written in-place to `out`.
"""
function derivative_z!(
    out::Array{Float64,3},
    f::Array{Float64,3},
    D::CentralDifference,
    grid::Grid
)
    L = length(D.coeffs)
    N = grid.N

    D.zwork[:, :, L+1:L+N] .= f
    D.zwork[:, :, 1:L] .= f[:, :, end-L+1:end]
    D.zwork[:, :, L+N+1:end] .= f[:, :, 1:L]

    inv_dz = 1 / grid.dx

    @inbounds for k in 1:N
        for j in 1:N
            for i in 1:N
                kdx = k + L
                acc = 0.0

                for l in 1:L
                    acc += D.coeffs[l] *
                        (D.zwork[i, j, kdx+l] -
                         D.zwork[i, j, kdx-l])
                end

                out[i, j, k] = acc * inv_dz
            end
        end
    end

    return out
end