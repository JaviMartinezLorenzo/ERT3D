"""
Basic verification and performance check of the central finite-difference
operators.

Tests the derivative of

    f(x,y,z) = sin(x)

whose exact x-derivative is

    df/dx = cos(x).

The verification checks the accuracy of the fourth-, sixth-, and
eighth-order central finite-difference operators.

The performance section measures the cost of one x-derivative evaluation
for increasing grid resolutions.
"""

using ERT3D
using BenchmarkTools


# ================================================================
# Accuracy verification
# ================================================================

N = 16
grid = Grid(N)

x = grid.x

# Test field: f(x,y,z) = sin(x)
f = reshape(sin.(x), N, 1, 1) .* ones(N, N, N)

# Exact derivative: df/dx = cos(x)
exact = reshape(cos.(x), N, 1, 1) .* ones(N, N, N)

# Central 4th order
D4 = Central4(grid)

# Central 6th order
D6 = Central6(grid)

# Central 8th order
D8 = Central8(grid)

out = similar(f)

derivative_x!(out, f, D4, grid)
error = maximum(abs.(out .- exact))
println("Central4 maximum error = ", error)

derivative_x!(out, f, D6, grid)
error = maximum(abs.(out .- exact))
println("Central6 maximum error = ", error)

derivative_x!(out, f, D8, grid)
error = maximum(abs.(out .- exact))
println("Central8 maximum error = ", error)


# ================================================================
# Performance / scaling
# ================================================================

function benchmark_derivative()

    println()
    println("Derivative performance")
    println("=======================")

    for N in (32, 64, 128)

        grid = Grid(N)
        x = grid.x

        f =
            reshape(sin.(x), N, 1, 1) .*
            ones(N, N, N)

        out = similar(f)
        D = Central4(grid)

        println()
        println("N = $N")

        # Warm-up
        derivative_x!(out, f, D, grid)

        @btime derivative_x!(
            $out,
            $f,
            $D,
            $grid,
        )
    end
end

benchmark_derivative()