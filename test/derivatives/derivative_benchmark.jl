using ERT3D
using BenchmarkTools

# ================================================================
# Setup
# ================================================================

N = 32

grid = Grid(N)

derivative = Central4(grid)

f = rand(N, N, N)
out = zeros(N, N, N)

println("Central4 derivative benchmark")
println("N = $N")
println()

# ================================================================
# Warm-up
# ================================================================

derivative_x!(out, f, derivative, grid)
derivative_y!(out, f, derivative, grid)
derivative_z!(out, f, derivative, grid)

# ================================================================
# Benchmarks
# ================================================================

println("derivative_x!:")
@btime derivative_x!($out, $f, $derivative, $grid)

println()

println("derivative_y!:")
@btime derivative_y!($out, $f, $derivative, $grid)

println()

println("derivative_z!:")
@btime derivative_z!($out, $f, $derivative, $grid)

# ================================================================
# Allocation counts
# ================================================================

println()
println("Allocation counts:")


println(
    "x: ",
    @allocated derivative_x!(out, f, derivative, grid)
)

println(
    "y: ",
    @allocated derivative_y!(out, f, derivative, grid)
)

println(
    "z: ",
    @allocated derivative_z!(out, f, derivative, grid)
)

println()
println("Individual operations:")

L = length(derivative.coeffs)
xwork = derivative.xwork

println("Full copy:")
@btime $xwork[$(L+1):$(L+N), :, :] .= $f

println("Left ghost:")
@btime $xwork[1:$L, :, :] .= $f[end-$L+1:end, :, :]

println("Left ghost with view:")
@btime $xwork[1:$L, :, :] .= @view $f[end-$L+1:end, :, :]

println("Right ghost:")
@btime $xwork[$(L+N+1):end, :, :] .= $f[1:$L, :, :]

println("Right ghost with view:")
@btime $xwork[$(L+N+1):end, :, :] .= @view $f[1:$L, :, :]