using ERT3D
using Printf


# ================================================================
# Helpers
# ================================================================

function max_state_error(a::State, b::State)

    return maximum([
        maximum(abs.(a.rho  .- b.rho)),
        maximum(abs.(a.rhou .- b.rhou)),
        maximum(abs.(a.rhov .- b.rhov)),
        maximum(abs.(a.rhow .- b.rhow)),
        maximum(abs.(a.rhoE .- b.rhoE)),
    ])
end


# ================================================================
# 1. copy_state!
# ================================================================

function test_copy_state()

    println()
    println("1. copy_state!")
    println("==============")

    N = 32

    grid = Grid(N)
    params = Parameters(1.4, 0.07)

    state = State(grid)
    initialize!(state, TaylorGreen(), grid, params)

    dest  = State(grid)

    copy_state!(dest, state)

    error = max_state_error(dest, state)

    @printf(
        "N = %d    max copy error = %.6e\n",
        N,
        error,
    )

    @assert error == 0.0
end


# ================================================================
# 2. axpy!
#
#     dest <- dest + α * src
# ================================================================

function test_axpy()

    println()
    println("2. axpy!")
    println("========")

    N = 32

    grid = Grid(N)
    params = Parameters(1.4, 0.07)

    a = State(grid)
    initialize!(a, TaylorGreen(), grid, params)
    b = State(grid)

    @. b.rho  = 0.5
    @. b.rhou = -0.25
    @. b.rhov = 0.75
    @. b.rhow = -1.0
    @. b.rhoE = 2.0

    α = 0.3

    expected = a

    @. expected.rho  += α * b.rho
    @. expected.rhou += α * b.rhou
    @. expected.rhov += α * b.rhov
    @. expected.rhow += α * b.rhow
    @. expected.rhoE += α * b.rhoE

    axpy!(a, α, b)

    error = max_state_error(a, expected)

    @printf(
        "N = %d    max AXPY error = %.6e\n",
        N,
        error,
    )

    @assert error == 0.0
end


# ================================================================
# 3. linear_combination!
#
#     dest <- α*a + β*b
# ================================================================

function test_linear_combination()

    println()
    println("3. linear_combination!")
    println("======================")

    N = 32

    grid = Grid(N)

    a    = State(grid)
    b    = State(grid)
    dest = State(grid)

    @. a.rho  = 1.0
    @. a.rhou = 2.0
    @. a.rhov = 3.0
    @. a.rhow = 4.0
    @. a.rhoE = 5.0

    @. b.rho  = -1.0
    @. b.rhou = -2.0
    @. b.rhov = -3.0
    @. b.rhow = -4.0
    @. b.rhoE = -5.0

    α = 0.25
    β = 0.75

    expected = State(grid)

    @. expected.rho  = α * a.rho  + β * b.rho
    @. expected.rhou = α * a.rhou + β * b.rhou
    @. expected.rhov = α * a.rhov + β * b.rhov
    @. expected.rhow = α * a.rhow + β * b.rhow
    @. expected.rhoE = α * a.rhoE + β * b.rhoE

    linear_combination!(
        dest,
        α,
        a,
        β,
        b,
    )

    error = max_state_error(dest, expected)

    @printf(
        "N = %d    max linear-combination error = %.6e\n",
        N,
        error,
    )

    @assert error == 0.0
end


# ================================================================
# 4. Workspace construction
# ================================================================

function test_rk3_workspace()

    println()
    println("4. RK3 workspace")
    println("================")

    N = 32

    grid = Grid(N)

    workspace = ERT3D.RK3Workspace(grid)

    @assert size(workspace.stage1.rho) == (N, N, N)
    @assert size(workspace.stage2.rho) == (N, N, N)
    @assert size(workspace.rhs.rho)    == (N, N, N)

    @assert size(workspace.spatial.tmp)  == (N, N, N)
    @assert size(workspace.spatial.flux) == (N, N, N)
    @assert size(workspace.spatial.H)    == (N, N, N)

    println("N = $N    workspace dimensions correct")
end


# ================================================================
# Run tests
# ================================================================

test_copy_state()
test_axpy()
test_linear_combination()
test_rk3_workspace()

println()
println("All RK3 workspace tests passed.")