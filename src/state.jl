"""
    State

Container for the conserved variables (rho, rhou, rhov, rhow, rhoE) on a
3D grid.

A struct-of-arrays layout (one full 3D field per conserved variable) is
used for cache efficiency and FFT-friendly operations.
"""
struct State
    rho::Array{Float64,3}
    rhou::Array{Float64,3}
    rhov::Array{Float64,3}
    rhow::Array{Float64,3}
    rhoE::Array{Float64,3}
end

"""
    State(grid::Grid)

Construct a zero-initialized state sized to `grid`.
"""
function State(grid::Grid)
    dims = (grid.N, grid.N, grid.N)
    return State(
        zeros(dims),
        zeros(dims),
        zeros(dims),
        zeros(dims),
        zeros(dims),
    )
end

"""
    PrimitiveState

Container for the primitive variables (rho, u, v, w, p) on a
3D grid.

A struct-of-arrays layout (one full 3D field per primitive variable) is
used for cache efficiency and FFT-friendly operations.
"""

struct PrimitiveState
    rho::Array{Float64,3}
    u::Array{Float64,3}
    v::Array{Float64,3}
    w::Array{Float64,3}
    p::Array{Float64,3}
end

"""
    PrimitiveState(grid::Grid)

Construct a zero-initialized primitive state sized to `grid`.
"""
function PrimitiveState(grid::Grid)
    dims = (grid.N, grid.N, grid.N)
    return PrimitiveState(
        zeros(dims),
        zeros(dims),
        zeros(dims),
        zeros(dims),
        zeros(dims),
    )
end


"""
    copy_state!(dest, src)

Copy all conserved variables from `src` into `dest` in-place.
"""
function copy_state!(
    dest::State,
    src::State,
)
    copyto!(dest.rho,  src.rho)
    copyto!(dest.rhou, src.rhou)
    copyto!(dest.rhov, src.rhov)
    copyto!(dest.rhow, src.rhow)
    copyto!(dest.rhoE, src.rhoE)

    return dest
end


"""
    axpy!(dest, α, src)

Compute

    dest = dest + α * src

for all conserved variables in-place.
"""
function axpy!(
    dest::State,
    α::Float64,
    src::State,
)
    @. dest.rho  += α * src.rho
    @. dest.rhou += α * src.rhou
    @. dest.rhov += α * src.rhov
    @. dest.rhow += α * src.rhow
    @. dest.rhoE += α * src.rhoE

    return dest
end


"""
    linear_combination!(dest, α, a, β, b)

Compute

    dest = α*a + β*b

for all conserved variables in-place.
"""
function linear_combination!(
    dest::State,
    α::Float64,
    a::State,
    β::Float64,
    b::State,
)
    @. dest.rho  = α * a.rho  + β * b.rho
    @. dest.rhou = α * a.rhou + β * b.rhou
    @. dest.rhov = α * a.rhov + β * b.rhov
    @. dest.rhow = α * a.rhow + β * b.rhow
    @. dest.rhoE = α * a.rhoE + β * b.rhoE

    return dest
end
