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
    copy(state::State)

Deep copy of every conserved field.
"""
Base.copy(s::State) = State(
    copy(s.rho),
    copy(s.rhou),
    copy(s.rhov),
    copy(s.rhow),
    copy(s.rhoE),
)

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
    copy(state::PrimitiveState)

Deep copy of every primitive field.
"""
Base.copy(s::PrimitiveState) = PrimitiveState(
    copy(s.rho),
    copy(s.u),
    copy(s.v),
    copy(s.w),
    copy(s.p),
)
