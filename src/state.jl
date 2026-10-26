
"""
    State{T}

Container for the conserved variables (rho, rhou, rhov, rhow, rhoE) on a
3D grid.

A struct-of-arrays layout (one full 3D field per conserved variable) is
used for cache efficiency and FFT-friendly operations.
"""
struct State{T<:AbstractFloat}
    rho::Array{T,3}
    rhou::Array{T,3}
    rhov::Array{T,3}
    rhow::Array{T,3}
    rhoE::Array{T,3}
end

"""
    State(grid::Grid{T})

Construct a zero-initialized state sized to `grid`, using its precision.
"""
function State(grid::Grid{T}) where {T<:AbstractFloat}
    dims = (grid.N, grid.N, grid.N)

    return State{T}(
        zeros(T, dims),
        zeros(T, dims),
        zeros(T, dims),
        zeros(T, dims),
        zeros(T, dims),
    )
end


"""
    PrimitiveState{T}

Container for the primitive variables (rho, u, v, w, p) on a
3D grid.

A struct-of-arrays layout (one full 3D field per primitive variable) is
used for cache efficiency and FFT-friendly operations.
"""
struct PrimitiveState{T<:AbstractFloat}
    rho::Array{T,3}
    u::Array{T,3}
    v::Array{T,3}
    w::Array{T,3}
    p::Array{T,3}
end

"""
    PrimitiveState(grid::Grid{T})

Construct a zero-initialized primitive state sized to `grid`, using its
precision.
"""
function PrimitiveState(grid::Grid{T}) where {T<:AbstractFloat}
    dims = (grid.N, grid.N, grid.N)

    return PrimitiveState{T}(
        zeros(T, dims),
        zeros(T, dims),
        zeros(T, dims),
        zeros(T, dims),
        zeros(T, dims),
    )
end


"""
    copy_state!(dest, src)

Copy all conserved variables from `src` into `dest` in-place.
"""
function copy_state!(
    dest::State{T},
    src::State{T},
) where {T<:AbstractFloat}
    copyto!(dest.rho,  src.rho)
    copyto!(dest.rhou, src.rhou)
    copyto!(dest.rhov, src.rhov)
    copyto!(dest.rhow, src.rhow)
    copyto!(dest.rhoE, src.rhoE)

    return dest
end


"""
    axpy!(dest, α, src)

Compute `dest = dest + α * src` for all conserved variables in-place.
"""
function axpy!(
    dest::State{T},
    α::S,
    src::State{T},
) where {T<:AbstractFloat, S<:AbstractFloat}
    @. dest.rho  += α * src.rho
    @. dest.rhou += α * src.rhou
    @. dest.rhov += α * src.rhov
    @. dest.rhow += α * src.rhow
    @. dest.rhoE += α * src.rhoE

    return dest
end


"""
    linear_combination!(dest, α, a, β, b)

Compute `dest = α*a + β*b` for all conserved variables in-place.
"""
function linear_combination!(
    dest::State{T},
    α::T,
    a::State{T},
    β::T,
    b::State{T},
) where {T<:AbstractFloat}
    @. dest.rho  = α * a.rho  + β * b.rho
    @. dest.rhou = α * a.rhou + β * b.rhou
    @. dest.rhov = α * a.rhov + β * b.rhov
    @. dest.rhow = α * a.rhow + β * b.rhow
    @. dest.rhoE = α * a.rhoE + β * b.rhoE

    return dest
end
