abstract type TimeIntegratorWorkspace end

"""
    SpatialWorkspace{T}

Reusable working memory for the spatial discretization of the
compressible Euler equations.
"""
struct SpatialWorkspace{T<:AbstractFloat}
    primitive::PrimitiveState{T}
    H::Array{T,3}
    tmp::Array{T,3}
    flux::Array{T,3}
end

function SpatialWorkspace(grid::Grid{T}) where {T<:AbstractFloat}
    primitive = PrimitiveState(grid)
    H = similar(primitive.rho)
    tmp = similar(primitive.rho)
    flux = similar(primitive.rho)
    return SpatialWorkspace{T}(primitive, H, tmp, flux)
end

"""
    RK3Workspace{T}

Reusable memory for ExplicitRK3. `rhs1` and `rhs2` store the first two
residuals for the increment formulation; `increment` stores the weighted
RK increment; `compensation` stores Kahan compensation for the final
state update. These extra states are allocated once and reused.
"""
struct RK3Workspace{T<:AbstractFloat} <: TimeIntegratorWorkspace
    stage1::State{T}
    stage2::State{T}
    rhs::State{T}
    rhs1::State{T}
    rhs2::State{T}
    increment::State{T}
    compensation::State{T}
    spatial::SpatialWorkspace{T}
end

function RK3Workspace(grid::Grid{T}) where {T<:AbstractFloat}
    return RK3Workspace{T}(
        State(grid),
        State(grid),
        State(grid),
        State(grid),
        State(grid),
        State(grid),
        State(grid),
        SpatialWorkspace(grid),
    )
end

"""
    ImplicitMidpointWorkspace{T}

Reusable working memory for the implicit midpoint time integrator.
"""
struct ImplicitMidpointWorkspace{T<:AbstractFloat} <: TimeIntegratorWorkspace
    Qk::State{T}
    Qnew::State{T}
    Qmid::State{T}
    R::State{T}
    spatial::SpatialWorkspace{T}
    maxiter::Int
end

function ImplicitMidpointWorkspace(
    grid::Grid{T};
    maxiter::Int = 200,
) where {T<:AbstractFloat}
    maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    return ImplicitMidpointWorkspace{T}(
        State(grid),
        State(grid),
        State(grid),
        State(grid),
        SpatialWorkspace(grid),
        maxiter,
    )
end
