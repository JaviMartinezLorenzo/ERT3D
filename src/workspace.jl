"""
    SpatialWorkspace

Reusable working memory for the spatial discretization of the
compressible Euler equations.

The workspace is allocated once for a given grid and reused across
spatial-operator evaluations to avoid repeated allocation of temporary
arrays.

Fields:

- `primitive` : reusable primitive-variable storage
- `H`         : total specific enthalpy
- `tmp`       : temporary derivative/result storage
- `flux`      : temporary flux/product storage
"""
struct SpatialWorkspace
    primitive::PrimitiveState
    H::Array{Float64,3}
    tmp::Array{Float64,3}
    flux::Array{Float64,3}
end


"""
    SpatialWorkspace(grid)

Allocate a spatial workspace compatible with `grid`.
"""
function SpatialWorkspace(grid::Grid)

    primitive = PrimitiveState(grid)

    H    = similar(primitive.rho)
    tmp  = similar(primitive.rho)
    flux = similar(primitive.rho)

    return SpatialWorkspace(
        primitive,
        H,
        tmp,
        flux,
    )
end


"""
    RK3Workspace

Reusable working memory for the Shu-Osher third-order Runge-Kutta
time integrator.

Fields:

- `stage1`  : first intermediate solution
- `stage2`  : second intermediate solution
- `rhs`     : reusable spatial right-hand-side storage
- `spatial` : reusable workspace for the spatial operator
"""
struct RK3Workspace
    stage1::State
    stage2::State
    rhs::State
    spatial::SpatialWorkspace
end


"""
    RK3Workspace(grid)

Allocate a workspace compatible with `grid` for `ExplicitRK3`.
"""
function RK3Workspace(grid::Grid)

    return RK3Workspace(
        State(grid),
        State(grid),
        State(grid),
        SpatialWorkspace(grid),
    )
end