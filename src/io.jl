"""
IO utilities for ERT3D.
"""

using WriteVTK
using WriteVTK

# -------------------------------------------------------------------------
# Internal helper
# -------------------------------------------------------------------------

function build_vtk(state::State,
                   grid::Grid,
                   params::Parameters,
                   path::String)

    primitive = primitive_variables(state, params)

    X, Y, Z = meshgrid(grid)

    vtk = vtk_grid(path, X, Y, Z)
    velocity = velocity_magnitude(state)

    vtk["density"]  = primitive.rho
    vtk["pressure"] = primitive.p
    vtk["velocity"] = (primitive.u,
                       primitive.v,
                       primitive.w)
    vtk["velocity magnitude"]    = velocity               

    return vtk
end

# -------------------------------------------------------------------------
# Export a single snapshot
# -------------------------------------------------------------------------

"""
    export_vtk(state, grid, params, path)

Export a single solution snapshot as a `.vts` file.
"""
function export_vtk(state::State,
                    grid::Grid,
                    params::Parameters,
                    path::String)

    vtk = build_vtk(state, grid, params, path)

    vtk_save(vtk)

    return nothing
end

# -------------------------------------------------------------------------
# Time series (.pvd)
# -------------------------------------------------------------------------

mutable struct VTKCollection
    pvd::WriteVTK.CollectionFile
    basepath::String
    counter::Int
end

"""
    VTKCollection(basepath)

Create a ParaView collection.

The resulting `.pvd` file references all exported snapshots and enables
time-dependent visualization.
"""
function VTKCollection(basepath::String)

    return VTKCollection(
        paraview_collection(basepath),
        basepath,
        0
    )

end

"""
    add_snapshot!(collection, state, grid, params, time)

Export one frame and register it in the ParaView collection.
"""
function add_snapshot!(collection::VTKCollection,
                       state::State,
                       grid::Grid,
                       params::Parameters,
                       time::Float64)

    collection.counter += 1

    filename = "$(collection.basepath)_$(collection.counter)"

    vtk = build_vtk(state,
                    grid,
                    params,
                    filename)

    collection.pvd[time] = vtk

    vtk_save(vtk)

    return nothing
end

"""
    close_collection!(collection)

Finalize the ParaView collection.
"""
function close_collection!(collection::VTKCollection)

    vtk_save(collection.pvd)

    return nothing

end