"""
IO helpers:

- save_checkpoint(state, path::String)   -> JLD2, for resumable/long runs
- load_checkpoint(path::String) -> State
- export_vtk(state, grid, path::String)  -> WriteVTK, for ParaView
                                             Q-criterion isosurfaces
                                             (do NOT attempt 3D isosurface
                                             rendering inside Julia itself)
"""

# TODO: function save_checkpoint(state, path) end
# TODO: function load_checkpoint(path) end
# TODO: function export_vtk(state, grid, path) end
