"""
IO utilities for ERT3D.
"""

using WriteVTK


# -------------------------------------------------------------------------
# Internal helper
# -------------------------------------------------------------------------

function build_vtk(state::State, grid::Grid, params::Parameters, path::String)
    primitive = primitive_variables(state, params)
    X, Y, Z = meshgrid(grid)

    vtk = vtk_grid(path, X, Y, Z)
    vtk["density"]  = primitive.rho
    vtk["pressure"] = primitive.p
    vtk["velocity"] = (primitive.u, primitive.v, primitive.w)
    vtk["velocity magnitude"] = sqrt.(primitive.u.^2 .+ primitive.v.^2 .+ primitive.w.^2)

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

# -------------------------------------------------------------------------
# Simulation Diagnostics
# -------------------------------------------------------------------------

"""
    Diagnostics

Stores scalar time histories generated during a simulation.

The diagnostics include:

- `kinetic_energy`: total kinetic energy.
- `density_rms`: RMS density fluctuation.
- `total_mass`: total mass.
- `total_momentum_x/y/z`: total momentum components.
- `total_energy`: total energy.
- `max_local_mach`: maximum local Mach number.
"""
mutable struct Diagnostics
    t::Vector{Float64}
    kinetic_energy::Vector{Float64}
    density_rms::Vector{Float64}
    total_mass::Vector{Float64}
    total_momentum_x::Vector{Float64}
    total_momentum_y::Vector{Float64}
    total_momentum_z::Vector{Float64}
    total_energy::Vector{Float64}
    max_local_mach::Vector{Float64}
end

Diagnostics() = Diagnostics(
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
)


"""
    record!(diag, sim)

Compute and store scalar diagnostics for the current simulation state.

All quantities are evaluated in a single grid traversal to avoid creating
temporary three-dimensional arrays.
"""
function record!(
    diag::Diagnostics,
    sim::Simulation,
)

    state = sim.state
    gamma = sim.params.gamma
    dV = sim.grid.Δx^3

    mass = 0.0
    momentum_x = 0.0
    momentum_y = 0.0
    momentum_z = 0.0
    energy = 0.0
    kinetic = 0.0
    density_sum = 0.0
    density_squared_sum = 0.0
    max_mach = 0.0

    @inbounds for k in axes(state.rho, 3)
        for j in axes(state.rho, 2)
            for i in axes(state.rho, 1)

                rho  = state.rho[i, j, k]
                rhou = state.rhou[i, j, k]
                rhov = state.rhov[i, j, k]
                rhow = state.rhow[i, j, k]
                rhoE = state.rhoE[i, j, k]

                u = rhou / rho
                v = rhov / rho
                w = rhow / rho

                velocity_squared = u^2 + v^2 + w^2

                pressure = (gamma - 1.0) * (
                    rhoE -
                    0.5 * rho * velocity_squared
                )

                sound_speed = sqrt(gamma * pressure / rho)

                mach = sqrt(velocity_squared) / sound_speed

                mass += rho
                momentum_x += rhou
                momentum_y += rhov
                momentum_z += rhow
                energy += rhoE
                kinetic += 0.5 * rho * velocity_squared

                density_sum += rho
                density_squared_sum += rho^2

                max_mach = max(max_mach, mach)
            end
        end
    end

    n_cells = length(state.rho)
    mean_density = density_sum / n_cells

    rho_rms = sqrt(max(0.0, density_squared_sum / n_cells - mean_density^2))

    push!(diag.t, sim.t)
    push!(diag.kinetic_energy, kinetic * dV)
    push!(diag.density_rms, rho_rms)
    push!(diag.total_mass, mass * dV)
    push!(diag.total_momentum_x, momentum_x * dV)
    push!(diag.total_momentum_y, momentum_y * dV)
    push!(diag.total_momentum_z, momentum_z * dV)
    push!(diag.total_energy, energy * dV)
    push!(diag.max_local_mach, max_mach)

    return diag
end


"""
    save_diagnostics(diag, path)

Save diagnostic histories to a JLD2 file.
"""
function save_diagnostics(
    diag::Diagnostics,
    path::String,
)
    jldsave(
        path;
        t = diag.t,
        kinetic_energy = diag.kinetic_energy,
        density_rms = diag.density_rms,
        total_mass = diag.total_mass,
        total_momentum_x = diag.total_momentum_x,
        total_momentum_y = diag.total_momentum_y,
        total_momentum_z = diag.total_momentum_z,
        total_energy = diag.total_energy,
        max_local_mach = diag.max_local_mach,
    )

    return nothing
end


"""
    load_diagnostics(path)

Load diagnostic histories from a JLD2 file.
"""
function load_diagnostics(
    path::String,
)
    d = load(path)

    return Diagnostics(
        d["t"],
        d["kinetic_energy"],
        d["density_rms"],
        d["total_mass"],
        d["total_momentum_x"],
        d["total_momentum_y"],
        d["total_momentum_z"],
        d["total_energy"],
        d["max_local_mach"],
    )
end
# -------------------------------------------------------------------------
# Simulation checkpoint
# -------------------------------------------------------------------------



"""
    save_checkpoint(sim::Simulation, path::String)

Full-state checkpoint for crash recovery. Unlike VTK export, this must
be exact/lossless — stores conserved variables directly (not primitives),
plus current simulation time, so a run can resume exactly where it left off.
"""
function save_checkpoint(sim::Simulation, path::String)
    jldsave(path;
        rho=sim.state.rho, rhou=sim.state.rhou, rhov=sim.state.rhov,
        rhow=sim.state.rhow, rhoE=sim.state.rhoE, t=sim.t)
    return nothing
end

"""
    load_checkpoint!(sim::Simulation, path::String)

Restore state and time from a checkpoint, in place.
"""
function load_checkpoint!(sim::Simulation, path::String)
    d = load(path)
    sim.state.rho  .= d["rho"]
    sim.state.rhou .= d["rhou"]
    sim.state.rhov .= d["rhov"]
    sim.state.rhow .= d["rhow"]
    sim.state.rhoE .= d["rhoE"]
    sim.t = d["t"]
    return sim
end