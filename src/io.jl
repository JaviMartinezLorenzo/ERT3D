"""
IO utilities for ERT3D.
"""

using WriteVTK


# -------------------------------------------------------------------------
# Internal helper
# -------------------------------------------------------------------------

function build_vtk(
    state::State{T},
    grid::Grid{T},
    params::Parameters{T},
    D::DerivativeOperator,
    path::String,
) where {T<:AbstractFloat}

    primitive = primitive_variables(state, params)
    X, Y, Z = meshgrid(grid)

    vtk = vtk_grid(path, X, Y, Z)
    vtk["density"] = primitive.rho
    vtk["pressure"] = primitive.p
    vtk["velocity"] = (primitive.u, primitive.v, primitive.w)
    vtk["velocity magnitude"] = sqrt.(primitive.u.^2 .+ primitive.v.^2 .+ primitive.w.^2)
    vtk["Q"] = q_criterion(state, params, D, grid)

    return vtk
end


# -------------------------------------------------------------------------
# Export a single snapshot
# -------------------------------------------------------------------------

function export_vtk(
    state::State{T},
    grid::Grid{T},
    params::Parameters{T},
    D::DerivativeOperator,
    path::String,
) where {T<:AbstractFloat}

    vtk = build_vtk(state, grid, params, D, path)
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

function VTKCollection(basepath::String)
    return VTKCollection(paraview_collection(basepath), basepath, 0)
end

function add_snapshot!(
    collection::VTKCollection,
    state::State{T},
    grid::Grid{T},
    params::Parameters{T},
    D::DerivativeOperator,
    time::Float64,
) where {T<:AbstractFloat}

    collection.counter += 1
    filename = "$(collection.basepath)_$(collection.counter)"

    vtk = build_vtk(state, grid, params, D, filename)
    collection.pvd[time] = vtk
    vtk_save(vtk)

    return nothing
end

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

All diagnostic histories use Float64, independently of the simulation
precision.
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

Compute and store scalar diagnostics for the current simulation state,
using Float64 arithmetic.
"""
function record!(diag::Diagnostics, sim::Simulation)
    state = sim.state
    gamma = Float64(sim.params.gamma)
    dV = Float64(sim.grid.dx)^3

    mass = 0.0
    momentum_x = 0.0
    momentum_y = 0.0
    momentum_z = 0.0
    energy = 0.0
    kinetic = 0.0
    density_sum = 0.0
    density_squared_sum = 0.0
    max_mach = 0.0

    @inbounds for i in eachindex(state.rho)
        rho = Float64(state.rho[i])
        rhou = Float64(state.rhou[i])
        rhov = Float64(state.rhov[i])
        rhow = Float64(state.rhow[i])
        rhoE = Float64(state.rhoE[i])

        u = rhou / rho
        v = rhov / rho
        w = rhow / rho

        velocity_squared = u^2 + v^2 + w^2

        pressure = (gamma - 1.0) * (
            rhoE - 0.5 * rho * velocity_squared
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

    n_cells = length(state.rho)
    mean_density = density_sum / n_cells
    rho_rms = sqrt(max(
        0.0,
        density_squared_sum / n_cells - mean_density^2,
    ))

    push!(diag.t, Float64(sim.t))
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
function save_diagnostics(diag::Diagnostics, path::String)
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
function load_diagnostics(path::String)
    d = load(path)

    return Diagnostics(
        Float64.(d["t"]),
        Float64.(d["kinetic_energy"]),
        Float64.(d["density_rms"]),
        Float64.(d["total_mass"]),
        Float64.(d["total_momentum_x"]),
        Float64.(d["total_momentum_y"]),
        Float64.(d["total_momentum_z"]),
        Float64.(d["total_energy"]),
        Float64.(d["max_local_mach"]),
    )
end


# -------------------------------------------------------------------------
# Simulation checkpoint
# -------------------------------------------------------------------------

"""
    save_checkpoint(sim::Simulation, path::String)

Save conserved variables and simulation time for restarting a run.
"""
function save_checkpoint(sim::Simulation, path::String)
    jldsave(
        path;
        rho = sim.state.rho,
        rhou = sim.state.rhou,
        rhov = sim.state.rhov,
        rhow = sim.state.rhow,
        rhoE = sim.state.rhoE,
        t = Float64(sim.t),
    )

    return nothing
end


"""
    load_checkpoint!(sim::Simulation, path::String)

Restore state and time from a checkpoint in place.
"""
function load_checkpoint!(sim::Simulation, path::String)
    d = load(path)

    sim.state.rho  .= d["rho"]
    sim.state.rhou .= d["rhou"]
    sim.state.rhov .= d["rhov"]
    sim.state.rhow .= d["rhow"]
    sim.state.rhoE .= d["rhoE"]
    sim.t = Float64(d["t"])

    return sim
end