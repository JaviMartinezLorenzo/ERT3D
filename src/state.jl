"""
    State

Container for the conserved variables (ρ, ρu, ρv, ρw, ρE) on a 3D grid.

Design notes:
- Keep this a plain struct of arrays (not array of structs) for FFT/vectorization friendliness.
- Should support a `copy(state)` and `state .= other` style interface for
  the forward/reverse bookkeeping in run_experiment().
- Consider whether primitive variables (ρ, u, p) are computed on-the-fly
  via accessor functions, or cached — decide once profiling shows which
  is the bottleneck.
"""
struct State
    # TODO: ρ::Array{Float64,3}
    # TODO: ρu::Array{Float64,3}, ρv::Array{Float64,3}, ρw::Array{Float64,3}
    # TODO: ρE::Array{Float64,3}
end

# TODO: primitive_variables(state, grid) -> (ρ, u, v, w, p, T)  [uses EOS]
# TODO: reverse_velocity!(state)  -> flips ρu, ρv, ρw in place (the u → -u operator)
