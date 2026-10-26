
"""
    Parameters(; gamma, Mt0, k0=6.0, T=Float64)

Container for the physical and simulation parameters.
Defaults to Float64 precision.
"""
struct Parameters{T<:AbstractFloat}
    gamma::T
    Mt0::T
    k0::T
end

function Parameters(;
    gamma,
    Mt0,
    k0=6.0,
    T::Type{<:AbstractFloat}=Float64,
)
    return Parameters{T}(
        T(gamma),
        T(Mt0),
        T(k0),
    )
end