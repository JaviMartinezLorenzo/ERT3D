"""
Parameters

Container for the physical and simulation parameters.
"""
struct Parameters
    gamma::Float64
    Mt0::Float64
    k0::Float64
end


Parameters(gamma::Float64, Mt0::Float64) = Parameters(gamma, Mt0, 6.0)