"""
    DerivativeOperator

Abstract supertype for spatial derivative operators.
"""
abstract type DerivativeOperator end

"""
    CentralDifference

Abstract supertype for centered finite-difference operators.
"""
abstract type CentralDifference <: DerivativeOperator end

"""
    Central4

Fourth-order central finite-difference operator.
"""
struct Central4 <: CentralDifference
    coeffs::NTuple{2,Float64}
    xwork::Array{Float64,3}
    ywork::Array{Float64,3}
    zwork::Array{Float64,3}
end

"""
    Central6

Sixth-order central finite-difference operator.
"""
struct Central6 <: CentralDifference
    coeffs::NTuple{3,Float64}
    xwork::Array{Float64,3}
    ywork::Array{Float64,3}
    zwork::Array{Float64,3}
end

"""
    Central8

Eighth-order central finite-difference operator.
"""
struct Central8 <: CentralDifference
    coeffs::NTuple{4,Float64}
    xwork::Array{Float64,3}
    ywork::Array{Float64,3}
    zwork::Array{Float64,3}
end

"""
    Spectral

Fourier spectral derivative operator for periodic domains.
"""
struct Spectral <: DerivativeOperator end