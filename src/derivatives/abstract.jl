abstract type DerivativeOperator end
abstract type CentralDifference <: DerivativeOperator end

struct Central4{T<:AbstractFloat} <: CentralDifference
    coeffs::NTuple{2,T}
    xwork::Array{T,3}
    ywork::Array{T,3}
    zwork::Array{T,3}
end

struct Central6{T<:AbstractFloat} <: CentralDifference
    coeffs::NTuple{3,T}
    xwork::Array{T,3}
    ywork::Array{T,3}
    zwork::Array{T,3}
end

struct Central8{T<:AbstractFloat} <: CentralDifference
    coeffs::NTuple{4,T}
    xwork::Array{T,3}
    ywork::Array{T,3}
    zwork::Array{T,3}
end

struct Spectral <: DerivativeOperator end