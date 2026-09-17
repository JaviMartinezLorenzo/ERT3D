abstract type InitialCondition end

struct TaylorGreen <: InitialCondition end

struct SyntheticTurbulence <: InitialCondition
    seed::Int
end

SyntheticTurbulence() = SyntheticTurbulence(12345)