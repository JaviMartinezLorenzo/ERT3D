abstract type InitialCondition end

struct TaylorGreen <: InitialCondition end

struct SyntheticTurbulence <: InitialCondition
    # parameters specific to this IC, if needed
end