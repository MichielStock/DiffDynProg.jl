module DiffDynProg

export MaxOperator, Max, LeakyMax, EntropyMax, SquaredMax, max_argmax, min_argmin

export DTW, dynamic_time_warping, ∂DPW, getD, getE


include("dynamictimewarping.jl")
include("maxoperators.jl")
include("utils.jl")

end