module DiffDynProg

export MaxOperator, Max, LeakyMax, EntropyMax, SquaredMax, max_argmax, min_argmin
export gap_cost_matrix, gumbel_softmax
export DTW, dynamic_time_warping, ∂DPW, getD, getE

using ChainRulesCore
import ChainRulesCore: rrule

include("maxoperators.jl")
include("utils.jl")
include("dynamictimewarping.jl")

end