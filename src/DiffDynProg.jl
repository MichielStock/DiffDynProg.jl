module DiffDynProg

export MaxOperator, Max, LeakyMax, EntropyMax, SquaredMax, max_argmax, min_argmin
export gap_cost_matrix, gumbel_softmax
export DP, getD, getE, getQ
export dynamic_time_warping, ∂DTW
export needleman_wunsch, ∂NW

using ChainRulesCore
import ChainRulesCore: rrule

include("maxoperators.jl")
include("dynprog.jl")
include("utils.jl")
include("dynamictimewarping.jl")
include("needlemanwunsch.jl")

end