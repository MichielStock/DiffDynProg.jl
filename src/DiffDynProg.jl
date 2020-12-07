module DiffDynProg

export MaxOperator, Max, LeakyMax, EntropyMax, SquaredMax, max_argmax, min_argmin

#include("dynamictimewarping.jl")
include("maxoperators.jl")
include("utils.jl")

end