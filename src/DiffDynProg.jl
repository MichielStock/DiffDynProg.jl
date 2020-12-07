module DiffDynProg

export MaxOperator, Max, LeakyMax, EntropyMax, SquaredMax

#include("dynamictimewarping.jl")
include("maxoperators.jl")
include("utils.jl")

end