module SingularVectors

# Data types
export SCMat, LieElement, EnvElement, LieAlgebra

# Constructions
export commutative_liealgebra, sl2

include("universalenv.jl")
using .UniversalEnvelope

end
