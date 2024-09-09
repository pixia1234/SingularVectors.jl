# Construction of Lie algebras

"""
    commutative_liealgebra(::Type{T}, n::Int) where T<:Number

Construct a commutative Lie algebra of dimension `n`.
"""
function commutative_liealgebra(::Type{T}, n::Int) where T<:Number
    scmat = [spzeros(T, n) for _ in 1:n, _ in 1:n]
    return AlgebraBySC(SCMat(scmat))
end
commutative_liealgebra(n::Int) = commutative_liealgebra(Int, n)

"""Structure constants of sl2 over Integer fields"""
const sl2scmat = SCMat([spzeros(Int, 3) for _ in 1:3, _ in 1:3])
sl2scmat[1, 2][1], sl2scmat[1, 3][2] = -2, 1
sl2scmat[2, 1][1], sl2scmat[2, 3][3] = 2, -2
sl2scmat[3, 1][2], sl2scmat[3, 2][3] = -1, 2
sl2scmat.syms .= @variables e h f

raw"""
Lie algebra $sl_2$ with chevalley basis $\{e, h, f\}$.
"""
const sl2 = AlgebraBySC(sl2scmat)

# Construction of Lie superalgebras

# Construction of Lie superalgebras

"""
    lie_superalgebra(::Type{T}, n_even::Int, n_odd::Int) where T<:Number

Construct a Lie superalgebra of dimension `n_even + n_odd`, where `n_even` is the dimension of the even part and `n_odd` is the dimension of the odd part.
"""
function lie_superalgebra(::Type{T}, n_even::Int, n_odd::Int) where T<:Number
    # Even-even, even-odd, and odd-odd structure constants
    scmat_even = [spzeros(T, n_even) for _ in 1:n_even, _ in 1:n_even]
    scmat_odd = [spzeros(T, n_odd) for _ in 1:n_odd, _ in 1:n_odd]
    scmat_mixed = [spzeros(T, n_even, n_odd) for _ in 1:n_even, _ in 1:n_odd]
    
    # Combine the structure constants in a superalgebra format
    return SuperAlgebraBySC(SCMat(scmat_even), SCMat(scmat_odd), SCMat(scmat_mixed))
end
lie_superalgebra(n_even::Int, n_odd::Int) = lie_superalgebra(Int, n_even, n_odd)

"""Structure constants of a basic Lie superalgebra over Integer fields"""
const sl2_scmat_even = SCMat([spzeros(Int, 3) for _ in 1:3, _ in 1:3])
const sl2_scmat_odd = SCMat([spzeros(Int, 2) for _ in 1:2, _ in 1:2])
const sl2_scmat_mixed = SCMat([spzeros(Int, 3, 2) for _ in 1:3, _ in 1:2])

# Define the structure constants (these would depend on the specific superalgebra)
sl2_scmat_even[1, 2][1], sl2_scmat_even[1, 3][2] = -2, 1
sl2_scmat_even[2, 1][1], sl2_scmat_even[2, 3][3] = 2, -2
sl2_scmat_even[3, 1][2], sl2_scmat_even[3, 2][3] = -1, 2

sl2_scmat_odd[1, 2][1] = 1  # Odd-odd commutator

# Define the mixed commutators between even and odd parts
sl2_scmat_mixed[1, 1][1] = 1

# Define variables for the basis elements
sl2_scmat_even.syms .= @variables e h f
sl2_scmat_odd.syms .= @variables ψ θ

raw"""
Lie superalgebra $sl(2|1)$ with even basis $\{e, h, f\}$ and odd basis $\{\psi, \theta\}$.
"""
const sl2_super = SuperAlgebraBySC(sl2_scmat_even, sl2_scmat_odd, sl2_scmat_mixed)
