# Operations of algebraic structures

"""
    ObjMatchError <: Exception

Operations on objects of different algebras.
"""
struct ObjMatchError <: Exception
    msg::String
end

# Basic operations
import Base: +, -, *, ÷, show, getindex, ==, iszero, keys, zero

"""
    getindex(scmat::AbstractSCMat, i::Int, j::Int)

Return the entry at position (i, j) in `scmat.mat`.
"""
getindex(scmat::AbstractSCMat, i::Int, j::Int) = scmat.mat[i, j]

"""
    sparse2dict(x::SparseVector)

Convert a sparse vector to a dictionary.
"""
function sparse2dict(x::AbstractSparseVector{T}) where T<:Number
    return Dict{Tuple, T}((i,)=>x[i] for i in findall(!iszero, x))
end

## Operations for AbstractElem
## as vector space
"""
    +(x::AbstractElem, y::AbstractElem)

Addition of two elements of the Lie algebra.
"""
function +(x::AbstractElem{T}, y::AbstractElem{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Addition of elements of different algebras"))
    return typeof(x)(x.scmat, x.element + y.element)
end

"""
    -(x::AbstractElem, y::AbstractElem)

Subtraction of two elements of the Lie algebra.
"""
function -(x::AbstractElem{T}, y::AbstractElem{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Subtraction of elements of different algebras"))
    return typeof(x)(x.scmat, x.element - y.element)
end

"""
    -(x::AbstractElem)

Negation of an element of the Lie algebra.
"""
-(x::AbstractElem{T}) where T<:Number = typeof(x)(x.scmat, -x.element)

"""
    *(a::Number, x::AbstractElem)

Multiplication of an element of the Lie algebra by an integer.
"""
function *(a::T, x::AbstractElem{T}) where T<:Number
    return typeof(x)(x.scmat, a * x.element)
end
*(x::AbstractElem, a::T) where T<:Number = a * x

raw"""
    getindex(x::AbstractElem, i::Int)

Return the coefficient of $x_i$ in `x`.
"""
getindex(x::AbstractElem, i::Int) = x.element[i]

## as algebra
"""
    *(x::LieElem, y::LieElem)

Multiplication of two elements of the Lie algebra.
"""
function *(x::AbstractElem{T}, y::AbstractElem{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Multiplication of elements of different algebras"))
    scmat = x.scmat
    z = spzeros(T, scmat.d)
    for i in x.element.nzind, j in y.element.nzind
        val = x[i] * y[j]
        iszero(val) && continue
        z += val * scmat[i, j]
    end
    return typeof(x)(scmat, z)
end

==(x::AbstractElem, y::AbstractElem) = iszero(x - y)
iszero(x::AbstractElem) = iszero(x.element)

## Operations for LieEnvElem
"""
    keys(x::AbstractEnvElem)

Return the decomposition basis of `x`.
"""
keys(x::AbstractEnvElem) = keys(x.element)

"""
    getindex(x::AbstractEnvElem, key::Tuple)

Return the coefficient of `key` in `x`.
"""
getindex(x::AbstractEnvElem{T}, key::Tuple) where T = get(x.element, key, zero(T))

"""
    zero(x::AbstractEnvElem)

Return the zero element of the same algebra as `x`.
"""
zero(x::AbstractEnvElem{T}) where T<:Number = typeof(x)(x.scmat, Dict{Tuple, T}())

"""
    unit(x::AbstractEnvElem)

Return the identity element of the same algebra as `x`.
"""
unit(x::AbstractEnvElem{T}) where T<:Number = typeof(x)(x.scmat, Dict{Tuple, T}(()=>one(T)))

"""
    +(x::AbstractEnvElem, y::AbstractEnvElem)

Addition of two elements of the universal enveloping algebra.
"""
function +(x::AbstractEnvElem{T}, y::AbstractEnvElem{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Addition of elements of different enveloping algebras"))
    ele = Dict{Tuple, T}()
    for k in union(keys(x), keys(y))
        val = x[k] + y[k]
        if !iszero(val)
            ele[k] = val
        end
    end
    return typeof(x)(x.scmat, ele)
end

"""Addition without creating new memory"""
function add!(x::AbstractEnvElem, y::AbstractEnvElem)
    for k in keys(y)
        val = x[k] + y[k]
        if iszero(val)
            delete!(x.element, k)
        else
            x.element[k] = val
        end
    end
    return x
end

"""
    -(x::AbstractEnvElem, y::AbstractEnvElem)

Subtraction of two elements of the universal enveloping algebra.
"""
function -(x::AbstractEnvElem{T}, y::AbstractEnvElem{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Subtraction of elements of different enveloping algebras"))
    ele = Dict{Tuple, T}()
    for k in union(keys(x), keys(y))
        val = x[k] - y[k]
        if !iszero(val)
            ele[k] = val
        end
    end
    return typeof(x)(x.scmat, ele)
end

"""
    -(x::AbstractEnvElem)

Negation of an element of the universal enveloping algebra.
"""
-(x::AbstractEnvElem{T}) where T<:Number = typeof(x)(x.scmat, Dict{Tuple, T}(k => -v for (k, v) in x.element))

"""
    *(a::Int, x::AbstractEnvElem)

Multiplication of an element of the universal enveloping algebra by an integer.
"""
function *(a::T, x::AbstractEnvElem{T}) where T<:Number
    return typeof(x)(x.scmat, Dict{Tuple, T}(k => a * v for (k, v) in x.element))
end
*(x::AbstractEnvElem{T}, a::T) where T<:Number = a * x

"""multiplication without creating new memory"""
function times!(x::AbstractEnvElem{T}, a::T) where T<:Number
    for k in keys(x.element)
        x.element[k] *= a
    end
    return x
end    

"""
    ==(x::AbstractEnvElem, y::AbstractEnvElem)

Test if two elements of the universal enveloping algebra are equal.
"""
==(x::AbstractEnvElem, y::AbstractEnvElem) = iszero(x - y)

"""
    iszero(x::AbstractEnvElem)

Test if an element of the universal enveloping algebra is zero.
"""
iszero(x::AbstractEnvElem) = all(iszero, values(x.element))

## **PBW basis**

# @inline tuplejoin(x) = x
@inline tuplejoin(x, y) = (x..., y...)
# @inline tuplejoin(x, y, z...) = tuplejoin(tuplejoin(x, y), z...)

"""
    mult(scmat::AbstractSCMat, x::Tuple, y::Tuple)

Simplify the product of two basis elements of the universal enveloping algebra.
"""
mult(scmat::AbstractSCMat, x::Tuple, y::Tuple) = simplify(scmat, tuplejoin(x, y))

"""
    simplify(scmat::SCMat, x::Tuple)

Reduce a sequence of indexes to the standard `LieEnvElem`.

Basic rule: xᵢxⱼ = xⱼxᵢ + [xᵢ, xⱼ]

Example: (1, 3, 2, 1) => (1, 2, 3, 1) + (1, [3, 2], 1)
"""
function simplify(scmat::SCMat{T}, x::Tuple) where T<:Number
    ind = findfirst(i -> x[i] > x[i+1], 1:length(x)-1)
    isnothing(ind) && return LieEnvElem(scmat, Dict{Tuple, T}(x=>one(T)))
    # decrease the inversion number by 1
    reducecase = simplify(scmat, (x[1:ind-1]..., x[ind+1], x[ind], x[ind+2:end]...))
    # decrease the length by 1
    ele = scmat[x[ind], x[ind+1]]
    for i in ele.nzind
        add!(reducecase, times!(simplify(scmat, (x[1:ind-1]..., i, x[ind+2:end]...)), ele[i]))
    end
    return reducecase
end

"""
    *(x::AbstractEnvElem, y::AbstractEnvElem)

Multiplication of two elements of the universal enveloping algebra.
"""
function *(x::AbstractEnvElem{T}, y::AbstractEnvElem{T}) where T<:Number
    x.scmat === y.scmat || throw(ObjMatchError("Multiplication of elements of different enveloping algebras"))
    ele = zero(x)
    for k1 in keys(x), k2 in keys(y)
        coef = x[k1] * y[k2]
        if !iszero(coef)
            ele += coef * mult(x.scmat, k1, k2)
        end
    end
    return ele
end

"""
    issortedbypbw(x::AbstractEnvElem)

Test if the key values of the `AbstractEnvElem` is sorted by PBW basis.
"""
issortedbypbw(x::AbstractEnvElem) = all(issorted, keys(x.element))

# Treat Lie algebra as a subalgebra of the Enveloping algebra.
*(x::AbstractElem{T}, y::AbstractEnvElem{T}) where T<:Number = oftype(y, x) * y
*(x::AbstractEnvElem{T}, y::AbstractElem{T}) where T<:Number = x * oftype(x, y)
+(x::AbstractEnvElem{T}, y::AbstractElem{T}) where T<:Number = x + oftype(x, y)
+(x::AbstractElem{T}, y::AbstractEnvElem{T}) where T<:Number = oftype(y, x) + y
-(x::AbstractEnvElem{T}, y::AbstractElem{T}) where T<:Number = x - oftype(x, y)
-(x::AbstractElem{T}, y::AbstractEnvElem{T}) where T<:Number = oftype(y, x) - y

# element type
eltype(::AbstractEnvElem{T}) where T = T
eltype(::AbstractElem{T}) where T = T
eltype(::AbstractSCMat{T}) where T = T
eltype(::AlgebraBySC{T}) where T = T

# dimension
dim(x::AbstractEnvElem) = x.scmat.d
dim(x::AbstractElem) = x.scmat.d
dim(x::AbstractSCMat) = x.d
dim(x::AlgebraBySC) = x.scmat.d

# should be replaced by Symbolic.jl
function show(io::IO, alg::AlgebraBySC)
    txt = "Lie algebra of dimension $(dim(alg)) with basis: " * join(string.(alg.basis), ", ")
    print(io,  txt)
end

show(io::IO, x::AbstractElem) = print(io, sum(x.element .* x.scmat.syms))

show(io::IO, scmat::AbstractSCMat) = print(io, "Structure constants with dimension $(dim(scmat))")

"""
    show(io::IO, x::LieEnvElem)

Show an element of the universal enveloping algebra.

Note: the output is sorted by dictionary order, not by PBW basis. However, the 
user should know that it represents the element defined by PBW ordering.

Also, one can use `issortedbypbw` to test if the basis is a standard PBW basis.
"""
function show(io::IO, x::AbstractEnvElem)
    # warning: not standard ordering !
    dict, syms = x.element, x.scmat.syms
    key2sym(key) = prod(syms[i] for i in key)
    print(io, sum(val * key2sym(k) for (k, val) in dict))
end

function +(x::LieElem{T}, y::LieElem{T}) where T <: Number
    return LieElem(x.scmat, x.vector + y.vector)
end

function *(a::T, x::LieElem{T}) where T <: Number
    return LieElem(x.scmat, a * x.vector)
end

function is_closed_under_lie_bracket(algebra::AlgebraBySC{T}) where T <: Number
    for i in 1:length(algebra.basis)
        for j in i+1:length(algebra.basis)
            result = lie_bracket(algebra, algebra.basis[i], algebra.basis[j])
            if !all(isapprox(result.vector[k], 0.0; atol=1e-10) for k in 1:length(result.vector))
                return false
            end
        end
    end
    return true
end

function center_of_algebra(algebra::AlgebraBySC{T}) where T <: Number
    d = length(algebra.basis)
    center_elements = []
    for i in 1:d
        elem = algebra.basis[i]
        is_central = true
        for j in 1:d
            if !isapprox(norm(lie_bracket(algebra, elem, algebra.basis[j]).vector), 0.0; atol=1e-10)
                is_central = false
                break
            end
        end
        if is_central
            push!(center_elements, elem)
        end
    end
    return center_elements
end

function dim(algebra::AlgebraBySC{T}) where T <: Number
    return length(algebra.basis)
end

function representation(algebra::AlgebraBySC{T}, x::LieElem{T}) where T <: Number
    d = dim(algebra.scmat)
    rep_matrix = zeros(T, d, d)
    for i in 1:d
        for j in 1:d
            rep_matrix[i, j] = lie_bracket(algebra, x, algebra.basis[j]).vector[i]
        end
    end
    return rep_matrix
end

function cartesian_product(algebra1::AlgebraBySC{T}, algebra2::AlgebraBySC{T}) where T <: Number
    new_basis = vcat(algebra1.basis, algebra2.basis)
    new_scmat = blockdiag(algebra1.scmat, algebra2.scmat)
    return AlgebraBySC(new_scmat, new_basis)
end

function derived_algebra(algebra::AlgebraBySC{T}) where T <: Number
    derived_basis = []
    for i in 1:length(algebra.basis)
        for j in i+1:length(algebra.basis)
            bracket = lie_bracket(algebra, algebra.basis[i], algebra.basis[j])
            push!(derived_basis, bracket)
        end
    end
    return AlgebraBySC(algebra.scmat, derived_basis)
end

function *(x::UEAElem{T}, y::UEAElem{T}, uea::UEA{T}) where T <: Number
    terms = Dict{Tuple{Int, Int}, T}()
    for ((i1, j1), v1) in x.terms
        for ((i2, j2), v2) in y.terms
            # 计算乘积 (i1 * i2, j1 * j2)
            new_term = (i1, j2)
            if new_term in terms
                terms[new_term] += v1 * v2
            else
                terms[new_term] = v1 * v2
            end
            # 添加反向乘积，依据李代数关系
            new_term_reverse = (i2, j1)
            if new_term_reverse in terms
                terms[new_term_reverse] -= v1 * v2
            else
                terms[new_term_reverse] = -v1 * v2
            end
        end
    end
    return UEAElem{T}(terms)
end
