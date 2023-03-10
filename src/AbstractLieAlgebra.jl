@attributes mutable struct AbstractLieAlgebra{C <: RingElement} <: LieAlgebra{C}
    R::Ring
    dim::Int
    struct_consts::Matrix{SRow{C}}
    s::Vector{Symbol}

    function AbstractLieAlgebra{C}(
        R::Ring,
        struct_consts::Matrix{SRow{C}},
        s::Vector{Symbol},
        cached::Bool=true;
        check::Bool=true,
    ) where {C <: RingElement}
        return get_cached!(AbstractLieAlgebraDict, (R, struct_consts, s), cached) do
            (n1, n2) = size(struct_consts)
            n1 == n2 || error("Invalid structure constants dimensions.")
            dim = n1
            length(s) == dim || error("Invalid number of basis element names.")
            if check
                all(iszero, struct_consts[i, i][k] for i in 1:dim, k in 1:dim) ||
                    error("Not anti-symmetric.")
                all(iszero, struct_consts[i, j][k] + struct_consts[j, i][k] for i in 1:dim, j in 1:dim, k in 1:dim) ||
                    error("Not anti-symmetric.")
                all(
                    iszero,
                    sum(
                        struct_consts[i, j][k] * struct_consts[k, l][m] +
                        struct_consts[j, l][k] * struct_consts[k, i][m] +
                        struct_consts[l, i][k] * struct_consts[k, j][m] for k in 1:dim
                    ) for i in 1:dim, j in 1:dim, l in 1:dim, m in 1:dim
                ) || error("Jacobi identity does not hold.")
            end
            new{C}(R, dim, struct_consts, s)
        end::AbstractLieAlgebra{C}
    end
end

const AbstractLieAlgebraDict = CacheDictType{Tuple{Ring, Matrix{SRow}, Vector{Symbol}}, AbstractLieAlgebra}()

struct AbstractLieAlgebraElem{C <: RingElement} <: LieAlgebraElem{C}
    parent::AbstractLieAlgebra{C}
    mat::MatElem{C}
end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{AbstractLieAlgebraElem{C}}) where {C <: RingElement} = AbstractLieAlgebra{C}

elem_type(::Type{AbstractLieAlgebra{C}}) where {C <: RingElement} = AbstractLieAlgebraElem{C}

parent(x::AbstractLieAlgebraElem{C}) where {C <: RingElement} = x.parent

base_ring(L::AbstractLieAlgebra{C}) where {C <: RingElement} = L.R

ngens(L::AbstractLieAlgebra{C}) where {C <: RingElement} = L.dim

@inline function Generic._matrix(x::AbstractLieAlgebraElem{C}) where {C <: RingElement}
    return (x.mat)::dense_matrix_type(C)
end


###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::AbstractLieAlgebra{C}) where {C <: RingElement}
    print(io, "AbstractLieAlgebra over ")
    print(IOContext(io, :compact => true), base_ring(V))
end

function expressify(v::AbstractLieAlgebraElem{C}, s=symbols(parent(v)); context=nothing) where {C <: RingElement}
    sum = Expr(:call, :+)
    for (i, c) in enumerate(_matrix(v))
        push!(sum.args, Expr(:call, :*, expressify(c, context=context), s[i]))
    end
    return sum
end

@enable_all_show_via_expressify AbstractLieAlgebraElem

function symbols(L::AbstractLieAlgebra{C}) where {C <: RingElement}
    return L.s
end


###############################################################################
#
#   Parent object call overload
#
###############################################################################

# no special ones


###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function bracket(x::AbstractLieAlgebraElem{C}, y::AbstractLieAlgebraElem{C}) where {C <: RingElement}
    check_parent(x, y)
    L = parent(x)
    mat =
        sum(cxi * cyj * L.struct_consts[i, j] for (i, cxi) in enumerate(_matrix(x)), (j, cyj) in enumerate(_matrix(y)))
    return L(mat)
end


###############################################################################
#
#   Constructor
#
###############################################################################

function liealgebra(
    R::Ring,
    struct_consts::Matrix{SRow{C}},
    s::Vector{<:Union{AbstractString, Char, Symbol}},
    cached::Bool=true;
    check::Bool=true,
) where {C <: RingElement}
    return AbstractLieAlgebra{elem_type(R)}(R, struct_consts, Symbol.(s), cached; check)
end

function liealgebra(
    R::Ring,
    struct_consts::Array{C, 3},
    s::Vector{<:Union{AbstractString, Char, Symbol}},
    cached::Bool=true;
    check::Bool=true,
) where {C <: RingElement}
    struct_consts2 = Matrix{SRow{elem_type(R)}}(undef, size(struct_consts, 1), size(struct_consts, 2))
    for i in axes(struct_consts, 1), j in axes(struct_consts, 2)
        struct_consts2[i, j] = sparse_row(R, collect(axes(struct_consts, 3)), struct_consts[i, j, :])
    end

    return AbstractLieAlgebra{elem_type(R)}(R, struct_consts2, Symbol.(s), cached; check)
end
