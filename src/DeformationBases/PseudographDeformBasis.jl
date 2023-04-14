"""
Concrete subtype of [`DeformBasis`](@ref).
Each element of the basis is induced by pseudograph with two vertices and
certain properties, which gets transformed to an arc diagram and then handled as
in [`ArcDiagDeformBasis`](@ref).
This process is due to [FM22](@cite).
"""
struct PseudographDeformBasis{C <: RingElement} <: DeformBasis{C}
    len::Int
    iter
    extra_data::Dict{DeformationMap{C}, Set{Tuple{Pseudograph2, Generic.Partition{Int}}}}
    normalize

    function PseudographDeformBasis{C}(
        sp::SmashProductLie{C},
        degs::AbstractVector{Int};
        no_normalize::Bool=false,
    ) where {C <: RingElement}
        get_attribute(sp.L, :type, nothing) == :special_orthogonal || error("Only works for so_n.")
        is_exterior_power(sp.V) && is_standard_module(get_attribute(sp.V, :inner_module)) ||
            error("Only works for exterior powers of the standard module.")

        e = get_attribute(sp.V, :power)

        extra_data = Dict{DeformationMap{C}, Set{Tuple{Pseudograph2, Generic.Partition{Int}}}}()
        normalize = no_normalize ? identity : normalize_default

        lens = []
        iters = []
        debug_counter = 0
        for d in degs
            pg_iter = pseudographs_with_partitions__so_extpowers_stdmod(e, d)
            len = length(pg_iter)
            iter = (
                begin
                    @debug "Basis generation deg $(d), $(debug_counter = (debug_counter % len) + 1)/$(len), $(floor(Int, 100*debug_counter / len))%"
                    diag = to_arcdiag(pg, part)
                    basis_elem = arcdiag_to_deformationmap__so(diag, sp)
                    if !no_normalize
                        basis_elem = normalize(basis_elem)
                    end
                    if haskey(extra_data, basis_elem)
                        push!(extra_data[basis_elem], (pg, part))
                    else
                        extra_data[basis_elem] = Set([(pg, part)])
                    end
                    basis_elem
                end for (pg, part) in pg_iter
            )
            push!(lens, len)
            push!(iters, iter)
        end
        len = sum(lens)
        iter = Iterators.flatten(iters)
        if !no_normalize
            iter = unique(Iterators.filter(b -> !iszero(b), iter))
            len = length(iter)
        end
        return new{C}(len, iter, extra_data, normalize)
    end
end

function Base.iterate(i::PseudographDeformBasis)
    return iterate(i.iter)
end

function Base.iterate(i::PseudographDeformBasis, s)
    return iterate(i.iter, s)
end

Base.length(basis::PseudographDeformBasis) = basis.len


function pseudographs_with_partitions__so_extpowers_stdmod(reg::Int, sumtotal::Int)
    iter = (
        begin
            (pg, Partition(copy(part)))
        end for sumpg in 0:sumtotal for pg in all_pseudographs(reg, sumpg; upto_iso=true) for
        part in AllParts(sumtotal - sumpg) if all(iseven, part) &&
        all(isodd, pg.loops1) &&
        all(isodd, pg.loops2) &&
        (pg.loops1 != pg.loops2 || isodd(sum(pg.edges)))
    )
    return collect(iter)
end
