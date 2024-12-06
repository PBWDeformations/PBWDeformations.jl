const SO = Val{:special_orthogonal}
const GL = Val{:general_linear}

"""
Concrete subtype of [`DeformBasis`](@ref).
Each element of the basis is induced by an arc diagram of a suitable size,
which gets symmetrized and specialised to the given smash product.
This process is due to [FM22](@cite).
"""
struct ArcDiagDeformBasis{T <: SmashProductLieElem} <: DeformBasis{T}
    len::Int
    iter
    extra_data::Dict{DeformationMap{T}, Set{ArcDiagram}}
    no_normalize::Bool

    function ArcDiagDeformBasis(
        sp::SmashProductLie{C, LieC, LieT},
        degs::AbstractVector{Int};
        no_normalize::Bool=false,
        check_all_diagrams::Bool=false,
    ) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
        LieType = get_attribute(base_lie_algebra(sp), :type, nothing)::Union{Nothing,Symbol}
        @req LieType in [:special_orthogonal, :general_linear] "Only works for so_n and gl_n."
        if LieType == :special_orthogonal && has_attribute(base_lie_algebra(sp), :form)
            @req isone(get_attribute(base_lie_algebra(sp), :form)::dense_matrix_type(C)) "Only works for so_n represented as skew-symmetric matrices."
        end
        return ArcDiagDeformBasis(Val(LieType), sp, degs; no_normalize, check_all_diagrams)
    end

    function ArcDiagDeformBasis(
        LieType::Union{SO, GL},
        sp::SmashProductLie{C, LieC, LieT},
        degs::AbstractVector{Int};
        no_normalize::Bool=false,
        check_all_diagrams::Bool=false,
    ) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
        V = base_module(sp)

        V_nice, h = isomorphic_module_with_simple_structure(V)
        fl, V_nice_summands = _is_direct_sum(V_nice)
        if !fl
            temp = direct_sum(V_nice)
            h = compose(h, hom(V_nice, temp, identity_matrix(coefficient_ring(temp), dim(temp))))
            V_nice_summands = [V_nice]
            V_nice = temp
        end

        extra_data = Dict{DeformationMap{elem_type(sp)}, Set{ArcDiagram}}()

        n_sum_cases = div(length(V_nice_summands) * (length(V_nice_summands) + 1), 2)
        lens = Int[]
        iters = []
        for d in degs
            sum_case = 0
            for (i_l, V_nice_summand_i_l) in enumerate(V_nice_summands),
                (i_r, V_nice_summand_i_r) in enumerate(V_nice_summands)

                if i_l > i_r
                    continue
                end

                sum_case += 1

                proj_to_summand_l = compose(h, canonical_projection(V_nice, i_l))
                proj_to_summand_r = compose(h, canonical_projection(V_nice, i_r))

                if i_l == i_r
                    W = exterior_power_obj(V_nice_summand_i_l, 2)
                    case = :exterior_power
                else
                    W = tensor_product(V_nice_summand_i_l, V_nice_summand_i_r)
                    case = :tensor_product
                end

                diag_iter = pbw_arc_diagrams(LieType, W, d; case, check_all_diagrams)
                len = length(diag_iter)
                prog_meter = ProgressMeter.Progress(len; output=stderr, enabled=true, desc="Basis generation: deg $d, case $(sum_case)/$(n_sum_cases)")
                generate_showvalues(counter, diag) = () -> [("iteration", (counter, diag))]
                iter = (
                    begin
                        # @vprintln :PBWDeformations 2 "Basis generation deg $(lpad(d, maximum(ndigits, degs))), case $(lpad(sum_case, ndigits(n_sum_cases)))/$(n_sum_cases), $(lpad(floor(Int, 100*counter / len), 3))%, $(lpad(counter, ndigits(len)))/$(len)"
                        _basis_elem = arcdiag_to_deformationmap(LieType, diag, sp, W, case)
                        basis_elem = matrix(proj_to_summand_l) * _basis_elem * transpose(matrix(proj_to_summand_r))
                        if i_l != i_r
                            basis_elem -= transpose(basis_elem)
                        end
                        @assert is_skew_symmetric(basis_elem)

                        if !no_normalize
                            basis_elem = normalize(basis_elem)
                        end
                        if haskey(extra_data, basis_elem)
                            push!(extra_data[basis_elem], diag)
                        else
                            extra_data[basis_elem] = Set([diag])
                        end
                        ProgressMeter.update!(prog_meter, counter; showvalues = generate_showvalues(counter, diag))
                        basis_elem
                    end for (counter, diag) in enumerate(diag_iter) if check_all_diagrams || is_in_canonical_form_for_deformation(diag; case)
                )
                # push!(lens, len)
                # push!(iters, iter)
                collected = collect(iter)
                push!(lens, length(collected))
                push!(iters, collected)
                ProgressMeter.finish!(prog_meter)
            end
        end
        len = sum(lens)
        iter = Iterators.flatten(iters)
        if !no_normalize
            iter = unique(Iterators.filter(b -> !iszero(b), iter))
            collected = Vector{DeformationMap{elem_type(sp)}}(collect(iter))::Vector{DeformationMap{elem_type(sp)}}
            _, rels = is_linearly_independent_with_relations(coefficient_ring(sp), collected)
            inds = [findlast(!iszero, vec(rels[i, :]))::Int for i in 1:nrows(rels)]
            deleteat!(collected, inds)
            return new{elem_type(sp)}(length(collected), collected, extra_data, no_normalize)
        end
        return new{elem_type(sp)}(len, iter, extra_data, no_normalize)
    end
end

function Base.iterate(i::ArcDiagDeformBasis)
    return iterate(i.iter)::Union{Tuple{eltype(i), Any}, Nothing}
end

function Base.iterate(i::ArcDiagDeformBasis, s)
    return iterate(i.iter, s)::Union{Tuple{eltype(i), Any}, Nothing}
end

Base.length(basis::ArcDiagDeformBasis) = basis.len


function pbw_arc_diagrams(T::Union{SO, GL}, V::LieAlgebraModule, d::Int; case::Symbol=:unknown, check_all_diagrams::Bool=false)
    upper_verts = arc_diagram_upper_points(T, V)
    lower_verts = arc_diagram_lower_points(T, V, d)
    upper_iss = arc_diagram_upper_iss(T, V)
    lower_iss = arc_diagram_lower_iss(T, V, d)
    indep_sets = Vector{Int}[[(-1) .* is for is in upper_iss]; [is for is in lower_iss]]
    if !check_all_diagrams
        ### if first arc goes from upper to lower, it goes to the first lower vertex
        for i in 2:arc_diagram_num_lower_points(T, V, d)
            push!(indep_sets, [-1, i])
        end
        if case == :exterior_power && d > 0
            ### first arc from upper right part may not got to the first lower vertex -> swap upper parts instead
            push!(indep_sets, [-(div(arc_diagram_num_upper_points(T, V), 2) + 1), 1])
        end
    end
    return all_arc_diagrams(arc_diagram_type(T), upper_verts, lower_verts; indep_sets)
end


arc_diagram_type(::SO) = Undirected

arc_diagram_type(::GL) = Directed


function is_tensor_generator(::SO, V::LieAlgebraModule)
    return _is_standard_module(V)
end

function is_tensor_generator(::GL, V::LieAlgebraModule)
    if _is_standard_module(V)
        return true
    end
    fl, base = _is_dual(V)
    return fl && _is_standard_module(base)
end


function arc_diagram_upper_points(T::SO, V::LieAlgebraModule)
    if _is_standard_module(V)
        return 1
    elseif ((fl, Ws) = _is_tensor_product(V); fl)
        return sum(arc_diagram_upper_points(T, W) for W in Ws)
    elseif ((fl, W, k) = is_power_with_data(V); fl)
        return arc_diagram_upper_points(T, W) * k
    else
        error("Not implemented.")
    end
end

function arc_diagram_upper_points(T::GL, V::LieAlgebraModule)
    if _is_standard_module(V)
        return 1
    elseif ((fl, W) = _is_dual(V); fl) && _is_standard_module(W)
        return 0
    elseif ((fl, Ws) = _is_tensor_product(V); fl)
        return reduce(vcat, arc_diagram_upper_points(T, W) for W in Ws)
    elseif ((fl, W, k) = is_power_with_data(V); fl)
        return reduce(vcat, [arc_diagram_upper_points(T, W) for _ in 1:k])
    else
        error("Not implemented.")
    end
end

function arc_diagram_num_upper_points(T::SO, V::LieAlgebraModule)
    return arc_diagram_upper_points(T, V)
end

function arc_diagram_num_upper_points(T::GL, V::LieAlgebraModule)
    return length(arc_diagram_upper_points(T, V))
end


function arc_diagram_upper_iss(T::Union{SO, GL}, V::LieAlgebraModule)
    if is_tensor_generator(T, V)
        return Vector{Int}[]
    elseif ((fl, inner_mods) = _is_tensor_product(V); fl)
        offset = 0
        iss = Vector{Int}[]
        for mod in inner_mods
            append!(iss, [is .+ offset for is in arc_diagram_upper_iss(T, mod)])
            offset += arc_diagram_num_upper_points(T, mod)
        end
        return iss
    elseif ((fl, inner_mod, power) = is_power_with_data(V); fl)
        if is_tensor_generator(T, inner_mod)
            if _is_exterior_power(V)[1]
                return [collect(1:power)]
            else
                return Vector{Int}[]
            end
        else
            iss = arc_diagram_upper_iss(T, inner_mod)
            return [is .+ k * arc_diagram_num_upper_points(T, inner_mod) for k in 0:power-1 for is in iss]
        end
    else
        error("Not implemented.")
    end
end


function arc_diagram_lower_points(::SO, _::LieAlgebraModule, d::Int)
    # L ≅ Sᵈ ⋀² V
    return 2d
end

function arc_diagram_lower_points(::GL, _::LieAlgebraModule, d::Int)
    # L ≅ Sᵈ (V ⊗ V*)
    return reduce(vcat, ([1, 0] for _ in 1:d); init=Int[])
end

function arc_diagram_num_lower_points(T::SO, V::LieAlgebraModule, d::Int)
    return arc_diagram_lower_points(T, V, d)
end

function arc_diagram_num_lower_points(T::GL, V::LieAlgebraModule, d::Int)
    return length(arc_diagram_lower_points(T, V, d))
end


function arc_diagram_lower_iss(::SO, _::LieAlgebraModule, d::Int)
    # L ≅ Sᵈ ⋀² V
    return collect([2i - 1, 2i] for i in 1:d)
end

function arc_diagram_lower_iss(::GL, _::LieAlgebraModule, _::Int)
    # L ≅ Sᵈ (V ⊗ V*)
    return Vector{Int}[]
end

function arc_diagram_lower_pair_to_L(::SO, dim_stdmod_V::Int)
    # L ≅ Sᵈ ⋀² V
    iso_pair_to_L = Dict{Tuple{Int, Int}, Int}()
    for (i, bs) in enumerate(combinations(dim_stdmod_V, 2) .|> NTuple{2})
        iso_pair_to_L[bs] = i
    end
    return function (k1::Int, k2::Int)
        if k1 == k2
            return 0, 0
        elseif k1 < k2
            return iso_pair_to_L[(k1, k2)], 1
        else
            return iso_pair_to_L[(k2, k1)], -1
        end
    end
end

function arc_diagram_lower_pair_to_L(::GL, dim_stdmod_V::Int)
    # L ≅ Sᵈ (V ⊗ V*)
    iso_pair_to_L = Dict{Tuple{Int, Int}, Int}()
    for (i, bs) in enumerate(ProductIterator(1:dim_stdmod_V, 2) .|> reverse .|> NTuple{2})
        iso_pair_to_L[bs] = i
    end
    return function (k1::Int, k2::Int)
        return iso_pair_to_L[(k1, k2)], 1
    end
    return iso_pair_to_L
end

function arc_diagram_label_iterator(T::Union{SO, GL}, V::LieAlgebraModule, base_labels::AbstractVector{Int})
    if is_tensor_generator(T, V)
        return [[l] for l in base_labels]
    elseif ((fl, inner_mods) = _is_tensor_product(V); fl)
        return ProductIterator([
                   arc_diagram_label_iterator(T, inner_mod, base_labels) for inner_mod in reverse(inner_mods)
               ]) .|>
               reverse .|>
               Iterators.flatten .|>
               collect
    elseif ((fl, inner_mod, power) = _is_exterior_power(V); fl)
        return combinations(collect(arc_diagram_label_iterator(T, inner_mod, base_labels)), power) .|>
               Iterators.flatten .|>
               collect
    elseif ((fl, inner_mod, power) = _is_symmetric_power(V); fl)
        return multicombinations(collect(arc_diagram_label_iterator(T, inner_mod, base_labels)), power) .|>
               Iterators.flatten .|>
               collect
    elseif ((fl, inner_mod, power) = _is_tensor_power(V); fl)
        return ProductIterator(arc_diagram_label_iterator(T, inner_mod, base_labels), power) .|>
               reverse .|>
               Iterators.flatten .|>
               collect
    else
        error("Not implemented.")
    end
end


function basis_index_mapping(V::LieAlgebraModule)
    if ((fl, inner_mods) = _is_tensor_product(V); fl)
        return ProductIterator([1:dim(inner_mod) for inner_mod in reverse(inner_mods)]) .|> reverse .|> collect
    elseif ((fl, inner_mod, power) = _is_exterior_power(V); fl)
        return combinations(dim(inner_mod), power) .|> collect
    elseif ((fl, inner_mod, power) = _is_symmetric_power(V); fl)
        return multicombinations(dim(inner_mod), power) .|> collect
    elseif ((fl, inner_mod, power) = _is_tensor_power(V); fl)
        return ProductIterator(1:dim(inner_mod), power) .|> reverse .|> collect
    else
        error("Not implemented.")
    end
end


function arc_diagram_label_permutations(T::Union{SO, GL}, V::LieAlgebraModule, label::AbstractVector{Int})
    if is_tensor_generator(T, V)
        @req length(label) == 1 "Number of labels mismatch."
        return [(label, 1)]
    elseif ((fl, inner_mods) = _is_tensor_product(V); fl)
        @req length(label) == sum(mod -> arc_diagram_num_upper_points(T, mod), inner_mods) "Number of labels mismatch."
        return [
            begin
                inner_label = reduce(vcat, first.(inner_iter))
                inner_sign = prod(last.(inner_iter))
                (inner_label, inner_sign)
            end for inner_iter in ProductIterator([
                arc_diagram_label_permutations(
                    T,
                    inner_mod,
                    label[sum(mod -> arc_diagram_num_upper_points(T, mod), inner_mods[1:i-1]; init=0)+1:sum(
                        mod -> arc_diagram_num_upper_points(T, mod),
                        inner_mods[1:i];
                        init=0,
                    )],
                ) for (i, inner_mod) in enumerate(inner_mods)
            ])
        ]
    elseif ((fl, inner_mod, power) = _is_exterior_power(V); fl)
        m = arc_diagram_num_upper_points(T, inner_mod)
        @req length(label) == m * power "Number of labels mismatch."
        return [
            begin
                inner_label = reduce(vcat, first.(inner_iter))
                inner_sign = prod(last.(inner_iter))
                (inner_label, inner_sign * outer_sign)
            end for (outer_perm, outer_sign) in permutations_with_sign(1:power) for inner_iter in ProductIterator([
                arc_diagram_label_permutations(T, inner_mod, label[(outer_perm[i]-1)*m+1:outer_perm[i]*m]) for
                i in 1:power
            ])
        ]
    elseif ((fl, inner_mod, power) = _is_symmetric_power(V); fl)
        m = arc_diagram_num_upper_points(T, inner_mod)
        @req length(label) == m * power "Number of labels mismatch."
        return [
            begin
                inner_label = reduce(vcat, first.(inner_iter))
                inner_sign = prod(last.(inner_iter))
                (inner_label, inner_sign)
            end for outer_perm in permutations(1:power) for inner_iter in ProductIterator([
                arc_diagram_label_permutations(T, inner_mod, label[(outer_perm[i]-1)*m+1:outer_perm[i]*m]) for
                i in 1:power
            ])
        ]
    elseif ((fl, inner_mod, power) = _is_tensor_power(V); fl)
        m = arc_diagram_num_upper_points(T, inner_mod)
        @req length(label) == m * power "Number of labels mismatch."
        return [
            begin
                inner_label = reduce(vcat, first.(inner_iter))
                inner_sign = prod(last.(inner_iter))
                (inner_label, inner_sign)
            end for inner_iter in
            ProductIterator([arc_diagram_label_permutations(T, inner_mod, label[(i-1)*m+1:i*m]) for i in 1:power])
        ]
    else
        error("Not implemented.")
    end
end


function arcdiag_is_lower_pair_label_bad(::SO, labeled_diag::Vector{Int}, k::Int)
    left_k = k % 2 == 1 ? k : k - 1
    return labeled_diag[left_k] == labeled_diag[left_k+1]
end

function arcdiag_is_lower_pair_label_bad(::GL, labeled_diag::Vector{Int}, k::Int)
    return false
end

function arcdiag_to_deformationmap(
    T::GL,
    diag::ArcDiagramDirected,
    sp::SmashProductLie{C},
    W::LieAlgebraModule=exterior_power_obj(base_module(sp), 2),
    case::Symbol=:exterior_power
) where {C <: RingElem}
    return arcdiag_to_deformationmap(T, arc_diagram(Undirected, diag), sp, W, case)
end

function arcdiag_to_deformationmap(
    T::Union{SO, GL},
    diag::ArcDiagramUndirected,
    sp::SmashProductLie{C},
    W::LieAlgebraModule=exterior_power_obj(base_module(sp), 2),
    case::Symbol=:exterior_power
) where {C <: RingElem}
    @req !_is_direct_sum(W)[1] "Not permitted for direct sums."
    ind_map = basis_index_mapping(W)

    dim_stdmod_V = base_lie_algebra(sp).n

    iso_pair_to_L = arc_diagram_lower_pair_to_L(T, dim_stdmod_V)

    if case == :exterior_power
        fl, Wbase, k = _is_exterior_power(W)
        @assert fl
        @assert k == 2
        nrows_kappa = ncols_kappa = dim(Wbase)
    elseif case == :tensor_product
        fl, W_factors = _is_tensor_product(W)
        @assert fl
        nrows_kappa, ncols_kappa = dim.(W_factors)
    else
        error("Unknown case")
    end

    kappa = zero_matrix(sp, nrows_kappa, ncols_kappa)
    for (label_index, upper_labels) in enumerate(arc_diagram_label_iterator(T, W, 1:dim_stdmod_V))

        i, j = ind_map[label_index]

        entry = arcdiag_to_deformationmap_entry(T, diag, W, upper_labels, sp, iso_pair_to_L, dim_stdmod_V)
        simplify(entry)

        kappa[i, j] += entry
        if case == :exterior_power
            kappa[j, i] -= entry
        end
    end
    return kappa
end


function arcdiag_to_deformationmap_entry(
    T::Union{SO, GL},
    diag::ArcDiagramUndirected,
    W::LieAlgebraModule{C},
    upper_labels::AbstractVector{Int},
    sp::SmashProductLie{C},
    iso_pair_to_L::Function,
    max_label::Int,
) where {C <: RingElem}
    entry = zero(sp)

    for (upper_labels, sgn_upper_labels) in arc_diagram_label_permutations(T, W, upper_labels)
        zeroprod = false
        lower_labels = [0 for _ in 1:n_lower_vertices(diag)]
        frees = Int[]
        for v in upper_vertices(diag)
            nv = neighbor(diag, v)
            if is_upper_vertex(nv) && upper_labels[vertex_index(v)] != upper_labels[vertex_index(nv)]
                zeroprod = true
                break
            elseif is_lower_vertex(nv)
                lower_labels[vertex_index(nv)] = upper_labels[vertex_index(v)]
            end
        end
        if zeroprod
            continue
        end
        for v in lower_vertices(diag)
            nv = neighbor(diag, v)
            if is_lower_vertex(nv) && vertex_index(v) < vertex_index(nv)
                push!(frees, vertex_index(v))
            end
        end

        entry_summand = zero(sp)

        # iterate over lower point labelings
        nextindex = 1
        while true
            if nextindex > length(frees)
                # begin inner
                coeff_lower_labels = 1
                basiselem = Int[]
                for k in 1:2:length(lower_labels)
                    gen_ind, coeff = iso_pair_to_L(lower_labels[k], lower_labels[k+1])
                    coeff_lower_labels *= coeff
                    if iszero(coeff_lower_labels)
                        break
                    end
                    append!(basiselem, gen_ind)
                end
                if !iszero(coeff_lower_labels)
                    symm_basiselem = sp(
                        underlying_algebra(sp)(
                            fill(C(1 // factorial(length(basiselem))), factorial(length(basiselem))),
                            [ind for ind in permutations(basiselem)],
                        ),
                    ) # TODO: benchmark use of `symmetrize` here once it is implemented with mutable arithmetics
                    entry_summand += coeff_lower_labels * symm_basiselem
                end
                # end inner

                nextindex -= 1
            end

            while nextindex >= 1 && lower_labels[frees[nextindex]] == max_label
                lower_labels[frees[nextindex]] = 0
                lower_labels[vertex_index(_neighbor_of_lower_vertex(diag, frees[nextindex]))] = 0
                nextindex -= 1
            end
            if nextindex == 0
                break
            end
            lower_labels[frees[nextindex]] += 1
            lower_labels[vertex_index(_neighbor_of_lower_vertex(diag, frees[nextindex]))] += 1

            if arcdiag_is_lower_pair_label_bad(T, lower_labels, frees[nextindex]) || arcdiag_is_lower_pair_label_bad(
                T,
                lower_labels,
                vertex_index(_neighbor_of_lower_vertex(diag, frees[nextindex])),
            )
                continue
            end

            nextindex += 1
        end

        entry += sgn_upper_labels * entry_summand
    end
    return entry
end

function is_in_canonical_form_for_deformation(A::ArcDiagramDirected; case::Symbol=:unknown)
    ### if first arc goes from upper to lower, it goes to the first lower vertex
    first_neigh = outneighbor(A, upper_vertex(A, 1))
    if is_lower_vertex(first_neigh) && vertex_index(first_neigh) != 1
        return false
    end
    ###
    if case == :exterior_power
        # first arc from upper right part may not got to the first lower vertex -> swap upper parts instead
        first_neigh_right = outneighbor(A, upper_vertex(A, div(n_upper_vertices(A), 2) + 1))
        if is_lower_vertex(first_neigh_right) && vertex_index(first_neigh_right) == 1
            return false
        end
    end
    ###

    ### lower part cycles in the end and desc sorted by length
    n_lower_pairs = div(n_lower_vertices(A), 2)
    lower_loop_lengths = fill(typemin(Int), n_lower_pairs)
    for i in 1:n_lower_pairs
        lower_loop_lengths[i] != typemin(Int) && continue
        loop_inds = [i]
        start_v = lower_vertex(A, 2i-1)
        curr_v = outneighbor(A, lower_vertex(A, 2i))
        contains_upper = false
        while curr_v != start_v
            if is_upper_vertex(curr_v)
                contains_upper = true
                break
            end
            curr_v_ind = vertex_index(curr_v) + 1
            push!(loop_inds, div(curr_v_ind, 2))
            curr_v = outneighbor(A, lower_vertex(A, curr_v_ind))
        end
        issorted(loop_inds) || return false
        length(loop_inds) == maximum(loop_inds) - minimum(loop_inds) + 1 || return false
        if contains_upper
            for i in loop_inds
                lower_loop_lengths[i] = typemax(Int)
            end
        else
            for i in loop_inds
                lower_loop_lengths[i] = length(loop_inds)
            end
        end
    end
    issorted(lower_loop_lengths; rev=true) || return false
    ###

    return true
end

function is_in_canonical_form_for_deformation(A::ArcDiagramUndirected; case::Symbol=:unknown)
    ### if first arc goes from upper to lower, it goes to the first lower vertex
    first_neigh = neighbor(A, upper_vertex(A, 1))
    if is_lower_vertex(first_neigh) && vertex_index(first_neigh) != 1
        return false
    end
    ###
    if case == :exterior_power
        # first arc from upper right part may not got to the first lower vertex -> swap upper parts instead
        first_neigh_right = neighbor(A, upper_vertex(A, div(n_upper_vertices(A), 2) + 1))
        if is_lower_vertex(first_neigh_right) && vertex_index(first_neigh_right) == 1
            return false
        end
    end
    ###

    ### lower part cycles in the end and desc sorted by length
    n_lower_pairs = div(n_lower_vertices(A), 2)
    lower_loop_lengths = fill(typemin(Int), n_lower_pairs)
    for i in 1:n_lower_pairs
        lower_loop_lengths[i] != typemin(Int) && continue
        loop_inds = [i]
        start_v = lower_vertex(A, 2i-1)
        curr_v = neighbor(A, lower_vertex(A, 2i))
        contains_upper = false
        while curr_v != start_v
            if is_upper_vertex(curr_v)
                contains_upper = true
                break
            end
            curr_v_ind = vertex_index(curr_v) + 1
            isodd(curr_v_ind) && return false
            push!(loop_inds, div(curr_v_ind, 2))
            curr_v = neighbor(A, lower_vertex(A, curr_v_ind))
        end
        issorted(loop_inds) || return false
        length(loop_inds) == maximum(loop_inds) - minimum(loop_inds) + 1 || return false
        if contains_upper
            for i in loop_inds
                lower_loop_lengths[i] = typemax(Int)
            end
        else
            for i in loop_inds
                lower_loop_lengths[i] = length(loop_inds)
            end
        end
    end
    issorted(lower_loop_lengths; rev=true) || return false
    ###

    return true
end
