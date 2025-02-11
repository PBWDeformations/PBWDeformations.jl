const GlnGraphEdge = Tuple{Int, Int}

struct GlnGraph
    n_left_verts::Int
    n_right_verts::Int
    parity_verts::BitVector    # true is out
    edges::Vector{GlnGraphEdge}

    function GlnGraph(
        n_left_verts::Int,
        n_right_verts::Int,
        parity_verts::BitVector,
        edges::Vector{GlnGraphEdge};
        check::Bool=true,
        sort::Bool=true,
    )
        if check
            @req n_left_verts + n_right_verts == length(parity_verts) "number of vertices mismatch"
            @req 2 * length(edges) == length(parity_verts) "number of edges mismatch"
            out_verts = first.(edges)
            @req allunique(out_verts) "out vertex with wrong degree"
            @req all(parity_verts[out_verts]) "out vertex with wrong parity"
            in_verts = last.(edges)
            @req allunique(in_verts) "in vertex with wrong degree"
            @req all(!, parity_verts[in_verts]) "in vertex with wrong parity"
        end
        if sort
            sort!(edges)
        end
        return new(n_left_verts, n_right_verts, parity_verts, edges)
    end
end

function n_edges(g::GlnGraph)
    return length(g.edges)
end

function all_gln_graphs(n_left_verts::Int, n_right_verts::Int, parity_verts::BitVector)
    @req n_left_verts + n_right_verts == length(parity_verts) "number of vertices mismatch"

    result = GlnGraph[]
    !iszero(parity_diff(parity_verts)) && return result # no graph if parity is not balanced

    out_indices = findall(parity_verts)
    in_indices = findall(!, parity_verts)

    return (
        begin
            edges = collect(zip(out_indices, in_indices_permuted))
            GlnGraph(n_left_verts, n_right_verts, parity_verts, edges)
        end for in_indices_permuted in permutations(in_indices)
    )
end

function number_of_gln_graphs(n_left_verts::Int, n_right_verts::Int, parity_verts::BitVector)
    @req n_left_verts + n_right_verts == length(parity_verts) "number of vertices mismatch"
    !iszero(parity_diff(parity_verts)) && return 0 # no graph if parity is not balanced

    return factorial(div(length(parity_verts), 2))
end


function gln_graph_labelings_with_partitions(n_edges::Int, total_weight::Int)
    return (
        (wcomp, part) for part_weight in 0:total_weight for part in partitions(part_weight) for
        wcomp in weak_compositions(total_weight - part_weight, n_edges)
    )
end

function number_of_gln_graph_labelings_with_partitions(n_edges::Int, total_weight::Int)
    return sum(
        part_weight -> n_weak_compositions(total_weight - part_weight, n_edges) * n_partitions(part_weight),
        0:total_weight;
        init=zero(ZZ),
    )
end


function arc_diagram(g::GlnGraph, labeling::Vector{Int}, part::Partition{Int}=Partition(Int[]))
    @req n_edges(g) == length(labeling) "number of edge labels mismatch"

    n_upper_verts = 2 * n_edges(g)
    n_lower_verts = 2 * (sum(labeling) + sum(part))

    parity_upper_verts = g.parity_verts
    @assert length(parity_upper_verts) == n_upper_verts
    parity_lower_verts = BitVector(undef, n_lower_verts)
    parity_lower_verts[1:2:end] .= true

    upper_adj = zeros(Int, n_upper_verts)
    lower_adj = zeros(Int, n_lower_verts)
    addarc_uu(i, j) = (upper_adj[i] = -j; upper_adj[j] = -i)
    addarc_ll(i, j) = (lower_adj[i] = j; lower_adj[j] = i)
    addarc_ul(i, j) = (upper_adj[i] = j; lower_adj[j] = -i)
    addarc_lu(i, j) = (lower_adj[i] = -j; upper_adj[j] = i)
    # i = 1
    # j = div(n_upper_verts, 2) + 1
    k = 1
    for ((v, w), label) in zip(g.edges, labeling)
        if label == 0
            addarc_uu(v, w)
        else
            addarc_ul(v, k)
            k += 1
            for _ in 2:label
                addarc_ll(k, k + 1)
                k += 2
            end
            addarc_lu(k, w)
            k += 1
        end
    end

    for p in part
        start = k
        k += 1
        for _ in 2:p
            addarc_ll(k, k + 1)
            k += 2
        end
        addarc_ll(k, start)
        k += 1
    end

    return arc_diagram(Directed, parity_upper_verts, parity_lower_verts, upper_adj, lower_adj)
end
