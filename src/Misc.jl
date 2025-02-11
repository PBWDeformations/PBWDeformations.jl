function parity_diff(v::BitVector)
    return 2 * sum(v) - length(v)
end

function symmetrize(f::FreeAssociativeAlgebraElem)
    R = parent(f)
    g = zero(R)
    for (c, exp) in zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_words(f))
        g += R(
            fill(c * QQ(1, factorial(length(exp))), factorial(length(exp))),
            [ind for ind in permutations(exp)],
        )

    end
    return g
end

function set_show_colorful_html(value::Bool)
    old_value = get_show_colorful_html()
    @set_preferences!("show_colorful_html" => value)
    return old_value
end

function get_show_colorful_html()
    return @load_preference("show_colorful_html", default = false)::Bool
end

Base.show(io::IO, ::MIME"text/html", A::ArcDiagram) = print(io, svg_string(A; colorful=get_show_colorful_html()))

function svg_string_defs(::ArcDiagram)
    return ""
end

function svg_string_defs(::ArcDiagramDirected)
    # marker to be used as an arrowhead
    return """
        <marker id=\"arrow\" refX=\"2\" refY=\"2\" markerWidth=\"3\" markerHeight=\"4\" orient=\"auto\">
            <path d=\"M 0 0 V 4 L 2 2 Z\" stroke=\"context-stroke\" fill=\"context-stroke\"/>
        </marker>
        """
end

function svg_string_edge_iterator(A::ArcDiagramUndirected)
    return ((v, neighbor(A, v)) for v in vertices(A) if _vertex_lt(v, neighbor(A, v)))
end

function svg_string_edge_iterator(A::ArcDiagramDirected)
    return ((v, outneighbor(A, v)) for v in vertices(A) if !isnothing(outneighbor(A, v)))
end

function svg_string_path(v::ArcDiagramVertex, nv::ArcDiagramVertex, dims)
    (;margin, h, w) = dims
    x_v = (vertex_index(v) - 1) * w + margin
    x_nv = (vertex_index(nv) - 1) * w + margin
    y_v = is_upper_vertex(v) ? margin : h + margin
    y_nv = is_upper_vertex(nv) ? margin : h + margin

    function y_c(x1, x2)
        return (1-exp(-0.15*(1+abs(x1 - x2))))*0.85
    end

    if is_upper_vertex(v)
        if is_upper_vertex(nv)
            x_c1 = (8*x_v + 2*x_nv) / 10
            y_c1 = y_c(vertex_index(v), vertex_index(nv)) * h + margin
            x_c2 = (2*x_v + 8*x_nv) / 10
            y_c2 = y_c(vertex_index(v), vertex_index(nv)) * h + margin
            return "M $(x_v) $(y_v) C $(x_c1) $(y_c1), $(x_c2) $(y_c2), $(x_nv) $(y_nv)"
        else
            return "M $(x_v) $(y_v) L $(x_nv) $(y_nv)"
        end
    else
        if is_upper_vertex(nv)
            return "M $(x_v) $(y_v) L $(x_nv) $(y_nv)"
        else
            x_c1 = (8*x_v + 2*x_nv) / 10
            y_c1 = (1 - y_c(vertex_index(v), vertex_index(nv))) * h + margin
            x_c2 = (2*x_v + 8*x_nv) / 10
            y_c2 = (1 - y_c(vertex_index(v), vertex_index(nv))) * h + margin
            return "M $(x_v) $(y_v) C $(x_c1) $(y_c1), $(x_c2) $(y_c2), $(x_nv) $(y_nv)"
        end
    end
end

function svg_string(A::ArcDiagram; colorful=false)
    color_mono = "#008"
    size = 40
    dims = (
        rad=0.1*size,
        margin=0.2*size,
        h=1*size,
        w=0.5*size,
        stroke_width=2, #px
    )

    height = dims.h + 2 * dims.margin
    width = (max(n_upper_vertices(A), n_lower_vertices(A)) - 1) * dims.w + 2 * dims.margin

    svg = "<svg height=\"$(height)\" width=\"$(width)\" style=\"vertical-align:middle;\">"

    svg *= "<defs>$(svg_string_defs(A))</defs>"

    i = 0
    for (v, nv) in svg_string_edge_iterator(A)
        color = colorful ? "hsl($(105*i),100%,45%)" : color_mono
        svg *= "<path d=\"$(svg_string_path(v, nv, dims))\" fill=\"transparent\" stroke=\"$color\" stroke-width=\"$(dims.stroke_width)px\" $(A isa ArcDiagramDirected ? "marker-end=\"url(#arrow)\"" : "")/>"
        i += 1
    end

    svg *= "</svg>"
    return svg
end

function Base.show(io::IO, mime::MIME"text/html", x::Tuple{Int,ArcDiagram})
    print(io, "($(first(x)), ")
    show(io, mime, x[2])
    print(io, ")")
end

function Base.show(io::IO, mime::MIME"text/html", x::Tuple{ArcDiagram,Vararg{ArcDiagram}})
    print(io, "<div>(")
    show(io, mime, first(x))
    for A in Base.tail(x)
        print(io, ", ")
        show(io, mime, A)
    end
    print(io, ")</div>")
end
