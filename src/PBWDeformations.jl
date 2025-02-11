module PBWDeformations

using Preferences: Preferences, @load_preference, @set_preferences!

import ProgressMeter

using Oscar

using Oscar.AbstractAlgebra: ProductIterator

using Oscar: IntegerUnion

using Oscar: _is_dual
using Oscar: _is_direct_sum
using Oscar: _is_exterior_power
using Oscar: _is_symmetric_power
using Oscar: _is_tensor_power
using Oscar: _is_tensor_product
using Oscar: _is_standard_module

using Oscar.LieAlgebras: combinations
using Oscar.LieAlgebras: multicombinations
using Oscar.LieAlgebras: permutations
using Oscar.LieAlgebras: permutations_with_sign


import Oscar.AbstractAlgebra: change_base_ring
import Oscar.AbstractAlgebra: coefficient_ring
import Oscar.AbstractAlgebra: elem_type
import Oscar.AbstractAlgebra: gen
import Oscar.AbstractAlgebra: gens
import Oscar.AbstractAlgebra: ngens
import Oscar.AbstractAlgebra: parent_type

import Oscar: base_lie_algebra
import Oscar: comm
import Oscar: data
import Oscar: edges
import Oscar: n_edges
import Oscar: neighbors
import Oscar: n_vertices
import Oscar: simplify
import Oscar: vertices

import Base: deepcopy_internal
import Base: hash
import Base: isone
import Base: iszero
import Base: length
import Base: one
import Base: parent
import Base: show
import Base: sum
import Base: zero


export AbstractLieAlgebra, AbstractLieAlgebraElem
export ArcDiagDeformBasis
export ArcDiagram
export ArcDiagramDirected
export ArcDiagramUndirected
export ArcDiagramVertex
export DeformationMap
export DeformBasis
export LieAlgebra, LieAlgebraElem
export LieAlgebraModule, LieAlgebraModuleElem
export LinearLieAlgebra, LinearLieAlgebraElem
export PseudographDeformBasis
export PseudographLabelled
export SmashProductLie, SmashProductLieElem
export SmashProductLieDeform, SmashProductLieDeformElem
export StdDeformBasis

export abstract_module
export all_arc_diagrams
export all_pbwdeformations
export all_pseudographs
export arc_diagram
export base_lie_algebra
export base_module
export deform
export edge_labels
export edges
export exterior_power
export exterior_power_obj
export general_linear_lie_algebra
export get_show_colorful_html
export inneighbor
export inneighbors
export is_crossing_free
export is_in_span
export is_linearly_independent
export is_linearly_independent_with_relations
export is_pbwdeformation
export is_span_equal
export is_span_subset
export isomorphic_module_with_simple_structure
export lookup_data
export lower_vertex, is_lower_vertex
export lower_vertices
export matrix_repr_basis
export n_edges
export n_lower_vertices
export n_upper_vertices
export n_vertices
export neighbor
export neighbors
export outneighbor
export outneighbors
export pbwdeform_eqs
export set_show_colorful_html
export simple_module
export smash_product
export special_linear_lie_algebra
export special_orthogonal_lie_algebra
export standard_module
export symmetric_deformation
export symmetric_power
export symmetric_power_obj
export symmetrize
export tensor_power
export tensor_power_obj
export to_arcdiag
export underlying_algebra
export upper_vertex, is_upper_vertex
export upper_vertices
export vertex_index
export vertices

function __init__()
    add_verbosity_scope(:PBWDeformations)
end

include("OscarPatches.jl")

include("Types.jl")

include("LinearIndependence.jl")
include("ModuleSimpleStructure.jl")
include("Misc.jl")

include("DeformationBases/DeformBasis.jl")

include("SmashProductLie.jl")
include("SmashProductLieDeform.jl")
include("SmashProductPBWDeformLie.jl")
include("ArcDiagram.jl")
include("Pseudograph.jl")
include("GlnGraph.jl")

include("DeformationBases/ArcDiagDeformBasis.jl")
include("DeformationBases/PseudographDeformBasis.jl")
include("DeformationBases/StdDeformBasis.jl")

end
