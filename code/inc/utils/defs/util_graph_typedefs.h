#ifndef UTIL_GRAPH_TYPEDEFS_H
#define UTIL_GRAPH_TYPEDEFS_H

class Polybezier;
class Intersection;

// ---------------------------------
// QT incompatability with boost fix
// ---------------------------------
#ifndef Q_MOC_RUN
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/function.hpp>
#include <boost/graph/copy.hpp>
#endif

// =========================================
// =========  UTIL_GRAPH_TYPEDEFS  =========
// =========================================

/* -----------------------------------------
 *
 * There are two types of custom graphs used
 *
 * curve graph
 * - main graph type
 * - vertexes hold polybez pointers
 * - edges hold intersection pointers
 * - edges also have an index property to
 *   enable subgraph cross referencing
 * - subgraph type to enable efficient
 *   storing and cross referencing of
 *   connected components, or what we would
 *   see in a curve network as a patch.
 *
 * intersection graph
 * - helper graph type
 * - vertexes hold curve graph edge descriptors
 *   which translate to intersection pointers
 * - edges hold curve graph vertex descriptors
 *   which translate to polybez pointers
 * - edges also have an edge weight property to
 *   enable djikstra calculations.
 *
 * ----------------------------------------- */

// --------------------------------------------------------
// polybez shared pointers and intersection shared pointers
// are commonly used in graphs and canvas
typedef boost::shared_ptr<Polybezier> polybez_ptr_t;
typedef boost::shared_ptr<Intersection> isec_ptr_t;
// --------------------------------------------------------

namespace boost {
    // ----------------------------------------------------------------
    // custom internal properties for vertex and edge data pointers
    enum vertex_data_pointer_t { vertex_data_pointer };
    enum edge_data_pointer_t   { edge_data_pointer   };

    BOOST_INSTALL_PROPERTY( vertex, data_pointer );
    BOOST_INSTALL_PROPERTY( edge,   data_pointer );

    // ----------------------------------------------------------------
    // custom internal properties for vertex and edge descriptor copies
    enum vertex_copied_descriptor_t { vertex_copied_descriptor };
    enum edge_copied_descriptor_t   { edge_copied_descriptor   };

    BOOST_INSTALL_PROPERTY( vertex, copied_descriptor );
    BOOST_INSTALL_PROPERTY( edge,   copied_descriptor );
}
// --------------------------------------------------------------------

// -----------------------------------------------------------------------------------------------------------------
// main graph type is a (polybez) curve graph
typedef boost::subgraph< // enable subgraphs
        boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              boost::property<boost::vertex_data_pointer_t, polybez_ptr_t, // vertexes point to polybez curves
                                              boost::property<boost::vertex_index_t, int>>,
                              boost::property<boost::edge_data_pointer_t, isec_ptr_t,      // edges point to intersections
                                              boost::property<boost::edge_index_t, int>>>> curve_graph_t;

// curve graph vertex descriptors are used to index and vertex iterators are used for vertex traversal
// dereferencing a vertex iterator gives the vertex descriptor
typedef boost::graph_traits<curve_graph_t>::vertex_descriptor curve_graph_vertex_desc_t;
typedef boost::graph_traits<curve_graph_t>::vertex_iterator   curve_graph_vertex_iter_t;

// curve graph edge descriptors are used to index and edge iterators are used for edge traversal
// out edge iterators are used to traverse all out edges from a vertex
// dereferencing an edge iterator gives the edge descriptor
typedef boost::graph_traits<curve_graph_t>::edge_descriptor   curve_graph_edge_desc_t;
typedef boost::graph_traits<curve_graph_t>::edge_iterator     curve_graph_edge_iter_t;
typedef boost::graph_traits<curve_graph_t>::out_edge_iterator curve_graph_out_edge_iter_t;

// curve graph vertex data pointer property map is used to read and write polybez pointers to curve graph vertexes
// curve graph edge data pointer property map is used to read and write intersection pointers to curve graph edges
typedef boost::property_map<curve_graph_t, boost::vertex_data_pointer_t>::type curve_graph_vertex_data_map_t;
typedef boost::property_map<curve_graph_t, boost::edge_data_pointer_t>::type curve_graph_edge_data_map_t;

// curve graph vertex index property map is used to read and write curve graph vertex indices
// curve graph edge index property map is used to read and write curve graph edge indices
// both of these internal properties must be explicitly set and updated
typedef boost::property_map<curve_graph_t, boost::vertex_index_t>::type curve_graph_vertex_index_map_t;
typedef boost::property_map<curve_graph_t, boost::edge_index_t>::type curve_graph_edge_index_map_t;

// filter template for curve graph that can be customized with lambdas for the predicates
// filters can be used to make iterator ranges on that will respect the filter predicates
// filters are filtered views of graphs that retain original descriptors
typedef boost::filtered_graph<curve_graph_t, boost::function<bool(curve_graph_edge_desc_t)>,
boost::function<bool(curve_graph_vertex_desc_t)>> curve_graph_filter_t;
// -----------------------------------------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------------------------------------
// helper graph type is an (intersection) isec graph derived from curve graph
// isec graph vertexes hold curve graph edge descriptors (intersections)
// isec graph edges hold curve graph vertexes (polybez curves) and weights for graph traversal algorithms
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              boost::property<boost::vertex_copied_descriptor_t, curve_graph_edge_desc_t>, // vertexes index to given edges (intersection pointers)
                              boost::property<boost::edge_copied_descriptor_t, curve_graph_vertex_desc_t,  // edges index to given vertexes (polybez pointers)
                                              boost::property<boost::edge_weight_t, int>>> isec_graph_t;

// isec graph vertex and edge descriptors are used to index
// dereferencing an iterator gives the descriptor
typedef boost::graph_traits<isec_graph_t>::vertex_descriptor isec_graph_vertex_desc_t;
typedef boost::graph_traits<isec_graph_t>::edge_descriptor isec_graph_edge_desc_t;

// isec graph vertex copied descriptor property map is used to read and write curve graph edge descriptors to isec graph vertexes
// isec graph edge copied descriptor property map is used to read and write curve graph vertex descriptors to isec graph edges
typedef boost::property_map<isec_graph_t, boost::vertex_copied_descriptor_t>::type isec_graph_copied_edge_map_t;
typedef boost::property_map<isec_graph_t, boost::edge_copied_descriptor_t>::type isec_graph_copied_vertex_map_t;
typedef boost::property_map<isec_graph_t, boost::edge_weight_t>::type isec_graph_weight_map_t;
// -----------------------------------------------------------------------------------------------------------------

#endif // UTIL_GRAPH_TYPEDEFS_H
