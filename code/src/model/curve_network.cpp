#include "model/curve_network.h"
#include "utils/gui/util_obj_writer.h"

#include <cmath>
#include <iostream>
#include <string>
#include <sstream>

#include <QFile>
#include <QDomElement>
#include <QTextStream>
#include <QQueue>

#include <boost/graph/dijkstra_shortest_paths.hpp>

// ===========================================================================
// ==========================        HELPERS        ==========================
// ===========================================================================

boost::optional<curve_graph_vertex_desc_t> Curve_Network::find_vertex_by_uuid(boost::uuids::uuid uuid, curve_graph_t graph)
{
    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, graph);
    for (auto polybez_desc : boost::make_iterator_range(boost::vertices(graph))) {
        if (polybez_ptr_map[polybez_desc]->m_uuid == uuid) {
            return polybez_desc;
        }
    }
    return boost::none;
}

boost::optional<curve_graph_edge_desc_t> Curve_Network::find_edge_by_uuid(boost::uuids::uuid uuid, curve_graph_t graph)
{
    curve_graph_edge_data_map_t isec_ptr_map = boost::get(boost::edge_data_pointer, graph);
    for (auto isec_desc : boost::make_iterator_range(boost::edges(graph))) {
        if (isec_ptr_map[isec_desc]->m_uuid == uuid) {
            return isec_desc;
        }
    }
    return boost::none;
}

// ===========================================================================
// ===========================================================================

Curve_Network::Curve_Network() : m_curve_network_graph(0), m_curve_network_subgraphs(0) { }
Curve_Network::~Curve_Network() { }

void Curve_Network::add_polybezier(QVector<point_t> ctrl_points, curve_type_t curve_type, bool closed_curve)
{
    // -------------------------------------------------------------------------
    // build the bezier polycurve by first finding all of its bezier curves
    // using the ctrl points returned by the fitting algorithm in groups of four
    QVector<bezier_curve_s> bezier_curves;
    for (int p = 0; p < ctrl_points.size() - 1; p += 3) {
        bezier_curve_s bezier_curve;
        bezier_curve.ctrl_points[0] = ctrl_points[p + 0];
        bezier_curve.ctrl_points[1] = ctrl_points[p + 1];
        bezier_curve.ctrl_points[2] = ctrl_points[p + 2];
        bezier_curve.ctrl_points[3] = ctrl_points[p + 3];
        bezier_curves.push_back(bezier_curve);
    }
    polybez_ptr_t polybez_ptr(boost::make_shared<Polybezier>(closed_curve, curve_type, bezier_curves, boost::none, boost::none));
    // -------------------------------------------------------------------------

    // --------------------------------------------------------------------
    // add the polybezier to the graph network and update all intersections
    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, m_curve_network_graph);
    curve_graph_vertex_desc_t curve_descriptor = boost::add_vertex(m_curve_network_graph);
    boost::put(polybez_ptr_map, curve_descriptor, polybez_ptr);

    update_intersections(curve_descriptor);
    update_connected_components();

    for (int p = 0; p < m_curve_network_subgraphs.size(); ++p) {
        regularize_tangent_frames(m_curve_network_subgraphs[p]);
        update_intersection_orientations_from_silhouettes(m_curve_network_subgraphs[p]);
        //integrate_polybeziers_normals(m_curve_network_subgraphs[p]);
    }
    // --------------------------------------------------------------------
}

// decompose bezier curves into sequences of straight lines, testing for intersections between each of them.
QVector<isec_ptr_t> Curve_Network::intersect_polybeziers(curve_graph_vertex_desc_t curve_desc_a, curve_graph_vertex_desc_t curve_desc_b)
{
    QVector<isec_ptr_t> out_isecs(0);

    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, m_curve_network_graph);
    polybez_ptr_t polycurve_a = polybez_ptr_map[curve_desc_a];
    polybez_ptr_t polycurve_b = polybez_ptr_map[curve_desc_b];

    point_t P, P0, P1, P2, P3;

    for (int a = 0; a < polycurve_a->m_bezier_curves.size(); ++a) {
        for(double t = T_STEP; t <= 1.0; t += T_STEP) {
            P0 = polycurve_a->compute_point_at(a, t - T_STEP);
            P1 = polycurve_a->compute_point_at(a, t);
            for (int b = 0; b < polycurve_b->m_bezier_curves.size(); ++b) {
                for(double r = T_STEP; r <= 1.0; r += T_STEP) {
                    P2 = polycurve_b->compute_point_at(b, r - T_STEP);
                    P3 = polycurve_b->compute_point_at(b, r);

                    double denominator = (P3.y - P2.y) * (P1.x - P0.x) - (P3.x - P2.x) * (P1.y - P0.y);
                    if (denominator * denominator >= 0.00000000000001) {
                        double ua = ((P3.x - P2.x) * (P0.y - P2.y) - (P3.y - P2.y) * (P0.x - P2.x)) / denominator;
                        double ub = ((P1.x - P0.x) * (P0.y - P2.y) - (P1.y - P0.y) * (P0.x - P2.x)) / denominator;
                        if (ua > 0.0 && ua < 1.0 && ub > 0.0 && ub < 1.0) {
                            point_t isec_pos = P0 + (P1 - P0) * ua;
                            vector_t isec_tan_a = polycurve_a->compute_tangent_at(a, (t - T_STEP) + (ua * T_STEP));
                            vector_t isec_tan_b = polycurve_b->compute_tangent_at(b, (r - T_STEP) + (ub * T_STEP));
                            double isec_tval_a = (t - T_STEP) + (ua * T_STEP) + a;
                            double isec_tval_b = (r - T_STEP) + (ub * T_STEP) + b;

                            isec_ptr_t isec(boost::make_shared<Intersection>(polycurve_a->m_uuid, polycurve_b->m_uuid, boost::none, boost::none, boost::none,
                                                                             isec_pos, isec_tval_a, isec_tval_b, boost::none, isec_tan_a, isec_tan_b, boost::none, boost::none));
                            isec->calculate_tangent_frame();
                            out_isecs.push_back(isec);
                        }
                    }
                }
            }
        }
    }
    return out_isecs;
}

void Curve_Network::update_intersections(curve_graph_vertex_desc_t curve_descriptor)
{
    // TODO: deletion
    // --------------------------------------------------------
    // iterate across all polycurves and test for intersections
    // against the recently added polycurve
    for (auto other_curve_descriptor : boost::make_iterator_range(boost::vertices(m_curve_network_graph))) {
        if (curve_descriptor == other_curve_descriptor) {
            continue;
        }
        QVector<isec_ptr_t> isecs = intersect_polybeziers(other_curve_descriptor, curve_descriptor);
        // --------------------------------------------------------

        // -------------------------------------------------------------
        // add the new intersections to the curve network graph as edges
        // as a reminder, curve network graph vertices are bezier polycurves
        // and curve network graph edges are intersections between bezier polycurves
        curve_graph_edge_data_map_t isec_ptr_map = boost::get(boost::edge_data_pointer, m_curve_network_graph);
        for (int i = 0; i < isecs.size(); ++i) {
            curve_graph_edge_desc_t isec_desc = boost::add_edge(curve_descriptor, other_curve_descriptor, m_curve_network_graph).first;
            boost::put(isec_ptr_map, isec_desc, isecs[i]);
        }
        // -------------------------------------------------------------
    }
}

void Curve_Network::update_connected_components()
{
    m_curve_network_subgraphs.clear();
    std::vector<curve_graph_filter_t> filtered_connected = filter_connected_components(m_curve_network_graph);
    for (size_t c = 0; c < filtered_connected.size(); ++c) {
        curve_graph_filter_t component = filtered_connected[c];
        curve_graph_t component_subgraph = m_curve_network_graph.create_subgraph();
        for (auto polybez : boost::make_iterator_range(boost::vertices(component))) {
            boost::add_vertex(polybez, component_subgraph);
        }
        m_curve_network_subgraphs.push_back(component_subgraph);
    }
}

void Curve_Network::update_intersection_orientations_from_silhouettes(curve_graph_t& curve_graph)
{
    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, curve_graph);
    curve_graph_edge_data_map_t isec_ptr_map = boost::get(boost::edge_data_pointer, curve_graph);

    for (auto polybez_desc : boost::make_iterator_range(boost::vertices(curve_graph))) {
        if (polybez_ptr_map[polybez_desc]->m_curve_type == CURVE_SILHOUETTE) {
            curve_graph_vertex_desc_t silhouette_polybez_desc = polybez_desc;
            for (auto silhouette_isec_desc : boost::make_iterator_range(boost::out_edges(silhouette_polybez_desc, curve_graph))) {
                boost::shared_ptr<Intersection> silhouette_isec_ptr = isec_ptr_map[silhouette_isec_desc];

                // get polybezier graph descriptors in the local graph scope
                boost::optional<curve_graph_vertex_desc_t> opt_polybez_desc_a = find_vertex_by_uuid(silhouette_isec_ptr->m_polybez_a_uuid, m_curve_network_graph);
                boost::optional<curve_graph_vertex_desc_t> opt_polybez_desc_b = find_vertex_by_uuid(silhouette_isec_ptr->m_polybez_b_uuid, m_curve_network_graph);
                if (opt_polybez_desc_a == boost::none || opt_polybez_desc_b == boost::none) {
                    std::cout << "[ ERROR ] could not find polybezier descriptors!" << std::endl;
                    continue;
                }
                curve_graph_vertex_desc_t polybez_desc_a = curve_graph.global_to_local(opt_polybez_desc_a.value());
                curve_graph_vertex_desc_t polybez_desc_b = curve_graph.global_to_local(opt_polybez_desc_b.value());

                // get information about the intersection in question
                curve_graph_vertex_desc_t cross_polybez_desc = -100; // curve descriptor of cross curve intersecting with silhouette

                vector_t silhouette_cross_tangent = { -100.0, -100.0, -100.0 }; // tangent along silhouette at intersecting point with cross curve
                vector_t cross_silhouette_tangent = { -100.0, -100.0, -100.0 }; // tangent along cross curve at intersecting point with silhouette

                double cross_silhouette_tval = -100.0; // tval along cross curve at intersecting point with silhouette

                if (polybez_desc_a == silhouette_polybez_desc) { // polybezier b must be the intersecting non-silhouette
                    if (polybez_ptr_map[polybez_desc_b]->m_curve_type != CURVE_CROSS) { continue; } // we do not care about non-cross curve intersections here
                    cross_polybez_desc = polybez_desc_b;

                    silhouette_cross_tangent = silhouette_isec_ptr->m_tangent_a;
                    cross_silhouette_tangent = silhouette_isec_ptr->m_tangent_b;
                    cross_silhouette_tval = silhouette_isec_ptr->m_tval_b;
                }

                else if (polybez_desc_b == silhouette_polybez_desc) { // polybezier a must be the intersecting non-silhouette
                    if (polybez_ptr_map[polybez_desc_a]->m_curve_type != CURVE_CROSS) { continue; } // we do not care about non-cross curve intersections here
                    cross_polybez_desc = polybez_desc_a;

                    silhouette_cross_tangent = silhouette_isec_ptr->m_tangent_b;
                    cross_silhouette_tangent = silhouette_isec_ptr->m_tangent_a;
                    cross_silhouette_tval = silhouette_isec_ptr->m_tval_a;
                }

                // find the xisec in the immediate vicinity of the silhouette intersection, along the cross curve in question
                boost::optional<curve_graph_edge_desc_t> cross_isec_desc = boost::none; // intersection descriptor of xisec immediately prior, or after, silhouette intersection
                double tval_proximity = 1000000; // running comparator to find closest tval between all tvals along the cross curve we are searching on
                double tval_closest = -100; // closest tval to silhouette intersection along the cross curve we are searching on

                for (auto search_isec_desc : boost::make_iterator_range(boost::out_edges(cross_polybez_desc, curve_graph))) {
                    if (search_isec_desc == silhouette_isec_desc) { continue; } // do not test the silhouette intersection itself

                    boost::shared_ptr<Intersection> search_isec_ptr = isec_ptr_map[search_isec_desc];
                    // get polybezier graph descriptors in the local graph scope
                    boost::optional<curve_graph_vertex_desc_t> opt_search_polybez_desc_a = find_vertex_by_uuid(search_isec_ptr->m_polybez_a_uuid, m_curve_network_graph);
                    boost::optional<curve_graph_vertex_desc_t> opt_search_polybez_desc_b = find_vertex_by_uuid(search_isec_ptr->m_polybez_b_uuid, m_curve_network_graph);
                    if (opt_search_polybez_desc_a == boost::none || opt_search_polybez_desc_b == boost::none) {
                        std::cout << "[ ERROR ] could not find polybezier descriptors!" << std::endl;
                        continue;
                    }
                    curve_graph_vertex_desc_t search_polybez_desc_a = curve_graph.global_to_local(opt_search_polybez_desc_a.value());
                    curve_graph_vertex_desc_t search_polybez_desc_b = curve_graph.global_to_local(opt_search_polybez_desc_b.value());

                    double search_tval = -100.0;
                    if (search_polybez_desc_a == cross_polybez_desc && polybez_ptr_map[search_polybez_desc_b]->m_curve_type == CURVE_CROSS) {
                        search_tval = search_isec_ptr->m_tval_a;
                    }
                    else if (search_polybez_desc_b == cross_polybez_desc && polybez_ptr_map[search_polybez_desc_a]->m_curve_type == CURVE_CROSS) {
                        search_tval = search_isec_ptr->m_tval_b;
                    }

                    if (abs(search_tval - cross_silhouette_tval) < tval_proximity) {
                        tval_proximity = abs(search_tval - cross_silhouette_tval);
                        tval_closest = search_tval;
                        cross_isec_desc = search_isec_desc;
                    }
                }

                // -------------------------------------------------------------------------------------------------------------------
                // begin aligning orientations with respect to the silhouette curve
                // this amounts to pointing tangents at cross-silhouette intersections towards the interior of the object
                // then flipping the tangents along the silhouettes at these intersections such that both interior pointing tangents and
                // silhouette tangents have maximal angle (maximal visibility)
                // then calculating the normal at the cross-silhouette intersections using these tangents
                // and finally flipping cross-hair orientations along the cross-curve that disagree with the cross-silhouette normal
                // -------------------------------------------------------------------------------------------------------------------

                bool should_reverse_direction = false;
                if (cross_isec_desc.has_value()) { // safety check, did we find a xisec in proximity to the cross-silhouette intersection?
                    if (cross_silhouette_tval > tval_closest) { // cross-silhouette intersection is at the end of the drawn cross curve
                        cross_silhouette_tangent = -cross_silhouette_tangent; // therefor tangents will be pointing in the wrong direction
                        should_reverse_direction = true; // and reversal will be necessary
                    }

                    // compute the tangent orientation along the silhouette curve at the cross-silhouette intersection
                    // that maximizes the angle between it, and the cross curves tangent at the cross-silhouette intersection
                    // sufficient to do this in 2D (see bisector condition and foreshortening)
                    cross_silhouette_tangent.z = 0;
                    silhouette_cross_tangent.z = 0;
                    if (silhouette_cross_tangent * cross_silhouette_tangent > silhouette_cross_tangent * (-cross_silhouette_tangent)) { silhouette_cross_tangent = -silhouette_cross_tangent; }
                }

                // align all subsequent cross-hairs along the cross curve with the silhouettes orientation, flip where necessary
                for (auto flip_isec_desc : boost::make_iterator_range(boost::out_edges(cross_polybez_desc, curve_graph))) {
                    if (flip_isec_desc == silhouette_isec_desc) { continue; } // do not flip the silhouette intersection itself

                    boost::shared_ptr<Intersection> flip_isec_ptr = isec_ptr_map[flip_isec_desc];
                    // get polybezier graph descriptors in the local graph scope
                    boost::optional<curve_graph_vertex_desc_t> opt_flip_polybez_desc_a = find_vertex_by_uuid(flip_isec_ptr->m_polybez_a_uuid, m_curve_network_graph);
                    boost::optional<curve_graph_vertex_desc_t> opt_flip_polybez_desc_b = find_vertex_by_uuid(flip_isec_ptr->m_polybez_b_uuid, m_curve_network_graph);
                    if (opt_flip_polybez_desc_a == boost::none || opt_flip_polybez_desc_b == boost::none) {
                        std::cout << "[ ERROR ] could not find polybezier descriptors!" << std::endl;
                        continue;
                    }
                    curve_graph_vertex_desc_t flip_polybez_desc_a = curve_graph.global_to_local(opt_flip_polybez_desc_a.value());
                    curve_graph_vertex_desc_t flip_polybez_desc_b = curve_graph.global_to_local(opt_flip_polybez_desc_b.value());

                    vector_t cross_isec_tangent = { -100.0, -100.0, -100.0 };; // tangents along cross curve in question at xhairs
                    vector_t isec_cross_tangent = { -100.0, -100.0, -100.0 }; // tangents along incident cross curves intersecting cross curve in question at xhairs

                    // get tangent frames at each xisec along this cross curve
                    if (flip_polybez_desc_a == cross_polybez_desc) { // polybez a must be the cross curve in question
                        if (polybez_ptr_map[flip_polybez_desc_b]->m_curve_type == CURVE_CROSS) { // we only care about xisecs
                            isec_cross_tangent = flip_isec_ptr->m_tangent_b;
                            cross_isec_tangent = flip_isec_ptr->m_tangent_a;
                        }
                    }
                    else if (flip_polybez_desc_b == cross_polybez_desc) { // polybez b must be the cross curve in question
                        if (polybez_ptr_map[flip_polybez_desc_a]->m_curve_type == CURVE_CROSS) { // we only care about xisecs
                            isec_cross_tangent = flip_isec_ptr->m_tangent_a;
                            cross_isec_tangent = flip_isec_ptr->m_tangent_b;
                        }
                    }

                    if (isec_cross_tangent.z < 0) { isec_cross_tangent = -isec_cross_tangent; }
                    if (should_reverse_direction) { cross_isec_tangent = -cross_isec_tangent; }

                    bool same_sign = ((isec_cross_tangent^cross_isec_tangent).z * (silhouette_cross_tangent^cross_silhouette_tangent).z) > 0;
                    if (!same_sign)	{ flip_isec_ptr->flip_tangent_frame(); }
                    flip_isec_ptr->m_is_set_silhouette = true;
                }
            }
        }
    }
}

void Curve_Network::regularize_tangent_frames(curve_graph_t& curve_graph)
{
    curve_graph_edge_data_map_t isec_ptr_map = boost::get(boost::edge_data_pointer, curve_graph);

    for (auto polybez_desc : boost::make_iterator_range(boost::vertices(curve_graph))) {
        QVector<curve_graph_edge_desc_t> sorted_intersections = get_sorted_isecs_along_curve(curve_graph, polybez_desc);
        for (int i = 0; i < sorted_intersections.size() - 1; ++i) {
            boost::shared_ptr<Intersection> current_isec = isec_ptr_map[sorted_intersections[i]];
            boost::shared_ptr<Intersection> next_isec = isec_ptr_map[sorted_intersections[i + 1]];

            vector_t current_tangent_a;
            vector_t current_tangent_b;
            vector_t next_tangent_b;
            vector_t next_tangent_a;

            if (curve_graph.global_to_local(find_vertex_by_uuid(current_isec->m_polybez_a_uuid, curve_graph).value()) == polybez_desc) {
                current_tangent_a = current_isec->m_tangent_a;
                current_tangent_b = current_isec->m_tangent_b;
            }
            else {
                current_tangent_a = current_isec->m_tangent_b;
                current_tangent_b = current_isec->m_tangent_a;
            }
            if (curve_graph.global_to_local(find_vertex_by_uuid(next_isec->m_polybez_a_uuid, curve_graph).value()) == polybez_desc) {
                current_tangent_a = next_isec->m_tangent_a;
                current_tangent_b = next_isec->m_tangent_b;
            }
            else {
                current_tangent_a = next_isec->m_tangent_b;
                current_tangent_b = next_isec->m_tangent_a;
            }
            if ((current_tangent_a.similar_to(next_tangent_a) || current_tangent_a.similar_to(-next_tangent_a)) &&
                    (current_tangent_b.similar_to(next_tangent_b) || current_tangent_b.similar_to(-next_tangent_b))) {
                current_isec->m_tangent_a = next_isec->m_tangent_a;
                current_isec->m_tangent_b = next_isec->m_tangent_b;
                current_isec->calculate_tangent_frame();
            }
        }
    }
}

/*--------------------- INTEGRATE_POLYBEZIER_NORMALS ---
 |
 |  Integrate normals along cross section curves
 |  through stepwise incremental point calculations in
 |  conjunction with minimal tangent frame rotations
 |  between intersections.
 |
 *-----------------------------------------------------*/
void Curve_Network::integrate_polybeziers_normals(curve_graph_t& curve_graph)
{
    // data maps
    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, curve_graph);
    curve_graph_edge_data_map_t isec_ptr_map = boost::get(boost::edge_data_pointer, curve_graph);

    for (auto polybezier_desc : boost::make_iterator_range(boost::vertices(curve_graph))) {
        // ----------------------------------------------------
        // requirements guard clauses
        // curves with no intersections are trivially skipped
        if (boost::out_degree(polybezier_desc, curve_graph) == 0) { continue; }

        // we are only interested in well formed cross curves
        boost::shared_ptr<Polybezier> polybezier_ptr = polybez_ptr_map[polybezier_desc];
        if (polybezier_ptr->m_curve_type != CURVE_CROSS) { continue; }
        if (polybezier_ptr->m_plane_normal.L2_norm3D() <= 0.0) { continue; }

        // ----------------------------------------------------
        // cross curves with no xhairs are trivially skipped
        bool at_least_one_xhair = false;
        for (auto verify_isec_desc : boost::make_iterator_range(boost::out_edges(polybezier_desc, curve_graph))) {
            curve_graph_vertex_desc_t polybez_desc_a = boost::source(verify_isec_desc, curve_graph);
            curve_graph_vertex_desc_t polybez_desc_b = boost::target(verify_isec_desc, curve_graph);
            if (polybezier_desc != polybez_desc_a) {
                if (polybez_ptr_map[polybez_desc_a]->m_curve_type == CURVE_CROSS) {
                    at_least_one_xhair = true;
                }
            }
            else if (polybezier_desc != polybez_desc_b) {
                if (polybez_ptr_map[polybez_desc_b]->m_curve_type == CURVE_CROSS) {
                    at_least_one_xhair = true;
                }
            }
        }
        if (!at_least_one_xhair) { continue; }
        // ----------------------------------------------------

        // ----------------------------------------------------
        // recover global tangent information across segments
        double min_tval =  100.0;
        double max_tval = -100.0;
        double search_tval = 0.0;
        // to recover depth information
        point_t p0 = { -100.0, -100.0, -100.0 };
        for (auto search_isec_desc : boost::make_iterator_range(boost::out_edges(polybezier_desc, curve_graph))) {
            boost::shared_ptr<Intersection> search_isec_ptr = isec_ptr_map[search_isec_desc];
            // get polybezier graph descriptors in the local graph scope
            boost::optional<curve_graph_vertex_desc_t> opt_search_polybez_desc_a = find_vertex_by_uuid(search_isec_ptr->m_polybez_a_uuid, m_curve_network_graph);
            boost::optional<curve_graph_vertex_desc_t> opt_search_polybez_desc_b = find_vertex_by_uuid(search_isec_ptr->m_polybez_b_uuid, m_curve_network_graph);
            if (opt_search_polybez_desc_a == boost::none || opt_search_polybez_desc_b == boost::none) {
                std::cout << "[ ERROR ] could not find polybezier descriptors!" << std::endl;
                continue;
            }
            curve_graph_vertex_desc_t search_polybez_desc_a = curve_graph.global_to_local(opt_search_polybez_desc_a.value());
            curve_graph_vertex_desc_t search_polybez_desc_b = curve_graph.global_to_local(opt_search_polybez_desc_b.value());

            if (search_polybez_desc_a == polybezier_desc) {
                search_tval = search_isec_ptr->m_tval_a;
                if (polybez_ptr_map[search_polybez_desc_b]->m_curve_type == CURVE_CROSS) {
                    p0 = search_isec_ptr->m_position;
                }
            }
            else if (search_polybez_desc_b == polybezier_desc) {
                search_tval = search_isec_ptr->m_tval_b;
                if (polybez_ptr_map[search_polybez_desc_a]->m_curve_type == CURVE_CROSS) {
                    p0 = search_isec_ptr->m_position;
                }
            }

            if (search_tval < min_tval) { min_tval = search_tval; }
            if (search_tval > max_tval) { max_tval = search_tval; }
        }
        // ----------------------------------------------------
        // ----------------------------------------------------

        // clear currently stored normal integrations across this cross curve
        polybezier_ptr->m_integrated_normals.clear();
        polybezier_ptr->m_integrated_normals.resize(0);

        // ----------------------------------------------------
        // begin the process, first get the point where we will
        // calculate a normal marker (integrate)
        for (int b = 0; b < polybezier_ptr->m_bezier_curves.size(); ++b) {
            point_t p_prev = polybezier_ptr->compute_point_at((double)b, 0);
            for (double t = 0.001; t < 1.0; t += 0.001) {
                while ((b + t) < min_tval) {
                    t += 0.001;
                }
                point_t p_curr = polybezier_ptr->compute_point_at(b, t);
                double p_delta = (p_curr - p_prev).L2_norm2D();

                while (p_delta < 0.02 && t <= 1.0) {
                    t += 0.001;
                    p_curr = polybezier_ptr->compute_point_at(b, t);
                    p_delta = (p_curr - p_prev).L2_norm2D();
                }
                if ((b + t) >= max_tval) { continue; }

                p_prev = p_curr;

                p_curr = p_curr - p0;
                p_curr = polybezier_ptr->project_onto_plane(p_curr);
                p_curr = p_curr + p0;

                vector_t p_tan = polybezier_ptr->compute_tangent_at(b, t);
                p_tan = polybezier_ptr->project_onto_plane(p_tan);
                p_tan.normalize();

                double t_global = b + t;
                // ----------------------------------------------------

                // ----------------------------------------------------
                // which two intersections are we between?
                curve_graph_edge_desc_t isec_prev_desc;
                curve_graph_edge_desc_t isec_next_desc;
                vector_t isec_prev_tan;
                vector_t isec_next_tan;
                bool is_isec_prev_xhair;
                bool is_isec_next_xhair;
                double tval_prev = -100.0;
                double tval_next = -100.0;
                double tval_proximity_prev = 100.0;
                double tval_proximity_next = 100.0;
                for (auto between_isec_desc : boost::make_iterator_range(boost::out_edges(polybezier_desc, curve_graph))) {
                    boost::shared_ptr<Intersection> between_isec_ptr = isec_ptr_map[between_isec_desc];
                    // get polybezier graph descriptors in the local graph scope
                    boost::optional<curve_graph_vertex_desc_t> opt_between_polybez_desc_a = find_vertex_by_uuid(between_isec_ptr->m_polybez_a_uuid, m_curve_network_graph);
                    boost::optional<curve_graph_vertex_desc_t> opt_between_polybez_desc_b = find_vertex_by_uuid(between_isec_ptr->m_polybez_b_uuid, m_curve_network_graph);
                    if (opt_between_polybez_desc_a == boost::none || opt_between_polybez_desc_b == boost::none) {
                        std::cout << "[ ERROR ] could not find polybezier descriptors!" << std::endl;
                        continue;
                    }
                    curve_graph_vertex_desc_t between_polybez_desc_a = curve_graph.global_to_local(opt_between_polybez_desc_a.value());
                    curve_graph_vertex_desc_t between_polybez_desc_b = curve_graph.global_to_local(opt_between_polybez_desc_b.value());

                    if (between_polybez_desc_a == polybezier_desc) {
                        if (t_global - between_isec_ptr->m_tval_a < tval_proximity_prev && t_global - between_isec_ptr->m_tval_a > 0) {
                            isec_prev_desc = between_isec_desc;
                            tval_prev = between_isec_ptr->m_tval_a;
                            tval_proximity_prev = t_global - tval_prev;
                            if (polybez_ptr_map[between_polybez_desc_b]->m_curve_type == CURVE_CROSS) {
                                is_isec_prev_xhair = true;
                                isec_prev_tan = between_isec_ptr->m_tangent_a;
                            }
                            else { is_isec_prev_xhair = false; }
                        }
                        else if (between_isec_ptr->m_tval_a - t_global < tval_proximity_next && between_isec_ptr->m_tval_a - t_global > 0) {
                            isec_next_desc = between_isec_desc;
                            tval_next = between_isec_ptr->m_tval_a;
                            tval_proximity_next = tval_next - t_global;
                            if (polybez_ptr_map[between_polybez_desc_b]->m_curve_type == CURVE_CROSS) {
                                is_isec_next_xhair = true;
                                isec_next_tan = between_isec_ptr->m_tangent_a;
                            }
                            else { is_isec_next_xhair = false; }
                        }
                    }
                    else if (between_polybez_desc_b == polybezier_desc) {
                        if (t_global - between_isec_ptr->m_tval_b < tval_proximity_prev && t_global - between_isec_ptr->m_tval_b > 0) {
                            isec_prev_desc = between_isec_desc;
                            tval_prev = between_isec_ptr->m_tval_b;
                            tval_proximity_prev = t_global - tval_prev;
                            if (polybez_ptr_map[between_polybez_desc_a]->m_curve_type == CURVE_CROSS) {
                                is_isec_prev_xhair = true;
                                isec_prev_tan = between_isec_ptr->m_tangent_b;
                            }
                            else { is_isec_prev_xhair = false; }
                        }
                        else if (between_isec_ptr->m_tval_b - t_global < tval_proximity_next && between_isec_ptr->m_tval_b - t_global > 0) {
                            isec_next_desc = between_isec_desc;
                            tval_next = between_isec_ptr->m_tval_b;
                            tval_proximity_next = tval_next - t_global;
                            if (polybez_ptr_map[between_polybez_desc_a]->m_curve_type == CURVE_CROSS) {
                                is_isec_next_xhair = true;
                                isec_next_tan = between_isec_ptr->m_tangent_b;
                            }
                            else { is_isec_next_xhair = false; }
                        }
                    }
                }
                // ----------------------------------------------------

                // ----------------------------------------------------
                // now that we have our two adjoining intersections
                // there are three possible scenarios
                // (1) neither intersection is a xhair (sanity check)
                // (2) only the previous intersection is a xhair
                // (3) only the next intersection is a xhair
                // (4) both intersections are xhairs
                vector_t new_normal = { -100, -100, -100 };
                if (!is_isec_prev_xhair && !is_isec_next_xhair) {
                    continue;
                }
                else if (is_isec_prev_xhair && !is_isec_next_xhair) {
                    vector_t frame_tangent = polybezier_ptr->project_onto_plane(isec_prev_tan);
                    frame_tangent.normalize();
                    new_normal = rotate(frame_tangent, p_tan, isec_ptr_map[isec_prev_desc]->m_normal);
                    new_normal.normalize();
                }
                else if (!is_isec_prev_xhair && is_isec_next_xhair) {
                    vector_t frame_tangent = polybezier_ptr->project_onto_plane(isec_next_tan);
                    frame_tangent.normalize();
                    new_normal = rotate(frame_tangent, p_tan, isec_ptr_map[isec_next_desc]->m_normal);
                    new_normal.normalize();
                }
                else if (is_isec_prev_xhair && is_isec_next_xhair) {
                    // rotate from the previous normal, forward
                    vector_t prev_frame_tangent = polybezier_ptr->project_onto_plane(isec_prev_tan);
                    prev_frame_tangent.normalize();
                    vector_t new_prev_normal = rotate(prev_frame_tangent, p_tan, isec_ptr_map[isec_prev_desc]->m_normal);
                    new_prev_normal.normalize();

                    // rotate from the next normal, backward
                    vector_t next_frame_tangent = polybezier_ptr->project_onto_plane(isec_next_tan);
                    next_frame_tangent.normalize();
                    vector_t new_next_normal = rotate(next_frame_tangent, p_tan, isec_ptr_map[isec_next_desc]->m_normal);
                    new_next_normal.normalize();

                    // take the weighted average of the two rotations
                    double dist_to_prev = sqrt((t_global - tval_prev) * (t_global - tval_prev));
                    double dist_to_next = sqrt((t_global - tval_next) * (t_global - tval_next));
                    vector_t avg_new_normal = new_prev_normal * dist_to_next + new_next_normal * dist_to_prev;
                    new_normal = avg_new_normal;
                    new_normal.normalize();
                }
                Polybezier::normal_marker_s normal_marker = Polybezier::normal_marker_s(t_global, p_curr, new_normal);
                polybezier_ptr->m_integrated_normals.push_back(normal_marker);
            }
        }
    }
}

isec_graph_t Curve_Network::calculate_intersections_network_graph(curve_graph_t& curve_graph)
{
    isec_graph_t isec_graph(0);
    for (auto polybez : boost::make_iterator_range(boost::vertices(curve_graph))) {
        QVector<curve_graph_edge_desc_t> polybez_isecs = get_sorted_isecs_along_curve(curve_graph, polybez);
        if (polybez_isecs.size() < 2) {
            if (polybez_isecs.size() == 1) {
                add_isec_to_isec_network_graph(isec_graph, polybez_isecs.first());
            }
            continue;
        }
        isec_graph_copied_vertex_map_t polybez_desc_map = boost::get(boost::edge_copied_descriptor, isec_graph);
        isec_graph_weight_map_t polybez_weight_map = boost::get(boost::edge_weight, isec_graph);
        for (int i = 0; i < polybez_isecs.size() - 1; ++i) {
            isec_graph_vertex_desc_t isec_desc_a = add_isec_to_isec_network_graph(isec_graph, polybez_isecs[i + 0]);
            isec_graph_vertex_desc_t isec_desc_b = add_isec_to_isec_network_graph(isec_graph, polybez_isecs[i + 1]);
            std::pair<isec_graph_edge_desc_t, bool> curve_desc = boost::add_edge(isec_desc_a, isec_desc_b, isec_graph);
            polybez_desc_map[curve_desc.first] = polybez;
            polybez_weight_map[curve_desc.first] = 1;
        }
    }
    return isec_graph;
}

isec_graph_vertex_desc_t Curve_Network::add_isec_to_isec_network_graph(isec_graph_t& isec_graph, curve_graph_edge_desc_t curve_net_isec_desc)
{
    isec_graph_copied_edge_map_t isec_desc_map = boost::get(boost::vertex_copied_descriptor, isec_graph);
    for (auto isec_desc : boost::make_iterator_range(boost::vertices(isec_graph))) {
        if (isec_desc_map[isec_desc] == curve_net_isec_desc) {
            return isec_desc;
        }
    }
    isec_graph_vertex_desc_t isec_desc = boost::add_vertex(curve_net_isec_desc, isec_graph);
    isec_desc_map[isec_desc] = curve_net_isec_desc;
    return isec_desc;
}

// return a qvector of isecs along a single polycurve sorted by their tvals
QVector<curve_graph_edge_desc_t> Curve_Network::get_sorted_isecs_along_curve(curve_graph_t& curve_graph, curve_graph_vertex_desc_t curve_desc)
{
    curve_graph_out_edge_iter_t isec_iter_first;
    curve_graph_out_edge_iter_t isec_iter_last;

    boost::tie(isec_iter_first, isec_iter_last) = boost::out_edges(curve_desc, curve_graph);

    double curr_tval = -1.0;
    double next_tval = -1.0;
    curve_graph_edge_desc_t curr_isec_desc;
    curve_graph_edge_desc_t next_isec_desc;
    curve_graph_edge_desc_t temp_isec_desc;
    QVector<curve_graph_edge_desc_t> sorted_isec_descs(0);

    for (auto isec_iter = isec_iter_first; isec_iter != isec_iter_last; ++isec_iter) {
        sorted_isec_descs.push_back(*isec_iter);
    }

    curve_graph_edge_data_map_t isec_ptr_map = boost::get(boost::edge_data_pointer, curve_graph);

    for (int i = 0; i < sorted_isec_descs.size(); ++i) {
        for (int j = 0; j < sorted_isec_descs.size() - i - 1; ++j) {
            curr_isec_desc = sorted_isec_descs[j + 0];
            next_isec_desc = sorted_isec_descs[j + 1];
            if (find_vertex_by_uuid(isec_ptr_map[curr_isec_desc]->m_polybez_a_uuid, curve_graph) == curve_desc)
                curr_tval = isec_ptr_map[curr_isec_desc]->m_tval_a;
            else curr_tval = isec_ptr_map[curr_isec_desc]->m_tval_b;

            if (find_vertex_by_uuid(isec_ptr_map[next_isec_desc]->m_polybez_a_uuid, curve_graph) == curve_desc)
                next_tval = isec_ptr_map[next_isec_desc]->m_tval_a;
            else next_tval = isec_ptr_map[next_isec_desc]->m_tval_b;

            // TODO use chord length
            if (curr_tval > next_tval) {
                temp_isec_desc = sorted_isec_descs[j];
                sorted_isec_descs[j] = sorted_isec_descs[j+1];
                sorted_isec_descs[j+1] = temp_isec_desc;
            }
        }
    }
    return sorted_isec_descs;
}

std::vector<curve_graph_filter_t> Curve_Network::filter_connected_components(curve_graph_t& graph)
{
    typedef std::map<curve_graph_vertex_desc_t, unsigned long> curve_mapping_t;
    typedef boost::shared_ptr<curve_mapping_t> curve_component_map_t;

    curve_component_map_t mapping = boost::make_shared<curve_mapping_t>();
    size_t num_connected_components = boost::connected_components(graph, boost::associative_property_map<curve_mapping_t>(*mapping));
    std::vector<curve_graph_filter_t> connected_components_graphs;

    for (size_t i = 0; i < num_connected_components; i++) {
        connected_components_graphs.emplace_back(graph, [=](curve_graph_edge_desc_t e) {
            return mapping->at(source(e, graph)) == i
                    || mapping->at(target(e, graph)) == i;
        },
        [mapping, i](curve_graph_vertex_desc_t v) {
            return mapping->at(v) == i;
        });
    }

    return connected_components_graphs;
}

curve_graph_filter_t Curve_Network::filter_cross_components(curve_graph_t& graph)
{
    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, graph);

    curve_graph_filter_t cross_graph = curve_graph_filter_t(graph, [=](curve_graph_edge_desc_t e) {
        return polybez_ptr_map[source(e, graph)]->m_curve_type == CURVE_CROSS
                && polybez_ptr_map[target(e, graph)]->m_curve_type == CURVE_CROSS;
    },
    [graph, polybez_ptr_map] (curve_graph_vertex_desc_t v) {
        return polybez_ptr_map[v]->m_curve_type == CURVE_CROSS;
    });

    return cross_graph;
}

void Curve_Network::run_nlopt_solve(double &time_tracker, int &nxhairs, int &nxcurves, int &nxpatches, int &npatches)
{
    // time tracker
    double total_time = 0.0;
    int num_xhairs = 0;
    int num_xcurves = 0;
    int num_xpatches = 0;
    // ===========================================================
    // run nlopt on every connected component (patch) in the graph
    for (int c = 0; c < m_curve_network_subgraphs.size(); ++c) {
        // --------------------------------------------------------------------
        // extract the cross curves from the patches, store in cross component,
        // and calculate (derive) the cross intersection graph using the cross
        // filtered component graph
        curve_graph_t component = m_curve_network_subgraphs[c];
        curve_graph_filter_t filtered_cross = filter_cross_components(component);
        curve_graph_t cross_component = component.create_subgraph();
        for (auto polybez : boost::make_iterator_range(boost::vertices(filtered_cross))) {
            boost::add_vertex(component.local_to_global(polybez), cross_component);
        }

        num_xhairs += boost::num_edges(cross_component);
        num_xcurves += boost::num_vertices(cross_component);

        // -------------------------------------------------
        // data maps for later
        curve_graph_vertex_index_map_t cross_polybezier_index_map = boost::get(boost::vertex_index, cross_component);
        curve_graph_vertex_data_map_t cross_polybezier_data_map = boost::get(boost::vertex_data_pointer, cross_component);
        curve_graph_edge_data_map_t cross_intersection_data_map = boost::get(boost::edge_data_pointer, cross_component);
        // -------------------------------------------------

        // patches with only one cross intersection are trivially solved
        if (boost::num_edges(cross_component) == 1) {
            num_xpatches += 1;

            curve_graph_vertex_desc_t polybez_desc_a = boost::vertex(0, cross_component);
            curve_graph_vertex_desc_t polybez_desc_b = boost::vertex(1, cross_component);
            curve_graph_edge_desc_t isec_desc = boost::edge(0, 1, cross_component).first;

            boost::shared_ptr<Polybezier> polybezier_a_ptr = cross_polybezier_data_map[polybez_desc_a];
            boost::shared_ptr<Polybezier> polybezier_b_ptr = cross_polybezier_data_map[polybez_desc_b];
            boost::shared_ptr<Intersection> isec_ptr = cross_intersection_data_map[isec_desc];

            if (isec_ptr->m_polybez_a_uuid == polybezier_a_ptr->m_uuid) {
                polybezier_a_ptr->m_plane_normal = isec_ptr->m_local_plane_a;
                polybezier_b_ptr->m_plane_normal = isec_ptr->m_local_plane_b;
            }
            else {
                polybezier_a_ptr->m_plane_normal = isec_ptr->m_local_plane_b;
                polybezier_b_ptr->m_plane_normal = isec_ptr->m_local_plane_a;
            }

            integrate_polybeziers_normals(component);
            continue;
        }

        std::vector<curve_graph_filter_t> filtered_cross_sub_components = filter_connected_components(cross_component);
        std::vector<curve_graph_t> cross_sub_components;
        cross_sub_components.resize(0);
        for (size_t fsc = 0; fsc < filtered_cross_sub_components.size(); ++fsc) {
            curve_graph_filter_t filtered_cross_sub_component = filtered_cross_sub_components[fsc];
            curve_graph_t cross_component_subgraph = cross_component.create_subgraph();
            for (auto polybez : boost::make_iterator_range(boost::vertices(filtered_cross_sub_component))) {
                boost::add_vertex(cross_component.local_to_global(polybez), cross_component_subgraph);
            }
            cross_sub_components.push_back(cross_component_subgraph);
        }

        num_xpatches += (int)cross_sub_components.size();

        for (int sc = 0; sc < (int)cross_sub_components.size(); ++sc) {
            curve_graph_t cross_sub_component = cross_sub_components[sc];

            // -------------------------------------------------
            // data maps for later
            curve_graph_vertex_index_map_t cross_polybezier_index_map = boost::get(boost::vertex_index, cross_sub_component);
            curve_graph_vertex_data_map_t cross_polybezier_data_map = boost::get(boost::vertex_data_pointer, cross_sub_component);
            curve_graph_edge_data_map_t cross_intersection_data_map = boost::get(boost::edge_data_pointer, cross_sub_component);
            // -------------------------------------------------

            // patches with only one cross intersection are trivially solved
            if (boost::num_edges(cross_sub_component) == 1) {
                curve_graph_vertex_desc_t polybez_desc_a = boost::vertex(0, cross_sub_component);
                curve_graph_vertex_desc_t polybez_desc_b = boost::vertex(1, cross_sub_component);
                curve_graph_edge_desc_t isec_desc = boost::edge(0, 1, cross_sub_component).first;

                boost::shared_ptr<Polybezier> polybezier_a_ptr = cross_polybezier_data_map[polybez_desc_a];
                boost::shared_ptr<Polybezier> polybezier_b_ptr = cross_polybezier_data_map[polybez_desc_b];
                boost::shared_ptr<Intersection> isec_ptr = cross_intersection_data_map[isec_desc];

                if (isec_ptr->m_polybez_a_uuid == polybezier_a_ptr->m_uuid) {
                    polybezier_a_ptr->m_plane_normal = isec_ptr->m_local_plane_a;
                    polybezier_b_ptr->m_plane_normal = isec_ptr->m_local_plane_b;
                }
                else {
                    polybezier_a_ptr->m_plane_normal = isec_ptr->m_local_plane_b;
                    polybezier_b_ptr->m_plane_normal = isec_ptr->m_local_plane_a;
                }

                integrate_polybeziers_normals(component);
                continue;
            }

            // -------------------------------------------------
            // do we need to quad guess?
            bool should_guess = true;
            int num_sil_sets = 0;
            for (auto isec_desc : boost::make_iterator_range(boost::edges(cross_sub_component))) {
                if (cross_intersection_data_map[isec_desc]->m_is_set_silhouette == true) {
                    should_guess = false;
                    num_sil_sets++;
                }
            }
            // -------------------------------------------------

            if (should_guess) {
                Nlopt_Solve::optimization_data_s optimization_datas[4]; // four initializations
                QPair<double, std::vector<double>> solutions[4]; // four different solutions
                isec_graph_t isec_network_graph = calculate_intersections_network_graph(cross_sub_component);
                isec_graph_copied_edge_map_t isec_desc_map = boost::get(boost::vertex_copied_descriptor, isec_network_graph);
                // --------------------------------------------------------------------

                // -------------------------------------------------------
                // calculate graph diameter of the (cross component) patch
                // by using djikstra and finding the two furthest cross
                // intersections in the derived cross intersection graph
                int patch_max_eccentricity = 0;
                QPair<curve_graph_edge_desc_t, curve_graph_edge_desc_t> patch_diameter_curves;

                for (isec_graph_vertex_desc_t isec : boost::make_iterator_range(boost::vertices(isec_network_graph))) {
                    std::vector<isec_graph_vertex_desc_t> predecessors(boost::num_vertices(isec_network_graph));
                    std::vector<int> distances(boost::num_vertices(isec_network_graph));

                    boost::dijkstra_shortest_paths(isec_network_graph, isec, boost::distance_map(boost::make_iterator_property_map(distances.begin(), boost::get(boost::vertex_index, isec_network_graph))));
                    isec_graph_vertex_desc_t isec_b = std::max_element(distances.begin(), distances.end()) - distances.begin();
                    if (*std::max_element(distances.begin(), distances.end()) > patch_max_eccentricity) {
                        patch_max_eccentricity = *std::max_element(distances.begin(), distances.end());
                        patch_diameter_curves.first = isec_desc_map[isec];
                        patch_diameter_curves.second = isec_desc_map[isec_b];
                    }
                }
                // -------------------------------------------------------

                // ---------------------------------------------------
                // prepare the data for the nlopt solver, nlopt should
                // not have to know about the graph structure
                boost::shared_ptr<Intersection> isec_a = cross_intersection_data_map[patch_diameter_curves.first];
                boost::shared_ptr<Intersection> isec_b = cross_intersection_data_map[patch_diameter_curves.second];

                for (int g = 0; g < 4; ++g) {
                    switch (g) {
                    case 0: { break; } // no flips
                    case 1: { // flip only the first isec
                        isec_a->flip_tangent_frame();
                        break;
                    }
                    case 2: { // flip the first isec back and flip the second one
                        isec_a->flip_tangent_frame();
                        isec_b->flip_tangent_frame();
                        break;
                    }
                    case 3: { // flip both isecs
                        isec_a->flip_tangent_frame();
                        break;
                    }
                    };

                    curve_graph_vertex_desc_t polybez_desc_a = find_vertex_by_uuid(isec_a->m_polybez_a_uuid, cross_sub_component).value();
                    curve_graph_vertex_desc_t polybez_desc_b = find_vertex_by_uuid(isec_a->m_polybez_b_uuid, cross_sub_component).value();
                    curve_graph_vertex_desc_t polybez_desc_c = find_vertex_by_uuid(isec_b->m_polybez_a_uuid, cross_sub_component).value();
                    curve_graph_vertex_desc_t polybez_desc_d = find_vertex_by_uuid(isec_b->m_polybez_b_uuid, cross_sub_component).value();

                    vector_t isec_a_plane_normal_a = isec_a->m_local_plane_a;
                    vector_t isec_a_plane_normal_b = isec_a->m_local_plane_b;
                    vector_t isec_b_plane_normal_a = isec_b->m_local_plane_a;
                    vector_t isec_b_plane_normal_b = isec_b->m_local_plane_b;

                    if (isec_a_plane_normal_a.z < 0.0) {
                        isec_a_plane_normal_a = -isec_a_plane_normal_a;
                    }
                    if (isec_a_plane_normal_b.z < 0.0) {
                        isec_a_plane_normal_b = -isec_a_plane_normal_b;
                    }
                    if (isec_b_plane_normal_a.z < 0.0) {
                        isec_b_plane_normal_a = -isec_b_plane_normal_a;
                    }
                    if (isec_b_plane_normal_b.z < 0.0) {
                        isec_b_plane_normal_b = -isec_b_plane_normal_b;
                    }
                    isec_a_plane_normal_a = isec_a_plane_normal_a/isec_a_plane_normal_a.z;
                    isec_a_plane_normal_b = isec_a_plane_normal_b/isec_a_plane_normal_b.z;
                    isec_b_plane_normal_a = isec_b_plane_normal_a/isec_b_plane_normal_a.z;
                    isec_b_plane_normal_b = isec_b_plane_normal_b/isec_b_plane_normal_b.z;

                    optimization_datas[g].initial_guess.resize(boost::num_vertices(cross_sub_component) * 3, 0.0);
                    optimization_datas[g].initial_guess[cross_polybezier_index_map[polybez_desc_a] * 3 + 0] = isec_a_plane_normal_a.x;
                    optimization_datas[g].initial_guess[cross_polybezier_index_map[polybez_desc_a] * 3 + 1] = isec_a_plane_normal_a.y;
                    optimization_datas[g].initial_guess[cross_polybezier_index_map[polybez_desc_b] * 3 + 0] = isec_a_plane_normal_b.x;
                    optimization_datas[g].initial_guess[cross_polybezier_index_map[polybez_desc_b] * 3 + 1] = isec_a_plane_normal_b.y;
                    optimization_datas[g].initial_guess[cross_polybezier_index_map[polybez_desc_c] * 3 + 0] = isec_b_plane_normal_a.x;
                    optimization_datas[g].initial_guess[cross_polybezier_index_map[polybez_desc_c] * 3 + 1] = isec_b_plane_normal_a.y;
                    optimization_datas[g].initial_guess[cross_polybezier_index_map[polybez_desc_d] * 3 + 0] = isec_b_plane_normal_b.x;
                    optimization_datas[g].initial_guess[cross_polybezier_index_map[polybez_desc_d] * 3 + 1] = isec_b_plane_normal_b.y;

                    std::map<boost::uuids::uuid, int> map_polybezier_uuid_to_index;
                    for (auto cross_polybezier_description : boost::make_iterator_range(boost::vertices(cross_sub_component))) {
                        map_polybezier_uuid_to_index.insert({cross_polybezier_data_map[cross_polybezier_description]->m_uuid,
                                                             cross_polybezier_index_map[cross_polybezier_description]});
                    }
                    optimization_datas[g].intersection_datas.resize(0);
                    for (auto cross_intersection_descripter : boost::make_iterator_range(boost::edges(cross_sub_component))) {
                        boost::shared_ptr<Intersection> read_intersection_ptr = cross_intersection_data_map[cross_intersection_descripter];
                        Nlopt_Solve::intersection_data_s write_intersection_data(map_polybezier_uuid_to_index[read_intersection_ptr->m_polybez_a_uuid],
                                map_polybezier_uuid_to_index[read_intersection_ptr->m_polybez_b_uuid],
                                read_intersection_ptr->m_tangent_a,
                                read_intersection_ptr->m_tangent_b,
                                read_intersection_ptr->m_position);
                        optimization_datas[g].intersection_datas.push_back(write_intersection_data);
                    }
                    // record start time
                    auto start = std::chrono::high_resolution_clock::now();
                    solutions[g] = Nlopt_Solve::solve(optimization_datas[g]);
                    // record end time
                    auto finish = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(finish - start);
                    total_time += time_span.count();
                }
                // ===========================================================

                // throw out poorly fitted solutions
                // keep well fitted solutions
                double best_min = 1000000.0;
                for (int j = 0; j < 4; ++j) {
                    if (solutions[j].first < best_min) { best_min = solutions[j].first; }
                }

                // ===========================================================
                // load nlopt solution back into the curve network
                // updating subgraph cross component updates the primary graph automagically
                // select the best solution wrt y-integrals
                int best_solution = 0;
                double best_metric = -100.0;
                for (int s = 0; s < 4; ++s) {
                    //if (solutions[s].first > best_min) { continue; }
                    std::vector<double> solution = solutions[s].second;
                    for (auto polybez_desc : boost::make_iterator_range(boost::vertices(cross_sub_component))) {
                        vector_t plane_normal = vector_t(solution[boost::get(boost::vertex_index, cross_sub_component, polybez_desc) * 3 + 0],
                                solution[boost::get(boost::vertex_index, cross_sub_component, polybez_desc) * 3 + 1],
                                1.0);
                        plane_normal.normalize();
                        cross_polybezier_data_map[polybez_desc]->m_plane_normal = plane_normal;
                        cross_polybezier_data_map[polybez_desc]->m_plane_offset = solution[boost::get(boost::vertex_index, cross_sub_component, polybez_desc) * 3 + 2];
                    }

                    for (auto intersection_descriptor : boost::make_iterator_range(boost::edges(cross_sub_component))) {

                        vector_t tangent_a = cross_intersection_data_map[intersection_descriptor]->m_tangent_a;
                        vector_t tangent_b = cross_intersection_data_map[intersection_descriptor]->m_tangent_b;
                        vector_t plane_a = cross_polybezier_data_map[find_vertex_by_uuid(cross_intersection_data_map[intersection_descriptor]->m_polybez_a_uuid, cross_sub_component).value()]->m_plane_normal;
                        vector_t plane_b = cross_polybezier_data_map[find_vertex_by_uuid(cross_intersection_data_map[intersection_descriptor]->m_polybez_b_uuid, cross_sub_component).value()]->m_plane_normal;

                        cross_intersection_data_map[intersection_descriptor]->m_tangent_a = project_into_plane(tangent_a, plane_a);
                        cross_intersection_data_map[intersection_descriptor]->m_tangent_b = project_into_plane(tangent_b, plane_b);
                        cross_intersection_data_map[intersection_descriptor]->m_tangent_a.normalize();
                        cross_intersection_data_map[intersection_descriptor]->m_tangent_b.normalize();

                        vector_t normal = cross_intersection_data_map[intersection_descriptor]->m_tangent_a ^ cross_intersection_data_map[intersection_descriptor]->m_tangent_b;
                        normal.normalize();

                        if (normal.z < 0.0) {
                            normal = -normal;
                        }
                        cross_intersection_data_map[intersection_descriptor]->m_normal = normal;
                    }
                    integrate_polybeziers_normals(component);

                    // integrate
                    double metric_x = 0.0;
                    double metric_y = 0.0;
                    for (auto polybez_desc : boost::make_iterator_range(boost::vertices(cross_sub_component))) {
                        QVector<Polybezier::normal_marker_s> integrated_normals = cross_polybezier_data_map[polybez_desc]->m_integrated_normals;
                        for (int i = 0; i < integrated_normals.size(); ++i) {
                            metric_x += integrated_normals[i].normal.x;
                            metric_y += integrated_normals[i].normal.y;
                        }
                    }
                    if (metric_y > best_metric) {
                        best_metric = metric_y;
                        best_solution = s;
                    }
                }

                // ===========================================================
                // apply best solution found
                std::vector<double> solution = solutions[best_solution].second;
                for (auto polybez_desc : boost::make_iterator_range(boost::vertices(cross_sub_component))) {
                    vector_t plane_normal = vector_t(solution[boost::get(boost::vertex_index, cross_sub_component, polybez_desc) * 3 + 0],
                            solution[boost::get(boost::vertex_index, cross_sub_component, polybez_desc) * 3 + 1],
                            1.0);
                    plane_normal.normalize();
                    cross_polybezier_data_map[polybez_desc]->m_plane_normal = plane_normal;
                    cross_polybezier_data_map[polybez_desc]->m_plane_offset = solution[boost::get(boost::vertex_index, cross_sub_component, polybez_desc) * 3 + 2];
                }

                for (auto intersection_descriptor : boost::make_iterator_range(boost::edges(cross_sub_component))) {

                    vector_t tangent_a = cross_intersection_data_map[intersection_descriptor]->m_tangent_a;
                    vector_t tangent_b = cross_intersection_data_map[intersection_descriptor]->m_tangent_b;
                    vector_t plane_a = cross_polybezier_data_map[find_vertex_by_uuid(cross_intersection_data_map[intersection_descriptor]->m_polybez_a_uuid, cross_sub_component).value()]->m_plane_normal;
                    vector_t plane_b = cross_polybezier_data_map[find_vertex_by_uuid(cross_intersection_data_map[intersection_descriptor]->m_polybez_b_uuid, cross_sub_component).value()]->m_plane_normal;

                    cross_intersection_data_map[intersection_descriptor]->m_tangent_a = project_into_plane(tangent_a, plane_a);
                    cross_intersection_data_map[intersection_descriptor]->m_tangent_b = project_into_plane(tangent_b, plane_b);
                    cross_intersection_data_map[intersection_descriptor]->m_tangent_a.normalize();
                    cross_intersection_data_map[intersection_descriptor]->m_tangent_b.normalize();

                    vector_t normal = cross_intersection_data_map[intersection_descriptor]->m_tangent_a ^ cross_intersection_data_map[intersection_descriptor]->m_tangent_b;
                    normal.normalize();

                    if (normal.z < 0.0) {
                        normal = -normal;
                    }
                    cross_intersection_data_map[intersection_descriptor]->m_normal = normal;
                }
                integrate_polybeziers_normals(component);
                // ===========================================================
            }
            else { // without guessing
                Nlopt_Solve::optimization_data_s optimization_data; // one initialization
                std::vector<double> solution; // one solution

                optimization_data.initial_guess.resize(boost::num_vertices(cross_sub_component) * 3, 0.0);

                for (auto isec_desc : boost::make_iterator_range(boost::edges(cross_sub_component))) {
                    if (cross_intersection_data_map[isec_desc]->m_is_set_silhouette != true) { continue; }
                    boost::shared_ptr<Intersection> isec_ptr = cross_intersection_data_map[isec_desc];

                    curve_graph_vertex_desc_t polybez_desc_a = find_vertex_by_uuid(isec_ptr->m_polybez_a_uuid, cross_sub_component).value();
                    curve_graph_vertex_desc_t polybez_desc_b = find_vertex_by_uuid(isec_ptr->m_polybez_b_uuid, cross_sub_component).value();

                    vector_t plane_normal_a = isec_ptr->m_local_plane_a;
                    vector_t plane_normal_b = isec_ptr->m_local_plane_b;

                    if (plane_normal_a.z < 0.0) {
                        plane_normal_a = -plane_normal_a;
                    }
                    if (plane_normal_b.z < 0.0) {
                        plane_normal_b = -plane_normal_b;
                    }

                    plane_normal_a = plane_normal_a/plane_normal_a.z;
                    plane_normal_b = plane_normal_b/plane_normal_b.z;

                    optimization_data.initial_guess[cross_polybezier_index_map[polybez_desc_a] * 3 + 0] = plane_normal_a.x;
                    optimization_data.initial_guess[cross_polybezier_index_map[polybez_desc_a] * 3 + 1] = plane_normal_a.y;
                    optimization_data.initial_guess[cross_polybezier_index_map[polybez_desc_b] * 3 + 0] = plane_normal_b.x;
                    optimization_data.initial_guess[cross_polybezier_index_map[polybez_desc_b] * 3 + 1] = plane_normal_b.y;
                }

                std::map<boost::uuids::uuid, int> map_polybezier_uuid_to_index;
                for (auto cross_polybezier_description : boost::make_iterator_range(boost::vertices(cross_sub_component))) {
                    map_polybezier_uuid_to_index.insert({cross_polybezier_data_map[cross_polybezier_description]->m_uuid,
                                                         cross_polybezier_index_map[cross_polybezier_description]});
                }

                optimization_data.intersection_datas.resize(0);
                for (auto cross_intersection_descripter : boost::make_iterator_range(boost::edges(cross_sub_component))) {
                    boost::shared_ptr<Intersection> read_intersection_ptr = cross_intersection_data_map[cross_intersection_descripter];
                    Nlopt_Solve::intersection_data_s write_intersection_data(map_polybezier_uuid_to_index[read_intersection_ptr->m_polybez_a_uuid],
                            map_polybezier_uuid_to_index[read_intersection_ptr->m_polybez_b_uuid],
                            read_intersection_ptr->m_tangent_a,
                            read_intersection_ptr->m_tangent_b,
                            read_intersection_ptr->m_position);
                    optimization_data.intersection_datas.push_back(write_intersection_data);
                }

                optimization_data.silhouette_set = true;

                // record start time
                auto start = std::chrono::high_resolution_clock::now();
                solution = Nlopt_Solve::solve(optimization_data).second;
                // record end time
                auto finish = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(finish - start);
                total_time += time_span.count();

                // ===========================================================

                // ===========================================================
                // load nlopt solution back into the curve network
                // updating subgraph cross component updates the primary graph automagically
                for (auto polybez_desc : boost::make_iterator_range(boost::vertices(cross_sub_component))) {
                    vector_t plane_normal = vector_t(solution[boost::get(boost::vertex_index, cross_sub_component, polybez_desc) * 3 + 0],
                            solution[boost::get(boost::vertex_index, cross_sub_component, polybez_desc) * 3 + 1],
                            1.0);
                    plane_normal.normalize();
                    cross_polybezier_data_map[polybez_desc]->m_plane_normal = plane_normal;
                    cross_polybezier_data_map[polybez_desc]->m_plane_offset = solution[boost::get(boost::vertex_index, cross_sub_component, polybez_desc) * 3 + 2];
                }

                for (auto intersection_descriptor : boost::make_iterator_range(boost::edges(cross_sub_component))) {

                    vector_t tangent_a = cross_intersection_data_map[intersection_descriptor]->m_tangent_a;
                    vector_t tangent_b = cross_intersection_data_map[intersection_descriptor]->m_tangent_b;
                    vector_t plane_a = cross_polybezier_data_map[find_vertex_by_uuid(cross_intersection_data_map[intersection_descriptor]->m_polybez_a_uuid, cross_sub_component).value()]->m_plane_normal;
                    vector_t plane_b = cross_polybezier_data_map[find_vertex_by_uuid(cross_intersection_data_map[intersection_descriptor]->m_polybez_b_uuid, cross_sub_component).value()]->m_plane_normal;

                    cross_intersection_data_map[intersection_descriptor]->m_tangent_a = project_into_plane(tangent_a, plane_a);
                    cross_intersection_data_map[intersection_descriptor]->m_tangent_b = project_into_plane(tangent_b, plane_b);
                    cross_intersection_data_map[intersection_descriptor]->m_tangent_a.normalize();
                    cross_intersection_data_map[intersection_descriptor]->m_tangent_b.normalize();

                    vector_t normal = cross_intersection_data_map[intersection_descriptor]->m_tangent_a ^ cross_intersection_data_map[intersection_descriptor]->m_tangent_b;
                    normal.normalize();

                    if (normal.z < 0.0) {
                        normal = -normal;
                    }
                    cross_intersection_data_map[intersection_descriptor]->m_normal = normal;
                }
                integrate_polybeziers_normals(component);
            }
        }
    }
    if (time_tracker != -1.0) { time_tracker += total_time; }
    if (nxhairs != -1.0) { nxhairs += num_xhairs; }
    if (nxcurves != -1.0) { nxcurves += num_xcurves; }
    if (nxpatches != -1.0) { nxpatches += num_xpatches; }
    if (npatches != -1.0) { npatches += m_curve_network_subgraphs.size(); }
}
