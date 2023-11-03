#ifndef CURVE_NETWORK_H
#define CURVE_NETWORK_H

#include "model/polybezier.h"
#include "model/intersection.h"

#include "solver/nlopt_solve.h"

#include "utils/defs/util_data_defs.h"

#include "utils/defs/util_graph_typedefs.h"
#include "utils/math/util_math_defs.h"
#include "utils/util_filepath_defs.h"

#include <chrono>
#include <ratio>
#include <ctime>
#include <map>

#include <QWidget>
#include <QVector>

class Curve_Network  : public QObject
{
    Q_OBJECT

public slots:

signals:
    void request_snapshot(QString, bool, bool, bool);

public:

    Curve_Network();
    ~Curve_Network();

    // >> HELPERS
    // --------------------------------------------------------------------------------------
    // curve graph helper functions for finding both vertices and edges using their IDs
    static boost::optional<curve_graph_vertex_desc_t> find_vertex_by_uuid(boost::uuids::uuid uuid, curve_graph_t graph);
    static boost::optional<curve_graph_edge_desc_t> find_edge_by_uuid(boost::uuids::uuid uuid, curve_graph_t graph);
    // --------------------------------------------------------------------------------------

    void add_polybezier(QVector<point_t> ctrl_points, curve_type_t curve_type, bool closed_curve);
    void update_connected_components();

    QVector<isec_ptr_t> intersect_polybeziers(curve_graph_vertex_desc_t curve_desc_a, curve_graph_vertex_desc_t curve_desc_b);
    void update_intersections(curve_graph_vertex_desc_t curve_descriptor);
    void update_intersection_orientations_from_silhouettes(curve_graph_t& curve_graph);

    isec_graph_t calculate_intersections_network_graph(curve_graph_t &curve_graph);
    void regularize_tangent_frames(curve_graph_t &curve_graph);

    isec_graph_vertex_desc_t add_isec_to_isec_network_graph(isec_graph_t& isec_graph, curve_graph_edge_desc_t curve_net_isec_desc);
    QVector<curve_graph_edge_desc_t> get_sorted_isecs_along_curve(curve_graph_t &curve_graph, curve_graph_vertex_desc_t curve_desc);

    void integrate_polybeziers_normals(curve_graph_t &curve_graph);

    std::vector<curve_graph_filter_t> filter_connected_components(curve_graph_t& graph);
    curve_graph_filter_t filter_cross_components(curve_graph_t& graph);

    void run_nlopt_solve(double &time_tracker, int &nxhairs, int &nxcurves, int &nxpatches, int &npatches);
    void run_analytic_solve() {} // TODO
    void run_graph_cut() {} // TODO

    curve_graph_t m_curve_network_graph;
    QVector<curve_graph_t> m_curve_network_subgraphs;

private:


};

#endif
