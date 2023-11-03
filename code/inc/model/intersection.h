#ifndef INTERSECTION_H
#define INTERSECTION_H

#include "utils/defs/util_data_defs.h"
#include "utils/defs/util_graph_typedefs.h"
#include "utils/math/util_math_defs.h"

// ---------------------------------
// QT incompatability with boost fix
// ---------------------------------
#ifndef Q_MOC_RUN
#include <boost/optional.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#endif

#include <QVector>

class Intersection
{
public:
    // -----------------------------------
    // INTERSECTION CONSTRUCTOR/DESTRUCTOR
    // -----------------------------------
    Intersection(boost::uuids::uuid polybez_a_uuid,
                 boost::uuids::uuid polybez_b_uuid,
                 boost::optional<bool> is_xhair = boost::none,
                 boost::optional<bool> is_set_user = boost::none,
                 boost::optional<bool> is_set_silhouette = boost::none,
                 boost::optional<point_t>  position = boost::none,
                 boost::optional<double> tval_a = boost::none,
                 boost::optional<double> tval_b = boost::none,
                 boost::optional<vector_t> normal = boost::none,
                 boost::optional<vector_t> tangent_a = boost::none,
                 boost::optional<vector_t> tangent_b = boost::none,
                 boost::optional<vector_t> local_plane_a = boost::none,
                 boost::optional<vector_t> local_plane_b = boost::none);

    ~Intersection();

    // ---------------------------
    // INTERSECTION IDENTIFICATION
    // ---------------------------
    boost::uuids::uuid  m_uuid;
    boost::uuids::uuid  m_polybez_a_uuid;
    boost::uuids::uuid  m_polybez_b_uuid;

    // -----------------
    // INTERSECTION TAGS
    // -----------------
    bool m_is_xhair;
    bool m_is_set_user;
    bool m_is_set_silhouette;
    int  m_label;

    // -----------------
    // INTERSECTION DATA
    // -----------------
    point_t  m_position;
    double m_tval_a;
    double m_tval_b;
    vector_t m_tangent_a;
    vector_t m_tangent_b;
    vector_t m_local_plane_a;
    vector_t m_local_plane_b;
    vector_t m_normal;

    // --------------------------
    // TANGENT FRAME MANIPULATION
    // --------------------------
    void calculate_tangent_frame();
    void calculate_local_planes();
    void flip_tangent_frame();
};


#endif
