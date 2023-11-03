#ifndef BEZIER_POLYCURVE_H
#define BEZIER_POLYCURVE_H

#include "utils/math/util_math_defs.h"
#include "utils/defs/util_data_defs.h"

// ---------------------------------
// QT incompatability with boost fix
// ---------------------------------
#ifndef Q_MOC_RUN
#include <boost/optional.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#endif

#include <QVector>

// -----------------
// BEZIER CURVE DATA
// -----------------
struct bezier_curve_s {
    point_t ctrl_points[4];
};

class Polybezier
{
public:
    // ---------------------------------
    // POLYBEZIER CONSTRUCTOR/DESTRUCTOR
    // ---------------------------------
    Polybezier(boost::optional<bool> is_closed_curve = boost::none,
               boost::optional<curve_type_t> curve_type = boost::none,
               boost::optional<QVector<bezier_curve_s>> bezier_curves = boost::none,
               boost::optional<vector_t> plane_normal = boost::none,
               boost::optional<double> plane_offset = boost::none);

    ~Polybezier();

    // -------------------------
    // POLYBEZIER IDENTIFICATION
    // -------------------------
    boost::uuids::uuid  m_uuid;

    // ----------------
    // POLYBEZIER FLAGS
    // ----------------
    bool m_is_closed_curve; // TODO is this actually needed, it's trivial to determine closed curve from control points..
    curve_type_t m_curve_type;

    // ---------------------
    // POLYBEZIER CURVE DATA
    // ---------------------
    QVector<bezier_curve_s> m_bezier_curves;
    vector_t m_plane_normal;
    double m_plane_offset;

    // -------------------------------------
    // NORMAL MARKERS ALONG POLYBEZIER CURVE
    // -------------------------------------
    typedef struct normal_marker_s {
         double tval;
         point_t position;
         vector_t normal;

         normal_marker_s() {
             this->tval = 0.0;
             this->position = { -100.0, -100.0, -100.0 };
             this->normal = { -100.0, -100.0, -100.0 };
         }
         normal_marker_s(double tval, point_t position, vector_t normal) {
             this->tval = tval;
             this->position = position;
             this->normal = normal;
         }
    } normal_marker_s;

    QVector<normal_marker_s> m_integrated_normals;

    // -------------------------------------
    // POLYBEZIER CURVE ELEMENT COMPUTATIONS
    // -------------------------------------
    point_t compute_point_at(int segment, double t);
    vector_t compute_tangent_at(int segment, double t);
    vector_t compute_normal_at(int segment, double t);
    vector_t project_onto_plane(const vector_t& v);
    point_t project_onto_plane(const point_t& p);

};

#endif // BEZIER_POLYCURVE_H
