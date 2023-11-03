#include "model/polybezier.h"

/*--------------------------------- BEZIER_POLYCURVE ---
 |
 |  Constructor for a polybez curve using boost optional
 |  and boost none to capture member variables that
 |  we may not know at initialization.
 |
 *-----------------------------------------------------*/
Polybezier::Polybezier(boost::optional<bool> is_closed_curve,
                                   boost::optional<curve_type_t> curve_type,
                                   boost::optional<QVector<bezier_curve_s>> bezier_curves,
                                   boost::optional<vector_t> plane_normal,
                                   boost::optional<double> plane_offset)
    : m_uuid(boost::uuids::random_generator()())
{
    m_is_closed_curve = is_closed_curve.value_or(false);
    m_curve_type = curve_type.value_or(CURVE_CROSS);

    m_bezier_curves = bezier_curves.value_or(QVector<bezier_curve_s>(0));
    m_plane_normal = plane_normal.value_or(vector_t(0.0, 0.0, 0.0));
    m_plane_offset = plane_offset.value_or(0.0);
}

/*-------------------------------- ~BEZIER_POLYCURVE ---
 |
 |  Polybez curve destructor.
 |
 *-----------------------------------------------------*/
Polybezier::~Polybezier() { }

/*--------------------------------- COMPUTE_*_ON_AT ---
 |
 |  Compute and return a point, tangent, or normal on
 |  the specified segment at the specified t value.
 |
 *-----------------------------------------------------*/
point_t  Polybezier::compute_point_at(int segment, double t)
{
    bezier_curve_s bezier_curve = m_bezier_curves[segment];

    point_t P;
    double x, y;
    double A, B, C;
    point_t P0, P1, P2, P3;
    P0 = bezier_curve.ctrl_points[0];
    P1 = bezier_curve.ctrl_points[1];
    P2 = bezier_curve.ctrl_points[2];
    P3 = bezier_curve.ctrl_points[3];

    C = 3 * (P1.x - P0.x);
    B = 3 * (P2.x - P1.x) - C;
    A = P3.x - P0.x - C - B;
    x = (A * t * t * t) + (B * t * t) + (C * t) + P0.x;

    C = 3 * (P1.y - P0.y);
    B = 3 * (P2.y - P1.y) - C;
    A = P3.y - P0.y - C - B;
    y = ( A * t * t * t ) + ( B * t * t ) + ( C * t ) + P0.y;

    P.x = x;
    P.y = y;
    P.z = 0.0;

    return P;
}

vector_t Polybezier::compute_tangent_at(int segment, double t)
{
    bezier_curve_s bezier_curve = m_bezier_curves[segment];

    vector_t T;
    double x, y;
    double A, B, C;
    point_t P0, P1, P2, P3;
    P0 = bezier_curve.ctrl_points[0];
    P1 = bezier_curve.ctrl_points[1];
    P2 = bezier_curve.ctrl_points[2];
    P3 = bezier_curve.ctrl_points[3];

    C = 3 * (P1.x - P0.x);
    B = 3 * (P2.x - P1.x) - C;
    A = P3.x - P0.x - C - B;
    x = 3 * ( A * t * t ) + 2 * ( B * t ) + C;

    C = 3 * (P1.y - P0.y);
    B = 3 * (P2.y - P1.y) - C;
    A = P3.y - P0.y - C - B;
    y = 3 * ( A * t * t ) + 2 * ( B * t ) + C;

    T.x = x;
    T.y = y;

    return T;
}

vector_t Polybezier::compute_normal_at(int segment, double t)
{
    return vector_t(0.0, 0.0, 0.0);
} 

/*------------------------------- PROJECT_ONTO_PLANE ---
 |
 |  Project the vector, or point, onto the plane in 3D.
 |
 *-----------------------------------------------------*/
vector_t Polybezier::project_onto_plane(const vector_t& v)
{
    vector_t projection = v;
    projection.z = -(m_plane_normal.x * projection.x + m_plane_normal.y * projection.y) / m_plane_normal.z;
    return projection;
}

point_t Polybezier::project_onto_plane(const point_t& p)
{
    point_t projection = p;
    vector_t n = m_plane_normal;

    n = n / n.z;
    double c = m_plane_offset;

    projection.z = -(c + (n.x * p.x) + (n.y * p.y)) / n.z;
    return projection;
}
