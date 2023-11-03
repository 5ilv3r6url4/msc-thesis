#include "model/intersection.h"

/*------------------------------------- INTERSECTION ---
 |
 |  Constructor for an intersection using boost optional
 |  and boost none to capture member variables that
 |  we may not know at initialization.
 |
 *-----------------------------------------------------*/
Intersection::Intersection(boost::uuids::uuid polybez_a_uuid,
                           boost::uuids::uuid polybez_b_uuid,
                           boost::optional<bool> is_xhair,
                           boost::optional<bool> is_set_user,
                           boost::optional<bool> is_set_silhouette,
                           boost::optional<point_t>  position,
                           boost::optional<double> tval_a,
                           boost::optional<double> tval_b,
                           boost::optional<vector_t> normal,
                           boost::optional<vector_t> tangent_a,
                           boost::optional<vector_t> tangent_b,
                           boost::optional<vector_t> local_plane_a,
                           boost::optional<vector_t> local_plane_b)
    : m_uuid(boost::uuids::random_generator()())
{
    m_polybez_a_uuid = polybez_a_uuid;
    m_polybez_b_uuid = polybez_b_uuid;

    m_is_xhair = is_xhair.value_or(false);
    m_is_set_user = is_set_user.value_or(false);
    m_is_set_silhouette = is_set_silhouette.value_or(false);

    m_position = position.value_or(point_t(0.0, 0.0, 0.0));
    m_tval_a = tval_a.value_or(-1.0);
    m_tval_b = tval_b.value_or(-1.0);
    m_normal = normal.value_or(vector_t(0.0, 0.0, 0.0));
    m_tangent_a = tangent_a.value_or(vector_t(0.0, 0.0, 0.0));
    m_tangent_b = tangent_b.value_or(vector_t(0.0, 0.0, 0.0));
    m_local_plane_a = local_plane_a.value_or(vector_t(0.0, 0.0, 0.0));
    m_local_plane_b = local_plane_b.value_or(vector_t(0.0, 0.0, 0.0));

    m_label = -1;
}

/*------------------------------------ ~INTERSECTION ---
 |
 |  Intersection destructor.
 |
 *-----------------------------------------------------*/
Intersection::~Intersection() {}


/*-------------------------- CALCULATE_TANGENT_FRAME ---
 |
 |  Local, analytic solve of normal and tangent z
 |  components at an intersection.
 |
 *-----------------------------------------------------*/
void Intersection::calculate_tangent_frame()
{
    // ==============================================================================================
    // local solution for a single cross intersection using orthogonality constraint between tangents
    // ==============================================================================================
    // ( tangent_a / ||tangent_a|| ) dot ( tangent_b / ||tangent_b|| ) = 0.0
    // ( tangent_a ) dot ( tangent_b ) = 0.0
    // ( tangent_a.x * tangent_b.x ) + ( tangent_a.y * tangent_b.y ) + ( tangent_a.z * tangent_b.z ) = 0.0
    // ( tangent_b.z ) = -tangent_a.z substitution
    // ( tangent_a.x * tangent_b.x ) + ( tangent_a.y * tangent_b.y ) + ( tangent_a.z * -tangent_a.z ) = 0.0
    // ( tangent_a.x * tangent_b.x ) + ( tangent_a.y * tangent_b.y ) - ( tangent_a.z )^2 = 0.0
    // ( tangent_a.z )^2 = ( tangent_a.x * tangent_b.x ) + ( tangent_a.y * tangent_b.y )
    // ( tangent_a.z ) = sqrt( ( tangent_a.x * tangent_b.x ) + ( tangent_a.y * tangent_b.y ) )
    // ==============================================================================================

    m_tangent_a.normalize();
    m_tangent_b.normalize();

    double discriminant = 4 * ( m_tangent_a.x * m_tangent_b.x + m_tangent_a.y * m_tangent_b.y );

    bool flip = false;
    if (discriminant < 0) {
        m_tangent_a = -m_tangent_a;
        flip = true;

        discriminant = 4 * ( m_tangent_a.x * m_tangent_b.x + m_tangent_a.y * m_tangent_b.y );
    }

    m_tangent_a.z = sqrt(discriminant)/(-2);
    m_tangent_b.z = -m_tangent_a.z;

    m_tangent_a.normalize();
    m_tangent_b.normalize();

    // ------------------------------------------------------------------------------
    // normal is the cross product between both tangents at an intersection
    // it is important here to note that order of the cross product does matter
    // and all future cross products use the same order here of
    // ( tangent_a ) cross ( tangent_b )
    // so if the normal points into the screen when we take this cross product
    // we not only flip the normal, but also switch the data stored for both tangents
    m_normal = m_tangent_a ^ m_tangent_b;
    m_normal.normalize();
    if (m_normal.z < 0.0) { m_normal = -m_normal; }
    if (flip) { m_tangent_a = -m_tangent_a; }

    calculate_local_planes();
}
// ------------------------------------------------------------------------------

/*--------------------------- CALCULATE_LOCAL_PLANES ---
 |
 |  Calculate the local solutions at this intersection.
 |
 *-----------------------------------------------------*/
void Intersection::calculate_local_planes() {
    m_local_plane_a = m_tangent_b;
    m_tangent_a.z = -(m_tangent_a.x * m_local_plane_a.x + m_tangent_a.y * m_local_plane_a.y) / m_local_plane_a.z;
    m_tangent_b.z = -(m_tangent_a.x * m_tangent_b.x + m_tangent_a.y * m_tangent_b.y) / m_tangent_a.z;

    m_normal = m_tangent_a ^ m_tangent_b;
    if (m_normal.z < 0.0) { m_normal = -m_normal; }

    m_local_plane_b = m_tangent_a;
}

/*------------------------------- FLIP_TANGENT_FRAME ---
 |
 |  Flip tangent frame orientation from facing upwards
 |  positive to downwards and vice versa, with respect to
 |  the normal at the intersection, or from facing left
 |  positive to right positive and vice versa.
 |
 *-----------------------------------------------------*/
void Intersection::flip_tangent_frame()
{
    m_normal = -m_normal;
    m_normal.z = -m_normal.z;
    m_tangent_a.z = -m_tangent_a.z;
    m_tangent_b.z = -m_tangent_b.z;

    calculate_local_planes();
}
