#ifndef UTIL_GUI_FORMS_H
#define UTIL_GUI_FORMS_H

#include "utils/math/util_math_defs.h"
#include "utils/defs/util_data_defs.h"

// ---------------------------------
// QT incompatability with boost fix
// ---------------------------------
#ifndef Q_MOC_RUN
#include <boost/uuid/uuid.hpp>
#endif

struct polybezier_info_form_s {
    boost::uuids::uuid uuid;
    curve_type_t curve_type;
    vector_t plane_normal;
    double plane_offset;
    int num_segments;
    int curr_segment;
    point_t ctrl_points[4];
};

struct intersection_info_form_s {
    boost::uuids::uuid uuid;
    boost::uuids::uuid uuid_polybez_a;
    boost::uuids::uuid uuid_polybez_b;
    vector_t position;
    double tval_a;
    double tval_b;
    vector_t normal;
    vector_t tangent_a;
    vector_t tangent_b;
};


#endif // UTIL_GUI_FORMS_H
