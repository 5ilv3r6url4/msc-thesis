#ifndef UTIL_OBJ_WRITER_H
#define UTIL_OBJ_WRITER_H

#include "utils/math/util_math_defs.h"

#include <QVector>

using namespace std;

class Util_Obj_Writer
{
public:
	// Writes a list of polylines, each defined by a series of points
	static void writePolylines(QVector< QVector<point_t> > lines);
};

#endif // UTIL_OBJ_WRITER_H

