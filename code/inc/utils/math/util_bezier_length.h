#ifndef UTIL_BEZIER_H
#define UTIL_BEZIER_H

#include "utils/math/util_math_defs.h"

#include <iostream>

#include <QVector>

using namespace std;
/**
 *  Assign parameter values to digitized points using relative distances between points.
 *  use for data represented by one cubic Bezier
 */
void
chord_length_parameterize(QVector<point_t>* d, QVector<double>* u );

/**
 *  Estimates length of tangent vectors P1, P2, when direction is given
 *  fills in bezier with the correct values for P1, P2
 *  Point bezier[4]
 */
void
estimate_lengths(point_t bezier[],
                 QVector<point_t> data, QVector<double> uPrime,
                 point_t const &tHat1, point_t const &tHat2);
#endif // UTIL_BEZIER_H
