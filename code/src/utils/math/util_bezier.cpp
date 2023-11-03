#include "utils/math/util_bezier_length.h"

/***************************************************************************/
#define INFTY 10000000	/* it suffices that this is longer than any pixel chain; it need not be really infinite */

/* ---------------------------------------------------------------------- */
/*
 *  B0, B1, B2, B3 : Bezier multipliers
 */

#define B0(u) ( ( 1.0 - u )  *  ( 1.0 - u )  *  ( 1.0 - u ) )
#define B1(u) ( 3 * u  *  ( 1.0 - u )  *  ( 1.0 - u ) )
#define B2(u) ( 3 * u * u  *  ( 1.0 - u ) )
#define B3(u) ( u * u * u )
/* ---------------------------------------------------------------------- */

double dot(point_t P1, point_t P2) {return P1.x*P2.x + P1.y*P2.y;}
/* ------------------------------------------------------------------------------------------------------------------------ */
/*Functions used for computing the Bezier control points;*/


/**
 *  Assign parameter values to digitized points using relative distances between points.
 *  use for data represented by one cubic Bezier
 */
void
chord_length_parameterize(QVector<point_t>* d, QVector<double>* u )
{
	int len = d->size();
	//assert( 2 <= len );

	(*u).resize(len);
	/* First let u[i] equal the distance travelled along the path from d[0] to d[i]. */
	(*u)[0] = 0.0;
	for (int i = 1; i < len; i++) {
		point_t p = d->at(i) - d->at(i-1);
		(*u)[i] = (*u)[i-1] + p.L2_norm2D();
	}

	/* Then scale to [0.0 .. 1.0]. */
	double tot_len = (*u)[len - 1];
	//assert( tot_len != 0 );
	for (int i = 1; i < len; ++i) {
		(*u)[i] /= tot_len;
	}
	(*u)[len - 1] = 1;
}

/**
 *  Estimates length of tangent vectors P1, P2, when direction is given
 *  fills in bezier with the correct values for P1, P2
 *  Point bezier[4]
 */
void
estimate_lengths(point_t bezier[],
                 QVector<point_t> data, QVector<double> uPrime,
                 point_t const &tHat1, point_t const &tHat2)
{
    double C[2][2];   /* Matrix C. */
    double X[2];      /* Matrix X. */

    /* Create the C and X matrices. */
    C[0][0] = 0.0;
    C[0][1] = 0.0;
    C[1][0] = 0.0;
    C[1][1] = 0.0;
    X[0]    = 0.0;
    X[1]    = 0.0;

    /* First and last control points of the Bezier curve are positioned exactly at the first and
       last data points. */
	unsigned len = data.size();
    bezier[0] = data[0];
    bezier[3] = data[len - 1];

    for (unsigned i = 0; i < len; i++) {
        /* Bezier control point coefficients. */
        double const b0 = B0(uPrime[i]);
        double const b1 = B1(uPrime[i]);
        double const b2 = B2(uPrime[i]);
        double const b3 = B3(uPrime[i]);

        /* rhs for eqn */
       point_t const a1 = tHat1 * b1;
       point_t const a2 = tHat2 * b2;

        C[0][0] += dot(a1, a1);
        C[0][1] += dot(a1, a2);
        C[1][0] = C[0][1];
        C[1][1] += dot(a2, a2);

        /* Additional offset to the data point from the predicted point if we were to set bezier[1]
           to bezier[0] and bezier[2] to bezier[3]. */
        point_t const shortfall
            = ( data[i]
                - ( bezier[0] *(b0 + b1) )
                - ( bezier[3] * (b2 + b3) ) );
        X[0] += dot(a1, shortfall);
        X[1] += dot(a2, shortfall);
    }

    /* We've constructed a pair of equations in the form of a matrix product C * alpha = X.
       Now solve for alpha. */
    double alpha_l, alpha_r;

    /* Compute the determinants of C and X. */
    double const det_C0_C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1];
    if ( det_C0_C1 != 0 ) {
        /* Apparently Kramer's rule. */
        double const det_C0_X  = C[0][0] * X[1]    - C[0][1] * X[0];
        double const det_X_C1  = X[0]    * C[1][1] - X[1]    * C[0][1];
        alpha_l = det_X_C1 / det_C0_C1;
        alpha_r = det_C0_X / det_C0_C1;
    } else {
        /* The matrix is under-determined.  Try requiring alpha_l == alpha_r.
         *
         * One way of implementing the constraint alpha_l == alpha_r is to treat them as the same
         * variable in the equations.  We can do this by adding the columns of C to form a single
         * column, to be multiplied by alpha to give the column vector X.
         *
         * We try each row in turn.
         */
        double const c0 = C[0][0] + C[0][1];
        if (c0 != 0) {
            alpha_l = alpha_r = X[0] / c0;
        } else {
            double const c1 = C[1][0] + C[1][1];
            if (c1 != 0) {
                alpha_l = alpha_r = X[1] / c1;
            } else {
                /* Let the below code handle this. */
                alpha_l = alpha_r = 0.;
            }
        }
    }

    /* If alpha negative, use the Wu/Barsky heuristic (see text).  (If alpha is 0, you get
       coincident control points that lead to divide by zero in any subsequent
       NewtonRaphsonRootFind() call.) */
    /// \todo Check whether this special-casing is necessary now that 
    /// NewtonRaphsonRootFind handles non-positive denominator.
    if ( alpha_l < 1.0e-6 ||
         alpha_r < 1.0e-6   )
    {
		point_t p = data[len - 1]- data[0];
		alpha_l = alpha_r = ( p.L2_norm2D()/ 3.0 );
    }

    /* Control points 1 and 2 are positioned an alpha distance out on the tangent vectors, left and
       right, respectively. */
	//alpha min is 4 ..so the points are far enough away
    /*bezier[1] = tHat1 * max(alpha_l,4.) + bezier[0];
    bezier[2] = tHat2 * max(alpha_r,4.) + bezier[3];*/

    bezier[1] = tHat1 * alpha_l + bezier[0];
    bezier[2] = tHat2 * alpha_r + bezier[3];
	

    return;
}


