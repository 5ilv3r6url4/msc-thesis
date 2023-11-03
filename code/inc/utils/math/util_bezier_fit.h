#include "utils/math/util_math_defs.h"

#include <QVector>

#define INFTY 10000000	/* it suffices that this is longer than any pixel chain; it need not be really infinite */
#define PI 3.14159265
/* ---------------------------------------------------------------------- */
/*
*  B0, B1, B2, B3 : Bezier multipliers
*/

#define B0(u) ( ( 1.0 - u )  *  ( 1.0 - u )  *  ( 1.0 - u ) )
#define B1(u) ( 3 * u  *  ( 1.0 - u )  *  ( 1.0 - u ) )
#define B2(u) ( 3 * u * u  *  ( 1.0 - u ) )
#define B3(u) ( u * u * u )
/* ---------------------------------------------------------------------- */

class Bezier_Fit
{
public:
    Bezier_Fit(){closedCurve = false;};
	//fits a multi-segment Bezier curve to the set of points 
	//return true if the curve is closed and false otherwise
	bool fitBezier(QVector<point_t>& freeHandDrawnPath, QVector<point_t>& bezierPoints, double zoom = 1);
private:
	//COMPUTE POLYLINE
	void bresenham(point_t p0, point_t p1);

	void calc_lon();
	double penalty3(int i, int j);
	int bestpolygon();

	//FIT BEZIERS ON THE POLYLINE
	int smooth(double zoom, double alphamax = 2);
	void fillBezierPoints(QVector<point_t> continuousCurve, double zoom, bool isFirst=false, bool isLast=false);
	/**
	* Fit a multi-segment Bezier curve to a set of digitized points, without
	* possible weedout of identical points and NaNs.
	**/
	int bezier_fit_cubic(point_t bezier[], int split_points[], 
		QVector<point_t> data, int const data_len, point_t const &tHat1, point_t const &tHat2,
		double const error, unsigned const max_beziers);
	
	/**
	*  Assign parameter values to digitized points using relative distances between points.
	*  use for data represented by one cubic Bezier
	*/
	void chord_length_parameterize(QVector<point_t>* d, QVector<double>* u, unsigned const len );
	/**
	* Fill in bezier[] based on the given data and tangent requirements, using
	* a least-squares fit.
	**/
	void generate_bezier(point_t bezier[], QVector<point_t> data, QVector<double> u, unsigned const len,
		point_t const &tHat1, point_t const &tHat2, double const tolerance_sq);
/**
 * Given set of points and their parameterization, try to find a better assignment of parameter
 * values for the points.
 */
void reparameterize(QVector<point_t> data, unsigned const len, 
			   QVector<double> u, const point_t bezCurve[]);

double compute_max_error_ratio(QVector<point_t> d, QVector<double> u, unsigned const len,
                        const point_t bezCurve[4], double const tolerance,
                        unsigned *const splitPoint);
double compute_hook(point_t const &a, point_t const &b, double const u, const point_t bezCurve[4],
             double const tolerance);
/**
 *  Use Newton-Raphson iteration to find better root.
 */
double NewtonRaphsonRootFind(const point_t Q[], point_t const &P, double const u);	
/** 
 * Evaluate a Bezier curve at parameter value \a t.
 */
point_t bezier_pt(unsigned const degree, point_t const V[], double const t);
/**
	*  Estimates length of tangent vectors P1, P2, when direction is given
	*  fills in bezier with the correct values for P1, P2
	*  Point bezier[4]
	*/
	void estimate_lengths(point_t bezier[],
		QVector<point_t> data, QVector<double> uPrime, unsigned const len, 
		point_t const &tHat1, point_t const &tHat2);
	//tangents
	void estimate_bi(point_t bezier[4], unsigned const ei, QVector<point_t> data, QVector<double> u, unsigned const len);
	point_t estimate_left_tangent(QVector<point_t>* d, unsigned const len);
	point_t estimate_right_tangent(QVector<point_t>* d, unsigned const len);
	point_t estimate_left_tangent(QVector<point_t>* d, unsigned const len, double const tolerance_sq);
	point_t estimate_right_tangent(QVector<point_t>* d, unsigned const len, double const tolerance_sq);
	point_t estimate_center_tangent(QVector<point_t>* d , unsigned const center, unsigned const len);


	//auxiliar functions (common operations)
	void swap_val(int& i, int& j);
	/*for two pixels that are not neighbours: P0, P1
	determine a pixel neighbouring P0 that is on the line uniting P0 to P1
	returns direction with center in P0: (-1,0,1)
	*/
	point_t ddir(point_t a);
	/* return 1 if a <= b < c < a, in a cyclic sense (mod n) */
	int cyclic(int a, int b, int c);
	int floordiv(int a, int n);
	int sign(double x);
	/* integer arithmetic */
	int mod(const int a, const int n);
	double xprod(point_t p1, point_t p2);
	/* return (p1-p0)x(p2-p0), the area of the parallelogram */
	double dpara(point_t p0, point_t p1, point_t p2);

	/* ddenom/dpara have the property that the square of radius 1 centered
	at p1 intersects the line p0p2 iff |dpara(p0,p1,p2)| <= ddenom(p0,p2) */
	double ddenom(point_t p0, point_t p2);
	double dot(point_t P1, point_t P2);

	//MEMBERS
	QVector<point_t> pixelChain;
	QVector<point_t> bezierControlPoints;
	std::vector<int> polLine;
	std::vector<int> lon;        /* lon[len]: (i,lon[i]) = longest straight line from i */
	bool closedCurve;
	static point_t UNCONSTRAINED_TANGENT;
	static double TOLERANCE_CALLIGRAPHIC;

};
