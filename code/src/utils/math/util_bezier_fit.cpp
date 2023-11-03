#include "utils/math/util_bezier_fit.h"

#include <assert.h>
#include <iostream>
#include <algorithm>
#include <cfloat>

using namespace std;

point_t Bezier_Fit::UNCONSTRAINED_TANGENT = point_t(0,0);
double Bezier_Fit::TOLERANCE_CALLIGRAPHIC = 8; //3

bool Bezier_Fit::fitBezier(QVector<point_t>& freeHandDrawnPath, QVector<point_t>& bezierPoints, double zoom)
{
	point_t p0, p;
	QVector<point_t>::iterator it;
	bool closedCurve = false;
	double dist_x, dist_y;
	int len = freeHandDrawnPath.size();
	int i;
	pixelChain.clear();

	p0 = freeHandDrawnPath[0];
	for(i=1;i<len;i++)
	{
		p = freeHandDrawnPath[i];
        //pixelChain.push_back(p);
		bresenham(p0, p);
		p0 = p;
	}

	//check for duplicate points and eliminate them
	it = pixelChain.begin();
	it++;
	while(it!=pixelChain.end())
	{
		if(*(it-1) == *it)
		{
			it = pixelChain.erase(it);
		}
		else
			it++;
	}

	//check to see whether it's a closed chain


	calc_lon();
	bestpolygon();
//    for(int i=0;i<pixelChain.size();i++)
//        polLine.push_back(i);
	smooth(zoom);

	len = bezierControlPoints.size();
	for(i=0;i<len;i++)
	{
		bezierPoints.push_back(bezierControlPoints[i]);
	}

	return closedCurve;
}


//the freeHandDrawnPath is not continuous; we need to fill in pixels missing (so we can use potrace polyline fitting)
void Bezier_Fit::bresenham(point_t p0, point_t p1)
{
	int x0 = p0.x;
	int y0 = p0.y;
	int x1 = p1.x;
	int y1 = p1.y;
	point_t p;
	QVector<point_t> computed_pixels;
	int computed_pixels_size;

	bool steep = abs(y1 - y0) > abs(x1 - x0);
	bool swap_points = false;
	if (steep)
	{
		swap_val(x0, y0);
		swap_val(x1, y1);
	}
	if (x0 > x1)
	{
		swap_points = true;
		swap_val(x0, x1);
		swap_val(y0, y1);
	}
	int deltax = x1 - x0;
	int deltay = abs(y1 - y0);
	int error = -(deltax + 1) / 2;
	int ystep;
	int y = y0;
	if (y0 < y1)
        ystep = 1;
	else 
        ystep = -1;
	for (int x = x0; x <= x1; x++)
	{
		if (steep)
		{
            p.x = y;
            p.y = x;
			computed_pixels.push_back(p);
		}
		else 
		{
            p.x = x;
            p.y = y;
			computed_pixels.push_back(p);
		}
		error = error + deltay;
		if (error >= 0)
		{
			y = y + ystep;
			error = error - deltax;
		}
	}

	computed_pixels_size = computed_pixels.size();
	if(swap_points)
	{
		//invert the pixels
		for (int i=computed_pixels_size-1; i>=0;i--)
		{
			pixelChain.push_back(computed_pixels[i]);
		}
	}
	else {
		for (int i=0; i<computed_pixels_size;i++)
		{
			pixelChain.push_back(computed_pixels[i]);
		}
	}
	computed_pixels.clear();

}


//POTRACE
/* ---------------------------------------------------------------------- */
/* Stage 1: determine the straight subpaths (Sec. 2.2.1). Fill in the
"lon" component of a path object (based on pt/len).	For each i,
lon[i] is the furthest index such that a straight line can be drawn
from i to lon[i]. Return 1 on error with errno set, else 0. */

/* this algorithm depends on the fact that the existence of straight
subpaths is a triplewise property. I.e., there exists a straight
line through squares i0,...,in iff there exists a straight line
through i,j,k, for all i0<=i<j<k<=in. (Proof?) */

/* this implementation of calc_lon is O(n^2). It replaces an older
O(n^3) version. A "constraint" means that future points must
satisfy xprod(constraint[0], cur) >= 0 and xprod(constraint[1],
cur) <= 0. */

/* Remark for Potrace 1.1: the current implementation of calc_lon is
more complex than the implementation found in Potrace 1.0, but it
is considerably faster. The introduction of the "nc" data structure
means that we only have to test the constraints for "corner"
points. On a typical input file, this speeds up the calc_lon
function by a factor of 31.2, thereby decreasing its time share
within the overall Potrace algorithm from 72.6% to 7.82%, and
speeding up the overall algorithm by a factor of 3.36. On another
input file, calc_lon was sped up by a factor of 6.7, decreasing its
time share from 51.4% to 13.61%, and speeding up the overall
algorithm by a factor of 1.78. In any case, the savings are
substantial. */

void Bezier_Fit::calc_lon() 
{
	int len = (int)pixelChain.size();
	lon.resize(len);

	int i; 
	int j;
	int k, k1, dir, dir_index;
	int ct[4];
	point_t constraint[2];
	point_t cur;
	point_t off;
	point_t p_dir;
	std::vector<int> pivk;  pivk.resize(len);
	std::vector<int> nc; nc.resize(len);  /* nc[len]: next corner */
	point_t dk;  /* direction of k-k1 */
	int a, b, c, d;

	/* initialize the nc data structure. Point from each point to the
	furthest future point to which it is connected by a vertical or
	horizontal segment or is the next neighbour  */
	k = len-1;
	for (i=len-2; i>=0; i--) {
		if (pixelChain[i].x != pixelChain[k].x && pixelChain[i].y != pixelChain[k].y) {
			k = (int)(i+1);  /* necessarily i<n-1 in this case */
		}
		nc[i] = k;
	}
	nc[len-1] = len-1;




	pivk[len-1] = len-1;
	/* determine pivot points: for each i, let pivk[i] be the furthest k
	such that all j with i<j<k lie on a line connecting i,k. */
	int direction[9];
	direction[0] = 1; direction[1] = 0; direction[2] = 7;
	direction[3] = 2; direction[4] = 8; direction[5] = 6;
	direction[6] = 3; direction[7] = 4; direction[8] = 5;

	for (i=len-2; i>=0; i--) {

		ct[0] = ct[1] = ct[2] = ct[3] = 0;


		/* keep track of "directions" that have occurred */ 
		p_dir.x = pixelChain[i+1].x-pixelChain[i].x;
		p_dir.y = pixelChain[i+1].y-pixelChain[i].y;
		p_dir = ddir(p_dir);
		dir_index= (4+3*(p_dir.x)+p_dir.y);
		assert(dir_index<9);


		dir = (direction[dir_index]+0.5)/2 - 1;
		ct[dir]=ct[dir]+1;

		if(direction[dir_index]%2>0) 
		{
			ct[mod(dir+1,4)]=ct[mod(dir+1,4)]+1;
		}
		constraint[0].x = 0;
		constraint[0].y = 0;
		constraint[1].x = 0;
		constraint[1].y = 0;

		/* find the next k such that no straight line from i to k */

		k = nc[i];
		k1 = i;
		while (1) {
			p_dir.x = pixelChain[i+1].x-pixelChain[i].x;
			p_dir.y = pixelChain[i+1].y-pixelChain[i].y;
			p_dir = ddir(p_dir);
			dir_index= (4+3*p_dir.x+p_dir.y);
			assert(dir_index<9);
			dir = (direction[dir_index]+0.5)/2 - 1; 
			ct[dir]++;
			if(direction[dir_index]%2>0) 
			{
				ct[mod(dir+1,4)]++;
			}
			/* if all four "directions" have occurred, cut this path */
			if (ct[0] && ct[1] && ct[2] && ct[3]) {
				pivk[i] = k1;
				goto foundk;
			}

			cur.x = pixelChain[k].x - pixelChain[i].x;
			cur.y = pixelChain[k].y - pixelChain[i].y;

			/* see if current constraint is violated */
			if (xprod(constraint[0], cur) < 0 || xprod(constraint[1], cur) > 0) {
				goto constraint_viol;
			}

			/* else, update constraint */
			if (fabs(cur.x) <= 1 && fabs(cur.y) <= 1) {
				/* no constraint */
			} else {
				off.x = cur.x + ((cur.y>=0 && (cur.y>0 || cur.x<0)) ? 1 : -1);
				off.y = cur.y + ((cur.x<=0 && (cur.x<0 || cur.y<0)) ? 1 : -1);
				if (xprod(constraint[0], off) >= 0) {
					constraint[0] = off;
				}
				off.x = cur.x + ((cur.y<=0 && (cur.y<0 || cur.x<0)) ? 1 : -1);
				off.y = cur.y + ((cur.x>=0 && (cur.x>0 || cur.y<0)) ? 1 : -1);
				if (xprod(constraint[1], off) <= 0) {
					constraint[1] = off;
				}
			}	
			k1 = k;
			k = nc[k1];
			if (!cyclic(k,i,k1)) {
				break;
			}
		} //end while

constraint_viol:
		/* k1 was the last "corner" satisfying the current constraint, and
		k is the first one violating it. We now need to find the last
		point along k1..k which satisfied the constraint. */

		dk.x = sign(pixelChain[k].x-pixelChain[k1].x);
		dk.y = sign(pixelChain[k].y-pixelChain[k1].y);
		cur.x = pixelChain[k1].x - pixelChain[i].x;
		cur.y = pixelChain[k1].y - pixelChain[i].y;

		/* find largest integer j such that xprod(constraint[0], cur+j*dk)
		>= 0 and xprod(constraint[1], cur+j*dk) <= 0. Use bilinearity
		of xprod. */
		a = xprod(constraint[0], cur);
		b = xprod(constraint[0], dk);
		c = xprod(constraint[1], cur);
		d = xprod(constraint[1], dk);
		/* find largest integer j such that a+j*b>=0 and c+j*d<=0. This
		can be solved with integer arithmetic. */
		j = INFTY;
		if (b<0) {
			j = floordiv(a,-b);
		}
		if (d>0) {
			j = min(j, floordiv(-c,d));
		}
		pivk[i] = mod(k1+j,len);
foundk:
		;
	} 

	j=pivk[len-1];
	lon[len-1]=j;
	i=len-2;
	for (i=len-2; i>=0; i--) {
		if (cyclic(i+1,pivk[i],j)) {
			j=pivk[i];
		}
		lon[i]=j;
	}


	pivk.clear();
	nc.clear();
	return;
}

/*********************************************************************************************/
/* ---------------------------------------------------------------------- */
/* Stage 2: calculate the optimal polygon (Sec. 2.2.2-2.2.4). */ 

/* Auxiliary function: calculate the penalty of an edge from i to j in
the given path. This needs the "lon" and "sum*" data. */

double Bezier_Fit::penalty3(int i, int j) {

	/* assume 0<=i<j<=n  */
	int len = (int)pixelChain.size();
	double x, y, x0, y0, x1, y1;
	double AB, sum;


	if(i==j) return 0;

	if(i>j)
	{
		int tmp = j;
		j = i; i = tmp;
	}
	if(j>len-1)
		j = len-1;

	int A, B, C;
    x0 = pixelChain[i].x; y0 = pixelChain[i].y;
    x1 = pixelChain[j].x; y1 = pixelChain[j].y;
	A = y0 - y1;
	B = x1 - x0;
	C = x0*y1 - x1*y0;

	AB = sqrt((double)A*A+B*B);
	sum = 0;
	for(int k=i;k<=j;k++)
	{
        x = pixelChain[k].x; y = pixelChain[k].y;
		sum +=fabs(A*x+B*y+C)/AB;
	}
	sum = AB*sqrt((double)sum/(double)(j-i+1));
	return sum;
}

/* find the optimal polygon. Fill in the m and po components. Return 1
on failure with errno set, else 0. Non-cyclic version: assumes i=0
is in the polygon. Fixme: ### implement cyclic version. */

int Bezier_Fit::bestpolygon()
{
	int i, j, m, k;     
	int n = (int)pixelChain.size();
	double *pen = NULL; // pen[n+1]: penalty vector 
	int *prev = NULL;   // prev[n+1]: best path pointer vector 
	int *clip0 = NULL;  // clip0[n]: longest segment pointer, non-cyclic 
	int *clip1 = NULL;  // clip1[n+1]: backwards segment pointer, non-cyclic 
	int *seg0 = NULL;    // seg0[m+1]: forward segment bounds, m<=n 
	int *seg1 = NULL;   // seg1[m+1]: backward segment bounds, m<=n 
	double thispen;
	double best;
	int c;

	pen = new double[n+1];
	prev = new int[n+1];
	clip0 = new int[n+1];
	clip1 = new int[n+1];
	seg0 = new int[n+1];
	seg1 = new int[n+1];

	// calculate clipped paths 
	for (i=0; i<n; i++) {
		c = lon[i];
		if (c == i) {
			c = mod(i+1,n);
		}
		if (c < i) {
			clip0[i] = n;
		} else {
			clip0[i] = c;
		}
	}

	// calculate backwards path clipping, non-cyclic. j <= clip0[i] if
	//clip1[j] <= i, for i,j=0..n. 
	j = 1;
	for (i=0; i<n; i++) {
		while (j <= clip0[i]) {
			clip1[j] = i;
			j++;
		}
	}

	// calculate seg0[j] = longest path from 0 with j segments 
	i = 0;
	for (j=0; i<n; j++) {
		seg0[j] = i;
		i = clip0[i];
	}
	seg0[j] = n;
	m = j;
	// calculate seg1[j] = longest path to n with m-j segments 
	i = n;
	for (j=m; j>0; j--) {
		seg1[j] = i;
		i = clip1[i];
	}
	seg1[0] = 0;

	// now find the shortest path with m segments, based on penalty3 
	// note: the outer 2 loops jointly have at most n interations, thus
	//the worst-case behavior here is quadratic. In practice, it is
	//close to linear since the inner loop tends to be short. 
	pen[0]=0;
	for (j=1; j<=m; j++) {
		for (i=seg1[j]; i<=seg0[j]; i++) {
			best = -1;
			for (k=seg0[j-1]; k>=clip1[i]; k--) {
				thispen = penalty3(k, i) + pen[k];
				if (best < 0 || thispen < best) {
					prev[i] = k;
					best = thispen;
				}
			}
			pen[i] = best;
		}
	}

	polLine.resize(m);


	// read off shortest path 
	for (i=n, j=m-1; i>0; j--) {
		i = prev[i];
		polLine[j] = i;
	}

	delete [] pen;
	delete [] prev;
	delete [] clip0;
	delete [] clip1;
	delete [] seg0;
	delete [] seg1;
	return 0;
}


int Bezier_Fit::smooth(double zoom, double alphamax) {
	int m = polLine.size();

	int i, j, k;
	double dd, denom, alpha;

	QVector<point_t> continuousCurve;
	bool isFirst = true;
	bool isLast = false;
	//for (i=0; i<m; i++) {
	//	polLineCorners.push_back(CURVE);
	//}

	/* examine each vertex and find its best fit */
	continuousCurve.push_back(pixelChain[polLine[0]]);
	for (i=0; i<m-2; i++) {
		j = mod(i+1, m);
		k = mod(i+2, m);

		denom = ddenom(pixelChain[polLine[i]], pixelChain[polLine[k]]);
		if (denom != 0.0) {
			dd = dpara(pixelChain[polLine[i]], pixelChain[polLine[j]], pixelChain[polLine[k]]) / denom;
			dd = fabs(dd);
			alpha = dd>1 ? (1 - 1.0/dd) : 0;
			alpha = alpha / 0.75;
		} 
		else {
			alpha = 4/3.0;
		}
		//curve->alpha0[j] = alpha;	 /* remember "original" value of alpha */

		if (alpha >= alphamax) {  /* pointed corner */
			//polLineCorners[j] = LINE;
			continuousCurve.push_back(pixelChain[polLine[j]]);
			fillBezierPoints(continuousCurve, zoom, isFirst, isLast);
			isFirst = false;
			continuousCurve.clear();
			continuousCurve.push_back(pixelChain[polLine[j]]);
		} 
		else {
			//polLineCorners[j] = CURVE;
			continuousCurve.push_back(pixelChain[polLine[j]]);
		}
	}//end for(i<m-2)
	//add the last point
	if(continuousCurve.isEmpty()) //does this really happens?
	{
		continuousCurve.push_back(pixelChain[polLine[m-2]]);
		continuousCurve.push_back(pixelChain[polLine[m-1]]);
	}
	else
		continuousCurve.push_back(pixelChain[polLine[m-1]]);
	isLast = true;
	fillBezierPoints(continuousCurve, zoom, isFirst, isLast);
	return 0;
}

//fit Beziers that pass through all points of the polyline (twice as many points)
void Bezier_Fit::fillBezierPoints(QVector<point_t> continuousCurve, double zoom, bool isFirst, bool isLast)
{
	
	int curveLen = continuousCurve.size();
	//int Max_Beziers = (curveLen+3)/4;
	int Max_Beziers = curveLen; //one bezier for each couple of points
	point_t *bezier = new point_t[4*Max_Beziers];
	double const tolerance_sq = fabs(zoom)*Bezier_Fit::TOLERANCE_CALLIGRAPHIC;
	point_t tHat1 = Bezier_Fit::UNCONSTRAINED_TANGENT; 
	point_t tHat2 = Bezier_Fit::UNCONSTRAINED_TANGENT;
	//int split_points[];
	int nsegs;

	nsegs = bezier_fit_cubic(bezier, NULL, continuousCurve, curveLen,
		tHat1, tHat2,
		tolerance_sq, Max_Beziers);
	for(int i=0;i<nsegs;i++)
	{
		int start = i*4;
        /*std::cerr<<"bezier["<<start<<"] = ("<<bezier[start].x<<"; "<<bezier[start].y<<"); ";
        std::cerr<<"bezier["<<start+1<<"] = ("<<bezier[start+1].x<<"; "<<bezier[start+1].y<<"); ";
        std::cerr<<"bezier["<<start+2<<"] = ("<<bezier[start+2].x<<"; "<<bezier[start+2].y<<"); ";
        std::cerr<<"bezier["<<start+3<<"] = ("<<bezier[start+3].x<<"; "<<bezier[start+3].y<<"); \n";*/

		bezierControlPoints.push_back(bezier[start]);
		bezierControlPoints.push_back(bezier[start+1]);
		bezierControlPoints.push_back(bezier[start+2]);
	}
	if(isLast && nsegs>0)
	{
		bezierControlPoints.push_back(bezier[nsegs*4-1]);
	}

}

//void BezierFit::fillBezierPoints(QVector<Point> continuousCurve, bool isFirst, bool isLast)
//{
//	int cLen;
//	int curveLen = continuousCurve.size();
//	Point bezier[4];
//	QVector<Point> data;
//	QVector<double> u;
//	Point tHat1; Point tHat2;
//	for(cLen=0;cLen<curveLen-3;cLen+=3)
//	{
//		data.push_back(continuousCurve[cLen]);
//		data.push_back(continuousCurve[cLen+1]);
//		data.push_back(continuousCurve[cLen+2]);
//		data.push_back(continuousCurve[cLen+3]);
//
//
//		chord_length_parameterize(&data, &u);
//		//compute tangent direction.
//		//Tangent P0-P1
//		if(cLen>0)//this lot is a continuation of another
//			tHat1 = continuousCurve[cLen+1] - continuousCurve[cLen-1];
//		else if ((cLen==0)&&(isFirst)&&(isClosedCurve())) //first point in a closed curve
//		{
//			//index of Point P(n-2) (for polyline P(0)....P(n-1), with P(0)=P(n-1))
//			int prev = polLine[polLine.size()-2];
//
//			tHat1 = continuousCurve[cLen+1] - pixelChain[prev];
//		}
//		else
//			tHat1 = continuousCurve[cLen+1] - continuousCurve[cLen];
//
//		//Tangent P2-P3
//		if(cLen+4<curveLen) //there's a continuation after this lot of 4 points
//			tHat2 = continuousCurve[cLen+2] - continuousCurve[cLen+4];
//		else if ((cLen+4==curveLen)&&(isLast)&&(isClosedCurve())) //last point in a closed curve; we're at the last group
//		{
//			//index of Point P(1) (for polyline P(0)....P(n-1), with P(0)=P(n-1))
//			int next = polLine[polLine[1]];
//			tHat2 = continuousCurve[cLen+2] - pixelChain[next];
//		}
//		else
//			tHat2 = continuousCurve[cLen+2] - continuousCurve[cLen+3];
//		tHat1.normalize();
//		tHat2.normalize();
//
//		estimate_lengths(bezier, data, u, tHat1, tHat2);
//		bezierControlPoints.push_back(bezier[0]);
//		bezierControlPoints.push_back(bezier[1]);
//		bezierControlPoints.push_back(bezier[2]);
//		if(isLast&&(cLen+4==curveLen))
//		{
//			bezierControlPoints.push_back(bezier[3]);
//		}
//		data.clear();
//		u.clear();
//	}//end for (all groups of 4 points)
//
//	//check whether there're any points left
//	switch((curveLen-1)%3)
//	{
//	case 0:
//		//nothing left
//		break;
//	case 2:
//		data.push_back(continuousCurve[curveLen-3]);
//	case 1:
//		data.push_back(continuousCurve[curveLen-2]);
//		data.push_back(continuousCurve[curveLen-1]);
//		//we have 2 or 3 points; 
//		cLen = curveLen-data.size();
//		chord_length_parameterize(&data, &u);
//		//compute tangent direction.
//		//Tangent P0-P1
//		if(cLen>0)//this lot is a continuation of another
//			tHat1 = continuousCurve[cLen+1] - continuousCurve[cLen-1];
//		else if ((cLen==0)&&(isFirst)&&(isClosedCurve())) //first point in a closed curve
//		{
//			//index of Point P(n-2) (for polyline P(0)....P(n-1), with P(0)=P(n-1))
//			int prev = polLine[polLine.size()-2];
//
//			tHat1 = continuousCurve[cLen+1] - pixelChain[prev];
//		}
//		else
//			tHat1 = continuousCurve[cLen+1] - continuousCurve[cLen];
//
//		//Tangent P2-P3
//
//		if ((isLast)&&(isClosedCurve())) //last point in a closed curve; we're at the last group of points anyway
//		{
//			//index of Point P(1) (for polyline P(0)....P(n-1), with P(0)=P(n-1))
//			int next = polLine[polLine[1]];
//			tHat2 = continuousCurve[curveLen-2] - pixelChain[next];
//		}
//		else
//			tHat2 = continuousCurve[curveLen-2] - continuousCurve[curveLen-1];
//		tHat1.normalize();
//		tHat2.normalize();
//
//		estimate_lengths(bezier, data, u, tHat1, tHat2);
//		bezierControlPoints.push_back(bezier[0]);
//		bezierControlPoints.push_back(bezier[1]);
//		bezierControlPoints.push_back(bezier[2]);
//		if(isLast)
//		{
//			bezierControlPoints.push_back(bezier[3]);
//		}
//		data.clear();
//		u.clear();
//		break;
//	}
//}

/* ------------------------------------------------------------------------------------------------------------------------ */
/*Functions used for computing the Bezier control points;*/

/** 
* Estimate the (forward) tangent at point d[first + 0.5].
*
* Unlike the center and right versions, this calculates the tangent in 
* the way one might expect, i.e., wrt increasing index into d.
* \pre (2 \<= len) and (d[0] != d[1]).
**/
point_t Bezier_Fit::estimate_left_tangent(QVector<point_t>* d, unsigned const len)
{
	//compute the left tangent using apportioned chords
	//the tangent ensures that the Bezier passes throught the second data point, also
	point_t p;
	point_t x0, x1, x2, x3;
	point_t tHat1;
	point_t L1, L2, L3;
	double t1, t2, m1, n1, m2, n2;
	assert( len >= 2 );
	
	switch (len)
	{
	case 2:
		p = (d->at(1) - d->at(0)); p.normalize();
		break;
	case 3:
		x0 = d->at(0);
		x1 = x2 = d->at(1);
		x3 = d->at(2);
		L1 = x1-x0;
		L2 = x2-x1;
		t1 = ( L1.L2_norm2D() )/(L1.L2_norm2D() + L2.L2_norm2D() );
		tHat1 = (x1 - x0*pow((double)(1-t1),(double)3.f) - x3*pow((double)t1, (double)3.f))/(3*t1*pow((double)(1-t1),(double)2.f))*(1-t1);
		p = tHat1-x0;
		p.normalize();
	default:
		x0 = d->at(0);
		x1 = d->at(1);
		x2 = d->at(2);
		x3 = d->at(3);
		L1 = x1-x0;
		L2 = x2-x1;
		L3 = x3-x2;
		t1 = ( L1.L2_norm2D() )/(L1.L2_norm2D() + L2.L2_norm2D() + L3.L2_norm2D() );
		t2 = ( L1.L2_norm2D() + L2.L2_norm2D() )/(L1.L2_norm2D() + L2.L2_norm2D() + L3.L2_norm2D() );

		m1 = B1(t1); m2 = B1(t2);
		n1 = B2(t1); n2 = B2(t2);
		L1 = x1 - x0*B0(t1) - x3*B3(t1);
		L2 = x2 - x0*B0(t2) - x3*B3(t2);
		tHat1 = L2/(m2-m1/n1) - L1/(n1*m2-m1);
		p = tHat1-x0;
		p.normalize();
	}
	return p;

}

/** 
* Estimates the (backward) tangent at d[last - 0.5].
*
* \note The tangent is "backwards", i.e. it is with respect to 
* decreasing index rather than increasing index.
*
* \pre 2 \<= len.
* \pre d[len - 1] != d[len - 2].
* \pre all[p in d] in_svg_plane(p).
*/
point_t Bezier_Fit::estimate_right_tangent(QVector<point_t>* d, unsigned const len)
{
	//compute the right tangent using apportioned chords
	//the tangent ensures that the Bezier passes throught the third data point, also
	point_t p;
	point_t x0, x1, x2, x3;
	point_t tHat2;
	point_t L1, L2, L3;
	double t1, t2, m1, n1, m2, n2;
	unsigned const last = len - 1;
	assert( len >= 2 );
	
	switch (len)
	{
	case 2:
		p = (d->at(last-1) - d->at(last)); p.normalize();
		break;
	case 3:
		x0 = d->at(last-2);
		x1 = x2 = d->at(last-1);
		x3 = d->at(last);
		L1 = x1-x0;
		L2 = x2-x1;
		t1 = ( L1.L2_norm2D() )/(L1.L2_norm2D() + L2.L2_norm2D() );
		tHat2 = (x1 - x0*pow((double)(1-t1),(double)3.f) - x3*pow((double)t1, (double)3.f))/(3*t1*pow((double)(1-t1),(double)2.f))*(1-t1);
		p = tHat2-x3;
		p.normalize();
	default:
		x0 = d->at(last-3);
		x1 = d->at(last-2);
		x2 = d->at(last-1);
		x3 = d->at(last);
		L1 = x1-x0;
		L2 = x2-x1;
		L3 = x3-x2;
		t1 = ( L1.L2_norm2D() )/(L1.L2_norm2D() + L2.L2_norm2D() + L3.L2_norm2D() );
		t2 = ( L1.L2_norm2D() + L2.L2_norm2D() )/(L1.L2_norm2D() + L2.L2_norm2D() + L3.L2_norm2D() );

		m1 = B1(t1); m2 = B1(t2);
		n1 = B2(t1); n2 = B2(t2);
		L1 = x1 - x0*B0(t1) - x3*B3(t1);
		L2 = x2 - x0*B0(t2) - x3*B3(t2);
		tHat2 = (L1*m2 - L2*m1)/(m2*n1-m1*n2);
		p = tHat2-x3;
		p.normalize();
	}
	return p;


	//assert( 2 <= len );
	//unsigned const last = len - 1;
	//unsigned const prev = last - 1;
	//assert( d->at(last) != d->at(prev) );
	//Point p = (d->at(prev) - d->at(last)); p.normalize();
	//return p;
}

void Bezier_Fit::estimate_bi(point_t bezier[4], unsigned const ei,
							QVector<point_t> data, QVector<double> u, unsigned const len)
{
	if(!(1 <= ei && ei <= 2))
		return;
	unsigned const oi = 3 - ei;
	point_t num = point_t(0.,0.);
	double den = 0.;
	for (unsigned i = 0; i < len; ++i) {
		double const ui = u[i];
		double const b[4] = {
			B0(ui),
			B1(ui),
			B2(ui),
			B3(ui)
		};

		{
			num.x += b[ei] * (b[0]  * bezier[0].x +
				b[oi] * bezier[oi].x +
				b[3]  * bezier[3].x +
				- data[i].x);

			num.y += b[ei] * (b[0]  * bezier[0].y +
				b[oi] * bezier[oi].y +
				b[3]  * bezier[3].y +
				- data[i].y);
		}
		den -= b[ei] * b[ei];
	}

	if (den != 0.) {

		bezier[ei] = num / den;
	} 
	else {
		bezier[ei] = (  bezier[0] * oi + bezier[3] * ei) / 3.;
	}
}

/** 
* Estimates the (backward) tangent at d[center], by averaging the two 
* segments connected to d[center] (and then normalizing the result).
*
* \note The tangent is "backwards", i.e. it is with respect to 
* decreasing index rather than increasing index.
*
* \pre (0 \< center \< len - 1) and d is uniqued (at least in 
* the immediate vicinity of \a center).
*/
point_t Bezier_Fit::estimate_center_tangent(QVector<point_t>* d ,
										 unsigned const center, unsigned const len)
{
	assert( center != 0 );
	assert( center < len - 1 );
	point_t ret;
	point_t P1 = d->at(center + 1);
	point_t P2 = d->at(center - 1);

	if( P1 == P2)
	{
		/* Rotate 90 degrees in an arbitrary direction. */
		point_t const diff = d->at(center) - d->at(center - 1);
		ret.x = -diff.y;
		ret.y = diff.x;
	}
	else
	{
		ret = d->at(center - 1) - d->at(center + 1); 
	}
	ret.normalize();
	return ret;
}


/** 
* Estimate the (forward) tangent at point d[0].
*
* Unlike the center and right versions, this calculates the tangent in 
* the way one might expect, i.e., wrt increasing index into d.
*
* \pre 2 \<= len.
* \pre d[0] != d[1].
* \pre all[p in d] in_svg_plane(p).
* \post is_unit_vector(ret).
**/

point_t Bezier_Fit::estimate_left_tangent(QVector<point_t>* d, unsigned const len, double const tolerance_sq)
{
	assert( 2 <= len );
	assert( 0 <= tolerance_sq );

	//tolerance_sq is always 0!!!)
	//return estimate_left_tangent(d, len);
	//std::cerr<<"estimate_left_tangent COMPUTE TANGENT LEFT: len = "<<len<<"; tolerance_sq = "<<tolerance_sq<<"\n";
	point_t p, pi, t;
	double distsq;
	for (unsigned i = 1;;) {
		pi = d->at(i);
		t = pi - d->at(0);
		distsq = dot(t, t);
		if ( tolerance_sq < distsq ) {
			p = t; p.normalize();
			return p;
		}
		++i;
		if (i == len) {
			p = t; p.normalize();
			return ( distsq == 0
				? estimate_left_tangent(d, len)
				: p );
		}
	}
}

/** 
* Estimates the (backward) tangent at d[last].
*
* \note The tangent is "backwards", i.e. it is with respect to 
* decreasing index rather than increasing index.
*
* \pre 2 \<= len.
* \pre d[len - 1] != d[len - 2].
* \pre all[p in d] in_svg_plane(p).
*/
point_t Bezier_Fit::estimate_right_tangent(QVector<point_t>* d, unsigned const len, double const tolerance_sq)

{
	//return estimate_right_tangent(d, len);
	//tolerance_sq =0 !
	assert( 2 <= len );
	assert( 0 <= tolerance_sq );
	unsigned const last = len - 1;
	point_t p, pi, t;
	double distsq;

	for (unsigned i = last - 1;; i--) {
		pi = d->at(i);
		t = pi - d->at(last);
		distsq = dot(t, t);
		if ( tolerance_sq < distsq ) {
			p = t; p.normalize();
			return p;
		}
		if (i == 0) {
			p = t; p.normalize();
			return ( distsq == 0
				? estimate_right_tangent(d, len)
				: p );
		}
	}
}

/**
*  Assign parameter values to digitized points using relative distances between points.
*  use for data represented by one cubic Bezier
*/
void Bezier_Fit::chord_length_parameterize(QVector<point_t>* d, QVector<double>* u, unsigned const len )
{
	//int len = d->size();
	assert( 2 <= len );

	(*u).resize(len);
	/* First let u[i] equal the distance travelled along the path from d[0] to d[i]. */
	(*u)[0] = 0.0;
	for (int i = 1; i < len; i++) {
		point_t p = d->at(i) - d->at(i-1);
		(*u)[i] = (*u)[i-1] + p.L2_norm2D();
	}

	/* Then scale to [0.0 .. 1.0]. */
	double tot_len = (*u)[len - 1];
	assert( tot_len != 0 );
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
void Bezier_Fit::estimate_lengths(point_t bezier[],
								 QVector<point_t> data, QVector<double> uPrime, unsigned const len, 
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


//common operations
void Bezier_Fit::swap_val(int& i, int& j) {
	int t = i;
	i = j;
	j = t;
}
/* integer arithmetic */
int Bezier_Fit::mod(const int a, const int n) {
	return a>=n ? a%n : a>=0 ? a : n-1-(-1-a)%n;
}

/*for two pixels that are not neighbours: P0, P1
determine a pixel neighbouring P0 that is on the line uniting P0 to P1
returns direction with center in P0: (-1,0,1)
*/
point_t Bezier_Fit::ddir(point_t a)
{
	point_t p;
	p.x = 0; p.y = 0;
	if(a.y==a.x==0)
		return p;
	//atan2 gives an error if both arguments are 0;
	double theta = atan2(a.y,a.x)* 180 / PI;
	if(theta<0)
	{
		if((a.x>=0)&&(a.y>=0))
			theta +=180; //
		else if((a.x<=0)&&(a.y<=0))
			theta +=360; //
		else
			std::cerr<<"oops!\n";
	}
	else
	{
		if((a.x>=0)&&(a.y<=0)&&theta<270)
			theta = 360-theta; //
	}
	//0/360
	if(theta<22.5 || theta>=337.5)
	{
		p.x = 1; p.y = 0;
		return p;
	}
	//45
	if(theta>=22.5 && theta<67.5)
	{
		p.x = 1; p.y = 1;
		return p;
	}
	//90
	if(theta>=67.5 && theta<112.5)
	{
		p.x = 0; p.y = 1;
		return p;
	}	
	//135
	if(theta>=112.5 && theta<157.5)
	{
		p.x = -1; p.y = 1;
		return p;
	}
	//180
	if(theta>=157.5 && theta<202.5)
	{
		p.x = -1; p.y = 0;
		return p;
	}
	//225
	if(theta>=202.5 && theta<247.5)
	{
		p.x = -1; p.y = -1;
		return p;
	}
	//270
	if(theta>=247.5 && theta<292.5)
	{
		p.x = 0; p.y = -1;
		return p;
	}
	//315
	if(theta>=292.5 && theta<337.5)
	{
		p.x = 1; p.y = -1;
		return p;
	}
	return p;
}


double Bezier_Fit::xprod(point_t p1, point_t p2)
{
	return p1.x*p2.y - p1.y*p2.x;
}
int Bezier_Fit::floordiv(int a, int n) {
	return (a>=0 ? a/n : -1-(-1-a)/n);
}

/* return 1 if a <= b < c < a, in a cyclic sense (mod n) */
int Bezier_Fit::cyclic(int a, int b, int c) {
	if (a<=c) {
		return (a<=b && b<c);
	} else {
		return (a<=b || b<c);
	}
}
int Bezier_Fit::sign(double x)
{
	return ((x)>0 ? 1 : (x)<0 ? -1 : 0);
}

/* return (p1-p0)x(p2-p0), the area of the parallelogram */
double  Bezier_Fit::dpara(point_t p0, point_t p1, point_t p2) {
	double x1, y1, x2, y2;

	x1 = p1.x-p0.x;
	y1 = p1.y-p0.y;
	x2 = p2.x-p0.x;
	y2 = p2.y-p0.y;

	return x1*y2 - x2*y1;
}

/* ddenom/dpara have the property that the square of radius 1 centered
at p1 intersects the line p0p2 iff |dpara(p0,p1,p2)| <= ddenom(p0,p2) */
double  Bezier_Fit::ddenom(point_t p0, point_t p2) {
	point_t r;
	r.y = sign(p2.x-p0.x);
	r.x = -sign(p2.y-p0.y);

	return r.y*(p2.x-p0.x) - r.x*(p2.y-p0.y);
}

double  Bezier_Fit::dot(point_t P1, point_t P2) {return P1.x*P2.x + P1.y*P2.y;}


//////////////////////////////////////////////////////////////////////////////////
/**
* Fill in bezier[] based on the given data and tangent requirements, using
* a least-squares fit.
*
* Each of tHat1 and tHat2 should be either a zero vector or a unit vector.
* If it is zero, then bezier[1 or 2] is estimated without constraint; otherwise,
* it bezier[1 or 2] is placed in the specified direction from bezier[0 or 3].
*
* \param tolerance_sq Used only for an initial guess as to tangent directions
*   when tHat1 or tHat2 is zero.
*/
void Bezier_Fit::generate_bezier(point_t bezier[],
								QVector<point_t> data, QVector<double> u, unsigned const len,
								const point_t &tHat1, const point_t &tHat2,
								double const tolerance_sq)
{
	bool const est1 = tHat1.is_zero();
	bool const est2 = tHat2.is_zero();
	point_t est_tHat1, est_tHat2;
	
	est_tHat1 =  est1
		? estimate_left_tangent(&data, len, tolerance_sq)
		: tHat1;
	est_tHat2 =  est2
		? estimate_right_tangent(&data, len, tolerance_sq)
		: tHat2;
	estimate_lengths(bezier, data, u, len, est_tHat1, est_tHat2);

	/* We find that sp_darray_right_tangent tends to produce better results
	for our current freehand tool than full estimation. */
	if (est1) {
		estimate_bi(bezier, 1, data, u, len);
		if (bezier[1] != bezier[0]) {
			est_tHat1 = (bezier[1] - bezier[0]);
			est_tHat1.normalize();
		}
		estimate_lengths(bezier, data, u,  len, est_tHat1, est_tHat2);
	}
}

/**
* Fit a multi-segment Bezier curve to a set of digitized points, without
* possible weedout of identical points and NaNs.
* 
* \pre data is uniqued, i.e. not exist i: data[i] == data[i + 1].
* \param max_beziers Maximum number of generated segments
* \param Result array, must be large enough for n. segments * 4 elements.
*/
int Bezier_Fit::bezier_fit_cubic(point_t bezier[], int split_points[], 
								QVector<point_t> data, int const len,
								point_t const &tHat1, point_t const &tHat2,
								double const error, unsigned const max_beziers)
{

	//Point bezier[4];
	//QVector<Point> data;
	QVector<double> u;
	//Point tHat1 = unconstrained_tangent; Point tHat2 = unconstrained_tangent;
	double dist;
	point_t difPoint;
	int const maxIterations = 4;   /* Max times to try iterating */
	int nsegs1, nsegs2;
	
	if ( len < 2 ) { 
		return 0;}

	if ( len == 2 ) {
		/* We have 2 points, which can be fitted by a line segment. */
		bezier[0] = data[0];
		bezier[3] = data[len - 1];
		/* Straight line segment. */
		difPoint = bezier[3] - bezier[0];
				
			bezier[1].x = bezier[2].x = (bezier[0].x + bezier[3].x)/2;
			bezier[1].y = bezier[2].y = bezier[0].y + (difPoint.y/(difPoint.x != 0 ? difPoint.x : 1))*(bezier[1].x - bezier[0].x);

		//curved Bezier
		/*	
		dist = (difPoint.L2_norm2D() / 3.0 );
		bezier[1] = ( tHat1.is_zero()
		? ( bezier[0] * 2 + bezier[3] ) / 3.
		: bezier[0] + tHat1 * dist);
		bezier[2] = ( tHat2.is_zero()
		? ( bezier[0] + bezier[3] * 2 ) / 3.
		: bezier[3] + tHat2 * dist );*/

		return 1;
	}

	/*  Parameterize points, and attempt to fit curve */
	unsigned splitPoint;   /* Point to split point set at. */
	bool is_corner;
	{
		chord_length_parameterize(&data, &u, len);
		if ( u[len - 1] == 0.0 ) {
			/* Zero-length path: every point in data[] is the same.
			*
			* (Clients aren't allowed to pass such data; handling the case is defensive
			* programming.)
			*/
			u.clear();
			return 0;
		}


		//std::cerr<<"\n        tHat1 = ("<<tHat1.x<<", "<<tHat1.y<<");\n";
		//std::cerr<<"        tHat2 = ("<<tHat2.x<<", "<<tHat2.y<<");\n\n";
		generate_bezier(bezier, data, u, len, tHat1, tHat2, error);
		reparameterize(data, len, u, bezier);

		/* Find max deviation of points to fitted curve. */
		double const tolerance = 1;//sqrt(error + 1e-9);
		double maxErrorRatio = compute_max_error_ratio(data, u, len, bezier, tolerance, &splitPoint);

		if ( fabs(maxErrorRatio) <= error ) {
			u.clear();
			return 1;
		}

		/* If error not too large, then try some reparameterization and iteration. */
		if ( 0.0 <= maxErrorRatio && maxErrorRatio <= 3.0*tolerance ) {
			for (int i = 0; i < maxIterations; i++) {
				generate_bezier(bezier, data, u, len, tHat1, tHat2, error);
				reparameterize(data, len, u, bezier);
				maxErrorRatio = compute_max_error_ratio(data, u, len, bezier, tolerance, &splitPoint);
				if ( fabs(maxErrorRatio) <= error ) {
					u.clear();
					return 1;
				}
			}
		}

		u.clear();
		is_corner = (maxErrorRatio < 0);
	}

	if (is_corner) {
		//std::cerr<<"\nis corner\n";
		assert(splitPoint < unsigned(len));
		if (splitPoint == 0) {
			if (tHat1.is_zero()) {
				/* Got spike even with unconstrained initial tangent. */
				++splitPoint;
			} else {
				int val = bezier_fit_cubic(bezier, split_points, data, len, Bezier_Fit::UNCONSTRAINED_TANGENT, tHat2,
					error, max_beziers);
				return val;
			}
		} else if (splitPoint == unsigned(len - 1)) {
			if (tHat2.is_zero()) {
				/* Got spike even with unconstrained final tangent. */
				--splitPoint;
			} else {
				int val;
				val = bezier_fit_cubic(bezier, split_points, data, len, tHat1, Bezier_Fit::UNCONSTRAINED_TANGENT,
					error, max_beziers);
				return val;
			}
		}
	}
	//return 1;
	if ( 1 < max_beziers ) {
		//std::cerr<<"RECURSIVE\n";
		/*
		*  Fitting failed -- split at max error point and fit recursively
		*/
		unsigned const rec_max_beziers1 = (splitPoint);

		point_t recTHat2, recTHat1;
		if (is_corner) {
			if(!(0 < splitPoint && splitPoint < unsigned(len - 1)))
			{
				return -1;
			}
			recTHat1 = recTHat2 = Bezier_Fit::UNCONSTRAINED_TANGENT;
		} 
		else {
			/* Unit tangent vector at splitPoint. */
			recTHat2 = estimate_center_tangent(&data, splitPoint, len);
			recTHat1 = -recTHat2;
		}
		////for(int i=0;i<split_points
		//std::cerr<<"        recTHat1 = ("<<recTHat1.x<<", "<<recTHat1.y<<");\n";
		//std::cerr<<"        recTHat2 = ("<<recTHat2.x<<", "<<recTHat2.y<<");\n\n";

		nsegs1 = bezier_fit_cubic(bezier, split_points, data, splitPoint + 1,
			tHat1, recTHat2, error, rec_max_beziers1);
		/*std::cerr<<"        nsegs1 = "<<nsegs1<<"; rec_max_beziers1 = "<<rec_max_beziers1<<"\n";*/
		if ( nsegs1 < 0 ) {
			return -1;
		}
		assert( nsegs1 != 0 );
		if (split_points != NULL) {
			split_points[nsegs1 - 1] = splitPoint;
		}
		unsigned const rec_max_beziers2 = max_beziers - splitPoint -1;

		QVector<point_t> data_temp = QVector<point_t>(data);
		data_temp.erase(data_temp.begin(), data_temp.begin()+splitPoint);
		int const nsegs2 = bezier_fit_cubic(bezier + nsegs1*4,
			( split_points == NULL ? NULL : split_points + nsegs1 ),
			data_temp, len - splitPoint,
			recTHat1, tHat2, error, rec_max_beziers2);
		//std::cerr<<"nsegs2 = "<<nsegs2<<"; rec_max_beziers2 = "<<rec_max_beziers2<<"\n";
		if ( nsegs2 < 0 ) {
			return -1;
		}
		return nsegs1 + nsegs2;

	}
	else {
		return 1;
	}
}




/**
* Given set of points and their parameterization, try to find a better assignment of parameter
* values for the points.
*
*  \param d  Array of digitized points.
*  \param u  Current parameter values.
*  \param bezCurve  Current fitted curve.
*  \param len  Number of values in both d and u arrays.
*              Also the size of the array that is allocated for return.
*/
void Bezier_Fit::reparameterize(QVector<point_t> d, unsigned const len, 
							   QVector<double> u, const point_t bezCurve[])
{
	assert( 2 <= len );

	unsigned const last = len - 1;
	assert( u[0] == 0.0 );
	assert( u[last] == 1.0 );
	/* Otherwise, consider including 0 and last in the below loop. */

	for (unsigned i = 1; i < last; i++) {
		u[i] = NewtonRaphsonRootFind(bezCurve, d[i], u[i]);
	}
}

/**
*  Use Newton-Raphson iteration to find better root.
*  
*  \param Q  Current fitted curve
*  \param P  Digitized point
*  \param u  Parameter value for "P"
*  
*  \return Improved u
*/
double Bezier_Fit::NewtonRaphsonRootFind(const point_t Q[], point_t const &P, double const u)
{
	assert( 0.0 <= u );
	assert( u <= 1.0 );

	/* Generate control vertices for Q'. */
	point_t Q1[3];
	for (unsigned i = 0; i < 3; i++) {
		Q1[i] = ( Q[i+1] - Q[i] )*3.0;
	}

	/* Generate control vertices for Q''. */
	point_t Q2[2];
	for (unsigned i = 0; i < 2; i++) {
		Q2[i] = ( Q1[i+1] - Q1[i] )*2.0;
	}

	/* Compute Q(u), Q'(u) and Q''(u). */
	point_t const Q_u  = bezier_pt(3, Q, u);
	point_t const Q1_u = bezier_pt(2, Q1, u);
	point_t const Q2_u = bezier_pt(1, Q2, u);

	/* Compute f(u)/f'(u), where f is the derivative wrt u of distsq(u) = 0.5 * the square of the
	distance from P to Q(u).  Here we're using Newton-Raphson to find a stationary point in the
	distsq(u), hopefully corresponding to a local minimum in distsq (and hence a local minimum
	distance from P to Q(u)). */
	point_t const diff = Q_u - P;
	double numerator = dot(diff, Q1_u);
	double denominator = dot(Q1_u, Q1_u) + dot(diff, Q2_u);

	double improved_u;
	if ( denominator > 0. ) {
		/* One iteration of Newton-Raphson:
		improved_u = u - f(u)/f'(u) */
		improved_u = u - ( numerator / denominator );
	} else {
		/* Using Newton-Raphson would move in the wrong direction (towards a local maximum rather
		than local minimum), so we move an arbitrary amount in the right direction. */
		if ( numerator > 0. ) {
			improved_u = u * .98 - .01;
		} else if ( numerator < 0. ) {
			/* Deliberately asymmetrical, to reduce the chance of cycling. */
			improved_u = .031 + u * .98;
		} else {
			improved_u = u;
		}
	}

	if (isnan(improved_u) || improved_u>=INFTY) {
		improved_u = u;
	} else if ( improved_u < 0.0 ) {
		improved_u = 0.0;
	} else if ( improved_u > 1.0 ) {
		improved_u = 1.0;
	}

	/* Ensure that improved_u isn't actually worse. */
	{
		double const diff_lensq = dot(diff,diff);
		for (double proportion = .125; ; proportion += .125) {
			point_t temp = bezier_pt(3, Q, improved_u) - P ;
			if ( dot(temp, temp)  > diff_lensq ) {
				if ( proportion > 1.0 ) {
					//g_warning("found proportion %g", proportion);
					improved_u = u;
					break;
				}
				improved_u = ( ( 1 - proportion ) * improved_u  +
					proportion         * u            );
			} else {
				break;
			}
		}
	}


	return improved_u;
}

/** 
* Evaluate a Bezier curve at parameter value \a t.
* 
* \param degree The degree of the Bezier curve: 3 for cubic, 2 for quadratic etc.
* \param V The control points for the Bezier curve.  Must have (\a degree+1)
*    elements.
* \param t The "parameter" value, specifying whereabouts along the curve to
*    evaluate.  Typically in the range [0.0, 1.0].
*
* Let s = 1 - t.
* BezierII(1, V) gives (s, t) * V, i.e. t of the way
* from V[0] to V[1].
* BezierII(2, V) gives (s**2, 2*s*t, t**2) * V.
* BezierII(3, V) gives (s**3, 3 s**2 t, 3s t**2, t**3) * V.
*
* The derivative of BezierII(i, V) with respect to t
* is i * BezierII(i-1, V'), where for all j, V'[j] =
* V[j + 1] - V[j].
*/
point_t Bezier_Fit::bezier_pt(unsigned const degree, point_t const V[], double const t)
{
	/** Pascal's triangle. */
	static int const pascal[4][4] = {{1},
	{1, 1},
	{1, 2, 1},
	{1, 3, 3, 1}};
	assert( degree < 4 );
	double const s = 1.0 - t;

	/* Calculate powers of t and s. */
	double spow[4];
	double tpow[4];
	spow[0] = 1.0; spow[1] = s;
	tpow[0] = 1.0; tpow[1] = t;
	for (unsigned i = 1; i < degree; ++i) {
		spow[i + 1] = spow[i] * s;
		tpow[i + 1] = tpow[i] * t;
	}

	point_t ret = V[0]*spow[degree];
	for (unsigned i = 1; i <= degree; ++i) {
		ret = ret + V[i] * pascal[degree][i] * spow[degree - i] * tpow[i];
	}
	return ret;
}


/**
* Find the maximum squared distance of digitized points to fitted curve, and (if this maximum
* error is non-zero) set \a *splitPoint to the corresponding index.
*
* \pre 2 \<= len.
* \pre u[0] == 0.
* \pre u[len - 1] == 1.0.
* \post ((ret == 0.0)
*        || ((*splitPoint \< len - 1)
*            \&\& (*splitPoint != 0 || ret \< 0.0))).
*/
double Bezier_Fit::compute_max_error_ratio(QVector<point_t> d, QVector<double> u, unsigned const len,
										  const point_t bezCurve[4], double const tolerance,
										  unsigned *const splitPoint)

{
	assert( 2 <= len );
	unsigned const last = len - 1;
	assert( u[0] == 0.0 );
	assert( u[last] == 1.0 );
	/* I.e. assert that the error for the first & last points is zero.
	* Otherwise we should include those points in the below loop.
	* The assertion is also necessary to ensure 0 < splitPoint < last.
	*/

	double maxDistsq = 0.0; /* Maximum error */
	double max_hook_ratio = 0.0;
	unsigned snap_end = 0;
	point_t prev = bezCurve[0];
	for (unsigned i = 1; i <= last; i++) {
		point_t const curr = bezier_pt(3, bezCurve, u[i]);
		double const distsq = dot( curr - d[i],  curr - d[i] );
		if ( distsq > maxDistsq ) {
			maxDistsq = distsq;
			*splitPoint = i;
		}
		double const hook_ratio = compute_hook(prev, curr, .5 * (u[i - 1] + u[i]), bezCurve, tolerance);
		
		if (max_hook_ratio < hook_ratio) {
			max_hook_ratio = hook_ratio;
			snap_end = i;
		}
		prev = curr;
	}
	double const dist_ratio = sqrt(maxDistsq) / tolerance;
	double ret;
	if (max_hook_ratio <= dist_ratio) {
		ret = dist_ratio;
	} else {
		assert(0 < snap_end);
		ret = -max_hook_ratio;
		*splitPoint = snap_end - 1;
	}
	assert( ret == 0.0
		|| ( ( *splitPoint < last )
		&& ( *splitPoint != 0 || ret < 0. ) ) );
	return ret;
}


/** 
* Whereas compute_max_error_ratio() checks for itself that each data point 
* is near some point on the curve, this function checks that each point on 
* the curve is near some data point (or near some point on the polyline 
* defined by the data points, or something like that: we allow for a
* "reasonable curviness" from such a polyline).  "Reasonable curviness" 
* means we draw a circle centred at the midpoint of a..b, of radius 
* proportional to the length |a - b|, and require that each point on the 
* segment of bezCurve between the parameters of a and b be within that circle.
* If any point P on the bezCurve segment is outside of that allowable 
* region (circle), then we return some metric that increases with the 
* distance from P to the circle.
*
*  Given that this is a fairly arbitrary criterion for finding appropriate 
*  places for sharp corners, we test only one point on bezCurve, namely 
*  the point on bezCurve with parameter halfway between our estimated 
*  parameters for a and b.  (Alternatives are taking the farthest of a
*  few parameters between those of a and b, or even using a variant of 
*  NewtonRaphsonFindRoot() for finding the maximum rather than minimum 
*  distance.)
*/
double Bezier_Fit::compute_hook(point_t const &a, point_t const &b, double const u, const point_t bezCurve[4],
							   double const tolerance)
{
	point_t P = bezier_pt(3, bezCurve, u);
	point_t const diff = (a + b)* 0.5 - P;
	double const dist = diff.L2_norm2D();
	if (dist < tolerance) {
		return 0;
	}
	P = b-a;
	double const allowed = P.L2_norm2D() + tolerance;
	return dist / allowed;
	/** \todo 
	* effic: Hooks are very rare.  We could start by comparing 
	* distsq, only resorting to the more expensive L2 in cases of 
	* uncertainty.
	*/
}
