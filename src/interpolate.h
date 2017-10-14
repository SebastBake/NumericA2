/***************************************************************************
 *
 *   File        : interpolate.h
 *   Student Id  : 757931
 *   Name        : Sebastian Baker
 *
 ***************************************************************************/

#ifndef INTERPOLATE_H

#define INTERP_INIT_ARRLEN 20
#define CUB_SPLINE_COMPUTE_SUCCESS 1
#define CUB_SPLINE_COMPUTE_FAIL -1

typedef struct interp_pt {

	double x;
	double fx;

} interp_pt_t;

typedef struct interp_set {

	interp_pt_t** pts;
	int N;
	int arrLen;

} interp_set_t;

typedef struct lagrange_term {

	int index;
	int order;
	double* root;
	double fx_i;
	double denominator;

} lagrange_term_t;

typedef struct lagrange_eqn {

	int num_terms;
	lagrange_term_t** terms;

} lagrange_eqn_t;

typedef struct cub_spline_segment {

	int index;
	double x_lo;
	double x_hi;
	double a;
	double b;
	double c;
	double d;

} cub_spline_seg_t;

typedef struct cub_spline {

	cub_spline_seg_t** segs;
	int num_segs;

} cub_spline_t;

interp_set_t* newInterpSet();
interp_pt_t* newInterpPt(double x, double fx);
void freeInterpSet(interp_set_t* set);
void freeInterpPt(interp_pt_t* pt);
void appendPtToSet(interp_set_t* set, interp_pt_t* pt);

lagrange_term_t* newLagrangeTerm(interp_set_t* set, int index);
lagrange_eqn_t* newLagrangeEqn(interp_set_t* set);
void freeLagrangeEqn(lagrange_eqn_t* eqn);
void freeLagrangeTerm(lagrange_term_t* term);
interp_pt_t* evaluateLagrangeEqn(lagrange_eqn_t* eqn, double x);
double evaluateLagrangeTerm(lagrange_term_t* term, double x);

cub_spline_t* newCubSpline(interp_set_t* set);
cub_spline_seg_t* newEmptyCubSplineSegment(int index, double x_lo, double x_hi);
void freeCubSpline(cub_spline_t* spline);
void freeCubSplineSegment(cub_spline_seg_t* seg);
int computeCubSplineConstants(interp_set_t* set, cub_spline_t* spline);
int computeCubSplineAs(interp_set_t* set, cub_spline_t* spline);
int computeCubSplineBs(interp_set_t* set, cub_spline_t* spline);
int computeCubSplineCs(interp_set_t* set, cub_spline_t* spline);
int computeCubSplineDs(interp_set_t* set, cub_spline_t* spline);
interp_pt_t* evaluateCubSpline(cub_spline_t* spline, double x);
double evaluateCubSplineSegment(cub_spline_seg_t* seg, double x);

#endif