/***************************************************************************
 *
 *   File        : tasks.h
 *   Student Id  : 757931
 *   Name        : Sebastian Baker
 *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "thomas_alg.h"
#include "interpolate.h"

#ifndef TASKS_H

#define FILE_READONLY "r"
#define FILE_OVERWRITE "w"

#define DEG2RAD(deg) deg * M_PI / 180.0
#define RAD2DEG(rad) rad * 180.0 / M_PI

#define NUM_PARAMS_2 5
#define F_M_2 0
#define F_T_2 1
#define F_G_2 2
#define F_2(b,M,t,g) (2.0/tan(b)) * (pow(M*sin(b), 2.0) - 1.0) / (M*M*(g + cos(2.0*b)) + 2.0) - tan(t)
#define B_l_GUESS_2 20.0
#define B_u_GUESS_2 90.0
#define FILENAME_2B "my_2b.csv"
#define FILENAME_2C "out_shock.csv"
#define FILE_HEADER_2BC "M,theta,beta_lower,beta_upper\n"
#define RADIAN_START_2B 0
#define RADIAN_INCREMENT_2B DEG2RAD(1)
#define M_LEN_2C 20

#define NUM_PARAMS_3 4
#define FILENAME_3 "out_linalsys.csv"
#define FILE_HEADER_3 "x\n"

#define NUM_PARAMS_5 2
#define PLOT_FILENAME_5 "my_5.csv"
#define FILENAME_5 "out_interp.csv"
#define PLOT_START_5 0
#define PLOT_END_5 8
#define PLOT_INTERVAL_5 0.01
#define PLOT_HEADER_LAGRANGE_5 "x,lagrange,spline\n"
#define HEADER_LAGRANGE_5 "lagrange\n"
#define HEADER_CUBIC_5 "cubic\n"

// Input type for shockwave question (Question 2)
typedef struct input_2 {

	double M_a;
	double t_a;
	double g_a;
	double b_l_a;
	double b_u_a;
	double *M_c;
	int num_M_c;

} input_2_t;

/* * * * * * * * * * * * * * * * * * * * * * * * * *  INTERP HELPER FUNCTIONS */

interp_set_t* parseInput_5(const char* filename);
void printInterp_5(double lag, double spline);
void plotInterp_5(lagrange_eqn_t* lagEqn, cub_spline_t* splineEqn);

/* * * * * * * * * * * * * * * * * * * * * * * * LINALGBYSYS HELPER FUNCTIONS */

tridiag_t *parseInput_3(const char* filename);
void printTridiag_3(tridiag_t *m);


/* * * * * * * * * * * * * * * * * * * * * * * * * SHOCKWAVE HELPER FUNCTIONS */

input_2_t *parseInput_2(const char* filename);
void parseInput2ndLine_2(FILE* fp, input_2_t* parsed);
void parseInputMvals_2(FILE* fp, input_2_t* parsed);
void freeInput_2(input_2_t *parsed);
double f_2(double beta, double* params);
void shockwave_2a(const input_2_t *parsed);
void shockwave_2b(const input_2_t *parsed);
void shockwave_2c(const input_2_t *parsed);


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * TASK FUNCTIONS */

void shockwave(const char* q2_file);
void linalgbsys(const char* q4_file);
void interp(const char* q5_file, const double xo);
void heateqn(const char* q6_file);

#endif
