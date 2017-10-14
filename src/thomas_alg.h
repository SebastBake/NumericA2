/***************************************************************************
 *
 *   File        : thomas_alg.h
 *   Student Id  : 757931
 *   Name        : Sebastian Baker
 *
 ***************************************************************************/

#include <math.h>

#ifndef THOMAS_ALG_H
#define THOMAS_ALG_H

#define COMPUTE_A_S_START 1
#define COMPUTE_Q_S_START 1
#define COMPUTE_x_END 1

#define THOMAS_INIT_ARR_LEN 20
#define TINY(x) fabs(x) < 1e-11
#define SOLVER_FAIL -1
#define SOLVER_SUCCESS 1

typedef struct tridiag_row_struct {

    double a;
    double b;
    double c;
    double Q;
    double a_s;
    double Q_s;
    double x;

} tridiag_row_t;

typedef struct tridiag_struct {

    tridiag_row_t *rows;
    int N;
    int arrLen;

} tridiag_t;

tridiag_t *newTridiag();
void freeTridiag(tridiag_t *m);
void appendTridiagRow(tridiag_t *m, double a, double b, double c, double Q);
tridiag_row_t *getTridiagRow(tridiag_t *m, int i); //retrieve the ith row
int solveTridiag(tridiag_t *m);

int computeTridiag_a_s(tridiag_t *m);
int computeTridiag_Q_s(tridiag_t *m);
int computeTrigiag_x(tridiag_t *m);

#endif