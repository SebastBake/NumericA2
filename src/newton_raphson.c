/***************************************************************************
 *
 *   File        : newton_raphson.c
 *   Student Id  : 757931
 *   Name        : Sebastian Baker
 *
 ***************************************************************************/

#include "newton_raphson.h"

// simple newton raphson implementation
int newtonRaphson(double* xi, double* params, double (*f)(double, double*)) {

    double fx = f((*xi), params);

    int i=0;
    while( (i < MAX_ITERATIONS) && (fabs(fx) > EPS) ) {

        fx = f((*xi), params);
        (*xi) = (*xi) - fx / dfdx((*xi), params, f);

        i++;
    }

    if (fabs(fx) > EPS) { 
        return ROOTFIND_FAIL;
    }

    return i;
}

// simple derivative approximation
double dfdx(double x, double* params, double (*f)(double, double*)) {

    double df = f(x+DX/2, params) - f(x-DX/2, params);
    return df/(DX);
}