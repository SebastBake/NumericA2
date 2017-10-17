/***************************************************************************
 *
 *   File        : newton_raphson.h
 *   Student Id  : 757931
 *   Name        : Sebastian Baker
 *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef NEWTON_RAPHSON_H
#define NEWTON_RAPHSON_H

#define MAX_ITERATIONS 40
#define EPS 1e-9
#define DX EPS/1000
#define ROOTFIND_FAIL -1

int newtonRaphson(double* xi, double* params, double (*f)(double, double*));
double dfdx(double x, double* params, double (*f)(double, double*));

#endif