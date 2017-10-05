/***************************************************************************
 *
 *   File        : newton_raphson.h
 *   Student Id  : 757931
 *   Name        : Sebastian Baker
 *
 ***************************************************************************/

#ifndef NEWTON_RAPHSON_H
#define NEWTON_RAPHSON_H

double newtonRaphson(double start, double* params, double (*f)(double, double*));

#endif