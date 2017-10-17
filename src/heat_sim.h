/***************************************************************************
 *
 *   File        : heat_sim.h
 *   Student Id  : 757931
 *   Name        : Sebastian Baker
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "thomas_alg.h"

#ifndef HEAT_SIM_H

#define EXPLICIT_RHS(mu, dx, f_0, f_1, f_2) mu*(f_2 - 2.0*f_1 + f_0)/(dx*dx)

// heat simulation cell
typedef struct heat_cell {

    double f;
    double x;
    double t;

} sim_cell_t;

// heat simulation struct
typedef struct heat_sim {

    sim_cell_t** cells; // matrix of cells, in the form f(x,t) = cells[i_x][j_t]
    double mu;
    double x_lo;
    double x_hi;
    double t_lo;
    double t_hi;
    int Nx;  // not including 0th x value
    int Nt; // not including 0th t value

} heat_sim_t;

heat_sim_t* newUnsolvedHeatSim(
    int Nx, 
    int Nt, 
    double x_lo, 
    double x_hi,
    double t_lo,
    double t_hi,
    double mu, 
    double (*init)(double x) 
);

void freeHeatSim(heat_sim_t* sim);
void eulerExFe(heat_sim_t* sim);
void eulerExVe(heat_sim_t* sim);
void eulerImFe(heat_sim_t* sim);


#endif 