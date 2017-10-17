/***************************************************************************
 *
 *   File        : heat_sim.c
 *   Student Id  : 757931
 *   Name        : Sebastian Baker
 *
 ***************************************************************************/

#include "heat_sim.h"

// creates a new unsolved heat simulation struct
heat_sim_t* newUnsolvedHeatSim(
	int Nx, 
	int Nt, 
	double x_lo, 
	double x_hi,
	double t_lo,
	double t_hi,
	double mu, 
	double (*init)(double x) 
) {

	assert(init!=NULL);

	heat_sim_t* sim = (heat_sim_t*)malloc(sizeof(heat_sim_t));
	assert(sim!=NULL);

	sim->mu = mu;
	sim->Nx = Nx;
	sim->Nt = Nt;
	sim->x_lo = x_lo;
	sim->x_hi = x_hi;
	sim->t_lo = t_lo;
	sim->t_hi = t_hi;

	sim->cells = (sim_cell_t**)calloc( sim->Nx+1, sizeof(sim_cell_t*));
	assert(sim->cells!=NULL);
	
	int i=0, j=0;
	double x_i=0, t_j=0;
	double dx = (x_hi-x_lo)/Nx;
	double dt = (t_hi-t_lo)/Nt;

	// initialise the matrix at t = 0
	for(i=0; i<=sim->Nx; i++) {
		sim->cells[i] = (sim_cell_t*)calloc( sim->Nt+1, sizeof(sim_cell_t));
		assert(sim->cells[i]!=NULL);

		x_i = x_lo + i*dx;
		(sim->cells[i][0]).f = init(x_i);

		for(j=1; j<=sim->Nt; j++) {

			t_j = t_lo + j*dt;
			(sim->cells[i][j]).x = x_i;
			(sim->cells[i][j]).t = t_j;
		}
	}
	return sim;
}

// free heat simulation struct
void freeHeatSim(heat_sim_t* sim) {

	assert(sim!=NULL);
	
	int i=0;
	for(i=0; i<=sim->Nx; i++) {
		free(sim->cells[i]);
	}
	free(sim->cells);
	free(sim);
}

// solve heat eqn simulation using explicit euler with fixed end points
void eulerExFe(heat_sim_t* sim) {

	assert(sim!=NULL);

	// let i be x index, j is t index
	int j=0,i=0;
	double dx = (sim->x_hi-sim->x_lo)/sim->Nx;
	double dt = (sim->t_hi-sim->t_lo)/sim->Nt;
	double mu = sim->mu;
	double rhs;
	sim_cell_t** m = sim->cells; // cell matrix

	for(j=0; j<sim->Nt; j++) { // each timestep
		for(i=0; i<=sim->Nx; i++) { // each x value

			// get RHS (fixed end points)
			if(i==0 || i == sim->Nx) {
				rhs = 0.0;
			} else {
				rhs = EXPLICIT_RHS(mu,dx, m[i-1][j].f, m[i][j].f, m[i+1][j].f);
			}

			// use RHS for explicit euler calculation
			m[i][j+1].f = m[i][j].f + dt * rhs;
		}
	}
}

// solve heat eqn simulation using explicit euler with variable end points
void eulerExVe(heat_sim_t* sim) {

	assert(sim!=NULL);

	// let i be x index, j is t index
	int j=0,i=0;
	double dx = (sim->x_hi-sim->x_lo)/sim->Nx;
	double dt = (sim->t_hi-sim->t_lo)/sim->Nt;
	double rhs;
	sim_cell_t** m = sim->cells; // cell matrix

	for(j=0; j<sim->Nt; j++) { // each timestep
		for(i=0; i<=sim->Nx; i++) { // each x value

			// get RHS (fixed end points)
			if(i==0) {
				rhs = EXPLICIT_RHS(sim->mu,dx, m[i][j].f, m[i+1][j].f, m[i+2][j].f);
			} else if(i == sim->Nx) {
				rhs = EXPLICIT_RHS(sim->mu,dx, m[i-2][j].f, m[i-1][j].f, m[i][j].f);
			} else {
				rhs = EXPLICIT_RHS(sim->mu,dx, m[i-1][j].f, m[i][j].f, m[i+1][j].f);
			}

			// use RHS for explicit euler calculation
			m[i][j+1].f = m[i][j].f + dt * rhs;
		}
	}
}

// solve heat eqn simulation using implicit euler with fixed end points
void eulerImFe(heat_sim_t* sim) {

	assert(sim!=NULL);

	// let i be x index, j is t index
	int j=0,i=0;

	double dx = (sim->x_hi-sim->x_lo)/sim->Nx;
	double dt = (sim->t_hi-sim->t_lo)/sim->Nt;
	double d = (sim->mu * dt)/(dx*dx);

	sim_cell_t** m = sim->cells; // cell matrix
	tridiag_t* tri;

	for(j=1; j<=sim->Nt; j++) { // each timestep

		tri = newTridiag();

		// add rows to tridiagonal matrix
		appendTridiagRow(tri, 1, 0, 0, 0);
		for(i=1; i<sim->Nx; i++) {
			appendTridiagRow(tri, 2*d+1, -d, -d, m[i][j-1].f);
		}
		appendTridiagRow(tri, 1, 0, 0, 0);

		assert(solveTridiag(tri) == SOLVER_SUCCESS);

		for(i=0; i<=sim->Nx; i++) {
			tridiag_row_t *r = getTridiagRow(tri, i+1);
			m[i][j].f = r->x;
		}

		freeTridiag(tri);
	}
}