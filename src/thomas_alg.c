/***************************************************************************
 *
 *   File        : thomas_alg.c
 *   Student Id  : 757931
 *   Name        : Sebastian Baker
 *
 ***************************************************************************/

#include "thomas_alg.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

tridiag_t *newTridiag() {
	
	tridiag_t *m = (tridiag_t*)malloc(sizeof(tridiag_t));
	assert(m != NULL);

	m->arrLen = INIT_ARR_LEN;
	m->N = 0;

	m->rows = (tridiag_row_t*)calloc( m->arrLen, sizeof(tridiag_row_t));
	assert(m->rows != NULL);

	return m;
}

void freeTridiag(tridiag_t *m) {
	
	assert(m!=NULL);
	free(m->rows);
	free(m);
}

void appendTridiagRow(tridiag_t *m, double a, double b, double c, double Q) {

	assert(m!=NULL);	

	m->N++;
	if(m->N == m->arrLen) {
		m->arrLen += INIT_ARR_LEN;
		m->rows = (tridiag_row_t*)realloc( m->rows, m->arrLen*sizeof(tridiag_row_t));
		assert(m->rows!=NULL);
	}

	tridiag_row_t *r = getTridiagRow(m, m->N);
	r->a = a;
	r->b = b;
	r->c = c;
	r->Q = Q;
}

tridiag_row_t *getTridiagRow(tridiag_t *m, int i) {
	assert(m!=NULL);
	assert(i<=m->N);
	return &(m->rows[i-1]);
}

int solveTridiag(tridiag_t *m) {

	assert(m!=NULL);
	
	if(
		computeTridiag_a_s(m) == SOLVER_SUCCESS &&
		computeTridiag_Q_s(m) == SOLVER_SUCCESS &&
		computeTrigiag_x(m) == SOLVER_SUCCESS
	) {
		return SOLVER_SUCCESS;
	}

	return SOLVER_FAIL;
}

int computeTridiag_a_s(tridiag_t *m) {

	assert(m!=NULL);

	// case for row 1 : a_s = a
	int i=COMPUTE_A_S_START;
	tridiag_row_t *r = getTridiagRow(m, i);	
	tridiag_row_t *rPrev = r;
	r->a_s = r->a;

	// case for row 2,3,4...N : a_s = a − c * bPrev / a_sPrev
	while(++i <= m->N) {
		rPrev = r;
		r = getTridiagRow(m, i);
		if(TINY(rPrev->a_s)) { return SOLVER_FAIL; }
		r->a_s = r->a - r->c * rPrev->b / rPrev->a_s;
	}

	return SOLVER_SUCCESS;
}

int computeTridiag_Q_s(tridiag_t *m) {
	
	assert(m!=NULL);

	// case for row 1 : Q_s = Q
	int i=COMPUTE_Q_S_START;
	tridiag_row_t *r = getTridiagRow(m, i);	
	tridiag_row_t *rPrev = r;
	r->Q_s = r->Q;

	// case for row 2,3,4...N : Q − c * Q_sPrev /a_sPrev
	while(++i <= m->N) {
		rPrev = r;
		r = getTridiagRow(m, i);
		if(TINY(rPrev->a_s)) { return SOLVER_FAIL; }
		r->Q_s = r->Q - r->c * rPrev->Q_s / rPrev->a_s;
	}

	return SOLVER_SUCCESS;
}

int computeTrigiag_x(tridiag_t *m) {
	
	assert(m!=NULL);
	
	// case for row N : x = Q_s/a_s
	int i=m->N;
	tridiag_row_t *r = getTridiagRow(m, i);	
	tridiag_row_t *rPrev = r;
	r->x = r->Q_s/r->a_s;

	// case for row N-1, N-2, N-3...1 : x = ( Q_s − b * xPrev ) / a_s
	while(--i >= COMPUTE_x_END) {
		rPrev = r;
		r = getTridiagRow(m, i);
		if(TINY(rPrev->a_s)) { return SOLVER_FAIL; }
		r->x = (r->Q_s - r->b * rPrev->x) / r->a_s;
	}

	return SOLVER_SUCCESS;
}