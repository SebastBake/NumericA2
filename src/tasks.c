/***************************************************************************
 *
 *   File        : tasks.c
 *   Student Id  : 757931
 *   Name        : Sebastian Baker
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <assert.h>
#include "tasks.h"
#include "newton_raphson.h"
#include "thomas_alg.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * TASK FUNCTIONS */

void shockwave(const char* filename) {

	input_2_t *parsed = parseInput_2(filename);
	shockwave_2a(parsed);
	shockwave_2b(parsed);
	shockwave_2c(parsed);
	freeInput_2(parsed);
}

void linalgbsys(const char* filename) {

	tridiag_t *m = parseInput_3(filename);
	assert(solveTridiag(m) == SOLVER_SUCCESS);
	printTridiag_3(m);
	freeTridiag(m);
}

void interp(const char* filename, const double xo) {
	
	interp_set_t* set = parseInput_5(filename);
	lagrange_eqn_t* lagEqn = newLagrangeEqn(set);

	interp_pt_t* lagEval;

	double i=0;
	for(i=0; i<20; i++) {
		lagEval = evaluateLagrangeEqn(lagEqn, i);
		printf("evaluated: %f, %f\n", lagEval->x, lagEval->fx);
		freeInterpPt(lagEval);
	}
	
	freeLagrangeEqn(lagEqn);
	freeInterpSet(set);
}

void heateqn(const char* filename) {
	printf("heateqn() - IMPLEMENT ME!\n");
}

/* * * * * * * * * * * * * * * * * * * * * * * * INTERP HELPER FUNCTIONS */

interp_set_t* parseInput_5(const char* filename) {

	assert(filename!=NULL);

	interp_set_t* newSet = newInterpSet(); 

	// open file
	FILE* fp = fopen(filename, FILE_READONLY);
	assert(fp != NULL);

	double tmp_x, tmp_fx;
	fscanf(fp, "x,f(x)\n");
	int read = 0;
	int i=0;
	while(1) {
		read = fscanf(fp,"%lf,%lf\n", &tmp_x, &tmp_fx);
		if(read != NUM_PARAMS_5) { break; }
		appendPtToSet(newSet, newInterpPt(tmp_x, tmp_fx) );
		printf("set: %f, %f\n", (newSet->pts[i])->x, (newSet->pts[i])->fx);
		i++;
	}

	fclose(fp);

	return newSet;
}


/* * * * * * * * * * * * * * * * * * * * * * * * LINALGBYSYS HELPER FUNCTIONS */

tridiag_t *parseInput_3(const char* filename) {
	
	assert(filename!=NULL);

	tridiag_t *m = newTridiag();

	// open file
	FILE* fp = fopen(filename, FILE_READONLY);
	assert(fp != NULL);

	double tmp_a, tmp_b, tmp_c, tmp_q;
	fscanf(fp, "a,b,c,q\n");
	int read = 0;
	while(1) {
		read = fscanf(fp,"%lf,%lf,%lf,%lf\n", &tmp_a, &tmp_b, &tmp_c, &tmp_q);
		if(read != NUM_PARAMS_3) { break; }
		appendTridiagRow(m,tmp_a, tmp_b, tmp_c, tmp_q);
	}

	fclose(fp);
	return m;
}

void printTridiag_3(tridiag_t *m) {

	assert(m!=NULL);

	FILE* fp = fopen(FILENAME_3, FILE_OVERWRITE);
	assert(fp != NULL);

	fprintf(fp, FILE_HEADER_3);
	int i=1;
	for(i=1; i <= m->N; i++) {
		tridiag_row_t *r = getTridiagRow(m, i);
		fprintf(fp,"%f\n",r->x);
	}

	fclose(fp);
}


/* * * * * * * * * * * * * * * * * * * * * * * * * SHOCKWAVE HELPER FUNCTIONS */

input_2_t *parseInput_2(const char* filename) {
	
	assert(filename!=NULL);
	
	// allocate memory
	input_2_t *parsed = (input_2_t*)malloc(sizeof(input_2_t));
	assert(parsed != NULL);
	parsed->M_c = (double*)calloc(M_LEN_2C, sizeof(double));
	assert(parsed->M_c != NULL);

	// open file
	FILE* fp = fopen(filename, FILE_READONLY);
	assert(fp != NULL);

	// parse file
	parseInput2ndLine_2(fp, parsed);
	parseInputMvals_2(fp, parsed);

	fclose(fp);
	return parsed;
}

void parseInput2ndLine_2(FILE* fp, input_2_t* parsed) {
	
	assert(parsed!=NULL && fp!= NULL);

	assert(NUM_PARAMS_2 == fscanf(
		fp,
		"M,theta,beta_l,beta_u,gamma\n%lf,%lf,%lf,%lf,%lf\nM\n",
		&(parsed->M_a),
		&(parsed->t_a),
		&(parsed->b_l_a),
		&(parsed->b_u_a),
		&(parsed->g_a)
	) );
}

void parseInputMvals_2(FILE* fp, input_2_t* parsed) {
	
	assert(parsed!=NULL && fp!= NULL);

	int i=0;
	int arrayLength = M_LEN_2C;
	while( 1 == fscanf(fp,"%lf", &(parsed->M_c[i])) ) {

		// Extend array if necessary
		if(i >= arrayLength) {
			arrayLength += M_LEN_2C;
			parsed->M_c = (double*)realloc(parsed->M_c, arrayLength*sizeof(double));
			assert(parsed->M_c != NULL);
		}
		i++;
	}
	parsed->num_M_c = i;
}

void freeInput_2(input_2_t *parsed) {

	assert(parsed!=NULL);
	free(parsed->M_c);
	free(parsed);
}

double f_2(double beta, double* params) {

	assert(params!=NULL);
	return F_2( beta, params[F_M_2], params[F_T_2], params[F_G_2]);
}

void shockwave_2a(const input_2_t *parsed) {

	assert(parsed!=NULL);

	double b_l = DEG2RAD(parsed->b_l_a);
	double b_u = DEG2RAD(parsed->b_u_a);
	double params[] = {
		parsed->M_a,
		DEG2RAD(parsed->t_a),
		parsed->g_a
	};

	int b_l_itr = newtonRaphson(&b_l, params, f_2);
	int b_u_itr = newtonRaphson(&b_u, params, f_2);
	assert( b_l_itr != ROOTFIND_FAIL && b_u_itr != ROOTFIND_FAIL);

	//printf("shockwave_2a: b_l=%.5f (Should be 29.80092), after %d iterations\n", RAD2DEG(b_l), b_l_itr);
	//printf("shockwave_2a: b_u=%.5f (Should be 84.55625), after %d iterations\n", RAD2DEG(b_u), b_u_itr);
}

void shockwavePrintThetaRange_2bc(FILE* fp, double* params) {

	int b_u_itr=0, b_l_itr=0;
	double b_l = DEG2RAD(B_l_GUESS_2), b_u = DEG2RAD(B_u_GUESS_2);

	while(1) {
		b_l_itr = newtonRaphson(&b_l, params, f_2);
		b_u_itr = newtonRaphson(&b_u, params, f_2);
		if( b_l_itr == ROOTFIND_FAIL || b_u_itr == ROOTFIND_FAIL) { break; }
		
		fprintf(fp, "%.1f,%.0f,%.6f,%.6f\n",
			params[F_M_2],
			RAD2DEG(params[F_T_2]),
			RAD2DEG(b_l),
			RAD2DEG(b_u)
		);

		params[F_T_2] += RADIAN_INCREMENT_2B;
	}
}

void shockwave_2b(const input_2_t *parsed) {

	assert(parsed!=NULL);

	double params[] = { parsed->M_a, RADIAN_START_2B, parsed->g_a };

	FILE* fp = fopen(FILENAME_2B, FILE_OVERWRITE);
	assert(fp != NULL);

	fprintf(fp,FILE_HEADER_2BC);
	shockwavePrintThetaRange_2bc( fp, params );

	fclose(fp);
}

void shockwave_2c(const input_2_t *parsed) {

	assert(parsed!=NULL);	

	FILE* fp = fopen(FILENAME_2C, FILE_OVERWRITE);
	assert(fp != NULL);

	fprintf(fp,FILE_HEADER_2BC);
	
	int i=0;
	for(i=0; i<parsed->num_M_c; i++) {
		double params[] = { parsed->M_c[i], RADIAN_START_2B, parsed->g_a };
		shockwavePrintThetaRange_2bc( fp, params);
	}

	fclose(fp);
}

