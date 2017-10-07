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

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * TASK FUNCTIONS */

void shockwave(const char* q2_filename) {

	input_2_t *parsed = parseInput_2(q2_filename);
	shockwave_2a(parsed);
	shockwave_2b(parsed);
	shockwave_2c(parsed);
	freeInput_2(parsed);
}

void linalgbsys(const char* q4_filename) {
	
}

void interp(const char* q5_filename, const double xo) {
	printf("interp() - IMPLEMENT ME!\n");
}

void heateqn(const char* q6_filename) {
	printf("heateqn() - IMPLEMENT ME!\n");
}

/* * * * * * * * * * * * * * * * * * * * * * * * LINALGBYSYS HELPER FUNCTIONS */



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
	parsed->NumM_c = i;
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
	for(i=0; i<parsed->NumM_c; i++) {
		double params[] = { parsed->M_c[i], RADIAN_START_2B, parsed->g_a };
		shockwavePrintThetaRange_2bc( fp, params);
	}

	fclose(fp);
}

