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

#define Q2_P2_M_ARR_LEN 20
#define Q2_HEADER_ITEMS 5
#define FILE_READONLY "r"
#define CSV_NEWLINE '\n'
#define CSV_NULLBYTE '\0'

typedef struct q2_params {
	double M;
	double theta;
	double beta_l;
	double beta_u;
	double gamma;
	double *p2_M;
	int M_n;
} q2_params_t;

q2_params_t *initQ2Params(const char* q2_filename);
void freeQ2Params(q2_params_t *q2Params);

void shockwave(const char* q2_filename) {

	q2_params_t *q2Params = initQ2Params(q2_filename);
	
	freeQ2Params(q2Params);
}

void linalgbsys(const char* q4_filename) {
	printf("linalgbsys() - IMPLEMENT ME!\n");
}

void interp(const char* q5_filename, const double xo) {
	printf("interp() - IMPLEMENT ME!\n");
}

void heateqn(const char* q6_filename) {
	printf("heateqn() - IMPLEMENT ME!\n");
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  */
q2_params_t *initQ2Params(const char* q2_filename) {
	
	// initialise params struct
	q2_params_t *q2Params = (q2_params_t*)malloc(sizeof(q2_params_t));
	assert(q2Params != NULL);
	q2Params->p2_M = (double*)calloc(Q2_P2_M_ARR_LEN, sizeof(double));
	assert(q2Params->p2_M != NULL);

	// open file
	FILE* q2File = fopen(q2_filename, FILE_READONLY);
	assert(q2File != NULL);

	// parse 2nd line (for part 1 parameters)
	char skip = CSV_NULLBYTE;
	while(skip != CSV_NEWLINE) { skip = fgetc(q2File); } // skip line
	assert(Q2_HEADER_ITEMS == fscanf(
		q2File,
		"%lf,%lf,%lf,%lf,%lf",
		&(q2Params->M),
		&(q2Params->theta),
		&(q2Params->beta_l),
		&(q2Params->beta_u),
		&(q2Params->gamma)
	) );
	// printf("%f,%f,%f,%f,%f\n",
	// (q2Params->M),
	// (q2Params->theta),
	// (q2Params->beta_l),
	// (q2Params->beta_u),
	// (q2Params->gamma));

	// parse M values for part 2
	skip = CSV_NULLBYTE;
	while(skip != CSV_NEWLINE) { skip = fgetc(q2File); } 
	fgetc(q2File);

	q2Params->M_n=0;
	int arrayLength = Q2_P2_M_ARR_LEN;
	while(1 == fscanf(q2File,"%lf\n", &(q2Params->p2_M[q2Params->M_n++])) ) {

		//printf("M: %f\n", q2Params->p2_M[q2Params->M_n -1]);
		if(q2Params->M_n >= arrayLength-1) {
			arrayLength += Q2_P2_M_ARR_LEN;
			q2Params->p2_M = (double*)realloc(q2Params->p2_M, arrayLength*sizeof(double));
			assert(q2Params->p2_M != NULL);
		}
	}

	// cleanup
	fclose(q2File);

	return q2Params;
}

void freeQ2Params(q2_params_t *q2Params) {
	
	free(q2Params->p2_M);
	free(q2Params);
}
