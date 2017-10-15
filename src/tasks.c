/***************************************************************************
 *
 *   File        : tasks.c
 *   Student Id  : 757931
 *   Name        : Sebastian Baker
 *
 ***************************************************************************/

#include "tasks.h"

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
	cub_spline_t* spline = newCubSpline(set);

	// generates data for plot
	plotInterp_5(lagEqn, spline);

	// calculate the interpolated values at xo
	interp_pt_t* lagEval = evaluateLagrangeEqn(lagEqn, xo);
	interp_pt_t* splineEval = evaluateCubSpline(spline, xo);
	printInterp_5(lagEval->fx, splineEval->fx);
	freeInterpPt(lagEval);
	freeInterpPt(splineEval);

	freeCubSpline(spline);
	freeLagrangeEqn(lagEqn);
	freeInterpSet(set);
}

void heateqn(const char* filename) {

	heat_sim_t* exFe = parseInput_6(filename);
	eulerExFe(exFe);
	printSim_6(exFe, FILENAME_EX_FE_6);
	freeHeatSim(exFe);

	heat_sim_t* exVe = parseInput_6(filename);
	eulerExFe(exVe);
	printSim_6(exVe, FILENAME_EX_VE_6);
	freeHeatSim(exVe);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * *  HEATEQN HELPER FUNCTIONS */

heat_sim_t* parseInput_6(const char* filename) {

	assert(filename!=NULL);

	// open file
	FILE* fp = fopen(filename, FILE_READONLY);
	assert(fp != NULL);

	double mu;
	int Nx, Nt;

	fscanf(fp, IN_HEADER_6);
	assert(fscanf(fp,"%lf,%d,%d\n", &mu, &Nx, &Nt) == NUM_PARAMS_6);
	fclose(fp);
	
	return newUnsolvedHeatSim(
		Nx, Nt, X_LO_6, X_HI_6, T_LO_6, T_HI_6, mu, myInitialCondition_6 );
}

double myInitialCondition_6(double x) {

	if(x>=IC_LO_6 && x<=IC_HI_6) {
		return IC_6(x);
	}
	return 0.0;
}

void printSim_6(heat_sim_t* sim, const char* filename) {
	
	assert(sim!=NULL);
	assert(filename!=NULL);

	sim_cell_t** m = sim->cells; // cell matrix

	FILE* fp = fopen(filename, FILE_OVERWRITE);
	assert(fp!=NULL);

	fprintf(fp, OUT_HEADER_6);

	int i=0;
	int j=100;
	for(i=0; i<=sim->Nt; i++) {
		fprintf(fp, "%.6f,%.6f\n", m[i][j].x, m[i][j].f);
    }
}


/* * * * * * * * * * * * * * * * * * * * * * * * * *  INTERP HELPER FUNCTIONS */

interp_set_t* parseInput_5(const char* filename) {

	assert(filename!=NULL);

	interp_set_t* newSet = newInterpSet(); 

	// open file
	FILE* fp = fopen(filename, FILE_READONLY);
	assert(fp != NULL);

	double tmp_x, tmp_fx;
	fscanf(fp, IN_HEADER_5);
	int read = 0;
	int i=0;
	while(1) {
		read = fscanf(fp,"%lf,%lf\n", &tmp_x, &tmp_fx);
		if(read != NUM_PARAMS_5) { break; }
		appendPtToSet(newSet, newInterpPt(tmp_x, tmp_fx) );
		i++;
	}

	fclose(fp);

	return newSet;
}

void printInterp_5(double lag, double spline) {

	// open file
	FILE* fp = fopen(FILENAME_5, FILE_OVERWRITE);
	assert(fp != NULL);

	fprintf(fp, HEADER_LAGRANGE_5);
	fprintf(fp, "%.6f\n", lag);
	fprintf(fp, HEADER_CUBIC_5);
	fprintf(fp, "%.6f\n", spline);

	fclose(fp);
}

void plotInterp_5(lagrange_eqn_t* lagEqn, cub_spline_t* spline) {

	// open file
	FILE* fp = fopen(PLOT_FILENAME_5, FILE_OVERWRITE);
	assert(fp != NULL);

	fprintf(fp, PLOT_HEADER_LAGRANGE_5);

	double x = PLOT_START_5;
	interp_pt_t* tmp_pt;
	while(x < PLOT_END_5) {

		fprintf(fp, "%f,", x);
		tmp_pt = evaluateLagrangeEqn(lagEqn, x);
		fprintf(fp, "%f,", tmp_pt->fx);
		freeInterpPt(tmp_pt);
		tmp_pt = evaluateCubSpline(spline, x);
		fprintf(fp, "%f\n", tmp_pt->fx);
		freeInterpPt(tmp_pt);
		x = x + PLOT_INTERVAL_5;
	}

	fclose(fp);
}


/* * * * * * * * * * * * * * * * * * * * * * * * LINALGBYSYS HELPER FUNCTIONS */

tridiag_t *parseInput_3(const char* filename) {
	
	assert(filename!=NULL);

	tridiag_t *m = newTridiag();

	// open file
	FILE* fp = fopen(filename, FILE_READONLY);
	assert(fp != NULL);

	double tmp_a, tmp_b, tmp_c, tmp_q;
	fscanf(fp, IN_HEADER_3);
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

	fscanf(fp, IN_1ST_HEADER_2);

	assert(NUM_PARAMS_2 == fscanf(
		fp,
		"%lf,%lf,%lf,%lf,%lf\n",
		&(parsed->M_a),
		&(parsed->t_a),
		&(parsed->b_l_a),
		&(parsed->b_u_a),
		&(parsed->g_a)
	) );
}

void parseInputMvals_2(FILE* fp, input_2_t* parsed) {
	
	assert(parsed!=NULL && fp!= NULL);

	fscanf(fp, IN_2ND_HEADER_2);

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

