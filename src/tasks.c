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
	cub_spline_t* spline = newCubSpline(set);
	lagrange_eqn_t* lagEqn = newLagrangeEqn(set);

	// generates data for plots
	//plotInterpAnyOrder_5(lagEqn, spline);
	//plotInterpQuadratic_5(set, spline);

	// calculate/output the interpolated values at xo
	interp_pt_t* lagEval = evaluateQuadLagrangeEqn(set, xo);
	interp_pt_t* splineEval = evaluateCubSpline(spline, xo);
	printInterp_5(lagEval->fx, splineEval->fx);
	freeInterpPt(lagEval);
	freeInterpPt(splineEval);

	freeLagrangeEqn(lagEqn);
	freeCubSpline(spline);
	freeInterpSet(set);
}

void heateqn(const char* filename) {

	heat_sim_t* exFe = parseInput_6(filename);
	eulerExFe(exFe);
	printSim_6(exFe, FILENAME_EX_FE_6);
	freeHeatSim(exFe);

	heat_sim_t* exVe = parseInput_6(filename);
	eulerExVe(exVe);
	printSim_6(exVe, FILENAME_EX_VE_6);
	freeHeatSim(exVe);

	heat_sim_t* imFe = parseInput_6(filename);
	eulerImFe(imFe);
	printSim_6(imFe, FILENAME_IM_FE_6);
	freeHeatSim(imFe);

	//runDxTest_6();
	//runDtTest_6();
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * HEATEQN HELPER FUNCTIONS */

// parses heat simulation file and returns a new unsolved heat simulation
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

// initial condition for f(x,t=0)
double myInitialCondition_6(double x) {

	if(x>=IC_LO_6 && x<=IC_HI_6) {
		return IC_6(x);
	}
	return 0.0;
}

// print x's at 100th timestep
void printSim_6(heat_sim_t* sim, const char* filename) {
	
	assert(sim!=NULL);
	assert(filename!=NULL);

	sim_cell_t** m = sim->cells; // cell matrix

	FILE* fp = fopen(filename, FILE_OVERWRITE);
	assert(fp!=NULL);

	fprintf(fp, OUT_HEADER_6);

	int i=0;
	int j=100;
	for(i=0; i<=sim->Nx; i++) {
		fprintf(fp, "%.4f,%.4f\n", m[i][j].x, m[i][j].f);
    }
}

heat_sim_t* generateRefSoln_6(void (*solve)(heat_sim_t*)) {

	double dx = (X_HI_6 - X_LO_6)/REF_Nx_6;
	double dt = Dn2Dt(REF_Dn_6, dx, REF_MU_6); // diff num to dt
	int Nt = round((T_HI_6 - T_LO_6)/dt);

	heat_sim_t* sim = newUnsolvedHeatSim(
		REF_Nx_6, Nt, X_LO_6, X_HI_6, T_LO_6, T_HI_6, REF_MU_6, myInitialCondition_6 );

	solve(sim);
	//printSim_6(sim, "reference.csv");
	return sim;
}

double calculateAbsAvgError_6(heat_sim_t* ref, heat_sim_t* sim) {

	double eSum = 0;
	double res = 100;
	double Rcx = ref->Nx / res;
	double Scx = sim->Nx / res;
	double Rct = ref->Nt / res; 
	double Sct = sim->Nt / res;

	int i=0, j=0;
	for(i=0;i<=res;i++) {
		for(j=0;j<=res;j++) {
			
			eSum += fabs(
				ref->cells[(int)round(Rcx*i)][(int)round(Rct*j)].f - 
				sim->cells[(int)round(Scx*i)][(int)round(Sct*j)].f
			);
		}
	}

	return eSum/((res+1)*(res+1));
}

void runDxTest_6() {

	FILE* fp = fopen(DX_TEST_FILE, FILE_OVERWRITE);

	fprintf(fp, "dx,Dn,exfe,exve,imfe\n");

	//printf("start dx tests\n");
	//printf("generate reference solutions tests\n");
	heat_sim_t* refExFe = generateRefSoln_6(eulerExFe);
	heat_sim_t* refExVe = generateRefSoln_6(eulerExVe);
	heat_sim_t* refImFe = generateRefSoln_6(eulerImFe);
	heat_sim_t* tmpSim;
	double tmpErr = 0.0;
	double dx = 0.0;
	double dt = (T_HI_6 - T_LO_6)/TEST_NX_NT;
	int nxStart = TEST_NX_START;
	int nxEnd = TEST_NX_END;

	//printf("begin error comparisons \n");

	int nx=0;
	for(nx=nxStart; nx<=nxEnd; nx++) {

		dx = (X_HI_6 - X_LO_6)/nx;
		fprintf(fp,"%f,%f,", dx, Dt2Dn(dt, dx, REF_MU_6));

		// explicit fe
		tmpSim = newUnsolvedHeatSim(
			nx, TEST_NX_NT, X_LO_6, X_HI_6, T_LO_6, T_HI_6, REF_MU_6, myInitialCondition_6 );

		eulerExFe(tmpSim);
		tmpErr = calculateAbsAvgError_6(refExFe, tmpSim);
		fprintf(fp,"%f,", tmpErr);
		freeHeatSim(tmpSim);

		// explicit ve
		tmpSim = newUnsolvedHeatSim(
			nx, TEST_NX_NT, X_LO_6, X_HI_6, T_LO_6, T_HI_6, REF_MU_6, myInitialCondition_6 );

		eulerExVe(tmpSim);
		tmpErr = calculateAbsAvgError_6(refExVe, tmpSim);
		fprintf(fp,"%f,", tmpErr);
		freeHeatSim(tmpSim);

		// implicit fe
		tmpSim = newUnsolvedHeatSim(
			nx, TEST_NX_NT, X_LO_6, X_HI_6, T_LO_6, T_HI_6, REF_MU_6, myInitialCondition_6 );

		eulerImFe(tmpSim);
		tmpErr = calculateAbsAvgError_6(refImFe, tmpSim);
		fprintf(fp,"%f\n", tmpErr);
		freeHeatSim(tmpSim);
	}

	freeHeatSim(refExFe);
	freeHeatSim(refExVe);
	freeHeatSim(refImFe);
	fclose(fp);
}

void runDtTest_6() {

	FILE* fp = fopen(DT_TEST_FILE, FILE_OVERWRITE);
	fprintf(fp, "dt,Dn,exfe,exve,imfe\n");

	//printf("start dt tests\n");
	//printf("generate reference solutions tests\n");
	heat_sim_t* refExFe = generateRefSoln_6(eulerExFe);
	heat_sim_t* refExVe = generateRefSoln_6(eulerExVe);
	heat_sim_t* refImFe = generateRefSoln_6(eulerImFe);
	heat_sim_t* tmpSim;
	double tmpErr = 0.0;
	double dx = (X_HI_6 - X_LO_6)/TEST_NT_NX;
	double dt = 0.0;
	int ntStart = TEST_NT_START;
	int ntEnd = TEST_NT_END;

	//printf("begin error comparisons \n");

	int nt=0;
	for(nt=ntStart; nt<=ntEnd; nt+=NT_STEP) {

		dt = (T_HI_6 - T_LO_6)/nt;
		fprintf(fp,"%f,%f,", dt, Dt2Dn(dt, dx, REF_MU_6));

		// explicit fe
		tmpSim = newUnsolvedHeatSim(TEST_NT_NX, nt, X_LO_6, X_HI_6, T_LO_6, T_HI_6, REF_MU_6, myInitialCondition_6 );
		eulerExFe(tmpSim);
		tmpErr = calculateAbsAvgError_6(refExFe, tmpSim);
		fprintf(fp,"%f,", tmpErr);
		freeHeatSim(tmpSim);

		// explicit ve
		tmpSim = newUnsolvedHeatSim(TEST_NT_NX, nt, X_LO_6, X_HI_6, T_LO_6, T_HI_6, REF_MU_6, myInitialCondition_6 );
		eulerExVe(tmpSim);
		tmpErr = calculateAbsAvgError_6(refExVe, tmpSim);
		fprintf(fp,"%f,", tmpErr);
		freeHeatSim(tmpSim);

		// implicit fe
		tmpSim = newUnsolvedHeatSim(TEST_NT_NX, nt, X_LO_6, X_HI_6, T_LO_6, T_HI_6, REF_MU_6, myInitialCondition_6 );
		eulerImFe(tmpSim);
		tmpErr = calculateAbsAvgError_6(refImFe, tmpSim);
		fprintf(fp,"%f\n", tmpErr);
		freeHeatSim(tmpSim);
	}

	freeHeatSim(refExFe);
	freeHeatSim(refExVe);
	freeHeatSim(refImFe);
	fclose(fp);
}




/* * * * * * * * * * * * * * * * * * * * * * * * * *  INTERP HELPER FUNCTIONS */

// parse input for interpolation task, return the set
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
		appendPtToSet(newSet, newInterpPt(tmp_x, tmp_fx));
		i++;
	}

	fclose(fp);

	return newSet;
}

// generate a quadratic lagrange eqn from a set, evaluate it at xo
interp_pt_t* evaluateQuadLagrangeEqn(interp_set_t* set, double xo) {

	// require at least 3 points to make a quadratic
	assert(set->N >= QUADRATIC_NUM_PTS);

	// insert the 3 closest points to xo into a new set to make the quadratic
	// note that my implementation does not require points to be in any order
	interp_set_t* newSet = newInterpSet();
	interp_pt_t* p1 = set->pts[0];
	double p1_diff = DBL_MAX;
	interp_pt_t* p2 = set->pts[0];
	double p2_diff = DBL_MAX;
	interp_pt_t* p3 = set->pts[0];
	double p3_diff = DBL_MAX;

	// find closest point <x
	int i = 0;
	for(i=0; i < set->N ; i++) {
		if(
			fabs(set->pts[i]->x-xo) <= p1_diff &&
			set->pts[i]->x <= xo
		) {
			p1 = set->pts[i];
			p1_diff = fabs(p1->x-xo);
		}
	}
	appendPtToSet(newSet, newInterpPt(p1->x, p1->fx));

	// find closest point >x
	for(i=0; i < set->N ; i++) {
		if(
			fabs(set->pts[i]->x-xo) <= p2_diff &&
			set->pts[i]->x > xo
		) {
			p2 = set->pts[i];
			p2_diff = fabs(p2->x-xo);
		}
	}
	appendPtToSet(newSet, newInterpPt(p2->x, p2->fx));

	// find last closest point which isn't p1 or p2
	for(i=0; i < set->N ; i++) {
		if(
			p3_diff >= fabs(set->pts[i]->x-xo) && 
			set->pts[i] != p1 &&
			set->pts[i] != p2
		) {
			p3 = set->pts[i];
			p3_diff = fabs(p3->x-xo);
		}
	}
	appendPtToSet(newSet, newInterpPt(p3->x, p3->fx));
	//printf("%f,%f,%f\n", p1->x, p2->x, p3->x);

	// generate & evaluate the equation
	lagrange_eqn_t* lagEqn = newLagrangeEqn(newSet);
	interp_pt_t* lagEval = evaluateLagrangeEqn(lagEqn, xo);
	freeLagrangeEqn(lagEqn);
	freeInterpSet(newSet);

	return lagEval;
}

// simple function to print the output for the interpolation output
void printInterp_5(double lag, double spline) {

	// open file
	FILE* fp = fopen(FILENAME_5, FILE_OVERWRITE);
	assert(fp != NULL);

	fprintf(fp, HEADER_LAGRANGE_5);
	fprintf(fp, "%.4f\n", lag);
	fprintf(fp, HEADER_CUBIC_5);
	fprintf(fp, "%.4f\n", spline);

	fclose(fp);
}

// generates data comparing lagrange (any order), spline interpolations
void plotInterpAnyOrder_5(lagrange_eqn_t* lagEqn, cub_spline_t* spline) {

	// open file
	FILE* fp = fopen(PLOT_FILENAME_ANY_ORDER_5, FILE_OVERWRITE);
	assert(fp != NULL);

	fprintf(fp, PLOT_HEADER_LAGRANGE_5);

	double x = PLOT_START_5;
	interp_pt_t* tmp_pt;
	while(x < PLOT_END_5) {

		fprintf(fp, "%.4f,", x);
		tmp_pt = evaluateLagrangeEqn(lagEqn, x);
		fprintf(fp, "%.4f,", tmp_pt->fx);
		freeInterpPt(tmp_pt);
		tmp_pt = evaluateCubSpline(spline, x);
		fprintf(fp, "%.4f\n", tmp_pt->fx);
		freeInterpPt(tmp_pt);
		x = x + PLOT_INTERVAL_5;
	}

	fclose(fp);
}

// generates data comparing lagrange (quadratic), spline interpolations
void plotInterpQuadratic_5(interp_set_t* set, cub_spline_t* spline) {

	// open file
	FILE* fp = fopen(PLOT_FILENAME_QUADRATIC_5, FILE_OVERWRITE);
	assert(fp != NULL);

	fprintf(fp, PLOT_HEADER_LAGRANGE_5);

	double x = PLOT_START_5;
	interp_pt_t* tmp_pt;
	while(x < PLOT_END_5) {

		fprintf(fp, "%.4f,", x);
		tmp_pt = evaluateQuadLagrangeEqn(set, x);
		fprintf(fp, "%.4f,", tmp_pt->fx);
		freeInterpPt(tmp_pt);
		tmp_pt = evaluateCubSpline(spline, x);
		fprintf(fp, "%.4f\n", tmp_pt->fx);
		freeInterpPt(tmp_pt);
		x = x + PLOT_INTERVAL_5;
	}

	fclose(fp);
}


/* * * * * * * * * * * * * * * * * * * * * * * * LINALGBYSYS HELPER FUNCTIONS */

// parses input for task 3, returns a tridiagonal matrix equation
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

// prints the solved tridiagonal linear system
void printTridiag_3(tridiag_t *m) {

	assert(m!=NULL);

	FILE* fp = fopen(FILENAME_3, FILE_OVERWRITE);
	assert(fp != NULL);

	fprintf(fp, FILE_HEADER_3);
	int i=1;
	for(i=1; i <= m->N; i++) {
		tridiag_row_t *r = getTridiagRow(m, i);
		fprintf(fp,"%.4f\n",r->x);
	}

	fclose(fp);
}


/* * * * * * * * * * * * * * * * * * * * * * * * * SHOCKWAVE HELPER FUNCTIONS */

// parses input from task 2, puts it into a nice struct
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

// parses the first part of the shockwave input file
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

// parses the 2nd part of the shockwave input file
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

// free's the input struct for shockwave qn
void freeInput_2(input_2_t *parsed) {
	assert(parsed!=NULL);
	free(parsed->M_c);
	free(parsed);
}

// returns the result of the shockwave function for a beta value
double f_2(double beta, double* params) {
	assert(params!=NULL);
	return F_2( beta, params[F_M_2], params[F_T_2], params[F_G_2]);
}

// shockwave part a
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
	//printf("%f,%f\n", RAD2DEG(b_l), RAD2DEG(b_u));
}

// prints the shockwave roots for a range of thetas until failure
void shockwavePrintThetaRange_2bc(FILE* fp, double* params) {

	int b_u_itr=0, b_l_itr=0;
	double b_l = DEG2RAD(B_l_GUESS_2), b_u = DEG2RAD(B_u_GUESS_2);

	while(1) {
		b_l_itr = newtonRaphson(&b_l, params, f_2);
		b_u_itr = newtonRaphson(&b_u, params, f_2);
		if( b_l_itr == ROOTFIND_FAIL || b_u_itr == ROOTFIND_FAIL) { break; }
		
		fprintf(fp, "%.4f,%.0f,%.4f,%.4f\n",
			params[F_M_2],
			RAD2DEG(params[F_T_2]),
			RAD2DEG(b_l),
			RAD2DEG(b_u)
		);

		params[F_T_2] += RADIAN_INCREMENT_2B;
	}
}

// shockwave part b
void shockwave_2b(const input_2_t *parsed) {

	assert(parsed!=NULL);

	double params[] = { parsed->M_a, RADIAN_START_2B, parsed->g_a };

	FILE* fp = fopen(FILENAME_2B, FILE_OVERWRITE);
	assert(fp != NULL);

	fprintf(fp,FILE_HEADER_2BC);
	shockwavePrintThetaRange_2bc( fp, params );

	fclose(fp);
}

// shockwave part c
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
