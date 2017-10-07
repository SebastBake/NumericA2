/***************************************************************************
 *
 *   File        : main.c
 *   Student Id  : 757931
 *   Name        : Sebastian Baker
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "tasks.h"

#define NUM_ARGS 6
#define INPUT_INSTRUCTIONS "Usage example: ./exec in_shock.csv in_linalsys.csv in_interp.csv 5 in_heateqn.csv"
#define ARGV_shockwave 1
#define ARGV_linalgbsys 2
#define ARGV_interp 3
#define ARGV_xo 4
#define ARGV_heateqn 5
#define MS_PER_SEC 1000.0

typedef struct timeval myTime_t;
myTime_t timer_start();
double timer_stop(myTime_t start);

int main(int argc, char *argv[]) {
	
	// Parse command line arguments
	char* shockwave_file_name = NULL;
	char* linalgbsys_file_name = NULL;
	char* interp_file_name = NULL;
	double xo = 0;
	char* heateqn_file_name = NULL;

	if(argc != NUM_ARGS) {
		printf(INPUT_INSTRUCTIONS);
		exit(EXIT_FAILURE);
	}

	shockwave_file_name = argv[ARGV_shockwave];
	linalgbsys_file_name = argv[ARGV_linalgbsys];
	interp_file_name = argv[ARGV_interp];
	xo = atof(argv[ARGV_xo]);
	heateqn_file_name = argv[ARGV_heateqn];
	//printf("shockwave %s, linalgbsys %s, interp %s, xo %f, heateqn %s\n", shockwave_file_name, linalgbsys_file_name, interp_file_name, xo, heateqn_file_name);

	/* Question 2 */
	myTime_t shockwave_time = timer_start();
	shockwave(shockwave_file_name);
	printf("shockwave: %.2lf milliseconds\n", timer_stop(shockwave_time));
	
	/* Question 4 */
	myTime_t linalgbsys_time = timer_start();
	linalgbsys(linalgbsys_file_name);
	printf("linalgbsys: %.2lf milliseconds\n", timer_stop(linalgbsys_time));
	
	/* Question 5 */
	myTime_t interp_time = timer_start();
	interp(interp_file_name, xo);
	printf("interp: %.2lf milliseconds\n", timer_stop(interp_time));
	
	/* Question 6 */
	myTime_t heateqn_time = timer_start();
	heateqn(heateqn_file_name);
	printf("heateqn: %.2lf milliseconds\n", timer_stop(heateqn_time));
    
	return EXIT_SUCCESS;
}

myTime_t timer_start() {
	myTime_t start;
	gettimeofday(&start, NULL);
	return start;
}

double timer_stop(myTime_t start) {
	myTime_t stop;
	gettimeofday(&stop, NULL);
	double elapsed_ms = (stop.tv_sec - start.tv_sec) * MS_PER_SEC;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / MS_PER_SEC;
	return elapsed_ms;
}
