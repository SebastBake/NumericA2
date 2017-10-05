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
#define ARGV_Q2 1
#define ARGV_Q4 2
#define ARGV_Q5 3
#define ARGV_xo 4
#define ARGV_Q6 5
#define MS_PER_SEC 1000.0

typedef struct timeval myTime_t;
myTime_t timer_start();
double timer_stop(myTime_t start);

int main(int argc, char *argv[]) {
	
	// Parse command line arguments
	char* q2_file_name = NULL;
	char* q4_file_name = NULL;
	char* q5_file_name = NULL;
	double xo = 0;
	char* q6_file_name = NULL;

	if(argc != NUM_ARGS) {
		printf(INPUT_INSTRUCTIONS);
		exit(EXIT_FAILURE);
	}

	q2_file_name = argv[ARGV_Q2];
	q4_file_name = argv[ARGV_Q4];
	q5_file_name = argv[ARGV_Q5];
	xo = atof(argv[ARGV_xo]);
	q6_file_name = argv[ARGV_Q6];
	//printf("q2 %s, q4 %s, q5 %s, xo %f, q6 %s\n", q2_file_name, q4_file_name, q5_file_name, xo, q6_file_name);

	/* Question 2 */
	myTime_t q2_time = timer_start();
	shockwave(q2_file_name);
	printf("q2: %.2lf milliseconds\n", timer_stop(q2_time));
	
	/* Question 4 */
	myTime_t q4_time = timer_start();
	linalgbsys(q4_file_name);
	printf("q4: %.2lf milliseconds\n", timer_stop(q4_time));
	
	/* Question 5 */
	myTime_t q5_time = timer_start();
	interp(q5_file_name, xo);
	printf("q5: %.2lf milliseconds\n", timer_stop(q5_time));
	
	/* Question 6 */
	myTime_t q6_time = timer_start();
	heateqn(q6_file_name);
	printf("q6: %.2lf milliseconds\n", timer_stop(q6_time));
    
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
