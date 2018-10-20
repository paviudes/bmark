#include <stdio.h>
#include <stdlib.h>
#include "arrays.h"
#include "mcresult.h"

void InitializeResult(struct mcresult *pmcr, long *breakpoints){
	pmcr->runtime = 0;
	pmcr->fails = malloc(sizeof(long) * 3);
	pmcr->trials = malloc(sizeof(long) * 3);
	pmcr->rates = malloc(sizeof(double) * 3);
	pmcr->bars = malloc(sizeof(double) * 3);
	pmcr->pass = malloc(sizeof(int) * 3);
	pmcr->running = malloc(sizeof(double *) * (breakpoints[0] + 1));
	(pmcr->running)[0] = malloc(sizeof(double) * 2);
	pmcr->running[0][0] = (double) breakpoints[0];
	pmcr->running[0][1] = 1;
	int i, j;
	for (i = 1; i <= (int) breakpoints[0]; i ++){
		(pmcr->running)[i] = malloc(sizeof(double) * 4);
		(pmcr->running)[i][0] = (double) breakpoints[i];
		for (j = 1; j <= 3; j ++)
			(pmcr->running)[i][j] = 0;
	}
}

void ResetResults(struct mcresult *pmcr){
	pmcr->runtime = 0;
	FillLongArrayWithConstants(pmcr->trials, 3, 0);
	FillDoubleArrayWithConstants(pmcr->rates, 3, 0);
	FillDoubleArrayWithConstants(pmcr->bars, 3, 0);
	FillLongArrayWithConstants(pmcr->fails, 3, 0);
	FillIntArrayWithConstants(pmcr->pass, 3, 0);
	int i;
	(pmcr->running)[0][1] = 1;
	for (i = 1; i <= (int) (pmcr->running)[0][0]; i ++)
		FillDoubleArrayWithConstants((pmcr->running)[i] + 1, 3, 0);
}

void FreeResult(struct mcresult *pmcr){
	free(pmcr->fails);
	free(pmcr->trials);
	free(pmcr->rates);
	free(pmcr->bars);
	int i;
	for (i = 0; i <= (int) (pmcr->running)[0][0]; i ++)
		free((pmcr->running)[i]);
	free(pmcr->running);
}

void PrintResult(struct mcresult *pmcr){
	fprintf(stderr, " 	X,Z		X		Z\n");
	fprintf(stderr, "Fails 	%ld 		%ld	 	%ld\n", (pmcr->fails)[0], (pmcr->fails)[1], (pmcr->fails)[2]);
	fprintf(stderr, "Trials  %ld 		%ld		%ld\n", (pmcr->trials)[0], (pmcr->trials)[1], (pmcr->trials)[2]);
	fprintf(stderr, "Rates 	%.3e 	%.3e	%.3e\n", (pmcr->rates)[0], (pmcr->rates)[1], (pmcr->rates)[2]);
	fprintf(stderr, "Bars 	%2.e		%.2e	%.2e\n", (pmcr->bars)[0], (pmcr->bars)[1], (pmcr->bars)[2]);
	fprintf(stderr, "Pass 	%d 		%d 		%d\n", (pmcr->pass)[0], (pmcr->pass)[1], (pmcr->pass)[2]);
}


void UpdateRunningAverage(long stat, struct mcresult *pmcr){
	// Compute the running estimated average at regular intervals in the Monte Carlo simulation.
	// The intervals (called breakpoints) are given in the first column of pmcr->running.
	// The next three columns contain the X & Z, X and Z failure rates.
	// The index of the next breakpoint for registering the running averge is given in pmcr->nextbp
	if ((stat + 1) == (long) (pmcr->running)[(int) ((pmcr->running)[0][1])][0]){
		// fprintf(stderr, "Checkpoint: %d\nrates: %.4e %.4e %.4e\n", (int) ((pmcr->running)[0][1]), (pmcr->rates)[0], (pmcr->rates)[1], (pmcr->rates)[2]);
		int i;
		for (i = 0; i < 3; i ++)	
			(pmcr->running)[(int) ((pmcr->running)[0][1])][i + 1] = (pmcr->rates)[i];
		(pmcr->running)[0][1] ++;
	}
}