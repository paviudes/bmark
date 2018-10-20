#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include "memory.h"
#include "notations.h"
#include "tiling.h"
#include "arrays.h"
#include "dual.h"
#include "save.h"
#include "rgens.h"
#include "rbim.h"
#include "simulation.h"
#include "mcresult.h"
#include "noise.h"
#include "correlations.h"
#include "homological.h"
#include "mcresult.h"
#include "noise.h"
#include "performance.h"

void FormDualErasure(struct tiling **pGs, struct simulation* pSim){
	// Compute the subgraph induced by the erasure pattern on the dual graph.
	FillIntArrayWithConstants((pSim->erasure)[2], pGs[1]->e, 0);
	int ei;
	for(ei = 0; ei < (pGs[0]->e); ei ++)
		if((pSim->erasure)[1][ei] == 1)
			if((pSim->edgeOrderings)[0][ei] != -1)
				(pSim->erasure)[2][(pSim->edgeOrderings)[0][ei]] = 1;
}


void GenerateErasure(struct tiling **pGs, struct simulation *pSim, struct noise *pn){
	// Generate an erasure pattern by following a particular error model
	int weight = 0;
	FillIntArrayWithConstants((pSim->erasure)[1], pGs[0]->e, 0);
	if(pn->model == 1){
		if (pn->dist == -1)
			BinomialBitString((pSim->erasure)[1], pGs[0]->e, pGs[0]->doE, (pn->params)[pn->current][0]);
		else{
			weight = Random("cumulative", pGs[0]->e, pn->impcumul);
			pn->currentbias = (pn->bias)[weight];
			RandBinWeight(pGs[0]->e, weight, (pSim->erasure)[1]);

			// printf("weight = %d, bias = %g\n", weight, pn->currentbias);
			// PrintIntArray((pSim->erasure)[1], pGs[0]->e);
		}
	}
	else if(pn->model == 2)
		BallErasure(pSim, pGs[0], (pn->params)[pn->current][0], pn->corrdist);
	else if(pn->model == 3)
		WalkErasure(pSim, pGs, (pn->params)[pn->current][0], pn->corrdist);
	else if(pn->model == 4)
		SpreadErasure(pSim, pGs, (pn->params)[pn->current][0], pn->corrdist);
	else if(pn->model == 5)
		RBIMErasure(pSim, pGs[0], pn->ising, (int) (pn->params)[pn->current][(int)(pn->params)[0][1] - 1]);
	else
		fprintf(stderr, "\033[91mUnable to comprehend the noise model (of type '%d'). No qubits will be erased.\033[0m\n", pn->model);
	// Compute the weight of the erasure and the induced erasure pattern on the dual graph.
	(pSim->erasure)[0][0] = Sum((pSim->erasure)[1], pGs[0]->e);
	(pSim->erasure)[0][1] = (pSim->erasure)[0][0];
	FormDualErasure(pGs, pSim);

	if (pn->dist == -1)
		pn->currentbias = 1;
}


void DecodingTest(struct tiling **pGs, struct simulation* pSim, struct mcresult *pmcr, int dist, double currentbias){
	// Perform one decoding test by generating an erasure and determining the type of logical operator covered by this erasure
	int ct, logType;
	/*
	// Increment the number of trials, only for the type of errors for which the numerics are performed
	for(ct = 0; ct < 3; ct ++)
		if((pSim->type == 0) || (pSim->type == ct))
			(pmcr->trials)[ct] ++;
	// Determine the type of logical operator supported by the erasure pattern
	logType = LogicalType(pGs, pSim);
	if(logType > 0){
		if(pSim->type == 0){
			(pmcr->fails)[0] ++;
			if (dist != -1)
				(pmcr->rates)[0] += currentbias/(double) pSim->stats;
		}
		if(logType == 1){
			(pmcr->fails)[1] ++;
			(pmcr->fails)[2] ++;
			if (dist != -1){
				(pmcr->rates)[1] += currentbias/(double) pSim->stats;
				(pmcr->rates)[2] += currentbias/(double) pSim->stats;
			}
		}
		else if(logType == 2){
			(pmcr->fails)[1] ++;
			if (dist != -1)
				(pmcr->rates)[1] += currentbias/(double) pSim->stats;
		}
		else if(logType == 3){
			(pmcr->fails)[2] ++;
			if (dist != -1)
				(pmcr->rates)[2] += currentbias/(double) pSim->stats;
		}
		else;
	}
	*/

	// Determine the type of logical operator supported by the erasure pattern
	logType = LogicalType(pGs, pSim);
	// Update the count of failures, failure rates and total trials accordingly.
	for(ct = 0; ct < 3; ct ++){
		if(((pSim->type) == 0) || ((pSim->type) == ct)){
			// Increment the number of trials, only for the type of errors for which the numerics are performed
			(pmcr->trials)[ct] ++;
			if (logType > 0){
				if ((ct == 0) || (logType == 1) || (logType == (ct + 1))){
					(pmcr->fails)[ct] ++;
					if (dist != -1)
						(pmcr->rates)[ct] += currentbias/(double) pSim->stats;
				}
			}
			if (dist == -1)
				(pmcr->rates)[ct] = (double) (pmcr->fails)[ct]/(double) (pmcr->trials)[ct];
		}
	}
	// fprintf(stderr, "(pmcr->fails)[0] = %ld, (pmcr->rates)[0] = %.4e\n", (pmcr->fails)[0], (pmcr->rates)[0]);
}



int IsBreakPoint(struct simulation *pSim, struct mcresult *pmcr){
	/*
	Decide if the Monte carlo rounds sampling over the erasures must be continued, based on the error bars on the current failure rate estimates.
	Stopping rules
		1. If the error bars (for X, Z failures) are below the maximum allowed size, stop numerics.
		2. Or else
			A. If both Z as well as X error bars are larger than the maximum, continue testing for both
			B. If the error bar corresponding to the estimation of the failure rates for Z-types errors is above the tolerated error, while those corresponding to X-type errors are within the tolerance, then continue sampling only Z-type errors.
			C. If the case is the opposite of the one described above, i.e, X-error bars are within tolerance while Z-error bars aren't, then sample only X-type errors.
	*/
	int ct, stop;
	for(ct = 0; ct < 3; ct ++){	
		if((pSim->type == 0) || (pSim->type == ct))
			(pmcr->bars)[pSim->type] = ErrorBar((pmcr->rates)[pSim->type], (pmcr->trials)[pSim->type]);
		if ((pmcr->bars)[ct] <= (pSim->accuracy) * (pmcr->rates)[ct])
			(pmcr->pass)[ct] = 1;
	}
	stop = 0;
	if(((pmcr->bars)[1] > (pSim->accuracy * (pmcr->rates)[1])) && ((pmcr->bars)[2] > (pSim->accuracy * (pmcr->rates)[2])))
		// Case 2A
		pSim->type = 0;
	else if((pmcr->bars)[1] > (pSim->accuracy * (pmcr->rates)[1]))
		// Case 2B
		pSim->type = 1;
	else if((pmcr->bars)[2] > (pSim->accuracy * (pmcr->rates)[2]))
		// Case 2C
		pSim->type = 2;
	else
		// Case 1
		stop = 1;

	// fprintf(stderr, "Currently checking for errors of type %d\n", pSim->type);
	// PrintResult(pmcr);
	return stop;
}


int MCStop(struct simulation* pSim, struct mcresult *pmcr, long running){
	// Determine if the estimates are good enough to test if any of the Montecarlo stopping rules in IsBreakPoint(...) apply.
	int isCheckPoint = 0;
	// Part 1.
	if(((pmcr->fails)[1] < (pSim->minfails)) && ((pmcr->fails)[2] < (pSim->minfails))){
		pSim->type = 0;
		isCheckPoint = 0;
	}
	else{
		// Part 2.
		if(((pmcr->trials)[pSim->type] > (int)((pSim->stats)/4)) && ((pmcr->fails)[pSim->type] == 0))
			isCheckPoint = 1;
		// Part 3
		else if((pmcr->fails)[1] < pSim->minfails){
			pSim->type = 1;
			isCheckPoint = 0;
		}
		else if((pmcr->fails)[2] < pSim->minfails){
			pSim->type = 2;
			isCheckPoint = 0;
		}
		else
			// Part 4.
			if(((pmcr->trials)[pSim->type] % ((long)((pSim->stats)/10))) == 0)
				isCheckPoint = 1;
	}
	if(running == (pSim->stats))
		isCheckPoint = 1;
	if (isCheckPoint == 1)
		return IsBreakPoint(pSim, pmcr);
	return 0;
}

void MCEstimate(struct tiling **pGs, struct simulation *pSim, struct noise *pn, struct mcresult *pmcr){
	// Run decoding trials to estimate the decoding failure probability
	long si;
	pSim->type = 0;
	double *clusterinfo = malloc(sizeof(double) * 3);
	for(si = 0; si < pSim->stats; si ++){
		pSim->totaltrials ++;
		GenerateErasure(pGs, pSim, pn);
		if ((pn->model == 2) || (pn->model == 3) || (pn->model == 4) || (pn->model == 5)){	
			(pn->properties)[1] += (double) ((pSim->erasure)[0][0]/(double) (pGs[0]->e));
			Cluster(pGs, pSim, clusterinfo);
			(pn->properties)[2] += clusterinfo[0];
			(pn->properties)[3] += clusterinfo[1];
			(pn->properties)[4] += clusterinfo[2];
			if (pn->model == 5){
				PhysicalProperties(pn->ising);
				(pn->properties)[5] += pn->ising->magnetization;
				(pn->properties)[6] += pn->ising->internalEnergy;
				(pn->properties)[7] += pn->ising->twopoint;
			}
		}
		DecodingTest(pGs, pSim, pmcr, pn->dist, pn->currentbias);

		// Get the running average of failure rates. This is just to see how the estimate average converges to the true average.
		UpdateRunningAverage(si, pmcr);

		if (pn->dist == -1)
			if(MCStop(pSim, pmcr, (si + 1)) == 1)
				break;
	}
	// Correction for the biased sampling over only errors whose weight was greater than distance.
	if (pn->model == 1){
		if (pn->dist != -1){
			int ct;
			for (ct = 0; ct < 3; ct ++)
				(pmcr->rates)[ct] = (1 - pn->trivial) * (pmcr->rates)[ct];
		}
	}
	free(clusterinfo);
	// Average the physical properties
	int i;
	for (i = 1; i <= (int) (pn->properties)[0]; i ++)
		(pn->properties)[i] = (pn->properties)[i]/(double) (si + 1);
	// fprintf(stderr, "\033[92m*********\033[0m\n");
	PrintNoise(pn);
	fprintf(stderr, "\033[92m......\033[0m\n");
	PrintResult(pmcr);
	fprintf(stderr, "\033[92m*********\033[0m\n");
}

void WriteResult(struct simulation *pSim, struct noise *pn, struct mcresult *pmcr, struct tiling *pG){
	// Write results of simulation with a noise rate into a file
	int pi, ct, nqubits = pG->e - Sum(pG->doE, pG->e);
	struct representation *repr = malloc(sizeof(struct representation));
	FixRepresentation(repr, pn->model, 0);
	char *datafile = malloc(sizeof(char) * 100);
	sprintf(datafile, "results/perf_%s_%s_%s.txt", pG->type, repr->mname, pSim->timestamp);
	char *selectedDataFile = malloc(sizeof(char) * 100);
	sprintf(selectedDataFile, "results/selected_perf_%s_%s_%s.txt", pG->type, repr->mname, pSim->timestamp);
	if (pn->current == 1){
		DeleteFile(datafile);
		DeleteFile(selectedDataFile);
	}
	FILE *res = fopen(datafile, "a");
	FILE *selectRes = fopen(selectedDataFile, "a");
	fprintf(res, "%lf", (pn->params)[pn->current][0]);
	fprintf(selectRes, "%lf", (pn->params)[pn->current][0]);
	for (pi = 1; pi < (pn->params)[0][1]; pi ++){
		fprintf(res, " %lf", (pn->params)[pn->current][pi]);
		fprintf(selectRes, " %lf", (pn->params)[pn->current][pi]);
	}
	
	for (pi = 1; pi <= (int) (pn->properties)[0]; pi ++){
		fprintf(res, " %.4e", (pn->properties)[pi]);
		fprintf(selectRes, " %.4e", (pn->properties)[pi]);
	}
	for(ct = 0; ct < 3; ct ++){
		fprintf(res, " %.4e %ld %ld %.4e", (pmcr->rates)[ct], (pmcr->fails)[ct], (pmcr->trials)[ct], (pmcr->bars)[ct]);
		if ((pmcr->pass)[ct])
			fprintf(selectRes, " %.4e %ld %ld %.4e", (pmcr->rates)[ct], (pmcr->fails)[ct], (pmcr->trials)[ct], (pmcr->bars)[ct]);
		else
			fprintf(selectRes, " 0 0 0 0 0");
	}
	fprintf(res, " %d %lf\n", nqubits, pmcr->runtime);
	fprintf(selectRes, " %d %lf\n", nqubits, pmcr->runtime);
	fclose(res);
	fclose(selectRes);
	free(datafile);
	free(selectedDataFile);

	// Write the running averages data
	char *runavgFile = malloc(sizeof(char) * 100);
	sprintf(runavgFile, "results/runavg_%s_%s_%s.txt", pG->type, repr->mname, pSim->timestamp);
	FILE *runavgfp = fopen(runavgFile, "a");
	fprintf(runavgfp, "%s = %g", (repr->varnames)[0], (pn->params)[pn->current][0]);
	for (pi = 1; pi < (pn->params)[0][1]; pi ++)
		fprintf(runavgfp, ", %s = %g", (repr->varnames)[pi], (pn->params)[pn->current][pi]);
	fprintf(runavgfp, ".\n");
	fprintf(runavgfp, "#N f_XZ f_X f_Z\n");
	for (pi = 1; pi <= (int) (pmcr->running)[0][0]; pi ++)
		fprintf(runavgfp, "%ld %.4e %.4e %.4e\n", (long) (pmcr->running)[pi][0], (pmcr->running)[pi][1], (pmcr->running)[pi][2], (pmcr->running)[pi][3]);
	fprintf(runavgfp, "\n\n");
	fclose(runavgfp);
	free(runavgFile);

	pSim->runtime += pmcr->runtime;
	if (pn->current == (pn->params)[0][0]){
		char *runfile = malloc(sizeof(char) * 100);
		sprintf(runfile, "results/runtimes_%s_%s_%s.txt", pG->type, repr->mname, pSim->timestamp);
		FILE *run = fopen(runfile, "w");
		fprintf(run, "%d %ld %g\n", nqubits, pSim->totaltrials, pSim->runtime);
		fclose(run);
		free(runfile);
		// Save the time stamp
		char *tsfname = malloc(sizeof(char) * 100);
		sprintf(tsfname, "results/timestamp_%s_%s.txt", pG->type, repr->mname);
		FILE *tsfp = fopen(tsfname, "w");
		fprintf(tsfp, "%s\n", pSim->timestamp);
		fclose(tsfp);
		free(tsfname);
	}
	FreeRep(repr);
}

void DecodingFailureRates(struct tiling **pGs, struct simulation *pSim, struct noise *pn, struct mcresult *pmcr){
	// Run montecarlo simulations to estimate the failure probabilty of the decoder for a specific noise range
	int pi, skip = 0;
	clock_t startTime, endTime;
	for (pi = 1; pi <= (pn->params)[0][0]; pi ++){
		skip = UpdateNoise(pn, pmcr, pGs[0]->type, pGs[0]->e);
		startTime = clock();
		ResetResults(pmcr);
		if (skip == 0)
			MCEstimate(pGs, pSim, pn, pmcr);
		else{
			if (skip == 1)
				FillDoubleArrayWithConstants(pmcr->rates, 3, 0);
			if (skip == 2)
				FillDoubleArrayWithConstants(pmcr->rates, 3, 1);
		}
		endTime = clock();
		pmcr->runtime = (double)(endTime - startTime)/CLOCKS_PER_SEC;
		WriteResult(pSim, pn, pmcr, pGs[0]);
	}
}

/* There are several forms of the Performance function.
	1. Performance: The number of Monte Carlo runs is an integer.
	2. PerformanceRunning: The number of Monte Carlo runs is speacified as an array. In this case, the running average is computed at at every array element.
*/
char* Performance(struct tiling G, int dist, int model, double **noiseParams, long stats){
	// Analyze the performance of a code (defined with a tiling) with montecarlo simulation of specific number of decoding trails
	int ti, nrates;
	struct tiling **pGs = malloc(sizeof(struct tiling*) * 3);
	for (ti = 0; ti < 3; ti ++)
		pGs[ti] = malloc(sizeof(struct tiling));
	pGs[0][0] = G;
	pGs[1][0] = DualTiling(pGs[0][0]);
	pGs[2][0] = DualTiling(pGs[1][0]);

	// printf("dist = %d\n", dist);
	
	struct mcresult *pmcr = malloc(sizeof(struct mcresult));
	long *breakpoints = malloc(sizeof(long) * 2);
	breakpoints[0] = 1;
	breakpoints[1] = stats;
	InitializeResult(pmcr, breakpoints);
	struct noise *pn = malloc(sizeof(struct noise));
	InitializeNoise(pn, pGs[0], dist, model, noiseParams);
	nrates = (int) ((pn->params)[0][0]);
	struct simulation *pSim = malloc(sizeof(struct simulation));
	InitializeSimulation(pSim, pGs, nrates, stats);

	DecodingFailureRates(pGs, pSim, pn, pmcr);
	char *summary = malloc(sizeof(char) * 500);
	sprintf(summary, "%s", pSim->timestamp);
	
	// Free memory
	for (ti = 1; ti < 3; ti ++)
		FreeTiling(pGs[ti]);
	free(pGs);
	FreeSimulation(pSim, nrates);
	FreeNoise(pn);
	FreeResult(pmcr);
	free(breakpoints);
	return summary;
}

char* PerformanceRunning(struct tiling G, int dist, int model, double **noiseParams, long *breakpoints){
	// Analyze the performance of a code (defined with a tiling) with montecarlo simulation of specific number of decoding trails
	int ti, nrates;
	struct tiling **pGs = malloc(sizeof(struct tiling*) * 3);
	for (ti = 0; ti < 3; ti ++)
		pGs[ti] = malloc(sizeof(struct tiling));
	pGs[0][0] = G;
	pGs[1][0] = DualTiling(pGs[0][0]);
	pGs[2][0] = DualTiling(pGs[1][0]);
	struct mcresult *pmcr = malloc(sizeof(struct mcresult));
	InitializeResult(pmcr, breakpoints);
	struct noise *pn = malloc(sizeof(struct noise));
	InitializeNoise(pn, pGs[0], dist, model, noiseParams);
	nrates = (int) ((pn->params)[0][0]);
	struct simulation *pSim = malloc(sizeof(struct simulation));
	InitializeSimulation(pSim, pGs, nrates, breakpoints[breakpoints[0]]);

	DecodingFailureRates(pGs, pSim, pn, pmcr);
	char *summary = malloc(sizeof(char) * 500);
	sprintf(summary, "%s", pSim->timestamp);
	
	// Free memory
	for (ti = 1; ti < 3; ti ++)
		FreeTiling(pGs[ti]);
	free(pGs);
	FreeSimulation(pSim, nrates);
	FreeNoise(pn);
	FreeResult(pmcr);
	return summary;
}
