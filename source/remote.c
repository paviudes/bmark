#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "tiling.h"
#include "memory.h"
#include "simulation.h"
#include "arrays.h"
#include "save.h"
#include "load.h"
#include "mcresult.h"
#include "noise.h"
#include "homological.h"
#include "correlations.h"
#include "performance.h"
#include "notations.h"
#include "plot.h"
#include "submission.h"
#include "remote.h"

int IndexOf(char *str, char **set, int size){
	// Find the index of a string in an array of strings. If the string doesn't exist, return -1.
	// if there are multiple occurences, we return the index of the last occurrence
	int i, index = -1;
	char *word1 = malloc(sizeof(char) * 10), *word2 = malloc(sizeof(char) * 10);
	for (i = 0; i < size; i ++){
		sprintf(word1, "%s", str);
		sprintf(word2, "%s", set[i]);
		if (strcmp(word1, word2) == 0)
			index = i;
	}
	free(word1);
	free(word2);
	return index;
}


void Remote(char *inputfile){
	// Analyze performance for a family of codes
	const int LINE_WIDTH = 100;
	FILE *infp = fopen(inputfile, "r");
	char line[LINE_WIDTH], *mode = malloc(sizeof(char) * 100), *option = malloc(sizeof(char) * 100), *queue = malloc(sizeof(char) * 100), *pname = malloc(sizeof(char) * 100), *jobname = malloc(sizeof(char) * 100), *tstamp = malloc(sizeof(char) * 100);
	// default values
	sprintf(mode, "plot");
	sprintf(jobname, "noname");
	sprintf(queue, "qwork");
	sprintf(tstamp, "empty_time_stamp");

	int i, sofar, read, family, ncodes, **lengths = NULL, model, wall = 10;
	long stats;
	double **params = NULL;
	struct representation *repr = NULL;
	struct axes *ax = NULL;

	while (feof(infp) == 0){
		sofar = 0;
		// ignore commented lines
		fgets(line, LINE_WIDTH, infp);
		if (line[0] == '#') continue;
		sscanf(line, "%s", option);
		if (strncmp(option, "mode", 4) == 0)
			sscanf(line, "%s %s", option, mode);
		if (strncmp(option, "family", 6) == 0)
			sscanf(line, "%s %d", option, &family);
		if (strncmp(option, "codes", 5) == 0){
			sscanf(line, "%s %d%n", option, &ncodes, &read);
			lengths = malloc(sizeof(int *) * (ncodes + 1));
			lengths[0] = malloc(sizeof(int));
			lengths[0][0] = ncodes;
			for (i = 0; i < ncodes; i ++){
				sofar += read;
				lengths[i + 1] = malloc(sizeof(int) * 3);
				sscanf(line + sofar, "%d%n", &lengths[i + 1][0], &read);
				lengths[i + 1][2] = -1;
			}
		}
		if (strncmp(option, "distances", 9) == 0){
			sscanf(line, "%s %d%n", option, &ncodes, &read);
			for (i = 0; i < ncodes; i ++){
				sofar += read;
				sscanf(line + sofar, "%d%n", &lengths[i + 1][2], &read);
			}
		}
		if (strncmp(line, "noise", 5) == 0){
			sscanf(line, "%s %d", option, &model);
			repr = malloc(sizeof(struct representation));
			FixRepresentation(repr, model, family);
			params = malloc(sizeof(double *) * (1 + repr->nvars));
			params[0] = malloc(sizeof(double) * 2);
			params[0][0] = (double) repr->nvars;
			params[0][1] = -1;
			int pidx;
			for (i = 0; i < repr->nvars; i ++){
				fgets(line, LINE_WIDTH, infp);
				sscanf(line, "%s", pname);
				pidx = IndexOf(pname, repr->varnames, repr->nvars);
				params[pidx + 1] = malloc(sizeof(double) * 3);
				sscanf(line, "%s %lf %lf %lf", pname, &params[pidx + 1][0], &params[pidx + 1][1], &params[pidx + 1][2]);
			}
		}
		if (strncmp(line, "implicit", 8) == 0){
			sscanf(line, "%s %s", option, pname);
			params[0][1] = (double) IndexOf(pname, repr->varnames, repr->nvars);
		}

		if (strncmp(line, "trials", 6) == 0)
			sscanf(line, "%s %ld", option, &stats);

		if (strncmp(line, "job", 3) == 0)
			sscanf(line, "%s %s %s %d", option, jobname, queue, &wall);

		if (strncmp(line, "axes", 4) == 0){
			ax = malloc(sizeof(struct axes));
			InitAxes(ax);
			sscanf(line, "%s %d %s %d %s %d %s %s", option, &ax->xcol, ax->xlabel, &ax->ycol, ax->ylabel, &ax->zcol, ax->zlabel, ax->scale);
		}

		if (strncmp(line, "mask", 4) == 0){
			if (ax != NULL){
				sscanf(line, "%s %f %f %f %f", option, &(ax->xlim)[0], &(ax->xlim)[1], &(ax->ylim)[0], &(ax->ylim)[1]);
			}
		}

		if (strncmp(line, "time", 4) == 0)
			sscanf(line, "%s %s", option, tstamp);
	}
	fclose(infp);

	if (strncmp(mode, "simulate", 8) == 0){
		printf("\033[92mInput done.\033[0m\n");
		printf("\033[92m******* SUMMARY *********\033[0m\n");
		printf("\033[92mFamily: %d\033[0m\n", family);
		printf("\033[92mCode lengths: \033[0m");
		for(i = 1; i <= lengths[0][0]; i ++)
			printf("\033[92m%d \033[0m", lengths[i][0]);
		printf("\n");
		printf("\033[92mNoise type: %s\n\033[0m", repr->mname);
		int nrates = 1;
		for (i = 1; i <= repr->nvars; i ++){
			nrates = nrates * (1 + (int) ((params[i][1] - params[i][0])/params[i][2]));
			printf("\033[92m%s : from %g to %g in steps of %g\033[0m\n", (repr->varnames)[i - 1], params[i][0], params[i][1], params[i][2]);
		}
		printf("\033[2m(In all, there are %d noise rates.)\033[0m\n", nrates);
		printf("\033[92mNumber of decoding trials: %ld\033[0m\n", stats);
		printf("\033[92m**************************\033[0m\n");

		// For each code in the family, we will load the tiling from respective files and launch Performance
		// One can parallelize the following for-loop.
		struct tiling T;
		char *tilingFile = malloc(sizeof(char) * 100), *resfile;
		clock_t startTime, endTime;
		double runtime;
		int ei, p;
		for(i = 1; i <= lengths[0][0]; i ++){
			TilingFile(tilingFile, family, lengths[i][0]);
			T = Load(tilingFile);
			// Number of qubits in the code defined by the above tiling
			lengths[i][1] = 0;
			for(ei = 0; ei < T.e; ei ++)
				if(T.doE[ei] == 0)
					lengths[i][1] ++;

			startTime = clock();
			printf("\033[93mEstimating performance of %s with noise rate\033[0m\n", tilingFile);
			for (p = 1; p <= repr->nvars ; p ++)
				printf("\033[93m\t%s -- from %g to %g in steps of %g.\033[0m\n", (repr->varnames)[p - 1], params[p][0], params[p][1], params[p][2]);
			resfile = Performance(T, lengths[i][2], model, params, stats);
			endTime = clock();
			runtime = (double)(endTime - startTime)/CLOCKS_PER_SEC;
			printf("\033[92m _/ Done %d qubits in %g seconds.\n\tOutput in %s.\033[0m\n", lengths[i][1], runtime, resfile);
		}
		printf("\033[92m_/ All simulations are done.\033[0m\n");
	}

	if (strncmp(mode, "plot", 4) == 0){
		// Gather all the performance data as well as runtime data and plot them
		time_t tme = time(NULL);
		struct tm tm = *localtime(&tme);
		char *current = malloc(sizeof(char) * 400);
		sprintf(current, "%d_%d_%d_%d_%d_%d", tm.tm_mday, tm.tm_mon, tm.tm_year + 1900, tm.tm_hour, tm.tm_min, tm.tm_sec);
		fprintf(stderr, "\033[2mSorting files ... \033[0m");
		if (strncmp(tstamp, "empty_time_stamp", 16) == 0){
			fprintf(stderr, "\033[2mError: no time stamp information provided. Plot cancelled.\033[0m\n");
		}
		else{
			Gather(family, lengths, model, tstamp, current);
			fprintf(stderr, "\033[2mdone\n\033[0m");

			if (ax == NULL){
				ax = malloc(sizeof(struct axes));
				InitAxes(ax);
				sprintf(ax->scale, "linear");
				ax->xcol = 0;
				ax->ycol = 1;
				ax->zcol = 2;
				sprintf(ax->xlabel, "column %d", ax->xcol);
				sprintf(ax->ylabel, "column %d", ax->ycol);
				sprintf(ax->zlabel, "column %d", ax->zcol);
			}
			PerformancePlot(family, model, lengths, ax, current);
			RuntimePlot(family, model, current);
			// fprintf(stderr, "Done plotting %s\n", current);
			printf("\033[92m_/ See the plot at reports/runtimes_%s_%s_%s.pdf\033[0m\n", repr->fname, repr->mname, current);

			free(current);
			FreeAxes(ax);
		}
	}

	if (strncmp(mode, "cluster", 7) == 0)
		Submit(inputfile, family, lengths, model, params, stats, jobname, queue, wall);

	// Free memory
	if (lengths != NULL)
		FreeIntArray(lengths, lengths[0][0] + 1);
	if (params != NULL)
		FreeDoubleArray(params, 1 + (int) params[0][0]);
	if (repr != NULL)
		FreeRep(repr);
	free(mode);
	free(option);
	free(jobname);
	free(queue);
	free(pname);
	free(tstamp);
}
