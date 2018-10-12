#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "memory.h"
#include "tiling.h"
#include "load.h"
#include "save.h"
#include "arrays.h"
#include "holes.h"
#include "regulartilings.h"
#include "performance.h"
#include "infer.h"

void ReadPerformance(char *timestamp, double *perfdata){
	// Read the performance data from a file
	int ct;
	char *fname = malloc(sizeof(char) * 100);
	sprintf(fname, "results/perf_%s.txt", timestamp);
	FILE *data = fopen(fname, "r");
	fscanf(data, "%*f");
	for (ct = 0; ct < 3; ct ++)
		fscanf(data, "%lf %*d %*d %*f", &perfdata[ct]);
	fscanf(data, "%*d %*f");
	fclose(data);
	// Delete the text files produced by Performance(...)
	DeleteFile(fname);
	char *selected = malloc(sizeof(char) * 100), *runtime = malloc(sizeof(char) * 100);
	sprintf(selected, "results/selected_perf_%s.txt", timestamp);
	DeleteFile(selected);
	sprintf(runtime, "results/runtimes_%s.txt", timestamp);
	DeleteFile(runtime);
	free(runtime);
	free(selected);
	free(fname);
}


double Validator(struct tiling G, double p, int errortype, double target, double ltol){
	// Test if the logical error rate for the given physical noise rate p, is equal to the target (up to tolerance).
	// If yes, return 1 and the associated error bar.
	// If not
	// 		if the logical error rate is either 0 or 1, it is not a good reference for predictions, so return 0.
	// 		else return the measured logical failure rate
	const double NUMERICAL_ZERO = 10E-20;
	const long MAX_STATS = 10E10, MIN_STATS = 1000;
	long stats = ceil(pow(10, 1 + abs(target)));
	stats = (stats > MAX_STATS) ? MAX_STATS : stats;
	stats = (stats < MIN_STATS) ? MIN_STATS : stats;

	double **noiseparams = malloc(sizeof(double) * 2);
    noiseparams[0] = malloc(sizeof(double));
	noiseparams[0][0] = 1;
    noiseparams[1] = malloc(sizeof(double) * 3);
    noiseparams[1][0] = p;
    noiseparams[1][1] = p;
	noiseparams[1][2] = 0.1;
	
	char *timestamp = malloc(sizeof(char) * 100);
	sprintf(timestamp, "%s", Performance(G, -1, 1, noiseparams, stats));
	FreeDoubleArray(noiseparams, 2);
	double *perfdata = malloc(sizeof(double) * 3);
	ReadPerformance(timestamp, perfdata);
	double validated;
	if ((abs(log10(perfdata[errortype])) >= abs(target) * (1 - ltol)) && (abs(log10(perfdata[errortype])) <= abs(target) * (1 + ltol)))
		validated = 1;
	else{
		if ((perfdata[errortype] >= 1 - NUMERICAL_ZERO) || (perfdata[errortype] <= NUMERICAL_ZERO))
			validated = 0;
		else
			validated = log10(perfdata[errortype]);
	}
	free(timestamp);
	free(perfdata);
	return validated;
}


double Predictor(double **trainingSet, double test){
	// Given the coordinates of a set of points, find the best fit line
	// Usig the best fit line and the y-coordinate of a new point, determine its corresponding x-coordinate.
	// The best fit line is given by y = a + bx
	// where b = sum((xi - xm)*(yi - ym))/sum((xi - xm)^2), where sum is over i = 1 to N and xm denotes the mean value.
	// and a = ym - b * xm.
	// So, the x-coordinate of the new point will have
	// xnew = (ynew - a)/b
	double xmean = 0, ymean = 0, bestslope, bestintercept;
	int i;
	for (i = 1; i <= (int) trainingSet[0][0]; i ++){
		xmean += trainingSet[i][0];
		ymean += trainingSet[i][1];
	}
	xmean = xmean/trainingSet[0][0];
	ymean = ymean/trainingSet[0][0];

	double correlation = 0, variance = 0;
	for (i = 1; i <= (int) trainingSet[0][0]; i ++){
		variance += pow((trainingSet[i][0] - xmean), 2);
		correlation += (trainingSet[i][0] - xmean) * (trainingSet[i][1] - ymean);
	}
	bestslope = correlation/variance;
	bestintercept = ymean - bestslope * xmean;
	double xnew = (test - bestintercept)/bestslope;

	// fprintf(stderr, "The best fit line is y = %g x + %g. So y = %g imples x = %g\n", bestslope, bestintercept, test, xnew);

	return xnew;
}


double Train(struct tiling G, double **trainingSet, int errortype, double target, double ltol){
	// Prepare a training set, containing a few physical noise rates and their correspondng failure rates.
	// This will help the Predictor(...) to do a simple linear regression.
	// Since we do regression on a log-log scale, it would be convenient to consider points that are euqally separated on a log scale.
	// The initial points in the training set we assume are of the form (2/3)^k for a few values of k, starting from 2.5, in steps of 0.1.
	// We populate the training set until we find 3 physical noise rates.
	// If the target failure rate is achieved using an erasure rate that is part of the training set itself,
	// return the noise rate value in the training set.
	// Else return -1.
	// Sometimes the points in the training set might be in a region where the curve has a small slope.
	// In this case, the interpolation in Predictor(...) will yield is very small physical noise rate.
	// To avoid this, after a training set has been built we will verify if the prediction is actually less than the target logical failure rate itself.
	// If so, then we redesign the training set by choosing points further down the performance curve, where holefully the curve is steeper.
	int startupSize = 10;
	double noiseBase = 2/(double) 3, noisesteps[2] = {2, 0.1}, noiseLimit = 7, measured, eraseRate, noiseExp = noisesteps[0], prediction = 0;

	do{
		trainingSet[0][0] = 0;
		while ((int) trainingSet[0][0] < startupSize){
			eraseRate = pow(noiseBase, noiseExp);
			measured = Validator(G, eraseRate, errortype, target, ltol);
			if (measured == 1)
				return eraseRate;
			else{
				if (measured < 0){
					trainingSet[0][0] ++;
					trainingSet[(int) trainingSet[0][0]][0] = log10(eraseRate);
					trainingSet[(int) trainingSet[0][0]][1] = measured;
				}
				noiseExp += noisesteps[1];
			}
			if (noiseExp >= noiseLimit)
				break;
		}
		prediction = Predictor(trainingSet, target);
		noiseExp ++;
	} while(prediction < target);
	return -1;
}

void RemoveInitData(double **trainingSet){
	// remove the initial data point in the training set and shift all the data points one cell to the left
	int i;
	for (i = 2; i <= trainingSet[0][0]; i ++){
		trainingSet[i - 1][0] = trainingSet[i][0];
		trainingSet[i - 1][1] = trainingSet[i][1];
	}
	trainingSet[0][0] --;
}

double Infer(struct tiling G, int errortype, double scaling){
	// Infer the physical noise rate that results in a given logical error rate
	const double ltol = 0.1;
	const int maxTrainingSize = 50, maxLearnSteps = 10;
	double **trainingSet = malloc(sizeof(double) * maxTrainingSize);
	int t;
	for (t = 0; t < (1 + maxTrainingSize); t ++)
		trainingSet[t] = malloc(sizeof(double) * 2);

	trainingSet[0][0] = 0;
	double physical = Train(G, trainingSet, errortype, scaling, ltol);
	if (physical > 0)
		return physical;
	else{
		double predicted, measured;
		for (t = 0; t < maxLearnSteps; t ++){
			predicted = Predictor(trainingSet, scaling);
			measured = Validator(G, pow(10, predicted), errortype, scaling, ltol);
			if (measured == 1)
				return pow(10, predicted);
			else{
				if (measured < 0){
					trainingSet[0][0] ++;
					trainingSet[(int) trainingSet[0][0]][0] = predicted;
					trainingSet[(int) trainingSet[0][0]][1] = measured;
				}
				else
					// the prediction results in a very small logical error rate.
					// So we will remove some initial points to steepen the best fit line.
					RemoveInitData(trainingSet);
			}
		}
	}
	FreeDoubleArray(trainingSet, 1 + maxTrainingSize);
	fprintf(stderr, "Inference failed after %d trails. Try increasing learing steps or the tolerance level.\n", t);
	return 0;
}
