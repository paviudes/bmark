#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "memory.h"
#include "arrays.h"
#include "notations.h"
#include "tiling.h"
#include "rgens.h"
#include "rbim.h"
#include "mcresult.h"
#include "simulation.h"
#include "correlations.h"
#include "noise.h"

void InitializeNoise(struct noise *pn, struct tiling *pG, int dist, int model, double **prange){
	// Initialize the parameters of a noise model
	int i, nvals, nparams, isexpanded;
	pn->model = model;
	// Part 2,3
	double **expanded = malloc(sizeof(double *) * (1 + (int) prange[0][0]));
	expanded[0] = malloc(sizeof(double));
	expanded[0][0] = prange[0][0];
	for(i = 1; i <= (int) prange[0][0]; i ++){
		nvals = IntervalSize(prange[i]);
		expanded[i] = malloc(sizeof(double) * (nvals + 1));
		expanded[i][0] = (double) nvals;
		isexpanded = 0;
		if (pn->model == 5){
			Expand(prange[i], expanded[i], 1);
			isexpanded = 1;
		}
		if ((pn->model == 2) || (pn->model == 3) || (pn->model == 4)){
			if (i > 1){
				Expand(prange[i], expanded[i], 1);
				isexpanded = 1;
			}
		}
		if (isexpanded == 0)
			Expand(prange[i], expanded[i], 0);
	}
	// Part 4,5
	nparams = (int)(1 + (int) Prod(Slice(expanded + 1, (int) expanded[0][0], 0), (int) expanded[0][0]));
	pn->params = malloc(sizeof(double *) * nparams);
	(pn->params)[0] = malloc(sizeof(double) * 2);
	for (i = 1; i < nparams; i ++)
		(pn->params)[i] = malloc(sizeof(double) * (int) ceil(prange[0][0]));
	CartesianProduct(expanded, pn->params);

	InitializeRBIMNoise(pn, pG);
	InitImportanceSampling(pn, pG);

	if (prange[0][1] != -1){
		if (pn->model == 2)
			for (i = 1; i < nparams; i ++)
				(pn->params)[i][(int) prange[0][1]] = BallInfer((pn->params)[i], (int) prange[0][1]);
		if ((pn->model == 3) || (pn->model == 4))
			for (i = 1; i < nparams; i ++)
				(pn->params)[i][(int) prange[0][1]] = WalkInfer((pn->params)[i], (int) prange[0][1]);
		if (pn->model == 5){
			int *sizes = malloc(sizeof(int) * ((int) (pn->params)[0][1]));
			for (i = 0; i < (int) (pn->params)[0][1]; i ++)
				sizes[i] = (int) expanded[i + 1][0];
			RBIMInfer(pn->ising, (pn->params), sizes, (int) prange[0][1]);
			free(sizes);
		}
	}

	pn->current = 0;
	pn->dist = dist;

	pn->corrdist = malloc(sizeof(double) * (1 + pG->e));
	// Free local memory
	FreeDoubleArray(expanded, (1 + (int) prange[0][0]));

	// physical properties of the noise model
	struct representation *repr = malloc(sizeof(struct representation));
	FixRepresentation(repr, model, 0);
	(pn->properties) = malloc(sizeof(double) * (1 + repr->nprops));
	(pn->properties)[0] = repr->nprops;
	FreeRep(repr);
}


void InitImportanceSampling(struct noise *pn, struct tiling *pG){
	// Initialize the distributions required for importance sampling. For now, it only applies to uncorrelated noise model (model 1).
	if (pn->model == 1){
		if (pn->dist == -1){
			pn->truedist = NULL;
			pn->impdist = NULL;
			pn->impcumul = NULL;
			pn->bias = NULL;
		}
		else{
			pn->truedist = malloc(sizeof(double) * pG->e);
			pn->impdist = malloc(sizeof(double) * pG->e);
			pn->impcumul = malloc(sizeof(double) * pG->e);
			pn->bias = malloc(sizeof(double) * pG->e);
		}
	}
}


void InitializeRBIMNoise(struct noise *pn, struct tiling *pG){
	// Initialize the Random bond Ising model according to a lattice, with arbitrary coupling strengths, field strengths and temperature
	if (pn->model == 5){
		pn->ising = malloc(sizeof(struct RandomBondIsingModel));
		InitRBIM(pn->ising);
		IntializeRBMWithCode(pn->ising, pG, 0, 0, 1);
	}
	else
		pn->ising = NULL;
}

double Binomial(int n, int k, double p){
	// Compute the term nCk p^k (1-p)^(n-k).
	// We will assume that k < n-k.
	// Note that the Binomial coefficient can be expressed as: (1-p)^n * prod_(i=0)^(n-k) [ (n-i)/(k-i) q ] , where q = p/(1-p).
	double q = p/(1-p), result = 1;
	int i;
	for (i = 0; i < k; i ++)
		result = result * (n - i)/((double) (k - i)) * q;
	result = result * pow(1 - p, (double) n);
	return result;
}

double BiasedBinomial(int n, int dist, double p, double *distribution){
	// Compute the probability distribution {Prob(a binary sequence of weight k) | k >= d}.
	double normalization = 0;
	int i;
	for (i = 0; i < n - dist; i ++){
		distribution[i] = Binomial(n, dist + i, p);
		normalization += distribution[i];
	}
	for (i = 0; i < n - dist; i ++)
		distribution[i] = distribution[i]/normalization;
	return normalization;
}

double GaussianPDF(double mean, double variance, int point){
	// Compute the value of the PDF of the normal distribution with given mean and variance at a given point.
	double pi = 3.1415, value = (1/sqrt(2 * variance * pi) * exp(-pow(point - mean, 2)/(2 * variance)));
	return value;
}

void Gaussian(double mean, double variance, int nelems, double *distribution){
	// Construct a gaussian distribution with given mean standard deviation.
	const double atol = 10E-40;
	double normalization = 0;
	int i;
	for (i = 0; i < nelems; i ++){
		distribution[i] = GaussianPDF(mean, variance, i);
		if (distribution[i] < atol)
			distribution[i] = 0;
		normalization += distribution[i];
	}
	for (i = 0; i < nelems; i ++)
		distribution[i] = distribution[i]/normalization;
	// printf("Gaussian distribution\n");
	// PrintDoubleArray(distribution, nelems);
}

void Poisson(double mean, double *cumulative, int cutoff){
	// Construct the CDF of a Poisson distribution with given mean
	int i, j;
	double normalization = 0;
	double *distribution = malloc(sizeof(double) * cutoff);
	for (i = 0; i < cutoff; i ++){
		distribution[i] = exp((-1) * mean);
		for (j = 1; j < i; j ++)
			distribution[i] = distribution[i] * mean/(double) j;
		normalization += distribution[i];
	}
	for (i = 0; i < cutoff; i ++)
		distribution[i] = distribution[i]/normalization;
	cumulative[0] = (double) cutoff;
	cumulative[1] = distribution[0];
	for (i = 1; i < cutoff; i ++)
		cumulative[i + 1] = cumulative[i] + distribution[i];
	free(distribution);
}

int UpdateNoise(struct noise *pn, struct mcresult *pmcr, char *tileName, int ne){
	// Update the parameters of a noise model and assign new filenames for the results of the various montecarlo runs.
	// const double atol = 10E-40;
	int skip = 0;
	pn->current ++;
	int cutoff = (int) Min(3 * (pn->params)[pn->current][1], (double) ne);
	if ((pn->model == 2) || (pn->model == 3) || (pn->model == 4)){
		Poisson((pn->params)[pn->current][1], pn->corrdist, cutoff);

	if (pn->current > 1){
		if (pn->model == 1)
			if (((pmcr->fails)[0] + (pmcr->fails)[1] + (pmcr->fails)[2]) == 0)
				skip = 1;
		}
		// if (pn->model == 4)
		// 	skip = RBIMSkipTest((pn->params)[pn->current - 1], pmcr->fails, pmcr->trials, (pn->params)[pn->current]);
	}
	if (skip == 0)
		FillDoubleArrayWithConstants(pn->properties + 1, (int) (pn->properties)[0], 0);
	
	if (pn->model == 1){
		if (pn->dist != -1){
			// The importance distribution PDF is a Gaussian whose mean is at d.
			// If np > d, then its the binomial distribution.

			// True weight distribution is given by the Binomial distribution
			int i;
			for (i = 0; i < ne; i ++){
				if (i < ne/2)
					(pn->truedist)[i] = Binomial(ne, i, (pn->params)[pn->current][0]);
				else
					(pn->truedist)[i] = 0;
			}

			// printf("True distribution for\n");
			// PrintNoise(pn);
			// PrintDoubleArray(pn->truedist, ne);

			// printf("n p = %d x %g = %g and distance = %g\n", ne, (pn->params)[pn->current][0], (pn->params)[pn->current][0] * (double) ne, (double) pn->dist);
			if ((pn->params)[pn->current][0] * (double) ne < (double) (pn->dist + 1)){
				// printf("Gaussian with mean = %g, variance = %g\n", (double) (pn->dist + 1), ceil(ne * (pn->params)[pn->current][0] * (1 - (pn->params)[pn->current][0])));
				Gaussian((double) pn->dist, ceil(ne * (pn->params)[pn->current][0]), ne, pn->impdist);
			}
			else{
				// printf("There is no need for importance sampling.\n");
				for (i = 0; i < ne; i ++){
					(pn->impdist)[i] = (pn->truedist)[i];
					(pn->bias)[i] = 1;
				}
			}

			// printf("Importance distribution for\n");
			// PrintNoise(pn);
			// PrintDoubleArray(pn->impdist, ne);

			// printf("atol = %g\n", atol);

			(pn->impcumul)[0] = (pn->impdist)[0];
			for (i = 0; i < ne; i ++){
				// printf("%d %g %g %g\n", i, (pn->truedist)[i], (pn->impdist)[i], atol);
				if ((pn->truedist)[i] <= 0){
					(pn->impdist)[i] = 0;
					(pn->bias)[i] = 0;
				}
				else
					(pn->bias)[i] = (pn->truedist)[i]/(pn->impdist)[i];
				if (i > 0)
					(pn->impcumul)[i] = (pn->impcumul)[i - 1] + (pn->impdist)[i];
			}

			// printf("Importance distribution for\n");
			// PrintNoise(pn);
			// PrintDoubleArray(pn->impdist, ne);
		}
	}
	if (pn->model == 5)
		UpdateRBIM(pn->ising, (pn->params)[pn->current][0], (pn->params)[pn->current][1], (pn->params)[pn->current][2]);
	ResetResults(pmcr);
	return skip;
}


void FreeNoise(struct noise *pn){
	// Free a noise structre
	FreeDoubleArray(pn->params, 1 + (pn->params)[0][0]);
	if (pn->ising != NULL)
		FreeRBIM(pn->ising);
	if (pn->truedist != NULL)
		free(pn->truedist);
	if (pn->impdist != NULL)
		free(pn->impdist);
	if (pn->bias != NULL)
		free(pn->bias);
	free(pn->properties);
	free(pn->corrdist);
}

void PrintNoise(struct noise *pn){
	// Print the values of the noise parameters, that correspond to the current simulation when the function is called.
	// "Current" noise parameters are in (pn->params)[pn->current]
	int i;
	struct representation *repr = malloc(sizeof(struct representation));
	FixRepresentation(repr, pn->model, 0);
	fprintf(stderr, "%s = %g", (repr->varnames)[0], (pn->params)[pn->current][0]);
	for (i = 1; i < (pn->params)[0][1]; i ++)
		fprintf(stderr, ", %s = %g", (repr->varnames)[i], (pn->params)[pn->current][i]);
	fprintf(stderr, ".\n");
	if (repr->nprops > 0)
		fprintf(stderr, "\033[92m%s = %g\033[0m", (repr->propnames)[0], (pn->properties)[1]);
	for (i = 1; i < repr->nprops; i ++)
		fprintf(stderr, "\033[92m, %s = %g\033[0m", (repr->propnames)[i], (pn->properties)[i + 1]);
	fprintf(stderr, "\033[92m.\n\033[0m");
	FreeRep(repr);
}
