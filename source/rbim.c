#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include "tiling.h"
#include "arrays.h"
#include "load.h"
#include "rgens.h"
#include "notations.h"
#include "memory.h"
#include "rbim.h"

void InitRBIM(struct RandomBondIsingModel *prbim){
	// Initialize a random bond Ising model
	prbim->nSites = 0;
	prbim->nCouplings = 0;
	prbim->range = 0;
	prbim->temperature = 0;
	prbim->sToC = NULL;
	prbim->cToS = NULL;
	prbim->config = NULL;
	prbim->interactions = NULL;
	prbim->fields = NULL;
	prbim->energy = 0;
	prbim->nclusters = 0;
	prbim->maxClusterSize = 0;
	prbim->avgClusterSize = 0;
	prbim->magnetization = 0;
	prbim->internalEnergy = 0;
	prbim->twopoint = 0;
}

void FreeRBIM(struct RandomBondIsingModel *prbim){
	// Free memory allocated to the structure
	if (prbim->interactions != NULL)
		free(prbim->interactions);
	if (prbim->fields != NULL)
		free(prbim->fields);
	if (prbim->config != NULL)
		free(prbim->config);
	if (prbim->sToC != NULL)
		FreeIntArray(prbim->sToC, prbim->nSites);
	if (prbim->cToS != NULL)
		FreeIntArray(prbim->cToS, prbim->nCouplings);
	free(prbim);
}


void PrintRBIM(struct RandomBondIsingModel *prbim, int isDetailed){
	// Print the details of a random bond Ising model object
	int si, ci;
	printf("Details of the Random bond Ising model\n");
	printf("Number of sites = %d\n", prbim->nSites);
	printf("Number of couplings = %ld\n", prbim->nCouplings);
	printf("Range of interaction = %d\n", prbim->range);
	printf("Temperature = %g\n", (double) (prbim->temperature));

	if (isDetailed == 1){
		printf("Hamiltonian = \n");
		for (ci = 0; ci < prbim->nCouplings; ci ++){
			if ((prbim->interactions)[ci] < 0)
				printf(" - ");
			else
				printf(" + ");
			printf("%g ( ", (double) fabs((prbim->interactions)[ci]));
			for (si = 1; si <= (prbim->cToS)[ci][0]; si ++)
				printf(" Z_%d", (prbim->cToS)[ci][si]);
			printf(" )\n");
		}
		for (si = 0; si < prbim->nSites; si ++){
			if ((prbim->fields)[si] < 0)
				printf(" - ");
			else
				printf(" + ");
			printf("%g ( X_%d ) ", (double) fabs((prbim->fields)[si]), si);
		}
		printf(".\n");
	}

	printf("The current configuration of spins is\n");
	if (isDetailed)
		for (si = 0; si < prbim->nSites; si ++)
			printf("%d ", (prbim->config)[si]);
	else
		printf("<A large array ... not printed here>");
	printf("\nIt has %d spin 'up's and its energy is = %g.\n", (int) Sum(prbim->config, prbim->nSites), prbim->energy);
}


long AddCouplings(struct RandomBondIsingModel *prbim, int *check, long coup){
	// Add the locally generated couplings from a check, to the set of all couplings recorded so far.
	long nlocal = Combination(check[0], prbim->range);
	int **localCouplings = malloc(sizeof(int *) * (1 + nlocal));
	localCouplings[0] = malloc(sizeof(int) * 2);
	int c;
	for (c = 1; c <= nlocal; c ++)
		localCouplings[c] = malloc(sizeof(int) * prbim->range);

	GenSubsets(check, prbim->range, localCouplings);

	int exists, l, ei;
 	for (l = 1; l <= nlocal; l ++){
 		// Test if this coupling already exists in the set of available couplings
 		exists = 0;
 		for (c = 0; c < coup; c ++)
 			exists += IsEqual(localCouplings[l], (prbim->cToS)[c] + 1, (prbim->cToS)[c][0]);
 		if (exists == 0){
 			(prbim->cToS)[coup][0] = prbim->range;
 			for (ei = 0; ei < prbim->range; ei ++){
 				(prbim->cToS)[coup][1 + ei] = localCouplings[l][ei];
 			}
	 		coup ++;
 		}
 	}
 	// Free memory of local couplings
	FreeIntArray(localCouplings, (int)(1 + nlocal));
 	return coup;
}

void AllocateCouplingsToSites(struct RandomBondIsingModel *prbim, struct tiling* pG){
	// Allocate values to a list describing the edges that are involved in a coupling, for every coupling
	// Generate all unique (up to permutation) subsets of size given by the range, out of qubits incident on the vertex vi
	int i;
	for (i = 0; i < pG->v; i ++)
		prbim->nCouplings = AddCouplings(prbim, (pG->VToE)[i], prbim->nCouplings);
	// Generate all unique (up to permutation) subsets of size given by the range, out of qubits incident on the face fi
	for (i = 0; i < pG->f; i ++)
		prbim->nCouplings = AddCouplings(prbim, (pG->F)[i], prbim->nCouplings);
}


void AllocateSiteToCouplings(struct RandomBondIsingModel *prbim){
	// Allocate values to a list describing the all couplings that involve a particular site, for every site.
	int c, s, site;
	int *degrees = malloc(sizeof(int) * prbim->nSites);
	for (s = 0; s < prbim->nSites; s ++)
		degrees[s] = 0;

	// compute the degree of every site
	for (c = 0; c < prbim->nCouplings; c ++)
		for (s = 1; s <= (prbim->cToS)[c][0]; s ++)
			degrees[(prbim->cToS)[c][s]] ++;

	prbim->sToC = malloc(sizeof(int *) * prbim->nSites);
	for (s = 0; s < prbim->nSites; s ++){
		(prbim->sToC)[s] = malloc(sizeof(int) * (1 + degrees[s]));
		(prbim->sToC)[s][0] = 0;
	}

	for (c = 0; c < prbim->nCouplings; c ++){
		for (s = 1; s <= (prbim->cToS)[c][0]; s ++){
			site = (prbim->cToS)[c][s];
			(prbim->sToC)[site][0] ++;
			(prbim->sToC)[site][(prbim->sToC)[site][0]] = c;		
		}
	}
	free(degrees);
}


void Energy(struct RandomBondIsingModel *prbim){
	// Compute the energy of a configuration of a random bond Ising model.
	int ci, si;
	double local;
	prbim->energy = 0;
	for (ci = 0; ci < prbim->nCouplings; ci ++){
		local = (-1) * (prbim->interactions)[ci];
		for (si = 1; si <= (prbim->cToS)[ci][0]; si ++){
			local *= pow(-1, (prbim->config)[(prbim->cToS)[ci][si]]);
		}
		prbim->energy += local;
	}
	// Contribution from the on-site terms
	for (si = 0; si < prbim->nSites; si ++)
		prbim->energy += (-1) * (prbim->fields)[si] * pow(-1, (prbim->config)[si]);
	// prbim->energy = (double) ((prbim->energy)/(prbim->temperature));
}

void UpdateRBIM(struct RandomBondIsingModel *prbim, double interact, double field, double temp){
	// Update a random bond Ising model by assigning new coupling strengths to interactions and few field strengths
	// The parameters are ordered as follows.
	// uniform interaction strength, uniform field strength, temperature
	int si, ci;
	for (ci = 0; ci < prbim->nCouplings; ci ++)
		(prbim->interactions)[ci] = interact;
	for (si = 0; si < prbim->nSites; si ++)
		(prbim->fields)[si] = field;
	prbim->temperature = temp;
	// if (field < 0)
	// 	FillIntArrayWithConstants(prbim->config, prbim->nSites, 1);
	// else
	// 	FillIntArrayWithConstants(prbim->config, prbim->nSites, 0);
	// Energy(prbim);
	// BinomialBitString(prbim->config, prbim->nSites, NULL, 0.5);
	FillIntArrayWithConstants(prbim->config, prbim->nSites, 0);
	prbim->nclusters = 0;
	prbim->maxClusterSize = 0;
	prbim->avgClusterSize = 0;
	prbim->magnetization = 0;
	prbim->internalEnergy = 0;
}


void IntializeRBMWithCode(struct RandomBondIsingModel *prbim, struct tiling* pG, double coupStrength, double fieldStrength, double temperature){
	// Initialize a random bond Ising model using the specifications of a tiling that defined an error correcting code
	prbim->nSites = pG->e;
	prbim->range = 2;
	
	long totalcouplings = 0;
	int i;
	for (i = 0; i < pG->v; i ++)
		totalcouplings += Combination((pG->VToE)[i][0], prbim->range);
	for (i = 0; i < pG->f; i ++)
		totalcouplings += Combination((pG->F)[i][0], prbim->range);

	prbim->cToS = malloc(sizeof(int *) * totalcouplings);
	int c;
	for (c = 0; c < totalcouplings; c ++)
		(prbim->cToS)[c] = malloc(sizeof(int) * (prbim->range + 1));

	AllocateCouplingsToSites(prbim, pG);
	// Reallocate the space reserved for couplings since extra space was allocated.
	for (c = prbim->nCouplings; c < totalcouplings; c ++)
		free((prbim->cToS)[c]);
	prbim->cToS = realloc(prbim->cToS, sizeof(int *) * prbim->nCouplings);

	AllocateSiteToCouplings(prbim);

	prbim->interactions = malloc(sizeof(double) * prbim->nCouplings);
	for (c = 0; c < prbim->nCouplings; c ++)
		(prbim->interactions)[c] = coupStrength; // Ferromagnetic interactions
	prbim->fields = malloc(sizeof(double) * prbim->nSites);
	int s;
	for (s = 0; s < prbim->nSites; s ++)
		(prbim->fields)[s] = fieldStrength;
	prbim->temperature = temperature;

	prbim->config = malloc(sizeof(int) * prbim->nSites);

	// compute the coordination number
	int *neighbours = malloc(sizeof(int) * prbim->nSites);
	FillIntArrayWithConstants(neighbours, prbim->nSites, 0);
	for (c = 1; c <= (prbim->sToC)[0][0]; c ++)
		for (s = 1; s <= (prbim->cToS)[(prbim->sToC)[0][c]][0]; s ++)
			neighbours[(prbim->cToS)[(prbim->sToC)[0][c]][s]] = 1;
	prbim->coordination = Sum(neighbours, prbim->nSites) - 1;
	free(neighbours);
}


double SiteContribution(struct RandomBondIsingModel *prbim, long site){
	// Compute the energy contribution coming from a particular site
	int c, s, coup, neigh;
	double energy = 0, local;
	// fprintf(stderr, "site contribution at %ld which has %d couplings\n", site, (prbim->sToC)[site][0]);
	for (c = 1; c <= (prbim->sToC)[site][0]; c ++){
		coup = (prbim->sToC)[site][c];
		local = 1;
		for (s = 1; s <= (prbim->cToS)[coup][0]; s ++){
			neigh = (prbim->cToS)[coup][s];
			local = local * pow(-1, (prbim->config)[neigh]);
		}
		energy = energy + local * (-1) * (prbim->interactions)[coup];
	}
	// Onsite interaction terms
	energy = energy - (prbim->fields)[site] * pow(-1, (prbim->config)[site]);
	// energy = (double) (energy/(prbim->temperature));
	return energy;
}


double FlipCost(struct RandomBondIsingModel *prbim, long site){
	// Compute the energy cost associated to flipping the spin at a particular site.
	// The spin is flipped and flipped back to it's original value.
	int ei;
	double cost;
	double *energies = malloc(sizeof(double) * 2);
	for (ei = 0; ei < 2; ei ++){
		energies[ei] = SiteContribution(prbim, site);
		(prbim->config)[site] = (1 + (prbim->config)[site]) % 2;
	}
	cost = energies[1] - energies[0];
	// fprintf(stderr, "cost = %g\n", cost);
	free(energies);
	return cost;
}


int ExploreSites(struct RandomBondIsingModel *prbim, int site, int *visited, int state){
	// Visit all the sites in a connected component of spin-ups(downs), starting from a site
	visited[site] = 1;
	int ci, si;
	for (ci = 1; ci <= (prbim->sToC)[site][0]; ci ++)
		for (si = 1; si <= (prbim->cToS)[(prbim->sToC)[site][ci]][0]; si ++)
			if (visited[(prbim->cToS)[(prbim->sToC)[site][ci]][si]] == 0)
				if ((prbim->config)[(prbim->cToS)[(prbim->sToC)[site][ci]][si]] == state)
					return (1 + ExploreSites(prbim, (prbim->cToS)[(prbim->sToC)[site][ci]][si], visited, state));
	return 1;
}


void ClusterSize(struct RandomBondIsingModel *prbim, int state, double *clusters){
	// Find the average and maximum sizes of clusters of spin-ups(downs)
	int *visited = malloc(sizeof(int) * prbim->nSites);
	FillIntArrayWithConstants(visited, prbim->nSites, 0);
	FillDoubleArrayWithConstants(clusters, 3, 0);
	int si, csize;
	for (si = 0; si < prbim->nSites; si ++){
		if ((prbim->config)[si] == state){
			if (visited[si] == 0){
				csize = ExploreSites(prbim, si, visited, state);
				clusters[0] ++;
				if (clusters[1] < csize)
					clusters[1] = csize;
				clusters[2] += (double) csize;
			}
		}
	}
	free(visited);
	// fprintf(stderr, "maximum cluster size = %g, number of cluster = %g and average cluster size = %g.\n", clusters[1], clusters[0], clusters[2]);
	if (clusters[0] > 0)
		clusters[2] = (double) (clusters[2]/((double) clusters[0]));
}

void Magnetization(struct RandomBondIsingModel *prbim){
	// compute the magnetization with the current configuration of the Ising model
	prbim->magnetization = 0;
	int s;
	for (s = 0; s < prbim->nSites; s ++)
		prbim->magnetization += pow(-1, (prbim->config)[s]);
	prbim->magnetization = prbim->magnetization/((double) (prbim->nSites));
}

void InternalEnergy(struct RandomBondIsingModel *prbim){
	// compute the internal energy with the current configuration of the Ising model
	prbim->internalEnergy = 0;
	int s;
	for (s = 0; s < prbim->nSites; s ++)
		prbim->internalEnergy += SiteContribution(prbim, s) + (prbim->fields)[s] * pow(-1, (prbim->config)[s]);
	prbim->internalEnergy = prbim->internalEnergy/((double) (prbim->nSites));
}

void TwoPointCorrelation(struct RandomBondIsingModel *prbim){
	// compute the two-point correlation fucntion with the current configuration of the Ising mdoel
	prbim->twopoint = 0;
	int s, os;
	for (s = 0; s < prbim->nSites; s ++)
		for (os = 0; os < s; os ++)
			prbim->twopoint += pow(-1, (prbim->config)[s]) * pow(-1, (prbim->config)[os]);
	prbim->twopoint = fabs(prbim->twopoint/(double) (prbim->nSites * (prbim->nSites - 1)/(double) 2));
}

void PhysicalProperties(struct RandomBondIsingModel *prbim){
	// Compute the average magnetization
	Magnetization(prbim);
	InternalEnergy(prbim);
	// TwoPointCorrelation(prbim);
	prbim->twopoint = 0;
	// compute the cluster sizes of spins pointing up.
	double *clusters = malloc(sizeof(double) * 3);
	ClusterSize(prbim, 1, clusters);
	prbim->nclusters = clusters[0];
	prbim->maxClusterSize = clusters[1];
	prbim->avgClusterSize = clusters[2];
	free(clusters);
}

void Cool(struct RandomBondIsingModel *prbim, int nsweeps){
	/*
	Cool a random bond Ising model, using the following procedure.
	1. At every site,
		A. Compute the energy cost associated to flipping a spin at this site.
		B. If the energy cost is negative (i.e., the spin flip lowers the energy), make the spin flip
		C. Else, make the spin flip, with probability exp(- energy cost).
	*/
	int m, s;
	double cost;
	for (m = 0; m < nsweeps; m ++){
		for (s = 0; s < prbim->nSites; s ++){
			cost = FlipCost(prbim, s);
			if (cost < 0)
				(prbim->config)[s] = (1 + (prbim->config)[s]) % 2;
			else
				if (Random("uniform", 0, NULL) < exp((-1) * cost/prbim->temperature))
					(prbim->config)[s] = (1 + (prbim->config)[s]) % 2;
		}
	}
}

double MeanFieldInfer(double *available, int coordination, int toinfer){
	// Infer the parameters of the Ising model using mean-field theory
	// the mean field equation is: m = tanh( (-B - J q m)/T )
	double infered = 0, magnetization = 1 - 2 * available[toinfer];
	if (toinfer == 0){
		// we need to infer the value of J. We have J = (- B - T atanh(m))/(q m)
		infered = ((-1) * available[1] - available[2] * atanh(magnetization))/((double) coordination * magnetization);
	}
	else if (toinfer == 1){
		// we need to infer the value of B. We have B = (- J q m - T atanh(m))
		infered = (-1) * available[0] * (double) coordination * magnetization - available[2] * atanh(magnetization);
	}
	else if (toinfer == 2){
		// we need to infer the value of T. We have T = (- B - J q m)/atanh(m)
		infered = ((-1) * available[1] - available[0] * (double) coordination * magnetization)/atanh(magnetization);
	}
	else;
	return infered;
}

void RBIMDryRun(struct RandomBondIsingModel *prbim, double *fixed, int tovarry, double **trainingset){
	// Perform a dry run of the Ising model to get an idea of how the magnetization behaves with one of the parameters, J or B.
	const int mcsteps = 100;
	const int coolsteps = 10;
	const int initialize = 1000;
	// const double tol = 10E-4;
	double empirical;
	int i, m;
	fprintf(stderr, "\033[2m[Training] With J = %g, T = %g and %d montecarlo steps, %d cooling steps.\033[0m\n", fixed[0], fixed[2], mcsteps, coolsteps);
	for (i = 1; i <= (int) trainingset[0][0]; i ++){
		fixed[tovarry] = trainingset[i][0];
		// Estimate the empirical noise rate for the chosen set of Ising model parameters
		UpdateRBIM(prbim, fixed[0], fixed[1], fixed[2]);
		// before collecting statistics, we will cool the Ising model to ensure that we sample the Gibbs distribution.
		Cool(prbim, initialize);
		empirical = 0;
		for (m = 0; m < mcsteps; m ++){
			Cool(prbim, coolsteps);
			Magnetization(prbim);
			empirical += (1 - prbim->magnetization)/(double) 2;
		}
		trainingset[i][1] = empirical/((double) mcsteps);
		// fprintf(stderr, "%d). B = %g, p = %g.\n", i, trainingset[i][0], trainingset[i][1]);
		fprintf(stderr, "\r\033[2m%.1f%% done.\033[0m", 100 * ((double) i)/trainingset[0][0]);
	}
	fprintf(stderr, "\n\033[2m[Training] Set has %d points.\nB           p\n\033[0m", (int) trainingset[0][0]);
	for (i = 1; i <= (int) trainingset[0][0]; i ++)
		fprintf(stderr, "\033[2m%.5f %.5f\033[0m\n", trainingset[i][0], trainingset[i][1]);
}

double Model(double *fitparams, double field){
	// Define a model to express the observed empirical noise rate as a function of the magnetic field strength.
	// Compute the magnetization predicted by this model.
	// The model is motivated from mean-field theory.
	double magnet = (1 - tanh(fitparams[0] * (field - fitparams[1])))/(double) 2;
	return magnet;
}

double ModelDerivative(double *fitparams, double field, int dcomp){
	// Compute the derivative of the model, with respect to the fit parameters.
	// The model is given by: p = H(a,b,x) where H = (1 - tanh(a (x - b)))/2.
	// We have dH/da = ( - sech^2(a (x - b)) * (x - b) )/2 and dH/db = ( sech^2(a (x - b)) * a )/2
	double dvalue = 0;
	if (dcomp == 0)
		dvalue = (-1) * pow(cosh(fitparams[0] * (field - fitparams[1])), -2) * (field - fitparams[1])/(double) 2;
	else if (dcomp == 1)
		dvalue = pow(cosh(fitparams[0] * (field - fitparams[1])), -2) * fitparams[0]/(double) 2;
	else{
		// The querry is for a second derivative
		// d^2H/dada = (b - x)^2 Sech[a (-b + x)]^2 Tanh[a (-b + x)]
		// d^2H/dadb = 1/2 Sech[a (-b + x)]^2 (1 + 2 a (b - x) Tanh[a (-b + x)])
		// d^2H/dbda = 1/2 Sech[a (-b + x)]^2 (1 + 2 a (b - x) Tanh[a (-b + x)])
		// d^2H/dbdb = a^2 Sech[a (-b + x)]^2 Tanh[a (-b + x)]
		if (dcomp == 2)
			dvalue = pow((field - fitparams[1]), 2) * pow(cosh(fitparams[0] * (field - fitparams[1])), -2) * tanh(fitparams[0] * (field - fitparams[1]));
		else if((dcomp == 3) || (dcomp == 4))
			dvalue = (1/(double) 2) * pow(cosh(fitparams[0] * (field - fitparams[1])), -2) * (1 + 2 * fitparams[0] * (field - fitparams[1]) * tanh(fitparams[0] * (field - fitparams[1])));
		else if(dcomp == 5)
			dvalue = pow(fitparams[0], 2) * pow(cosh(fitparams[0] * (field - fitparams[1])), -2) * tanh(fitparams[0] * (field - fitparams[1]));
		else;
	}
	return dvalue;
}


void NonLinearFit(double *fitparams, double **trainingset){
	// Non-linear regression using Newton Gauss method
	// https://en.wikipedia.org/wiki/Gauss–Newton_algorithm
	/*
		We need to do a curve fitting, with data points (xi, yi) using the model yi = F(a,b,xi) for some non-linear function F.
		We will use Gauss Newton method to minimize a function L(a,b) that is the sum of the least sequares of the residuals of this fit model
		L(a,b) = 1/2 \sum_i [ yi - F(a,b,xi) ]^2
		a and b can be found by iterating the following update rule.
		[a_k , b_k] = [a_(k-1) , b_(k-1)] - (gradient L)/(hessian L).
		where (gradient L) = [ \sum_i [yi - F(a,b,xi)] da(yi - F(a,b,xi)) ]
							 [ \sum_i [yi - F(a,b,xi)] db(yi - F(a,b,xi)) ]
			and J is the (2x2) Jacobian matrix and
				(hessian L) ~ (J^T*J) = \sum_i [ da(yi - F(a,b,xi))*da(yi - F(a,b,xi)) ]	\sum_i [ da(yi - F(a,b,xi))*db(yi - F(a,b,xi)) ]
						  				\sum_i [ db(yi - F(a,b,xi))*da(yi - F(a,b,xi)) ]	\sum_i [ db(yi - F(a,b,xi))*db(yi - F(a,b,xi)) ]
		Note that the inverse of a (2x2) matrix [[a , b] , [c , d]] is just 1/(ad - bc) * [[d , -b] , [-c , a]].
	*/
	const int steps = 1000;
	const double tol = 10E-20;
	double determinant, gradient[2], **hessian = malloc(sizeof(double *) * 2), invhessian[2][2];
	int i, k, ci, cj;
	for (i = 0; i < 2; i ++)
		hessian[i] = malloc(sizeof(double) * 2);

	for (k = 0; k < steps; k ++){
		// Compute the two components of the gradient
		for (ci = 0; ci < 2; ci ++){
			gradient[ci] = 0;
			for (i = 1; i <= (int) trainingset[0][0]; i ++)
				gradient[ci] = gradient[ci] + (trainingset[i][1] - Model(fitparams, trainingset[i][0])) * (-1) * ModelDerivative(fitparams, trainingset[i][0], ci);
		}
		// We will compute the hessian matrix first and then compute its inverse.
		for (ci = 0; ci < 2; ci ++){
			for (cj = 0; cj < 2; cj ++){
				hessian[ci][cj] = 0;
				for (i = 1; i <= (int) trainingset[0][0]; i ++)
					hessian[ci][cj] = hessian[ci][cj] + ModelDerivative(fitparams, trainingset[i][0], ci) * ModelDerivative(fitparams, trainingset[i][0], cj) ; // + (trainingset[i][1] - Model(fitparams, trainingset[i][0])) * ModelDerivative(fitparams, trainingset[i][0], 2 + 2 * ci + cj);
			}
		}
		// Compute the inverse of the hessian matrix
		determinant = hessian[0][0] * hessian[1][1] - hessian[0][1] * hessian[1][0];
		if (fabs(determinant) < tol){
			fprintf(stderr, "\033[2m[Gauss–Newton algorithm] At step %d, the hessian matrix is singular (determinant = %.2e)!\033[0m\n", k, determinant);
			break;
		}
		// else if(IsPositiveDefinite(hessian, 2) == 0){
		// 	fprintf(stderr, "\033[2m[Gauss–Newton algorithm] At step %d, the hessian matrix is not positive definite!\033[0m\n", k);
		// 	break;
		// }
		else{
			invhessian[0][0] = hessian[1][1]/determinant;
			invhessian[0][1] = (-1) * hessian[0][1]/determinant;
			invhessian[1][0] = (-1) * hessian[1][0]/determinant;
			invhessian[1][1] = hessian[0][0]/determinant;

			// Update rule
			fitparams[0] = fitparams[0] - (invhessian[0][0] * gradient[0] + invhessian[0][1] * gradient[1]);
			fitparams[1] = fitparams[1] - (invhessian[1][0] * gradient[0] + invhessian[1][1] * gradient[1]);
		}
	}
	// compute the residue
	double residue = 0;
	for (i = 1; i <= (int) trainingset[0][0]; i ++)
		residue = residue + pow(trainingset[i][1] - Model(fitparams, trainingset[i][0]), 2);
	residue = residue/(double) 2;

	fprintf(stderr, "\033[2m[Training] After %d iteations of the Gauss–Newton algorithm, the fitting function\n\tp = (1 - tanh(%g * (x - %g)))/2\nwas obtained with an uncertainity of %g.\033[0m\n", k, fitparams[0], fitparams[1], residue);
	free(hessian);
}


double NonLinearPredictor(double *fitparams, double empirical){
	// predict the next point in a data set, assuming that the next point lies exactly on a cosh function that is obtained by fitting the earlier points in the dataset.
	// We are given a, b and p where 1 - 2 p = tanh( a (x - b)), we must determine x.
	// x is given by: x = atanh(1 - 2 p)/a + b
	double new = atanh(1 - 2 * empirical)/fitparams[0] + fitparams[1];
	fprintf(stderr, "\033[2m%.5f %.5f.\033[0m\n", new, empirical);
	return new;
}


double LinearPredictor(double ynew, double fitslope, double fitintercept){
	// predict a point in the dataset assuming that it lies on the best fit line of the dataset
	double xnew = (ynew - fitintercept)/fitslope;
	fprintf(stderr, "%.5f %.5f.\n", xnew, ynew);
	return xnew;
}


void LinearFit(double **trainingSet, double* fitparams){
	// Given the coordinates of a set of points, find the best fit line
	// Usig the best fit line and the y-coordinate of a new point, determine its corresponding x-coordinate.
	// The best fit line is given by y = a + bx
	// where b = sum((xi - xm)*(yi - ym))/sum((xi - xm)^2), where sum is over i = 1 to N and xm denotes the mean value.
	// and a = ym - b * xm.
	// So, the x-coordinate of the new point will have
	// xnew = (ynew - a)/b
	double xmean = 0, ymean = 0;
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
	fitparams[0] = correlation/variance;
	fitparams[1] = ymean - fitparams[0] * xmean;

	fprintf(stderr, "\033[2m[Training] The linear fit produced y = %g x + %g where the error was %g.\033[0m\n", fitparams[0], fitparams[1], variance);
}

void SettingsForJ(double interaction, double *fitparams, double trainRange[3]){
	// Determine the training range and the initial values of the fit paramters based on the value of J, for which the model must be trained.
	// Read settings from a file called settings.txt.
	// The settings are presented as follows. Each line stands for a setting of fit parameters and training ranges for a particular window in J.
	// <from J> <to J> <fit param 1> <fit param 2> <train range from> <train range to> <train range steps>
	
	// Default values. These values are to be assigned to the training set ranges and the fit parameters if no appropriate setting is found.
	fitparams[0] = 1;
	fitparams[1] = 1;
	trainRange[0] = -1;
	trainRange[1] = 1;
	trainRange[2] = 0.1;

	const int LINE_WIDTH = 100;
	char line[LINE_WIDTH];
	double fromJ, toJ;
	if(access("source/settings.txt", R_OK) == -1){
		fprintf(stderr, "\033[2m[Training] Settings file not found. Assigning default values to the training set and initial fit parameters.\033[0m\n");
	}
	else{
		FILE *setfp = fopen("source/settings.txt", "r");
		while (!feof(setfp)){
			fgets(line, LINE_WIDTH, setfp);
			if (line[0] != '#'){
				sscanf(line, "%lf %lf", &fromJ, &toJ);
				if ((interaction >= fromJ) && (interaction < toJ))
					sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &fromJ, &toJ, &(fitparams[0]), &(fitparams[1]), &(trainRange[0]), &(trainRange[1]), &(trainRange[2]));
			}
		}
		fclose(setfp);
	}
}

void RBIMInfer(struct RandomBondIsingModel *prbim, double **params, int *sizes, int toinfer){
	// The mean field back-calculation does not work well if J > meanfieldValidity, which in this case is 0.2.
	// In this case, we resort to making a few dry-runs and learning from the accumulated data.
	// The learning set is simply drawing a best fit line.
	const int maxtrain = 100;
	double *available = malloc(sizeof(double) * (int) (params[0][1]));
	double trainRange[3];
	double **trainingset = malloc(sizeof(double) * (1 + maxtrain));
	trainingset[0] = malloc(sizeof(double));
	int p;
	for (p = 1; p <= maxtrain; p ++)
		trainingset[p] = malloc(sizeof(double) * 2);

	// training must be re-done for every new value of J. Once training is done, we can use the data to predict all the values of B that would correspond to the desired J.
	double *fitparams = malloc(sizeof(double) * 2);
	int j, b, npoints = 0;
	for (j = 0; j < sizes[0]; j ++){
		for (p = 0; p < (int) params[0][1]; p ++)
			available[p] = params[1 + j][p];
		
		// Design the range of the training set points based on the strength of correlation.
		SettingsForJ(available[0], fitparams, trainRange);

		npoints = 1 + (int) ((trainRange[1] - trainRange[0])/trainRange[2]);
		if (npoints > maxtrain){
			fprintf(stderr, "\033[2m[Training] The training set has too many points (%d). Taking only the first %d points.\033[0m\n", npoints, maxtrain);
			npoints = maxtrain;
		}
		trainingset[0][0] = (double) npoints;
		for (p = 1; p <= npoints; p ++)
			trainingset[p][0] = trainRange[0] + (p - 1) * trainRange[2];
		RBIMDryRun(prbim, available, toinfer, trainingset);
		
		// If J is negative, fit the data with a straight line and if J is positive, fit with 1/(1+e^B) function.
		if (available[0] < 0){
			LinearFit(trainingset, fitparams);
			fprintf(stderr, "\033[2mPredictions\nB           p\033[0m\n");
			for (b = 0; b < sizes[1]; b ++)
				params[1 + b * sizes[0] + j][1] = LinearPredictor(params[1 + b * sizes[0] + j][1], fitparams[0], fitparams[1]);
		}
		else{
			NonLinearFit(fitparams, trainingset);
			fprintf(stderr, "\033[2mPredictions\nB           p\033[0m\n");
			for (b = 0; b < sizes[1]; b ++)
				params[1 + b * sizes[0] + j][1] = NonLinearPredictor(fitparams, params[1 + b * sizes[0] + j][1]);
		}
	}
	FreeDoubleArray(trainingset, 1 + maxtrain);
	free(fitparams);
}

#if 0
int main(void)
{
	/* 
		Only for code testing purposes
		Run this code independently by
		gcc rbim.c tiling.c rgens.c arrays.c notations.c load.c memory.c -o rbim.o
	*/
	
	// See the random number generator
	srand(time(NULL));

	int side;
	printf("Enter the side-length of the square lattice defining the Toric code.\n");
	printf(">>");
	scanf("%d", &side);

	char *tilingFile = malloc(sizeof(char) * 100);
	TilingFile(tilingFile, 1, side);
	struct tiling T = Load(tilingFile);
	struct tiling *pT = &T;

	// Create an Ising model.
	struct RandomBondIsingModel *pIM = malloc(sizeof(struct RandomBondIsingModel));
	InitRBIM(pIM);
	// Set the specifications of an Ising model according to a given code.
	double *interactions = malloc(sizeof(double) * 3);
	printf("Enter the ranges for J and B as <low high step>, ignoring the sign.\n");
	printf(">>J: ");
	scanf("%lf %lf %lf", &interactions[0], &interactions[1], &interactions[2]);
	double *fields = malloc(sizeof(double) * 3);
	printf(">>B: ");
	scanf("%lf %lf %lf", &fields[0], &fields[1], &fields[2]);
	double temperature = 1;
	// printf(">>T = ");
	// scanf("%f", &temperature);
	
	// Metroplis algorithm -- Cooling the configuration
	int coolSteps, montecarlo;
	printf(">>M1, M2 = ");
	scanf("%d %d", &coolSteps, &montecarlo);
	
	int si, mi, comp = 1, nrates;
	nrates = (1 + (int)((interactions[1] - interactions[0])/interactions[2])) * (1 + (int)((fields[1] - fields[0])/fields[2]));
	double interac, field;
	IntializeRBMWithCode(pIM, pT, 0, 0, 0);
	char *datafileB = malloc(sizeof(char) * 100);
	char *datafileJ = malloc(sizeof(char) * 100);
	FILE *dataB, *dataJ;
	double magnet, internal, nclusters, maxcluster, avgcluster, *params = malloc(sizeof(double) * 3);
	for (interac = interactions[1]; interac >= interactions[0]; interac -= interactions[2]){
		sprintf(datafileJ, "./../test/data_J_%g.txt", interac);
		dataJ = fopen(datafileJ, "w");
		fclose(dataJ);
	}
	params[2] = temperature;
	fprintf(stderr, "\ncoordination number = %d\n", pIM->coordination);
	for (field = fields[1]; field >= fields[0]; field -= fields[2]){
		sprintf(datafileB, "./../test/data_B_%g.txt", field);
		dataB = fopen(datafileB, "w");
		// fprintf(data, "J,B,<m>,C,<u>\n");
		for (interac = interactions[1]; interac >= interactions[0]; interac -= interactions[2]){
			UpdateRBIM(pIM, interac, field, temperature);
			Cool(pIM, 1000);
			magnet = 0;
			internal = 0;
			avgcluster = 0;
			for (mi = 0; mi < montecarlo; mi ++){
				Cool(pIM, coolSteps);
				PhysicalProperties(pIM);
				magnet += pIM->magnetization;
				internal += pIM->internalEnergy;
				avgcluster += pIM->avgClusterSize;
			}
			avgcluster = avgcluster/((double) montecarlo);
			magnet = magnet/((double) montecarlo);
			internal = internal/((double) montecarlo);
			params[0] = interac;
			params[1] = (1 - magnet)/(double) 2;
			fprintf(stderr, "%g  %g  %g  %g  %g  %g\n", interac, field, magnet, (1 - magnet)/(double) 2, pIM->avgClusterSize, pIM->internalEnergy);

			fprintf(dataB, "%g  %g  %g  %g  %g %g\n", interac, field, pIM->magnetization, (1 - magnet)/(double) 2, pIM->avgClusterSize, pIM->internalEnergy);
			
			sprintf(datafileJ, "./../test/data_J_%g.txt", interac);
			dataJ = fopen(datafileJ, "a");
			fprintf(dataJ, "%g  %g  %g  %g  %g %g\n", interac, field, pIM->magnetization, (1 - magnet)/(double) 2, pIM->avgClusterSize, pIM->internalEnergy);
			fclose(dataJ);
			
			// fprintf(stderr, "\r%d%% done", (int)(100 * comp/(double) nrates));
			comp ++;
		}
		fclose(dataB);
	}
	fprintf(stderr, ".\n");
	free(datafileB);
	free(datafileJ);
	free(params);
	FreeRBIM(pIM);
	// Free the structure
	FreeTiling(pT);
	return 0;
}
#endif
