#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "tiling.h"
#include "simulation.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "RandomNumberGenerators.h"
#include "arrays.h"
#include "holes.h"
#include "draw.h"
#include "dual.h"
#include "components.h"
#include "weights.h"
#include "cell.h"
#include "grid.h"
#include "Save.h"
#include "Load.h"
#include "regulartilings.h"
#include "homological.h"
#include "correlations.h"
#include "Performance.h"

int main()
{
	/*
	// Testing the generation of correlated erasure

	struct tiling G = SquareTiling(5, 5);
	VToNAllocation(&G);

	// printf("Print the tiling in the terminal:\n");
	// PrintTiling(G);
	// Draw(G, "indexed", 'b');

	// breakpoint
	// printf(">>Press any key to continue");
	// scanf("%*d");

	double loss = 0.9;
	struct simulation S;
	S.erasure = malloc(sizeof(int) * G.e);
	FillIntArrayWithConstants(S.erasure, G.e, 0);
	S.randomGenerator = 0;
	SeedGenerator(S.randomGenerator);
	
	GenerateErasure(&G, loss, &S);
	// Highlight erased edges
	HighlightEdges(G, S.erasure, 'b');
	// Ripple error
	int rei;
	printf(">>Correlations around:");
	scanf("%d", &rei);

	S.correlationLength = 5;
	// ErrorBall(rei, &G, &S);
	ErrorCluster(rei, &G, &S);
	printf("done.\n");

	// Highlight erased edges
	HighlightEdges(G, S.erasure, 'b');

	// breakpoint
	// printf(">>Press any key to continue");
	// scanf("%*d");
	*/

	
	// testing the random number generator function
	// gsl_rng_default_seed = time(NULL);
	unsigned int randNum;
	int ri = 0, nRandNums = 100, mean = 10;
	printf("\033[93mRandom numbers:\n");
	for (ri = 0; ri < nRandNums; ri ++){
		randNum = PoissonRandomNumber(mean);
		printf(" %d", randNum);
	}
	printf(".\033[0m\n");
	
	/*float correlationLength;
	int randomGenerator = 0, toricLength;
	char* latticeName;

	printf(">>Toric code linear dimension: ");
	scanf("%d", &(toricLength));

	latticeName = malloc(sizeof(char) * 200);
	sprintf(latticeName, "%d_toric_square", toricLength);
	struct tiling G = Load(latticeName);
	// Count the number of qubits
	int nqubits = 0, ei;
	for(ei = 0; ei < G.e; ei ++)
		if(G.doE[ei] == 0)
			nqubits ++;

	// Testing Performance
	double *noiseRange = malloc(sizeof(double) * 3);
	printf("Analyze the Performance\n");
	printf(">>Mean correlation length: ");
	scanf("%f", &(correlationLength));
	printf(">>Lower limit: ");
	scanf("%lf", &(noiseRange[0]));
	printf(">>Upper limit: ");
	scanf("%lf", &(noiseRange[1]));
	printf(">>Step size: ");
	scanf("%lf", &(noiseRange[2]));
	long nStats = 0;
	printf(">>number of decoding trials = ");
	scanf("%ld", &nStats);

	// printf("Create a hole from 11 to 7.\n");
	// CreateRectangleHole(&G, 11, 7);
	clock_t startTime = clock();
	char* perfResultFile = malloc(sizeof(char) * 100);
	printf("Launching performance.\n");
	perfResultFile = Performance(G, noiseRange, nStats, correlationLength, randomGenerator);
	clock_t endTime = clock();
	double runtime = (double)(endTime - startTime)/CLOCKS_PER_SEC;
	printf("\033[92m _/ Done %d qubits with correlation length %.2f in %.3f seconds.\n\tOutput in %s.\033[0m\n", nqubits, correlationLength, runtime, perfResultFile);
	
	// Free memory
	free(latticeName);
	free(noiseRange);
	return 0;
	*/
}
