#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "memory.h"
#include "tiling.h"
#include "rgens.h"
#include "dual.h"
#include "arrays.h"
#include "save.h"
#include "simulation.h"
#include "homological.h"


void VToNAllocation(struct tiling* pG){
	// Compute the vertex edge incidence table of a Tiling, from the vertex edge incidence table
	if((pG->VToN) != NULL){
		int vi = 0;
		for(vi = 0; vi < pG->v; vi ++)
			free(pG->VToN[vi]);

		free(pG->VToN);
		pG->VToN = NULL;
	}

	(pG->VToN) = malloc(sizeof(int*) * (pG->v));
	int vi = 0, ei = 0, endPoint = 0, outEdge = 0;
	for(vi = 0; vi < pG->v; vi ++){
		(pG->VToN)[vi] = malloc(sizeof(int) * (1 + (pG->VToE)[vi][0]));
		(pG->VToN)[vi][0] = (pG->VToE)[vi][0];
		for(ei = 1; ei <= (pG->VToE)[vi][0]; ei ++){
			outEdge = (pG->VToE)[vi][ei];
			if(vi == ((pG->E)[outEdge][0]))
				endPoint = (pG->E)[outEdge][1];
			else
				endPoint = (pG->E)[outEdge][0];

			(pG->VToN)[vi][ei] = endPoint;
		}
	}
}

int NonOpenVertices(struct tiling* pG){
	// Count the number of non-open vertices of a Tiling.
	int vi = 0, noV = 0;
	for(vi = 0; vi < pG->v; vi ++)
		noV += (pG->doV)[vi];
	return (pG->v - noV);
}

double ErrorBar(double pfail, long nSamps){
	return (1.96 * sqrt(pfail * (1 - pfail)/((double) nSamps)));
}

/********************************** Prepare the simulation data structure **********************************/

void InitializeSimulation(struct simulation* pSim, struct tiling **pGs, int nRates, long stats){
	// Preprocessing for the round of Montecarlo simulations to estimate decoding failure probability
	int ti, ei;
	int *edgeOrd;
	pSim->nOV = malloc(sizeof(int) * 2);
	pSim->comps = malloc(sizeof(int) * 2);
	pSim->edgeOrderings = malloc(sizeof(int *) * 3);
	for (ti = 0; ti < 3; ti ++){
		VToNAllocation(pGs[ti]); // Initializing the vertex-vertex incidence matrix for the Tiling, Dual Tiling and the double Dual Tiling
		if (ti < 2){
			(pSim->nOV)[ti] = NonOpenVertices(pGs[ti]);
			(pSim->comps)[ti] = ConnectedComponents(pGs[ti], NULL, NULL, 1);
		}
		(pSim->edgeOrderings)[ti] = malloc(sizeof(int) * pGs[ti]->e);
		FillIntArrayWithConstants((pSim->edgeOrderings)[ti], pGs[ti]->e, -1);
		if (ti == 0){
			edgeOrd = EToDualE(pGs[ti][0]);
			for (ei = 0; ei < pGs[ti]->e; ei ++)
				(pSim->edgeOrderings)[ti][ei] = edgeOrd[ei];
		}
		else{
			edgeOrd = EToDualE(pGs[ti - 1][0]);
			for (ei = 0; ei < pGs[ti - 1]->e; ei ++)
				if (edgeOrd[ei] != -1)
					(pSim->edgeOrderings)[ti][edgeOrd[ei]] = ei;
		}

	}
	pSim->accuracy = 0.05;
	pSim->minfails = 1000;
	pSim->totaltrials = 0;

	pSim->stats = stats;
	pSim->type = 0;
	
	pSim->erasure = malloc(sizeof(int *) * 3);
	(pSim->erasure)[0] = malloc(sizeof(int) * 2);
	(pSim->erasure)[0][0] = 0;
	(pSim->erasure)[0][1] = 0;
	for (ti = 0; ti < 2; ti ++)
		(pSim->erasure)[ti + 1] = malloc(sizeof(int) * pGs[ti]->e);

	// When running on a cluster, timestamp is not a good variable since there are parallel processes. So we include the time stamp with the node ID.
	char* jobid = malloc(sizeof(char) * 100);
	FILE *jobp = popen("echo $PBS_JOBID", "r");
	if (jobp == NULL)
		sprintf(jobid, "ju_%lf", Random("uniform", 0, NULL));
	else
		fgets(jobid, 100, jobp);
	fclose(jobp);
	jobid[strcspn(jobid, "\n")] = 0;
	if (strlen(jobid) == 0)
		sprintf(jobid, "ju_%lf", Random("uniform", 0, NULL));
	time_t tme = time(NULL);
	struct tm tm = *localtime(&tme);
	pSim->timestamp = malloc(sizeof(char) * 400);
	sprintf(pSim->timestamp, "%d_%d_%d_%d_%d_%d_%s", tm.tm_mday, tm.tm_mon, tm.tm_year + 1900, tm.tm_hour, tm.tm_min, tm.tm_sec, jobid);
	free(jobid);
}


void FreeSimulation(struct simulation* pSim, int nRates){
	// Free memory allocated for simulation
	free(pSim->nOV);
	free(pSim->comps);
	FreeIntArray((pSim->edgeOrderings), 3);
	FreeIntArray((pSim->erasure), 3);
	free(pSim->timestamp);
}
