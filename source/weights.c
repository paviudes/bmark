/*
 * weights.c
 *
 *  Created on: Aug 22, 2015
 *      Author: ndelfosse
 */
#include <stdlib.h>
#include "tiling.h"


//
int MaxWeight(struct tiling G, struct tiling H)
{
	int i;
	int maxWeight = 0;
	for (i=0; i<G.f; i++)
		if (G.F[i][0] > maxWeight)
			maxWeight = G.F[i][0];
	for (i=0; i<H.f; i++)
		if (H.F[i][0] > maxWeight)
			maxWeight = H.F[i][0];
	return maxWeight;
}


//
int *FaceMeasurementWeights(struct tiling G, int maxWeight)
{
	int i, j, w;
	int *Weight = malloc((maxWeight+1)*sizeof(int));
	for (i=0; i<maxWeight+1; i++)
		Weight[i] = 0;
	for (i=0; i<G.f; i++)
	{
		w = 0;
		for (j=1; j<=G.F[i][0]; j++)
			if (G.doE[G.F[i][j]] == 0)
				w++;
		Weight[w]++;
	}
	return Weight;
}


