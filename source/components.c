#include <stdlib.h>
#include <stdio.h>
#include "tiling.h"


void Explore(struct tiling G, int v, int *visited)
{
	//mark the visited vertex
	visited[v] = 1;

	int i, ei, endpoint1, endpoint2;
	//recursive call
	for (i=1; i<=G.VToE[v][0]; i++)
	{
		ei = G.VToE[v][i];
		endpoint1 = G.E[ei][0];
		endpoint2 = G.E[ei][1];
		if (v == endpoint1)
		{
			if (!visited[endpoint2])
				Explore(G, endpoint2, visited);
		}
		else
		{
			if (!visited[endpoint1])
				Explore(G, endpoint1, visited);
		}
	}
}


int Components(struct tiling G)
{
	int i;
	//init visited
	int *visited = malloc(G.v*sizeof(int));
	for (i=0; i<G.v; i++) visited[i] = 0;

	int ncomp = 0;
	for (i=0; i<G.v; i++)
	{
		if (!visited[i])
		{
			Explore(G, i, visited);
			ncomp++;
		}
		// i++;
	}
		
	return ncomp;
}


int NoOpenVerticesComponents(struct tiling G)
{
	int i;
	//init visited
	int *visited = malloc(G.v*sizeof(int));
	for (i=0; i<G.v; i++) visited[i] = 0;

	//mark all components with an open vertex as visited
	for (i=0; i<G.v; i++)
		if (G.doV[i] == 1)
			Explore(G, i, visited);
			
	int ncomp = 0;
	for (i=0; i<G.v; i++)
	{
		if (!visited[i])
		{
			Explore(G, i, visited);
			ncomp++;
		}
		// i++;
	}
		
	return ncomp;
}


int NoClosedEdgesComponents(struct tiling G)
{
	int i;
	//init visited
	int *visited = malloc(G.v*sizeof(int));
	for (i=0; i<G.v; i++) visited[i] = 0;

	//mark all components with an closed edge as visited
	for (i=0; i<G.e; i++)
		if ((G.dE[i] == 1) && (G.doE[i] == 0))
			Explore(G, G.E[i][0], visited);
			
	int ncomp = 0;
	for (i=0; i<G.v; i++)
	{
		if (!visited[i])
		{
			Explore(G, i, visited);
			ncomp++;
		}
		// i++;
	}
	
	return ncomp;
}


//compute the dimension of the first homology group of the tiling
//This is also the number of logical qubits
//It is a direct appplication of Prop 1.6:
int HomDim(struct tiling G)
{
	//compute the nb of non-open edges and non-open vertices
	int i, nov = 0;
	for (i=0; i<G.v; i++)
		if (G.doV[i] != 1) nov++;

	int noe = 0;
	for (i=0; i<G.e; i++)
		if (G.doE[i] != 1) noe++;

	return (- nov + noe - G.f + NoOpenVerticesComponents(G) + NoClosedEdgesComponents(G));
}


