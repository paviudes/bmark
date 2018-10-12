#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stddef.h>
#include "memory.h"
#include "tiling.h"


//=======================================
//		Initialize tiling
//=======================================


//initialize length, v, e and f to 0
//initialize all pionters to NULL
//This function does not free the pointers in the structure
//Do not apply it to an allocated pointer.
void InitTiling(struct tiling *G)
{
	G->length = 0;
	G->type = NULL;
	G->v = 0;
	G->e = 0;
	G->f = 0;
	G->E = NULL;
	G->F = NULL;
	G->VToE = NULL;
	G->VToN = NULL;
	G->dV = NULL;
	G->dE = NULL;
	G->doE = NULL;
	G->doV = NULL;
	G->Coord = NULL;
}


void InitStorage(struct storePerformance *perf)
{
	perf->tilingName = NULL;
	perf->tilingName = NULL;
	perf->noiseRange = NULL;
	perf->n = 0;
	perf->k = 0;
	perf->maxWeight = 0;
	perf->ZWeight = NULL;
	perf->XWeight = NULL;
}



//=======================================
//		Free memory
//=======================================


//
void FreeTiling(struct tiling *pG)
{
	if (pG->type != NULL)
		free(pG->type);
	if (pG->E != NULL)
		FreeIntArray(pG->E, pG->e);
	if (pG->F != NULL)		
		FreeIntArray(pG->F, pG->f);
	if (pG->VToE != NULL)
		FreeIntArray(pG->VToE, pG->v);
	if (pG->VToN != NULL)
		FreeIntArray(pG->VToN, pG->v);
	if (pG->dV != NULL)
		free(pG->dV);
	if (pG->dE != NULL)
		free(pG->dE);
	if (pG->doV != NULL)
		free(pG->doV);
	if (pG->doE != NULL)
		free(pG->doE);
	if (pG->Coord != NULL)
		FreeFloatArray(pG->Coord, pG->v);
	InitTiling(pG);
}


void FreeStorage(struct storePerformance *perf)
{
	if (perf->tilingName != NULL)
		free(perf->tilingName);
	if (perf->dataFile != NULL)
		free(perf->dataFile);
	if (perf->noiseRange != NULL)
		free(perf->noiseRange);
	if (perf->ZWeight != NULL)
		free(perf->ZWeight);
	if (perf->XWeight != NULL)
		free(perf->XWeight);
	InitStorage(perf);
}

//=======================================
//		Print the tiling in the terminal
//=======================================


void PrintVertices(struct tiling G)
{
	int i;
	printf("     VERTICES\n");
	for (i=0; i<G.v; i++)
	{
		printf("v%d ", i);
		if (G.Coord != NULL)
			printf("= (%f, %f) ", G.Coord[i][0], G.Coord[i][1]);
		if (G.doV[i] == 1) 
			printf("open boundary");
		if (G.dV[i] == 1 && G.doV[i] == 0) 
			printf("closed boundary");
		printf("\n");
	}
	printf("\n");
}


void PrintEdges(int **E, int *dE, int *doE, int e)
{
	int i;
	printf("     EDGES\n");
	for (i=0; i<e; i++)
	{
		printf("e%d = (v%d,v%d) ", i, E[i][0], E[i][1]);
		if (doE[i] == 1) printf("open boundary");
		if (dE[i] == 1 && doE[i] == 0) printf("closed boundary");
		printf("\n");
	}
	printf("\n");
}


void PrintFaces(int **F, int f)
{
	int i, j;
	printf("     FACES\n");
	for (i=0; i<f; i++)
	{
		printf("face f%d of degree %d = { ", i, F[i][0]);
		for (j=1; j<F[i][0]+1; j++)
			printf("e%d ", F[i][j]);
		printf("}\n");
	}
	printf("\n");
}


void PrintVToE(int **VToE, int v)
{
	int i, j;
	printf("     VERTEX/EDGE INCIDENCE\n");
	for (i=0; i<v; i++)
	{
		printf("vertex v%d of degree %d is incident to { ", i, VToE[i][0]);
		for (j=1; j<VToE[i][0]+1; j++)
			printf("e%d ", VToE[i][j]);
		printf("}\n");
	}
	printf("\n");
}

/*
void PrintVToN(int **VToN, int v)
{
	int i, j;
	printf("     VERTEX/VERTEX INCIDENCE\n");

	for (i=0; i<v; i++)
	{
		printf("v%d of degre %d is neighbor to { ", i, VToN[i][0]);
		for (j=1; j<VToN[i][0]+1; j++)
			printf("v%d ", VToN[i][j]);
		printf("}\n");
	}
	printf("\n");
}
*/


void PrintTiling(struct tiling G)
{
	if (G.type == NULL)
		printf("No Tiling\n");
	else
	{
		printf("TILING PRINTING:\n");
		printf(">> type: %s\n", G.type);
		printf(">> size: v=%d, e=%d, f=%d\n", G.v, G.e, G.f);
		printf(">> vertices:\n");
		PrintVertices(G);
		printf(">> edges:\n");
		PrintEdges(G.E, G.dE, G.doE, G.e);
		printf(">> faces:\n");
		PrintFaces(G.F, G.f);
	}
}


//=======================================
//		Update a tiling
//=======================================



//return the array VToE that encodes the incident relations from vertices to edges
//VToE[i] represents the set of edges incident to the vertex i, for i=0,1,...,v-1
//VToE[i][0] is the degree of the vertex i
//VToE[i][j] is the index of the j-th edge indicent to the vertex i, for j=1,2,...,VToE[i][0].
int **VToEAllocation(int **E, int e, int v)
{
	int i, x, y;

	//computation of the degree max of the graph
	int *degree = malloc(v*sizeof(int));
	for (i=0; i<v; i++)
		degree[i] = 0;
	for (i=0; i<e; i++)
	{
		x = E[i][0];
		y = E[i][1];
		degree[x] = degree[x]+1;
		degree[y] = degree[y]+1;
	}
	int degreeMax = 0;
	for (i=0; i<v; i++)
	{
		if (degree[i] > degreeMax)
			degreeMax = degree[i];
	}
	free(degree);
	
	//construction of VToE
	int **VToE = malloc(v*sizeof(int *));
	for (i=0; i<v; i++)
	{
		VToE[i] = malloc((degreeMax+1)*sizeof(int));
		VToE[i][0] = 0;
	}
	for (i=0; i<e; i++)
	{
		x = E[i][0];
		y = E[i][1];
		VToE[x][0] = VToE[x][0]+1;
		VToE[x][VToE[x][0]] =  i;
		VToE[y][0] = VToE[y][0]+1;
		VToE[y][VToE[y][0]] =  i;
	}
	return VToE;
}

/*
//return the array VToN that encodes the incidence relation between vertices.
//VToN[i] represents the set of vertices incident to the vertex i, for i=0,1,...,v-1.
//VToN[i][0] is the degree of the vertex i.
//VToN[i][j] is the index of the j-th vertex indicent to the vertex i, for j=1,2,...,VToN[i][0].
int **VToNAllocation(int **E, int e, int v)
{
	int i, x, y;
	
	int *degree = malloc(v*sizeof(int));
	for (i=0; i<v; i++)
		degree[i] = 0;
	for (i=0; i<e; i++)
	{
		x = E[i][0];
		y = E[i][1];
		degree[x] = degree[x]+1;
		degree[y] = degree[y]+1;
	}
	int degreeMax = 0;
	for (i=0; i<v; i++)
	{
		if (degree[i] > degreeMax)
			degreeMax = degree[i];
	}
	free(degree);
	
	int **VToN = malloc(v*sizeof(int *));
	for (i=0; i<v; i++)
	{
		VToN[i] = malloc((degreeMax+1)*sizeof(int));
		VToN[i][0] = 0;
	}
	
	for (i=0; i<e; i++)
	{
		x = E[i][0];
		y = E[i][1];
		VToN[x][0]++;
		VToN[x][VToN[x][0]] = y;
		VToN[y][0]++;
		VToN[y][VToN[y][0]] = x;
	}
	
	return VToN;
}
*/


//we start by finding the edges that belong to single face
int *dEAllocation(int **F, int f, int e)
{
	int i, j;
	
	int *dE = malloc(e*sizeof(int));
	for (i=0; i<e; i++)
		dE[i] = 2;
	for (i=0; i<f; i++)
		for (j=1; j<=F[i][0]; j++)
		{
			dE[F[i][j]]--;
		}
	
	return dE;
}

//from the boundary edges, the boundary vertices are their endpoints
int *dVAllocation(int **E, int *dE, int e, int v)
{
	int i;
	
	int *dV = malloc(v*sizeof(int));	
	for (i=0; i<v; i++)
		dV[i] = 0;
	for (i=0; i<e; i++)
		if (dE[i] == 1)
		{
			dV[E[i][0]] = 1;
			dV[E[i][1]] = 1;		
		}
	return dV;
}



//Once v, e, f, E, F and doE are defined, use it to update the tiling
//update dE
//update dvE
//update dovE
//update VToE
void UpdateTiling(struct tiling *G)
{
	//update boundaries
	if ((*G).dE != NULL) free((*G).dE);
	(*G).dE = dEAllocation((*G).F, (*G).f, (*G).e);
	if ((*G).dV != NULL) free((*G).dV);
	(*G).dV = dVAllocation((*G).E, (*G).dE, (*G).e, (*G).v);

	//update open vertices
	int i;
	int *doV = malloc((*G).v*sizeof(int));
	for (i=0; i<(*G).v; i++)
		doV[i]=0;
	
	for (i=0; i<(*G).e; i++)
 		if ((*G).doE[i] == 1) 
		{
			doV[(*G).E[i][0]] = 1;
			doV[(*G).E[i][1]] = 1;	
		}
	if ((*G).doV != NULL) free((*G).doV);
	(*G).doV = doV;
	
	//update the incidence relations
	if ((*G).VToE != NULL) FreeIntArray((*G).VToE, (*G).v);
	(*G).VToE = VToEAllocation((*G).E, (*G).e, (*G).v);
}





//=======================================
//		Basic functions on tilings
//=======================================


//takes a set of k edges forming a cycle and return the
//corresponding k consecutive vertices
//c[0] is the length of the cycle
int *FaceToVertexList(struct tiling G, int *c)
{
	int i, j;
	int *v = malloc((c[0]+1)*sizeof(int));
	v[0] = c[0];
	v[1] = G.E[c[1]][0];
	v[2] = G.E[c[1]][1];
	
	int endpoint0;
	int endpoint1;
	
	i=2;
	while (i<c[0])
	{
		j=2;
		while (j<=c[0])
		{
			//non backtracking walk
			endpoint0 = G.E[c[j]][0];
			endpoint1 = G.E[c[j]][1];
			if ((endpoint0 == v[i]) && (endpoint1 != v[i-1]))
			{
				v[i+1] = endpoint1;
				j=c[0]+1;
			}
			if ((endpoint1 == v[i]) && (endpoint0 != v[i-1]))
			{
				v[i+1] = endpoint0;
				j=c[0]+1;
			}
			j++;
		}
		i++;
	}
	
	return v;
}





