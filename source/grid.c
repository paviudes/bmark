#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tiling.h"


//=======================================
//		Construction of a square grid of size lx x ly
//=======================================


//return the coordinates of the vertex v
int GridVToCoordX(int lx, int ly, int v)
{
	return v%(lx+1);
}

//return the coordinates of the vertex v
int GridVToCoordY(int lx, int ly, int v)
{
	return v/(lx+1);
}


// construct the array of coordinates of the planar square tiling
float **GridCoordAllocation(int v, int lx, int ly)
{
	int i;
	float **Coord = malloc(v*sizeof(float *));
	for (i=0; i<v; i++) Coord[i] = malloc(2*sizeof(float));
	
	for (i=0; i<v; i++)
	{
		Coord[i][0] = GridVToCoordX(lx, ly, i);
		Coord[i][1] = GridVToCoordY(lx, ly, i);
	}
	return Coord;
}



//Fill and return the tiling structure of the square lattice of the torus
struct tiling Grid(int lx, int ly)
{
	struct tiling G;
	InitTiling(&G);
	G.length = lx*ly;

	char *type = malloc(50*sizeof(char));
	sprintf(type, "grid_%dx%d", lx, ly);
	G.type = type;

	G.v = (lx+1)*(ly+1);
	G.e = 0;
	G.f = 0;

	//vertex coordinates
	G.Coord = GridCoordAllocation(G.v, lx, ly);

	//update incidence, boundaries, open vertices
	UpdateTiling(&G);

	return G;
}


//=======================================
//		Add a face to a tiling
//=======================================


//add the face newface to the tiling
//this face is given in the usual format of faces
//newFace[0] is the length the face
//The indices of its egdes are the newFace[i] for i=1,..,newFace[0].
void AddFaceFromEdges(struct tiling *G, int *newFace)
{
	int i, j;
	(*G).f++;

	//allocate the new face array
	int maxlength = 0;
	for (i=0; i<(*G).f-1; i++)
		if ((*G).F[i][0] > maxlength)
			maxlength = (*G).F[i][0];
	if (newFace[0] > maxlength)
		maxlength = newFace[0];
			
	int **nextF = malloc(((*G).f)*sizeof(int *));
	for (i=0; i<(*G).f; i++)
		nextF[i] = malloc((maxlength+1)*sizeof(int));

	//copy (*G).F in the new array	
	for (i=0; i<(*G).f-1; i++)
		for (j=0; j<=(*G).F[i][0]; j++)
			nextF[i][j] = (*G).F[i][j];
			
	for (j=0; j<=newFace[0]; j++)
		nextF[(*G).f-1][j] = newFace[j];
	
	free((*G).F);
	(*G).F = nextF;

	UpdateTiling(G);
}



//return the index of the edge added
int AddEdgeFromVertices(struct tiling *G, int *newEdge)
{	
	int i;
	int indexNewEdge = -1;
	
	//test if this edge already exists
	for (i=0; i<G->e-1; i++)
		if (((newEdge[0] == G->E[i][0]) && (newEdge[1] == G->E[i][1])) || ((newEdge[0] == G->E[i][1]) && (newEdge[1] == G->E[i][0])))
			indexNewEdge = i;
	
	if (indexNewEdge == -1)
	{
		indexNewEdge = G->e++;
		
		//update E
		int **nextE = malloc((G->e)*sizeof(int *));
		for (i=0; i<G->e; i++)
			nextE[i] = malloc(2*sizeof(int));

		for (i=0; i<G->e-1; i++)
			{
				nextE[i][0] = G->E[i][0];
				nextE[i][1] = G->E[i][1];
			}
	
		nextE[G->e-1][0] = newEdge[0];
		nextE[G->e-1][1] = newEdge[1];
		
		free(G->E);
		G->E = nextE;
	
		//update doE
		int *nextdoE = malloc(G->e*sizeof(int));
	
		for (i=0; i<G->e-1; i++)
			nextdoE[i] = G->doE[i];
		nextdoE[G->e-1] = 0;
	
		free(G->doE);
		G->doE = nextdoE;
	
		UpdateTiling(G);
	}
	return indexNewEdge;
}



//We combine the two previous functions to add a face from
//a vertex list
void AddFace(struct tiling *G, int *vertexList, int length)
{
	int i;
	
	int *newEdge = malloc(2*sizeof(int));
	int *newFace = malloc((length+1)*sizeof(int));
	newFace[0] = length;

	for (i=0; i<length; i++)
	{
		newEdge[0] = vertexList[i];
		newEdge[1] = vertexList[(i+1)%length];
		newFace[i+1] = AddEdgeFromVertices(G, newEdge);
	}
	
	AddFaceFromEdges(G, newFace);
}




//=======================================
//		Remove isolated vertices
//=======================================


//return the array nextIndV 
//nextIndV[i] gives the index of the vertex i after removing isolated vertices
//nextIndV[i] = -1 if vi is isolated
//update Coord
//update G.v
int *NonIsolatedIndex(struct tiling *G)
{
	int i;

	//determine the correspondance between indices of vertices 
	//before and after remving isolated vertices.
	//nextIndV[i] = -1 if vi is isolated.
	int *nextIndV = malloc((*G).v*sizeof(int));
	for (i=0; i<(*G).v; i++) nextIndV[i] = -1;
	
	//keep only vertices of degree>0
	int nextv = 0;
	for (i=0; i<(*G).v; i++)
		if ((*G).VToE[i][0] > 0) nextIndV[i] = nextv++;
		
	//update vertices coordinates
	float **nextCoord = malloc(nextv*sizeof(float *));
	for (i=0; i<nextv; i++)
		nextCoord[i] = malloc(2*sizeof(float));
	
	for (i=0; i<(*G).v; i++)
	{
		if (nextIndV[i] != -1)
		{
			nextCoord[nextIndV[i]][0] = (*G).Coord[i][0];
			nextCoord[nextIndV[i]][1] = (*G).Coord[i][1];
		}
	}
	
	for (i=0; i<(*G).v; i++) free((*G).Coord[i]);
	free((*G).Coord);
	(*G).Coord = nextCoord;
	(*G).v = nextv;
	
	return nextIndV;
}



void RemoveIsolatedVertices(struct tiling *G)
{
	int i;
	
	//nextIndV[i] is the index of vi after removing isolated vertices
	//nextIndV[i] = -1 if vi is isolated
	//G.v is updated
	//G.Coord is updated
	int *nextIndV = NonIsolatedIndex(G);

	//update G.E
	for (i=0; i<(*G).e; i++)
		{
			(*G).E[i][0] = nextIndV[(*G).E[i][0]];
			(*G).E[i][1] = nextIndV[(*G).E[i][1]];
		}
	
	//G.F is unchanged
	//G.doE is unchanged
		
	UpdateTiling(G);	
}








