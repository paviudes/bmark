#include <stdlib.h>
#include <stdio.h>
#include "tiling.h"
#include "memory.h"


//=======================================
//		Open some boundaries
//=======================================


//OpenE is the incidator vector of a set of edges
//change the edges ek such that OpenE[k]=1 to open
//When OpenE contains an edge that is not a boundary nothing happens
//to this edge. Only boundaries can be open
void OpenBoundary(struct tiling *G, int *OpenE)
{
	int i;

	//open edges
	int *doE = malloc((*G).e*sizeof(int));
	for (i=0; i<(*G).e; i++)
		doE[i] = (*G).doE[i];

	for (i=0; i<(*G).e; i++)
	{
		if ((OpenE[i] == 1) && ((*G).dE[i] == 1))
			doE[i] = 1;
	}
	free((*G).doE);
	(*G).doE = doE;
	
	//open vertices	
	int *doV = malloc((*G).v*sizeof(int));
	for (i=0; i<(*G).v; i++)
		doV[i]=0;
	
	for (i=0; i<(*G).e; i++)
 		if ((*G).doE[i] == 1) 
		{
			doV[(*G).E[i][0]] = 1;
			doV[(*G).E[i][1]] = 1;	
		}
	free((*G).doV);
	(*G).doV = doV;
}


void CloseBoundary(struct tiling *G, int *CloseE)
{
	int i;

	//open edges
	int *doE = malloc((*G).e*sizeof(int));
	for (i=0; i<(*G).e; i++)
		doE[i] = (*G).doE[i];

	for (i=0; i<(*G).e; i++)
	{
		if (CloseE[i] == 1) doE[i] = 0;
	}
	free((*G).doE);
	(*G).doE = doE;

	//open vertices
	int *doV = malloc((*G).v*sizeof(int));
	for (i=0; i<(*G).v; i++)
		doV[i]=0;
	
	for (i=0; i<(*G).e; i++)
 		if ((*G).doE[i] == 1) 
		{
			doV[(*G).E[i][0]] = 1;
			doV[(*G).E[i][1]] = 1;	
		}
	free((*G).doV);
	(*G).doV = doV;
}



//=======================================
//		Creation of closed holes
//=======================================


//return the indicator vector  HoleE of the edges inside the hole
int *HoleEAllocation(struct tiling *G, int *HoleF)
{
	int i, j;
	int *HoleE = malloc((*G).e*sizeof(int));
	for (i=0; i<(*G).e; i++) HoleE[i] = 0;
	
	//find the edges that lives inside the hole by computing
	//dualDegE[k]== nb of faces incident to ek after hole
	int *degreeEToF = malloc((*G).e*sizeof(int));
	for (i=0; i<(*G).e; i++) degreeEToF[i] = 0;
	
	for (i=0; i<(*G).f; i++)
		if (!HoleF[i])
			for (j=1; j<=(*G).F[i][0]; j++)
				degreeEToF[(*G).F[i][j]]++;
		
	for (i=0; i<(*G).e; i++)
		if (degreeEToF[i] == 0) HoleE[i] = 1;		
	
	free(degreeEToF);
	
	return HoleE; 
}


//return the array nextIndV 
//nextIndV[i] gives the index of the vertex i after removing the hole
//nextIndV[i] = -1 if vi lives inside the hole
//update Coord
//update G.v
int *RemoveHolesVertices(struct tiling *G, int *HoleE)
{

	int i, j;
	//determine the correspondance between indices of vertices 
	//before and after hole.
	//nextIndV[i] = -1 if vi lives inside the hole.
	int *nextIndV = malloc((*G).v*sizeof(int));
	for (i=0; i<(*G).v; i++) nextIndV[i] = -1;
	
	int deg;
	int nextv = 0;
	for (i=0; i<(*G).v; i++)
	{
		//compute the degree of vi after the hole
		deg = 0;
		for (j=1; j<=(*G).VToE[i][0]; j++)
			if (!HoleE[(*G).VToE[i][j]])
				deg++;
		//keep this vertex if is does not become isolated
		if (deg > 0) nextIndV[i] = nextv++;
	}
		
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


//return the array nextIndE
//nextIndE[i] gives the index of the edge i after removing the hole
//nextIndE[i] = -1 if ei lives inside the hole
//update G.E
//update G.e
int *RemoveHolesEdges(struct tiling *G, int *HoleE, int *nextIndV)
{
	int i;
	//contruct nextIndE and compute nexte
	int *nextIndE = malloc((*G).e*sizeof(int));
	for (i=0; i<(*G).e; i++) nextIndE[i] = -1;

	int nexte = 0;
	for (i=0; i<(*G).e; i++)
		if (!HoleE[i]) nextIndE[i] = nexte++;

	//remove the hole's edges and update edges endpoints and 
	//open edges
	int **nextE = malloc(nexte*sizeof(int *));
	for (i=0; i<nexte; i++)
		nextE[i] = malloc(2*sizeof(int));
	int *nextdoE = malloc(nexte*sizeof(int));

	int k=0;
	for (i=0; i<(*G).e; i++)
		if (!HoleE[i])
		{
			nextE[k][0] = nextIndV[(*G).E[i][0]];
			nextE[k][1] = nextIndV[(*G).E[i][1]];
			nextdoE[k] = (*G).doE[i];
			k++;
		}
		
	FreeIntArray((*G).E, (*G).e);
	free((*G).doE);
	(*G).E = nextE;
	(*G).doE = nextdoE;	
	(*G).e = nexte;

	return nextIndE;
}

//return the array nextIndF
//nextIndF[i] gives the index of the face i after removing the hole
//nextIndF[i] = -1 if fi lives inside the hole
//update G.F
//update G.f
int *RemoveHolesFaces(struct tiling *G, int *HoleF, int *nextIndE)
{
	int i, j;
	
	//contruct nextIndF and compute nextf
	int *nextIndF = malloc((*G).f*sizeof(int));
	int nextf = 0;
	for (i=0; i<(*G).f; i++)
		if (!HoleF[i])
			nextIndF[i] = nextf++;
	
	//allocate face set after hole nextF
	int maxlength = 0;
	for (i=0; i<(*G).f; i++)
		if ((*G).F[i][0]>maxlength) maxlength = (*G).F[i][0];
	
	int **nextF = malloc(nextf*sizeof(int *));
	for (i=0; i<nextf; i++)
	{
		nextF[i] = malloc((maxlength+1)*sizeof(int));
		for (j=0; j<=maxlength; j++) nextF[i][j] = 0;
	}

	//fill nextF using G.F and nextIndF
	for (i=0; i<(*G).f; i++)
		if (HoleF[i] == 0)
		{
			nextF[nextIndF[i]][0] = (*G).F[i][0];
			for (j=1; j<=(*G).F[i][0]; j++)
				nextF[nextIndF[i]][j] = nextIndE[(*G).F[i][j]];
		}

	//uddate faces and f
	FreeIntArray((*G).F, (*G).f);
	(*G).F = nextF;
	(*G).f = nextf;
	
	return nextIndF;
}



//create a hole by removing faces
//HoleF is the indicator vector of the hole faces: HoleF[k]=1 iff
//the k-th face belong to the hole
//To create a hole, we keep the same graph but we change the length 
//of the hole faces and edges to 0. Some vertices may then be isolated
//but these vertices are not removed.
void CreateHole(struct tiling *G, int *HoleF)
{
	//Compute the indicator vector of the hole's edges
	int *HoleE = HoleEAllocation(G, HoleF);

	//remove isolated vertices after hole
	//update G.v
	//udpate G.Coord
	int *nextIndV = RemoveHolesVertices(G, HoleE);
	
	//remove hole's edges
	//update G.e
	//update G.E and endpoints indices
	//udpate G.doE
	int *nextIndE = RemoveHolesEdges(G, HoleE, nextIndV);

	//remove hole's faces
	//update G.f
	//update G.F and edges indices
	int *nextIndF = RemoveHolesFaces(G, HoleF, nextIndE);
		
	//udpate tiling
	//update dE
	//update dV
	//update doV
	//update VToE
	UpdateTiling(G);
	
	//free memory
	free(nextIndV);
	free(nextIndE);
	free(nextIndF);
	free(HoleE);
}



//=======================================
//		Rectangle holes
//=======================================


//create a hole containing all the faces whose center is inside 
//the rectangle defined by any two oposite corners
void CreateRectangleHole(struct tiling *G, int a, int b)
{
	int i, j, ej;
	float ux, uy, vx, vy;

	float ax, ay, bx, by;
	ax = (*G).Coord[a][0];
	ay = (*G).Coord[a][1];
	bx = (*G).Coord[b][0];
	by = (*G).Coord[b][1];

	int *HoleF = malloc(((*G).f)*sizeof(int));
	for (i=0; i<(*G).f; i++) HoleF[i] = 0;
	
	//find the list of faces inside this rectangle
	//run over all the faces and add a face to the hole if all its
	//vertices lives inside the rectangle
	for (i=0; i<(*G).f; i++)
	{
		float centerx = 0;
		float centery = 0;
		for (j=1; j<=(*G).F[i][0]; j++)
		{
			ej = (*G).F[i][j];
			ux = (*G).Coord[(*G).E[ej][0]][0];
			vx = (*G).Coord[(*G).E[ej][1]][0];
			uy = (*G).Coord[(*G).E[ej][0]][1];
			vy = (*G).Coord[(*G).E[ej][1]][1];
			centerx += ux + vx;
			centery += uy + vy;
		}
		centerx = centerx/(2*(*G).F[i][0]);
		centery = centery/(2*(*G).F[i][0]);
		//test if the center belong to one of the 4 possible rectangles
		if (((ax <= centerx) && (centerx <= bx) && (ay <= centery) && (centery <= by))\
||((ax <= centerx) && (centerx <= bx) && (by <= centery) && (centery <= ay))\
||((bx <= centerx) && (centerx <= ax) && (ay <= centery) && (centery <= by))\
||((bx <= centerx) && (centerx <= ax) && (by <= centery) && (centery <= ay)) == 1)
			HoleF[i] = 1;
	}

	CreateHole(G, HoleF);
	free(HoleF);
}





//if p=0 return the x coordinate of the center of the edge ei
//if p=1 return the y coordinate of the center of the edge ei 
float MiddleEdgeCoord(struct tiling G, int i, int p)
{
	int endpoint1 = G.E[i][0];
	int endpoint2 = G.E[i][1];
	return ((G.Coord[endpoint1][p]+G.Coord[endpoint2][p])/2.0);
}

//
void CreateOpenRectangleHole(struct tiling *G, int a, int b)
{	
	float ax, ay, bx, by;
	ax = (*G).Coord[a][0];
	ay = (*G).Coord[a][1];
	bx = (*G).Coord[b][0];
	by = (*G).Coord[b][1];
		
	int i;
	CreateRectangleHole(G, a, b);
	
	int *OpenE = malloc(((*G).e)*sizeof(int));
	for (i=0; i<(*G).e; i++) OpenE[i] = 0;


	//the rectangle is
	//(ax, ay) -- (ax, by) -- (bx, by) -- (bx, ay)
	//run over all the edges and add those of the rectangle to OpenE
	float middlex;
	float middley;
	for (i=0; i<(*G).e; i++)
	{
		middlex = MiddleEdgeCoord((*G), i, 0);
		middley = MiddleEdgeCoord((*G), i, 1);
//		printf("middle = (%f, %f)\n", middlex, middley);
		//test if the middle belong ta a side of one of the 4 possible rectangles
		if (((ax <= middlex) && (middlex <= bx) && (ay <= middley) && (middley <= by))\
||((ax <= middlex) && (middlex <= bx) && (by <= middley) && (middley <= ay))\
||((bx <= middlex) && (middlex <= ax) && (ay <= middley) && (middley <= by))\
||((bx <= middlex) && (middlex <= ax) && (by <= middley) && (middley <= ay)))
			OpenE[i] = 1;
	}

	OpenBoundary(G, OpenE);
	free(OpenE);
}





