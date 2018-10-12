#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "tiling.h"
#include "cell.h"
#include "draw.h"
#include "grid.h"

//==================================
//	hexagonal Tiling
//==================================


//
struct tiling HexagonalTiling(int Lx, int Ly)
{
	int length = 6,i;
	float **vertexCoord = malloc(length*sizeof(float *));
	for (i=0; i<length; i++)
		vertexCoord[i]=malloc(2*sizeof(float));
		
	vertexCoord[0][0] = 0;
	vertexCoord[0][1] = 0;

	vertexCoord[1][0] = 1;
	vertexCoord[1][1] = 1;

	vertexCoord[2][0] = 1;
	vertexCoord[2][1] = 2;

	vertexCoord[3][0] = 0;
	vertexCoord[3][1] = 3;

	vertexCoord[4][0] = -1;
	vertexCoord[4][1] = 2;

	vertexCoord[5][0] = -1;
	vertexCoord[5][1] = 1;

	struct tiling G;
	InitTiling(&G);
	struct tiling Gtemp;
	InitTiling(&Gtemp);

	float *t1 = malloc(2*sizeof(float));
	t1[0] = 2;
	t1[1] = 0;
	float *t2 = malloc(2*sizeof(float));
	t2[0] = 1;
	t2[1] = 2;

	G = ConstructCell(vertexCoord, length);

	Gtemp = PeriodicTiling(G, t1, Lx);
	FreeTiling(&G);
	G = PeriodicTiling(Gtemp, t2, Ly);
	FreeTiling(&Gtemp);

	for (i=0; i<length; i++)
		free(vertexCoord[i]);
	free(vertexCoord);
	free(t1);
	free(t2);
	
	free(G.type);
	char *type = malloc(100*sizeof(char));
	sprintf(type, "hex_%dx%d", Lx, Ly);
	G.type = type;

	return G;
}


//==================================
//	triangular Tiling
//==================================


//
struct tiling TriangularTiling(int Lx, int Ly)
{
	struct tiling G;
	InitTiling(&G);
	struct tiling Gtemp;
	InitTiling(&Gtemp);	
	
	G = Grid(1, 1);

	int *vertexList = malloc(3*sizeof(int));

	vertexList[0] = 0;
	vertexList[1] = 1;
	vertexList[2] = 2;
	AddFace(&G, vertexList, 3);

	vertexList[0] = 3;
	vertexList[1] = 1;
	vertexList[2] = 2;
	AddFace(&G, vertexList, 3);

	float *t1 = malloc(2*sizeof(float));
	t1[0] = 1;
	t1[1] = 0;
	float *t2 = malloc(2*sizeof(float));
	t2[0] = 0;
	t2[1] = 1;

	Gtemp = PeriodicTiling(G, t1, Lx);
	FreeTiling(&G);
	G = PeriodicTiling(Gtemp, t2, Ly);
	FreeTiling(&Gtemp);

	free(vertexList);
	free(t1);
	free(t2);

	free(G.type);
	char *type = malloc(100*sizeof(char));
	sprintf(type, "triangle_%dx%d", Lx, Ly);
	G.type = type;	
	
	return G;
}



//==================================
//	square Tiling
//==================================


//
struct tiling SquareTiling(int Lx, int Ly)
{
	struct tiling G;
	InitTiling(&G);
	struct tiling Gtemp;
	InitTiling(&Gtemp);	
	
	G = Grid(1, 1);

	int *vertexList = malloc(4*sizeof(int));

	vertexList[0] = 0;
	vertexList[1] = 1;
	vertexList[2] = 3;
	vertexList[3] = 2;	
	AddFace(&G, vertexList, 4);

	float *t1 = malloc(2*sizeof(float));
	t1[0] = 1;
	t1[1] = 0;
	float *t2 = malloc(2*sizeof(float));
	t2[0] = 0;
	t2[1] = 1;

	Gtemp = PeriodicTiling(G, t1, Lx);
	FreeTiling(&G);
	G = PeriodicTiling(Gtemp, t2, Ly);
	FreeTiling(&Gtemp);

	free(vertexList);
	free(t1);
	free(t2);
	
	free(G.type);
	char *type = malloc(100*sizeof(char));
	sprintf(type, "square_%dx%d", Lx, Ly);
//	sprintf(type, "square_%dx%d_%ld", Lx, Ly, time(0));
	G.type = type;
	
	return G;
}



