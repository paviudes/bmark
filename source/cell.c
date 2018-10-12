#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tiling.h"

static const float epsilon = 0.001;




//==================================
//  SINGLE FACE CELL
//==================================


//
float **ConstructCoord(float **vertexCoord, int length)
{
	int i;
	float **Coord = malloc(length*sizeof(float *));
	for (i=0; i<length; i++) Coord[i] = malloc(2*sizeof(float));

	for (i=0; i<length; i++)
	{
		Coord[i][0] = vertexCoord[i][0];
		Coord[i][1] = vertexCoord[i][1];
//		printf("debug v%d = (%f,%f)\n", i, vertexCoord[i][0], vertexCoord[i][1]);
	}
	return Coord;
}


//
int **ConstructEdges(int length)
{
	int i;
	int **E = malloc(length*sizeof(int *));
	for (i=0; i<length; i++) E[i] = malloc(2*sizeof(int));
		
	for (i=0; i<length; i++)
	{
		E[i][0] = i;
		E[i][1] = (i+1)%length;
	}
	return E;
}


//
int **ConstructFace(int **E, int length)
{
	int i;
	int **F = malloc(1*sizeof(int *));
	F[0] = malloc((length+1)*sizeof(int));

	F[0][0] = length;
	for (i=1; i<=F[0][0]; i++)
		F[0][i] = i-1;
	return F;
}


//
int *ConstructOpenEgdes(int length)
{
	int i;
	int *doE = malloc(length*sizeof(int));
	for (i=0; i<length; i++) doE[i] = 0;
	return doE;
}


//return a tiling representing a single face from the ordered
//list of coordinates of its vertices
//vertexCoord[i][0] is the x-coord of the i-th vertex
//vertexCoord[i][1] is the y-coord of the i-th vertex
struct tiling ConstructCell(float **vertexCoord, int length)
{
	struct tiling G;
	InitTiling(&G);
	G.length = length;

	char *type = malloc(50*sizeof(char));
	sprintf(type, "face_%d", length);
	G.type = type;

	G.v = length;
	G.e = length;
	G.f = 1;

	G.E = ConstructEdges(length);
	G.F = ConstructFace(G.E, length);
	G.doE = ConstructOpenEgdes(length);
	G.Coord = ConstructCoord(vertexCoord, length);
	UpdateTiling(&G);

	return G;
}




//==================================
//	TRANSLATION OF A CELL
//==================================


//return the array identifyVertex such that
//the vertex i of the translated cell t(C)
//is identified with the vertex identifiedV[i] in the cell C
//IdentifiedV[i] = -1 if the vertex vi of t(C) is not
//identified to a vertex of C.
int *FindIdentifiedCellVertex(struct tiling C, float *t)
{
	int i, j;
	int *IdentifiedV = malloc(C.v*sizeof(int));
	for (i=0; i<C.v; i++)
		IdentifiedV[i] = -1;
	
	//compare the coordinates of t(vi) and vj
	for (i=0; i<C.v; i++)
		for (j=0; j<C.v; j++)
			if ((fabs(C.Coord[i][0]+t[0]-C.Coord[j][0]) < epsilon) && (fabs(C.Coord[i][1]+t[1]-C.Coord[j][1]) < epsilon))
				IdentifiedV[i] = j;
	
	return IdentifiedV;
}



//compare middle point of edges to find identified edges
//IdentifiedE[i] = j if the edge ei of t(C) is identified with ej in C
//IdentifiedE[i] = -1 if ei is not identified to an edge of C.
int *FindIdentifiedCellEdge(struct tiling C, float *t)
{
	int i, j;	
	int *identifiedE = malloc(C.e*sizeof(int));
	for (i=0; i<C.e; i++)
		identifiedE[i] = -1;

	int ui, vi, uj, vj;
	double middle_ei_x, middle_ei_y, middle_ej_x, middle_ej_y;
	for (i=0; i<C.e; i++)
		for (j=0; j<C.e; j++)
		{
			//compare the edges ei={ui,vi} and ej={uj,vj} using their middle
			ui = C.E[i][0];
			vi = C.E[i][1];

			uj = C.E[j][0];
			vj = C.E[j][1];

			middle_ei_x = (C.Coord[ui][0] + C.Coord[vi][0])/2;
			middle_ei_y = (C.Coord[ui][1] + C.Coord[vi][1])/2;

			middle_ej_x = (C.Coord[uj][0] + C.Coord[vj][0])/2;
			middle_ej_y = (C.Coord[uj][1] + C.Coord[vj][1])/2;
			
			if ((fabs(middle_ei_x+t[0]-middle_ej_x) < epsilon) && (fabs(middle_ei_y+t[1]-middle_ej_y) < epsilon))
				identifiedE[i] = j;
		}
		
	
	return identifiedE;
}


//return the barycenter of the face fi in the tiling G.
float *FaceBarycenter(struct tiling G, int i)
{
	int j;
	float *center = malloc(2*sizeof(float));
	center[0] = 0;
	center[1] = 0;
	
	int ej;
	float ux, uy, vx, vy;
	for (j=1; j<=G.F[i][0]; j++)
	{
		ej = G.F[i][j];
		ux = G.Coord[G.E[ej][0]][0];
		vx = G.Coord[G.E[ej][1]][0];
		uy = G.Coord[G.E[ej][0]][1];
		vy = G.Coord[G.E[ej][1]][1];
		center[0] += ux + vx;
		center[1] += uy + vy;
	}
	center[0] = center[0]/(2*G.F[i][0]);
	center[1] = center[1]/(2*G.F[i][0]);
	
	return center;
}


//Compare barycenter of faces to find identified faces.
//IdentifiedF[i] = j if the face fi of t(C) is identified with fj in C
//IdentifiedF[i] = -1 if fi is not identified to a face of C.
int *FindIdentifiedCellFace(struct tiling C, float *t)
{
	int i, j;	
	int *identifiedF = malloc(C.f*sizeof(int));
	for (i=0; i<C.f; i++)
		identifiedF[i] = -1;
	
	float *center_fi, *center_fj;
	for (i=0; i<C.f; i++)
		for (j=0; j<C.f; j++)
		{
			center_fi = FaceBarycenter(C, i);
			center_fj = FaceBarycenter(C, j);
			if ((fabs(center_fi[0]+t[0]-center_fj[0]) < epsilon) && (fabs(center_fi[1]+t[1]-center_fj[1]) < epsilon))
				identifiedF[i] = j;
		}
	
	return identifiedF;
}



//Before gluing the sides of the n copies of the cell C,
//there are n*C.v vertices.
//Denote by C(i) the cell translated by t1^i.
//The k-th vertex of the cell C(i) is indexed by i*C.v+k,
//with i=0,..,n-1, k=0,..,C.v-1
//This function returns an array gridV such that
//gridV[i*C.v+k] is the vertex after identification or
//gridV[i*C.v+k] = i*C.v+k if it is not identified
int *GridVertices(struct tiling C, float *t, int n)
{
	int i, k;
	//identifications by the translation
	int *identV = FindIdentifiedCellVertex(C, t);
	
	//initialize all the vertices to gridV[i*C.v+k] = i*C.v+k 
	int *gridV = malloc((n*C.v)*sizeof(int));
	for (i=0; i<n; i++)
		for (k=0; k<C.v; k++)
			gridV[i*C.v+k] = i*C.v+k;
	
	//replace the vertex i*C.v+k of the cell C(i) with i>0 by
	//(i-1)*C.v + identV[k] in the cell C(i-1) when identV[k] != -1
	for (i=0; i<n; i++)
		for (k=0; k<C.v; k++)
			if (identV[k] != -1 && i>0)
				gridV[i*C.v+k] = (i-1)*C.v+identV[k];
			
	return gridV;
}

//Return an array gridE such that
//gridE[i*C.e+k] is the vertex after identification or
//gridE[i*C.e+k] = i*C.e+k if it is not identified
int *GridEdges(struct tiling C, float *t, int n)
{
	int i, k;
	//identifications by the 2 translations
	int *identE = FindIdentifiedCellEdge(C, t);
	
	//initialize all the edges to gridE[i*C.e+k] = i*C.e+k
	int *gridE = malloc((n*C.e)*sizeof(int));
	for (i=0; i<n; i++)
		for (k=0; k<C.e; k++)
			gridE[i*C.e+k] = i*C.e+k;
	
	//replace the edge i*C.e+k of the cell C(i) with i>0 by
	//(i-1)*C.e + identE[k] in the cell C(i-1) when identE[k] != -1
	for (i=0; i<n; i++)
		for (k=0; k<C.e; k++)
			if (identE[k] != -1 && i>0)
				gridE[i*C.e+k] = (i-1)*C.e+identE[k];
	
	return gridE;
}


//Return an array gridF of size n*C.f such that
//gridF[i*C.f+k] is the vertex after identification or
//gridF[i*C.f+k] = i*C.f+k if it is not identified
int *GridFaces(struct tiling C, float *t, int n)
{
	int i, k;
	//identifications by the 2 translations
	int *identF = FindIdentifiedCellFace(C, t);
	
	//initialize all the edges to gridF[i*C.f+k] = i*C.f+k
	int *gridF = malloc((n*C.f)*sizeof(int));
	for (i=0; i<n; i++)
		for (k=0; k<C.f; k++)
			gridF[i*C.f+k] = i*C.f+k;
			
	//replace the edge i*C.f+k of the cell C(i) with i>0 by
	//(i-1)*C.v + identE[k] in the cell C(i-1) when identE[k] != -1
	for (i=0; i<n; i++)
		for (k=0; k<C.f; k++)
			if (identF[k] != -1 && i>0)
				gridF[i*C.f+k] = (i-1)*C.f+identF[k];
					
	return gridF;
}


//return a tiling obtained by translating a generating cell 
//n-1 time in the direction t
struct tiling PeriodicTiling(struct tiling C, float *t, int n)
{
	int i, k;
	
	struct tiling G;
	InitTiling(&G);
	G.length = C.length*n;

	char *type = malloc(50*sizeof(char));
	sprintf(type, "%s_x%d", C.type, n);
	G.type = type;
	
	//vertices:
	//the i-th vertex is identified with gridV[i]
	int *gridV = GridVertices(C, t, n);

	//Count the nb of distinct vertices in G	
	k = 0;
	for (i=0; i<n*C.v; i++)
		if (gridV[i] == i)
			k++;
	G.v = k;


	//Remove unused vertex indices in gridV obtain vertices of G
	//the vertice to conserve are those with gridV[i]=i.
	//indexV[k] contains the k-th vertex of gridV after identification
	//and removing unused indices
	int *indexV = malloc(n*C.v*sizeof(int));
	k = 0;
	for (i=0; i<n*C.v; i++)
		if (gridV[i] == i)
			indexV[i] = k++;

	int fixpoint;
	for (i=0; i<n*C.v; i++)
		if (gridV[i] != i)
		{
			fixpoint = i;
			while (gridV[fixpoint] != fixpoint)
				fixpoint = gridV[fixpoint];
			indexV[i] = indexV[fixpoint];
		}


	//construct Coord by translating Coord of the cell C by t
	float **Coord = malloc(G.v*sizeof(float*));
	for (i=0; i<G.v; i++)
		Coord[i] = malloc(2*sizeof(float));

	for (i=0; i<n; i++)
		for (k=0; k<C.v; k++)
			if (gridV[i*C.v+k] == i*C.v+k)
				{
					Coord[indexV[i*C.v+k]][0] = C.Coord[k][0]+i*t[0];
					Coord[indexV[i*C.v+k]][1] = C.Coord[k][1]+i*t[1];				}

	G.Coord = Coord;
	
	
	//Edges:
	//i*C.e+k is the index of the k-th edges of the cell C(i)
	//This edge will be identified with the edge gridE[i]
	//gridE[i] = i means that ei is not identified
	int *gridE = GridEdges(C, t, n);
	
	
	//Count the nb of distinct edges in G
	k = 0;
	for (i=0; i<n*C.e; i++)
		if (gridE[i] == i)
			k++;
	G.e = k;
	
	//Remove unused edges indices in G
	//The index of the i-th edge of the grid becomes
	//indexE[i] after removing unused indices
	int *indexE = malloc(n*C.e*sizeof(int));
	k = 0;
	for (i=0; i<n*C.e; i++)
		if (gridE[i] == i)
			indexE[i] = k++;

	for (i=0; i<n*C.e; i++)
		if (gridE[i] != i)
		{
			fixpoint = i;
			while (gridE[fixpoint] != fixpoint)
				fixpoint = gridE[fixpoint];
			indexE[i] = indexE[fixpoint];
		}
	
	
	//construct E
	int **E = malloc(G.e*sizeof(int*));
	for (i=0; i<G.e; i++)
		E[i] = malloc(2*sizeof(int));
	
	int endpoint1, endpoint2;
	for (i=0; i<n; i++)
		for (k=0; k<C.e; k++)
			if (gridE[i*C.e+k] == i*C.e+k)
			{
				endpoint1 = i*C.v + C.E[k][0];
				endpoint2 = i*C.v + C.E[k][1];
				E[indexE[i*C.e+k]][0] = indexV[endpoint1];
				E[indexE[i*C.e+k]][1] = indexV[endpoint2];
			}
 			
	G.E = E;


	//Faces:
	//i*C.f+k is the index of the k-th face of the cell C(i)
	//This edge will be identified with the edge gridF[i]
	//gridF[i] = i means that fi is not identified
	int *gridF = GridFaces(C, t, n);
	
	//Count the nb of distinct faces in G
	k = 0;
	for (i=0; i<n*C.f; i++)
		if (gridF[i] == i)
			k++;
	G.f = k;
	
	//Remove unused faces indices in G
	//The index of the i-th face of the grid becomes
	//indexF[i] after removing unused indices
	int *indexF = malloc(n*C.f*sizeof(int));
	k = 0;
	for (i=0; i<n*C.f; i++)
		if (gridF[i] == i)
			indexF[i] = k++;

	for (i=0; i<n*C.f; i++)
		if (gridF[i] != i)
		{
			fixpoint = i;
			while (gridF[fixpoint] != fixpoint)
				fixpoint = gridF[fixpoint];
			indexF[i] = indexF[fixpoint];
		}
	
	//construct F
	int maxlength = 0;
	for (i=0; i<C.f; i++)
		if (C.F[i][0] > maxlength) maxlength = C.F[i][0];
	
	int **F = malloc(G.f*sizeof(int*));
	for (i=0; i<G.f; i++)
		F[i] = malloc((maxlength+1)*sizeof(int));
	
	int j;
	int face, ej;
	for (i=0; i<n; i++)
		for (k=0; k<C.f; k++)
			if (gridF[i*C.f+k] == i*C.f+k)
			{
				face = indexF[i*C.f+k];
				F[face][0] = C.F[k][0];
				for (j=1; j<=C.F[k][0]; j++)
				{
					ej = i*C.e + C.F[k][j];
					F[face][j] = indexE[ej];
				}
			}
 			
	G.F = F;
	
	//Open edges
	int *doE = malloc(G.e*sizeof(int));
	for (i=0; i<G.e; i++)
		doE[i] = 0;

	for (i=0; i<n; i++)
		for (k=0; k<C.f; k++)
			doE[indexE[i*C.e+k]] = C.doE[k];

	G.doE = doE;
	
	//Update the rest of the tiling
	UpdateTiling(&G);

	free(gridV);
	free(indexV);
	free(gridE);
	free(indexE);
	free(gridF);
	free(indexF);

	return G;
}








