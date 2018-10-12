#include <stdlib.h>
#include <stdio.h>
#include "tiling.h"
#include "holes.h"


//============================
// dual size
//============================


//v dual = b of faces + nb of closed boundary edges
int Dualv(struct tiling G)
{
	int i;
	int v2 = G.f;
	for (i=0; i<G.e; i++)
		if ((G.dE[i] == 1) && (G.doE[i] != 1))
			v2++;
	return v2;
}

// e dual = nb of non open edges + nb of closed boundary vertices
int Duale(struct tiling G)
{
	int i;
	int e2 = 0;
	for (i=0; i<G.e; i++)
		if (G.doE[i] != 1) e2++;
		
	for (i=0; i<G.v; i++)
		if (G.dV[i] == 1 && G.doV[i] == 0) 	e2++;
	
	return e2;
}

//f dual = number of non-open vertices
int Dualf(struct tiling G)
{
	int i;
	int f2 = 0;
	for (i=0; i<G.v; i++)
		if (G.doV[i] == 0) f2++;
	return f2;
}


//nb of non-open faces of the dual = nb of interior vertices of G
int Dualnof(struct tiling G)
{
	int i;
	int nof2 = 0;
	for (i=0; i<G.v; i++)
		if (G.dV[i] == 0) nof2++;
	return nof2;
}


//Compute the number of non-open edges
int NonOpene(struct tiling G)
{
	int noe = 0,i;
	for (i=0; i<G.e; i++)
		if (G.doE[i] != 1) noe++;
	return noe;
}


//Compute the number of closed boundary edges
int ClosedBoundarye(struct tiling G)
{
	int dce = 0,i;
	for (i=0; i<G.e; i++)
		if (G.dE[i] == 1 && G.doE[i] == 0) dce++;
	return dce;
}



//============================
// dual tiling
//============================

//There are two types of vertices:
//(i) v_f corresponding to faces
//(ii) v_e corresponding to closed boundary edges
//In the dual graph we use the following indexation of vertices:
//for k<f: k is the vertex v_fk of type (i) associated with the k-th face fk
//then:f+k is the index v_ek, of type (ii) associated with the k-th closed boundary edge ek

//There are three types of dual edges:
//(i) v_f, v_{f'} corresponding with pairs of incident faces and 
//(ii) v_f v_e corresponding where e is a closed boundary edge incident to f
//(iii) v_e v_e' where e and e' are closed boundary edges and share a closed vertex v; that is e = {u,v} and e'={v,u'}.
//In the dual graph we use the following indexation of vertices:
//for k<noe: the k-th edge of type (i) or (ii) of the dual corresponds to the k-th non-open edge of the original graph.
//then: noe+k is the index of the k-th edge of type (iii) corresponding to the k-th closed vertex of G. It is the set of open edges.



//============================
// dual vertices index
//============================


//return the index of the dual vertex v_f of the face f
//It is the index of f
int *FToDualV(struct tiling G)
{
	int *FToDualV = malloc(G.f*sizeof(int));
	int i;
	for (i=0; i<G.f; i++)
		FToDualV[i] = i;
	return FToDualV;	
}


//return the index of the dual vertex v_e of the closed boundary edge e
//if e is the k-th closed boundary v_e is indexed by f+k
//if e is not a closed boundary returns -1.
int *EToDualV(struct tiling G)
{
	int *EToDualV = malloc(G.e*sizeof(int));
	int k=0,i;
	for (i=0; i<G.e; i++)
		if (G.dE[i] == 1 && G.doE[i] == 0)
		{
			EToDualV[i] = G.f + k;
			k++;
		}
		else EToDualV[i] = -1;
	return EToDualV;
}

//============================
// dual edges index
//============================


//return the index of the dual edge of the k-th edge
//the dual edge associated with the k-th non-open edge is k
//if e is not a non-open edge it has no dual, then return -1
int *EToDualE(struct tiling G)
{
	int *EToDualE = malloc(G.e*sizeof(int));
	int k=0,i;
	for (i=0; i<G.e; i++)
		if (G.doE[i] == 0)
		{
			EToDualE[i] = k;
			k++;
		}
		else EToDualE[i] = -1;
	return EToDualE;
}


//return the index of the dual edge of the k-th vertex
//the index of the dual edge of the k-th closed boundary vertex is noe+k
//if v is not a closed boundary vertex return -1
int *VToDualE(struct tiling G)
{
	int noe = NonOpene(G);
	int *VToDualE = malloc(G.v*sizeof(int));
	int k=0,i;
	for (i=0; i<G.v; i++)
		if (G.dV[i] == 1 && G.doV[i] == 0)
		{
			VToDualE[i] = noe + k;
			k++;
		}
		else VToDualE[i] = -1;
	return VToDualE;
}


//============================
// dual faces index
//============================


//return the index of the dual face of the k-th vertex
//dual face correspond to non-open vertices
//the dual face associated with the k-th non-open vertex is k
//if f is not a non-open face it has no dual, then return -1
int *VToDualF(struct tiling G)
{
	int *VToDualF = malloc(G.v*sizeof(int));
	int k=0,i;
	for (i=0; i<G.v; i++)
		if (G.doV[i] == 0)
		{
			VToDualF[i] = k;
			k++;
		}
		else VToDualF[i] = -1;
	return VToDualF;	
}



//============================
// 
//============================



//returns the array of indices of an edge among only non-open edges
//nonOpen[k] = index in E of the k-th non-open edge
int *NonOpenE(struct tiling G)
{	
	int i, k;
	int noe = NonOpene(G);
	int *nonOpenE = malloc(noe*sizeof(int));
	k=0;
	for (i=0; i<G.e; i++)
		if (G.doE[i] != 1)
		{
			nonOpenE[k] = i;
			k++;
		}
	return nonOpenE;
}



//construct the incidence edges/face
int **EToFConstruction(struct tiling G)
{
	int i, j, ej;
	int **EToF = malloc(G.e*sizeof(int *));
	for (i=0; i<G.e; i++)
	{
		EToF[i] = malloc(3*sizeof(int));
		EToF[i][0] = 0;
	}

	for (i=0; i<G.f; i++)
	{
		for (j=1; j<=G.F[i][0]; j++) 
		{
			ej = G.F[i][j];
			EToF[ej][0]++;	//EToF[ej][0] is the nb of faces incident to ej
			EToF[ej][EToF[ej][0]] = i; //add the face i to this list
		}
	}
	return EToF;
}



//============================
// Dual edges construction
//============================


//construct the dual edges
int **DualEdges(struct tiling G, int v2, int e2)
{
	int i, j;
	
	int **E2 = malloc(e2*sizeof(int *));
	for (i=0; i<e2; i++)
		E2[i] = malloc(2*sizeof(int));
		
	//construct the correspondances: 
	//from faces and closed boundary edges to dual vertices
	//from non-open edges and closed boundary vertices to dual edges
	int *FToDV = FToDualV(G);
	int *EToDV = EToDualV(G);
	int *EToDE = EToDualE(G);
	int *VToDE = VToDualE(G);

	//use EToF to construct dual edges of type (i) and (ii).
	int **EToF = EToFConstruction(G);
	int dualedge, endpoint1, endpoint2;
	for (i=0; i<G.e; i++)
	{
		if (G.doE[i] == 0)	//if ei is non-open
		{
			if (EToF[i][0] == 2)	//type (i) edge
			{
				dualedge = EToDE[i];
				endpoint1 = FToDV[EToF[i][1]];
				endpoint2 = FToDV[EToF[i][2]];
				E2[dualedge][0] = endpoint1;
				E2[dualedge][1] = endpoint2;
			}
			if (EToF[i][0] == 1)	//type (ii) edge
			{
				dualedge = EToDE[i];
				endpoint1 = FToDV[EToF[i][1]];
				endpoint2 = EToDV[i];
				E2[dualedge][0] = endpoint1;
				E2[dualedge][1] = endpoint2;
			}
		}
	}
	
	//add type (iii) open edges
	//add a closed edge for each closed boundary vertex v
	for (i=0; i<G.v; i++)
		if (G.dV[i] == 1 && G.doV[i] == 0) //if vi is a closed boundary
		{
			dualedge = VToDE[i];
			j=1;
			while (j<=G.VToE[i][0] && G.dE[G.VToE[i][j]]!=1) j++;
			endpoint1 = EToDV[G.VToE[i][j]];
			E2[dualedge][0] = endpoint1;
			j++;
			while (j<=G.VToE[i][0] && G.dE[G.VToE[i][j]]!=1) j++;
			endpoint2 = EToDV[G.VToE[i][j]];
			E2[dualedge][1] = endpoint2;
		}

	return E2;
}


//============================
// Dual faces construction
//============================


//
int **DualFaces(struct tiling G)
{
	int i, j;
	int edual, fdual;

	//the max degree of a face of the dual is at most the max degree of a
	//vertex +1 in G
	int lengthmax = 0;
	for (i=0; i<G.v; i++)
		if (G.VToE[i][0] + 1 > lengthmax) 
			lengthmax = G.VToE[i][0] + 1;
	
	int f2 = Dualf(G);
	int **F2 = malloc(f2*sizeof(int *));
	for (i=0; i<f2; i++)
	{
		F2[i] = malloc((lengthmax+1)*sizeof(int));
		F2[i][0] = 0;
	}

	//EToDE[k] is the index of ek of G or -1 if ek is open
	int *EToDE = EToDualE(G);
	int *VToDF = VToDualF(G);
	int *VToDE = VToDualE(G);

	//loop over non-open vertices 
	//transform edge set F_v incident to a vertex v in G 
	//into a face (or path when v is a boundary) of the dual
	for (i=0; i<G.v; i++)
		if (G.doV[i] == 0)
			for (j=1; j<=G.VToE[i][0]; j++)
				{
					edual = EToDE[G.VToE[i][j]];
					fdual = VToDF[i];
					F2[fdual][0]++;
					F2[fdual][F2[fdual][0]] = edual; 
					//add edual to the edges of the dual face fdual
				}

	//loop over closed boundary vertices to add the closed open to each open face
	//This edge edual is the edge dual to the vertex v.
	for (i=0; i<G.v; i++)
		if (G.dV[i] == 1 && G.doV[i] == 0)
		{
			edual = VToDE[i];
			fdual = VToDF[i];
			F2[fdual][0]++;
			F2[fdual][F2[fdual][0]] = edual; //add edual around fdual
		}
	return F2;
}



//============================
// construction of the dual coordinate array
//============================


//
float **DualCoordAllocate(struct tiling G, int dualv)
{
	int i, j;
	float **DualCoord = malloc(dualv*sizeof(float *));
	for (i=0; i<dualv; i++)
	{
		DualCoord[i] = malloc(2*sizeof(float));
		DualCoord[i][0] = 0;
		DualCoord[i][1] = 0;
	}

	//fill the coordinate of vertices of the dual coming from faces of G
	int *FToDV = FToDualV(G);
	int ej;
	float centerx, centery, ux, uy, vx, vy;	
	for (i=0; i<G.f; i++)
	{
		if (FToDV[i] != -1)
		{			
			//compute the coordinate of the center fo the face fi
			centerx = 0;
			centery = 0;
			for (j=1; j<=G.F[i][0]; j++)
			{
				ej = G.F[i][j];
				ux = G.Coord[G.E[ej][0]][0];
				vx = G.Coord[G.E[ej][1]][0];
				uy = G.Coord[G.E[ej][0]][1];
				vy = G.Coord[G.E[ej][1]][1];
				centerx += ux + vx;
				centery += uy + vy;
			}
			centerx = centerx/(2*G.F[i][0]);
			centery = centery/(2*G.F[i][0]);
			DualCoord[FToDV[i]][0] = centerx;
			DualCoord[FToDV[i]][1] = centery;
		}
	}
	free(FToDV);
	
	//fill the coordinate of vertices of the dual coming from closed edges of G
	int *EToDV = EToDualV(G);
	for (i=0; i<G.e; i++)
	{
		if (EToDV[i] != -1)
		{
			ux = G.Coord[G.E[i][0]][0];
			vx = G.Coord[G.E[i][1]][0];
			uy = G.Coord[G.E[i][0]][1];
			vy = G.Coord[G.E[i][1]][1];
			DualCoord[EToDV[i]][0] = (ux+vx)/2;
			DualCoord[EToDV[i]][1] = (uy+vy)/2;
		}
	}
	free(EToDV);
	
	return DualCoord;
}



//============================
// construction of the dual tiling
//============================


//
struct tiling DualTiling(struct tiling G)
{
	int i;
	struct tiling H;
	InitTiling(&H);
	
	H.length = G.length;
	char *type = malloc(50*sizeof(char));
	sprintf(type, "%s_%s", G.type, "dual");
	H.type = type;

	//dual size
	H.v = Dualv(G);
	H.e = Duale(G);
	H.f = Dualf(G);

	//dual edges	
	H.E = DualEdges(G, H.v, H.e);

	//dual faces
	H.F = DualFaces(G);

	//dual open edges
	H.doE = malloc(H.e*sizeof(int));
	for (i=0; i<H.e; i++) H.doE[i] = 0;

	int *VToDE = VToDualE(G);		
	for (i=0; i<G.v; i++)
		if (VToDE[i] != -1)
			H.doE[VToDE[i]] = 1;

	//dual doordinates
	if (G.Coord != NULL)
		H.Coord = DualCoordAllocate(G, H.v);
	
	//update incidence, boundaries and open vertices
	UpdateTiling(&H);

	return H;
}
