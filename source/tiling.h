#ifndef TILING_H
#define TILING_H


//length: the size of the lattices or simply an index of a lattice in a family
//type: name of the tiling
//v: number of vertices
//e: number of edges
//f: number of faces

//vertices, edges and faces are indexed by i=0,1,...,s-1 for s=v,e or f and are denoted vi, ei fi.

//E is the array representing the edges.
//E[i][0] and E[i][1] are the indices of the endpoints of i-th edge ei.

//F is the array representing the faces.
//F[i][0] is the length of the i-th face fi.
//F[i][j] is the index of the j-th edge incident to fi

//VToE encodes the vertex/edge incidence relations.
//VToE[i][0] is the degree of vi.
//VToE[i][j] is the index of the j-th edge indicent to vi for j=1,2,...,VToE[i][0].

//VToN is the array that encodes the incidence relation between vertices.
//VToN[i][0] is the degree of vi.
//VToN[i][j] is the index of the j-th vertex indicent to vi, for j=1,2,...,VToN[i][0].

//dV and dE define the boundaries. dF is not needed here.
//dE[k]=1 if the edge k is a boundary and 0 otherwise
//dV[k]=1 if the vertex k is a boundary and 0 otherwise

//doV and doE define open boundaries.
//dE[k]=1 if the edge k is an open boundary and 0 otherwise
//dV[k]=1 if the vertex k is an open boundary and 0 otherwise

//Coord stores the euclidean coordintates of the vertices
//It is required to draw an euclidean tiling
//We can leave it to NULL
//Coord[i][0] is the x coordinate of vi
//Coord[i][1] is the y coordinate of vi

struct tiling
{
	int length;
	char *type;
	int v;
	int e;
	int f;
	int **E;
	int **F;
	int **VToE;
	int **VToN;
	int *dV;
	int *dE;
	int *doE;
	int *doV;
	float **Coord;
};


//initialize length, v, e and f to 0.
//initialize all pionters to NULL.
//This function does not free the pointers in the structure.
//Do not apply it to an allocated pointer.
extern void InitTiling(struct tiling *G);

//free all the pointers of the tiling
//set length, v, e,  f, 0 to 0.
//set pointers to NULL.
extern void FreeTiling(struct tiling *G);

//print the properties of th tiling in the terminal:
//type, size, vertices, edges and faces
extern void PrintTiling(struct tiling G);

//Once v, e, f, E, F and doE are defined, use it to update the tiling
//update dE
//update dvE
//update dovE
//update VToE
extern void UpdateTiling(struct tiling *G);

//return an array v containing an ordered list of the vertices of a face
//The second input is a array c containing the edges of a face:
//c[0] is the length and c[1], ..., c[c[0]] are the indices of these edges 
//v[0] is the length and v[1], ..., v[v[0]] are the indices fo these vertices
extern int *FaceToVertexList(struct tiling, int *);



//tilingName contains the name G.graphType of the tiling.
//dataFile contains the name of the file containing the
//	performance results in the folder results/.
//[noiseRange[0], noiseRange[1]] is the interval over which the
//	simulation is performed.
//codeParameter[0] = n = nb of physical qubits.
//codeParameter[1] = k = nb of logical qubits.
//The overhead is n/k.
//maxWeight is the maximum weight of a measured operator
//ZWeight[i] contains the number of Z-measurements of weight i
//XWeight[i] contains the number of Z-measurements of weight i
//	for i=0, 1, ..., maxWeight.
struct storePerformance
{
	char *tilingName;
	char *dataFile;
	double *noiseRange;
	int n;
	int k;
	int maxWeight;
	int *ZWeight;
	int *XWeight;
};

extern void InitStorage(struct storePerformance *perf);
extern void FreeStorage(struct storePerformance *perf);


#endif	/* TILING_H */
