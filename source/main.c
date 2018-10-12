#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "memory.h"
#include "tiling.h"
// #include "holes.h"
// #include "draw.h"
// #include "dual.h"
// #include "weights.h"
// #include "cell.h"
// #include "grid.h"
// #include "regulartilings.h"
// #include "components.h"
#include "simulation.h"
#include "rgens.h"
#include "arrays.h"
#include "save.h"
#include "load.h"
#include "notations.h"
#include "mcresult.h"
#include "noise.h"
#include "homological.h"
#include "performance.h"
#include "plot.h"
#include "submission.h"
#include "remote.h"

// static const float maxFamilySize = 50;

int main(int argc, char **argv)
{
	srand(time(NULL));
	if (argc > 1){
		// If the user exclusively wants to estimate performances, plot already available data, try out different noise model including correlations.
		printf("\033[92mReading inputs from %s.\033[0m\n", argv[1]);
		Remote(argv[1]);
	}
	return 0;
}
#if 0
	// ***********************************************

	//initialize rand()
	srand((unsigned)(time(0)));

	char function[100];
	function[0] = 0;
	int *parameter = malloc(2*sizeof(int));
	parameter[0] = 0;
	parameter[1] = 0;
	int *size = malloc(sizeof(int));
	size[0] = 0;
	
	//tiling parameters
	struct tiling G; 
	InitTiling(&G);
	struct tiling Gtemp; 
	InitTiling(&Gtemp);

	int *H;
	int i;
	int *OpenE;
	int *CloseE;

	//perf parameters
	struct storePerformance *perf = malloc(maxFamilySize*sizeof(struct storePerformance));
	int familySize = 0;

	/*
	double *noiseRange = malloc(sizeof(double) * 3);

	char **DataFName = malloc(maxFamilySize*sizeof(char *));	//increase to have larger families
	for (int i=0; i<maxFamilySize; i++)
		DataFName[i] = NULL;
	char **tilingName = malloc(maxFamilySize*sizeof(char *));	//increase to have larger families
	for (int i=0; i<maxFamilySize; i++)
		tilingName[i] = NULL;

	int **codeParameters = malloc(maxFamilySize*sizeof(int *));
	for (int i=0; i<maxFamilySize; i++)
		codeParameters[i] = NULL;

	int **ZWeight = malloc(maxFamilySize*sizeof(int *));
	int **XWeight = malloc(maxFamilySize*sizeof(int *));
	for (int i=0; i<maxFamilySize; i++)
	{
		ZWeight[i] = NULL;
		XWeight[i] = NULL;
	}
*/

/*
//================================
//		test 
//================================	

	int Lx;
	int Ly;
	int Lmax = 10;
	int N = 1000;

//================================
//		test Grid
//================================

	printf("test Grid\n");
	for (int i=0; i<N; i++)
	{
		Lx = Lmax*((float)rand()/RAND_MAX);
		Ly = Lmax*((float)rand()/RAND_MAX);
//		printf("%d x %d\n", Lx, Ly);
		G = Grid(Lx, Ly);
		FreeTiling(&G);
	}
	printf("done\n");

//================================
//		test AddFace
//================================
		
	printf("test AddFace\n");
	for (int i=0; i<N; i++)
	{
		Lx = Lmax*((float)rand()/RAND_MAX);
		Ly = Lmax*((float)rand()/RAND_MAX);
//		printf("%d x %d", Lx, Ly);
		G = Grid(Lx, Ly);
		int *vertexList;
		int l;
		
		for (int j=0; j<2; j++)
		{
			l = (int)(3+5*((float)rand()/RAND_MAX));
			vertexList = malloc(l*sizeof(float *));
			for (int i=0; i<l; i++)
				vertexList[i] = (int)(G.v-1)*((float)rand()/RAND_MAX);
		}
		AddFace(&G, vertexList, l);
		FreeTiling(&G);
	}
	printf("done\n");

//================================
//		test HexagonalTiling
//================================

	printf("test Hexagon\n");
	for (int i=0; i<N; i++)
	{
		Lx = Lmax*((float)rand()/RAND_MAX);
		Ly = Lmax*((float)rand()/RAND_MAX);
//		printf("hex %d x %d\n", Lx, Ly);
		G = HexagonalTiling(Lx, Ly);
		FreeTiling(&G);
	}
	printf("done\n");

//================================
//		test TriangularTiling
//================================

	printf("test Triangular\n");
	for (int i=0; i<N; i++)
	{
		Lx = Lmax*((float)rand()/RAND_MAX);
		Ly = Lmax*((float)rand()/RAND_MAX);
//		printf("triangular: %d x %d\n", Lx, Ly);
		G = TriangularTiling(Lx, Ly);
		FreeTiling(&G);
	}
	printf("done\n");

//================================
//		test SquareTiling
//================================

	printf("test Square\n");
	for (int i=0; i<N; i++)
	{
		Lx = Lmax*((float)rand()/RAND_MAX);
		Ly = Lmax*((float)rand()/RAND_MAX);
//		printf("square: %d x %d\n", Lx, Ly);
		G = SquareTiling(Lx, Ly);
		FreeTiling(&G);
	}
	printf("done\n");

//================================
//		test dual
//================================

	N = 10;
	printf("test dual and dual dual\n");
	int u, v;
	struct tiling Dual;
	struct tiling DDual;
	for (int i=0; i<N; i++)
	{
		Lx = 1+Lmax*((float)rand()/RAND_MAX);
		Ly = 1+Lmax*((float)rand()/RAND_MAX);
//		printf("\nLx = %d, Ly = %d\n", Lx, Ly);
		G = SquareTiling(Lx, Ly);
		u = (int)((G.v-1)*((float)rand()/RAND_MAX));
		v = (int)((G.v-1)*((float)rand()/RAND_MAX));
//		printf("hole 1: u = %d, v = %d\n", u, v);
		CreateOpenRectangleHole(&G, u, v);
		if (G.v != 0)
		{
			u = (int)((G.v-1)*((float)rand()/RAND_MAX));
			v = (int)((G.v-1)*((float)rand()/RAND_MAX));
//			printf("hole 2: u = %d, v = %d\n", u, v);
			CreateOpenRectangleHole(&G, u, v);
		}
		UpdateOpenEdges(&G);
//		Draw(G, "figure", 'b' );
		Dual = DualTiling(G);
//		Draw(Dual, "figuredual", 'r' );
		DDual = DualTiling(Dual);
	}
	printf("done\n");

//================================
//		test perf
//================================


	printf("test Perf\n");
	perf[familySize].noiseRange = malloc(3*sizeof(double));
	perf[familySize].dataFile = malloc(200*sizeof(char));
	for (int i=0; i<N; i++)
	{
		Lx = 1+Lmax*((float)rand()/RAND_MAX);
		Ly = 1+Lmax*((float)rand()/RAND_MAX);
		//squaretiling with a hole

		G = SquareTiling(Lx, Ly);
		H = malloc(G.f*sizeof(int));
		for (int i=0; i<G.f; i++) H[i]=0;
		H[(int)((G.f-1)*((float)rand()/RAND_MAX))] = 1;
		CreateHole(&G, H);
		perf[familySize].noiseRange[0] = .5*((float)rand()/RAND_MAX);
		perf[familySize].noiseRange[1] = 1-.5*((float)rand()/RAND_MAX);
		perf[familySize].noiseRange[2] = .1;
		perf[familySize].dataFile = Performance(G, perf[familySize].noiseRange, 100);
		free(H);
		FreeTiling(&G);

		//squaretiling with two open rectangle hole
		G = SquareTiling(Lx, Ly);
		u = (int)((G.v-1)*((float)rand()/RAND_MAX));
		v = (int)((G.v-1)*((float)rand()/RAND_MAX));
		CreateOpenRectangleHole(&G, u, v);
		if (G.v != 0)
		{
			u = (int)((G.v-1)*((float)rand()/RAND_MAX));
			v = (int)((G.v-1)*((float)rand()/RAND_MAX));
			CreateOpenRectangleHole(&G, u, v);
		}
		UpdateOpenEdges(&G);
		perf[familySize].noiseRange[0] = .5*((float)rand()/RAND_MAX);
		perf[familySize].noiseRange[1] = 1-.5*((float)rand()/RAND_MAX);
		perf[familySize].noiseRange[2] = .1;
		perf[familySize].dataFile = Performance(G, perf[familySize].noiseRange, 100);
		FreeTiling(&G);

		//hex tiling with a hole
 		G = HexagonalTiling(Lx, Ly);
		H = malloc(G.f*sizeof(int));
		for (int i=0; i<G.f; i++) H[i]=0;
		H[(int)((G.f-1)*((float)rand()/RAND_MAX))] = 1;
		CreateHole(&G, H);
		perf[familySize].noiseRange[0] = .5*((float)rand()/RAND_MAX);
		perf[familySize].noiseRange[1] = 1-.5*((float)rand()/RAND_MAX);
		perf[familySize].noiseRange[2] = .1;
		perf[familySize].dataFile = Performance(G, perf[familySize].noiseRange, 100);
		free(H);
		FreeTiling(&G);

		//hex with two open rectangle hole
		G = HexagonalTiling(Lx, Ly);
		u = (int)((G.v-1)*((float)rand()/RAND_MAX));
		v = (int)((G.v-1)*((float)rand()/RAND_MAX));
		CreateOpenRectangleHole(&G, u, v);
		if (G.v != 0)
		{
			u = (int)((G.v-1)*((float)rand()/RAND_MAX));
			v = (int)((G.v-1)*((float)rand()/RAND_MAX));
			CreateOpenRectangleHole(&G, u, v);
		}
		UpdateOpenEdges(&G);
		perf[familySize].noiseRange[0] = .5*((float)rand()/RAND_MAX);
		perf[familySize].noiseRange[1] = 1-.5*((float)rand()/RAND_MAX);
		perf[familySize].noiseRange[2] = .1;
		perf[familySize].dataFile = Performance(G, perf[familySize].noiseRange, 100);
		FreeTiling(&G);

		//triangle tiling with a hole
		G = TriangularTiling(Lx, Ly);
		H = malloc(G.f*sizeof(int));
		for (int i=0; i<G.f; i++) H[i]=0;
		H[(int)((G.f-1)*((float)rand()/RAND_MAX))] = 1;
		CreateHole(&G, H);
		perf[familySize].noiseRange[0] = .5*((float)rand()/RAND_MAX);
		perf[familySize].noiseRange[1] = 1-.5*((float)rand()/RAND_MAX);
		perf[familySize].noiseRange[2] = .1;
		perf[familySize].dataFile = Performance(G, perf[familySize].noiseRange, 100);
		free(H);
		FreeTiling(&G);

		//triangle with two open rectangle hole
		G = TriangularTiling(Lx, Ly);
		u = (int)((G.v-1)*((float)rand()/RAND_MAX));
		v = (int)((G.v-1)*((float)rand()/RAND_MAX));
		CreateOpenRectangleHole(&G, u, v);
		if (G.v != 0)
		{
			u = (int)((G.v-1)*((float)rand()/RAND_MAX));
			v = (int)((G.v-1)*((float)rand()/RAND_MAX));
			CreateOpenRectangleHole(&G, u, v);
		}
		UpdateOpenEdges(&G);
		perf[familySize].noiseRange[0] = .5*((float)rand()/RAND_MAX);
		perf[familySize].noiseRange[1] = 1-.5*((float)rand()/RAND_MAX);
		perf[familySize].noiseRange[2] = .1;
		perf[familySize].dataFile = Performance(G, perf[familySize].noiseRange, 100);
		FreeTiling(&G);
	}

	system("cd results && rm *.txt");

//================================
//		test end
//================================
*/

	G = SquareTiling(7, 7);
	H = malloc(G.f*sizeof(int));
	for (int i=0; i<G.f; i++) H[i]=0;
	H[24] = 1;
	CreateHole(&G, H);
	


	while (strncmp(function, "Quit", 4) != 0)
	{
		//test that G is created
/*		if ((G.v == 0)
		&&(( strncmp(function, "AddFace", 7) == 0 )
		|| ( strncmp(function, "Repeat", 6) == 0 )
		|| ( strncmp(function, "CleanV", 6) == 0 )
		|| ( strncmp(function, "PrintDual", 9) == 0 )
		|| ( strncmp(function, "Print", 5) == 0 )
		|| ( strncmp(function, "DrawDual", 8) == 0 )		
		|| ( strncmp(function, "DrawMin", 7) == 0 )		
		|| ( strncmp(function, "Draw", 4) == 0 )		
		|| ( strncmp(function, "HRO", 3) == 0 )
		|| ( strncmp(function, "HR", 2) == 0 )
		|| ( strncmp(function, "HoleList", 8) == 0 )
		|| ( strncmp(function, "Hole", 4) == 0 )		
		|| ( strncmp(function, "OpenList", 8) == 0 )
		|| ( strncmp(function, "Open", 4) == 0 )
		|| ( strncmp(function, "CloseList", 9) == 0 )
		|| ( strncmp(function, "Close", 5) == 0 )
		|| ( strncmp(function, "DRE", 3) == 0 )
		))
			{
				printf("Create a tiling before\n");
				function[0] = 0;
			}
*/
//TILING CONSTRUCTION

		if (strncmp(function, "RegularTiling", 13) == 0)
		{
			scanf(" %d", size);
			scanf(" %d %d", parameter, &(parameter[1]));
			if (size[0] == 3)
			{
				printf("Creation of a hexagonal tiling\n");
				FreeTiling(&G);
				G = HexagonalTiling(parameter[0], parameter[1]);
			}
			else if (size[0] == 6)
			{
				printf("Creation of a triangular tiling\n");
				FreeTiling(&G);
				G = TriangularTiling(parameter[0], parameter[1]);
			}
			else if (size[0] == 4)
			{
				printf("Creation of a square tiling\n");				
				FreeTiling(&G);
				G = SquareTiling(parameter[0], parameter[1]);
			}
			else
				printf("The degree must be 3, 4 or 6\n");
			function[0] = 0;
		}

		if (strncmp(function, "Grid", 4) == 0)
		{
			printf("Create a grid of vertices:\n");
			scanf(" %d %d", parameter, &(parameter[1]));
			G = Grid(parameter[0], parameter[1]);

			function[0] = 0;
			printf("done\n");
		}

		if (strncmp(function, "AddFace", 7) == 0)
		{
			printf("Add a face from an ordered list of vertices:\n");

			scanf(" %d", size);
			int check = 0;
			if (size[0] < 3) check = 1;
			if (check == 1)
				printf("Failure: the minimum length of a face is 3\n");
			else
			{
				int *vertexList = malloc(size[0]*sizeof(int));
				for (i=0; i<size[0]; i++)
					scanf(" %d", &(vertexList[i]));
			
				//check that these numbers well vertex indices
				for (i=0; i<size[0]; i++)
					if ((vertexList[i] > G.v-1) || (vertexList[i] < 0))
						check = 1;
			
				if (check == 1)
					printf("Failure: These numbers are not vertex indices\n");
				else 
					AddFace(&G, vertexList, size[0]);			
				free(vertexList);
			}
			function[0] = 0;
			printf("done\n");
		}
		
		
		if (strncmp(function, "Repeat", 6) == 0)
		{
			printf("Translate cell:\n");

			float *t = malloc(2*sizeof(float *));
			int *n = malloc(1*sizeof(int *));

			scanf(" %d", n);
			scanf(" %f %f", &(t[0]), &(t[1]));

			printf("%d times t = (%f,%f)\n", n[0], t[0], t[1]);

			FreeTiling(&Gtemp);
			Gtemp = G;
			G = PeriodicTiling(G, t, n[0]);
			
			free(n);
			free(t);

			function[0] = 0;
			printf("done\n");
		}
		
		if (strncmp(function, "CleanV", 5) == 0)
		{
			printf("Remove isolated vertices:\n");
			RemoveIsolatedVertices(&G);
			function[0] = 0;
			printf("done\n");
		}


//PRINTING
		
		if (strncmp(function, "PrintDual", 9) == 0)
		{
			printf("Print the dual graph:\n");
			struct tiling H = DualTiling(G);
			PrintTiling(H);
			function[0] = 0;
			FreeTiling(&H);
			printf("done\n");
		}
		
		if (strncmp(function, "Print", 5) == 0)
		{
			printf("Print the tiling in the terminal:\n");
			PrintTiling(G);
			function[0] = 0;			
			parameter[0] = 0;
			printf("done\n");
		}
		
//DRAWING

		if (strncmp(function, "DrawDualMin", 11) == 0)
		{
			printf("Draw dual tiling:\n");
			struct tiling H = DualTiling(G);
			DrawMin(H, "figure", 'r');
			function[0] = 0;
			FreeTiling(&H);
			printf("done\n");
		}
		
		if (strncmp(function, "DrawDual", 8) == 0)
		{
			printf("Draw indexed dual tiling:\n");
			struct tiling H = DualTiling(G);
			Draw(H, "figure", 'r');
			function[0] = 0;
			FreeTiling(&H);
			printf("done\n");
		}

		if (strncmp(function, "DrawMin", 7) == 0)
		{
			printf("Draw tiling:\n");
			DrawMin(G, "figure", 'b');
			function[0] = 0;
			printf("done\n");
		}
		
		if (strncmp(function, "Draw", 4) == 0)
		{
			printf("Draw indexed tiling:\n");
			Draw(G, "figure", 'b');
			function[0] = 0;
			printf("done\n");
		}
		
//HOLES

		if (strncmp(function, "HoleRectOpen", 12) == 0)
		{
			printf("Create a rectangle hole with open boundaries:\n");
			scanf(" %d %d", parameter, &(parameter[1]));
			CreateOpenRectangleHole(&G, parameter[0], parameter[1]);
			function[0] = 0;
			parameter[0] = 0;
			parameter[1] = 0;
			printf("done\n");
		}
		
		if (strncmp(function, "HoleRect", 8) == 0)
		{
			printf("Create a rectangle hole:\n");
			scanf(" %d %d", parameter, &(parameter[1]));
			CreateRectangleHole(&G, parameter[0], parameter[1]);
			function[0] = 0;			
			parameter[0] = 0;
			parameter[1] = 0;
			printf("done\n");
		}
		
		if (strncmp(function, "HoleList", 8) == 0)
		{
			printf("Create a hole in a list of faces:\n");
			scanf(" %d", size);
	
			H = malloc(G.f*sizeof(int));
			for (int i=0; i<G.f; i++) H[i]=0;
			for (i=0; i<size[0]; i++)
			{
				scanf(" %d", parameter);
				H[parameter[0]] = 1;
			}
			CreateHole(&G, H);
			free(H);

			function[0] = 0;
			parameter[0] = 0;
			printf("done\n");
		}
		
		if (strncmp(function, "Hole", 4) == 0)
		{
			printf("Create a hole in a face:\n");
			scanf(" %d", parameter);
			H = malloc(G.f*sizeof(int));
			for (int i=0; i<G.f; i++) H[i]=0;
			H[parameter[0]] = 1;
			CreateHole(&G, H);
			free(H);

			function[0] = 0;			
			parameter[0] = 0;
			printf("done\n");
		}
		
		
//OPEN BOUNDARIES
		
		if (strncmp(function, "OpenList", 8) == 0)
		{
			scanf(" %d", size);
			printf("Open a list of edges:\n");
	
			OpenE = malloc(G.e*sizeof(int));
			for (int i=0; i<G.e; i++) OpenE[i]=0;
			for (i=0; i<size[0]; i++)
			{
				scanf(" %d", parameter);
				OpenE[parameter[0]] = 1;
			}
			OpenBoundary(&G, OpenE);
			free(OpenE);

			function[0] = 0;
			parameter[0] = 0;
			printf("done\n");
		}
		
		if (strncmp(function, "Open", 4) == 0)
		{
			printf("Open an edge:\n");
			scanf(" %d", parameter);

			OpenE = malloc(G.e*sizeof(int));
			for (i=0; i<G.e; i++) OpenE[i] = 0;
			OpenE[parameter[0]] = 1;
			OpenBoundary(&G, OpenE);
			free(OpenE);

			function[0] = 0;			
			parameter[0] = 0;
			printf("done\n");
		}
		
		if (strncmp(function, "CloseList", 9) == 0)
		{
			scanf(" %d", size);
			printf("Close a list of edges:\n");
	
			CloseE = malloc(G.e*sizeof(int));
			for (int i=0; i<G.e; i++) CloseE[i]=0;
			for (i=0; i<size[0]; i++)
			{
				scanf(" %d", parameter);
				CloseE[parameter[0]] = 1;
			}
			CloseBoundary(&G, CloseE);
			free(CloseE);

			function[0] = 0;
			parameter[0] = 0;
			printf("done\n");
		}
		
		if (strncmp(function, "Close", 5) == 0)
		{
			printf("Close an edge:\n");
			scanf(" %d", parameter);

			CloseE = malloc(G.e*sizeof(int));
			for (i=0; i<G.e; i++) CloseE[i] = 0;
			CloseE[parameter[0]] = 1;
			CloseBoundary(&G, CloseE);
			free(CloseE);

			function[0] = 0;			
			parameter[0] = 0;
			printf("done\n");
		}
		
//PERFORMANCE

		if (strncmp(function, "Code", 4) == 0)
		{
			printf("Code parameters:\n");
			int n = G.e;
			int k = HomDim(G);
			printf("number of physical qubits: n=%d\n", n);
			printf("number of logical qubits: k=%d\n", k);
			printf("overhead: n/k=%.2f\n", ((float) n)/((float) k));
			function[0] = 0;
			printf("done\n");
		}

		if (strncmp(function, "RandomErasure", 13) == 0)
		{
			printf("Draw a random erasure:\n");
			float *p = malloc(sizeof(float));
			scanf(" %f", p);

			//generate random subset of non-open edges
			int *subsetE = malloc(G.e*sizeof(int));
			for (i=0; i<G.e; i++)
				if (G.doE[i] != 1)
					subsetE[i] = ((float)rand()/RAND_MAX < p[0]);
				else subsetE[i] = 0;

			char color = 'b';
			HighlightEdges(G, subsetE, color);

			if (IsCorrectable(G, subsetE))
				printf("correctable\n");
			else printf("not correctable\n");

			free(subsetE);
			function[0] = 0;
			printf("done\n");
		}

		if (strncmp(function, "Report", 5) == 0)
		{
			InitStorage(&perf[familySize]);

			//store tiling name
			perf[familySize].tilingName = malloc(100*sizeof(char));
			strcpy(perf[familySize].tilingName, G.type);
			printf("tiling name: %s\n", perf[familySize].tilingName);

			//store code parameters
			perf[familySize].n = G.e;
			perf[familySize].k = HomDim(G);

			//store measurement weights
			struct tiling H = DualTiling(G);
			perf[familySize].maxWeight = MaxWeight(G, H);
			perf[familySize].ZWeight = FaceMeasurementWeights(G, perf[familySize].maxWeight);
			perf[familySize].XWeight = FaceMeasurementWeights(H, perf[familySize].maxWeight);

			//store performance file name
			perf[familySize].noiseRange = malloc(3*sizeof(double));
			printf("Analyze the performance of the code\n");
			printf("Enter the noise range as: <lower limit>, <upper limit>, <step size>\n>>>>");
			scanf("%lf %lf %lf", &(perf[familySize].noiseRange[0]), &(perf[familySize].noiseRange[1]), &(perf[familySize].noiseRange[2]));
			long stats;
			printf("Number of decoding trials\n>>>>");
			scanf("%ld", &stats);
			int rgen;
			printf("Type of random number generator\n[0] Built in rand() function\n[1] Mersenne primes.\n>>>>");
			scanf("%d", &rgen);
            
            // The noise model is assumed to be of type 1 (see noise.h), the uncorrelated erasure model.
            double **noiseparams = malloc(sizeof(double) * 2);
            noiseparams[0] = malloc(sizeof(double));
            noiseparams[0][0] = 1;
            noiseparams[1] = malloc(sizeof(double) * 3);
            noiseparams[1][0] = perf[familySize].noiseRange[0];
            noiseparams[1][1] = perf[familySize].noiseRange[1];
            noiseparams[1][2] = perf[familySize].noiseRange[2];
            perf[familySize].dataFile = Performance(G, noiseparams, 1, stats, rgen);
            FreeDoubleArray(noiseparams, 2);
			
			familySize++;
			function[0] = 0;
			
			//report
			Report(G, H, &perf[familySize-1]);
			
			FreeTiling(&H);
			function[0] = 0;
		}

		if (strncmp(function, "Compare", 7) == 0)
		{
			for (i=0; i<familySize; i++)
				printf("compare%d: %s\n", i, perf[i].tilingName);

			ComparisonReport(perf, familySize);

			function[0] = 0;
			printf("done\n");
		}

		if (strncmp(function, "Perf", 4) == 0)
		{
			InitStorage(&perf[familySize]);

			//store tiling name
			perf[familySize].tilingName = malloc(100*sizeof(char));
			strcpy(perf[familySize].tilingName, G.type);

			//store code parameters
			perf[familySize].n = G.e;
			perf[familySize].k = HomDim(G);

			//store measurement weights
			struct tiling H = DualTiling(G);
			perf[familySize].maxWeight = MaxWeight(G, H);
			perf[familySize].ZWeight = FaceMeasurementWeights(G, perf[familySize].maxWeight);
			perf[familySize].XWeight = FaceMeasurementWeights(H, perf[familySize].maxWeight);

			//store performance file name
			perf[familySize].noiseRange = malloc(3*sizeof(double));
			printf("Analyze the performance of the code\n");
			printf("Enter the noise range as: <lower limit>, <upper limit>, <step size>\n>>>>");
			scanf("%lf %lf %lf", &(perf[familySize].noiseRange[0]), &(perf[familySize].noiseRange[1]), &(perf[familySize].noiseRange[2]));
			long stats;
			printf("Number of decoding trials\n>>>>");
			scanf("%ld", &stats);
			int rgen;
			printf("Type of random number generator\n[0] Built in rand() function\n[1] Mersenne primes.\n>>>>");
			scanf("%d", &rgen);

			// The noise model is assumed to be of type 1 (see noise.h), the uncorrelated erasure model.
			double **noiseparams = malloc(sizeof(double) * 2);
			noiseparams[0] = malloc(sizeof(double));
			noiseparams[0][0] = 1;
			noiseparams[1] = malloc(sizeof(double) * 3);
			noiseparams[1][0] = perf[familySize].noiseRange[0];
			noiseparams[1][1] = perf[familySize].noiseRange[1];
			noiseparams[1][2] = perf[familySize].noiseRange[2];
			perf[familySize].dataFile = Performance(G, noiseparams, 1, stats, rgen);
			FreeDoubleArray(noiseparams, 2);
			familySize++;
			function[0] = 0;
		}

		if (strncmp(function, "PlotF", 5) == 0)
		{
			for (i=0; i<familySize; i++)
				printf("plot%d: %s\n", i, perf[i].tilingName);

			PlotPerformance(perf, familySize, 1);

			function[0] = 0;
			printf("done\n");
		}

		if (strncmp(function, "Plot", 4) == 0)
		{			
			if (familySize == 0)
				printf("Run Perf before\n");
			else
				PlotPerformance(&(perf[familySize-1]), 1, 1);

			function[0] = 0;
			printf("done\n");
		}

		if (strncmp(function, "ResetPlot", 9) == 0)
		{
			for (i=0; i<familySize; i++)
				FreeStorage(&perf[i]);
			familySize = 0;
			function[0] = 0;
			printf("done\n");
		}

//LOAD AND SAVE

		if (strncmp(function, "Save", 4) == 0)
		{
			printf("Write tiling files:\n");
			char *tilingName = malloc(sizeof(char) * 100);
			scanf("%s", tilingName);
			Save(G, tilingName);
			function[0] = 0;
			parameter[0] = 0;
			printf("done\n");
		}

		if (strncmp(function, "Load", 4) == 0)
		{
			printf("Load Tiling from files:\n");
			char *tilingName = malloc(sizeof(char) * 100);
			scanf("%s", tilingName);
			G = Load(tilingName);
			function[0] = 0;
			parameter[0] = 0;
			printf("done\n");
		}
				
		printf(">>");
		scanf("%s", function);
	}
}
#endif
