#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include "tiling.h"
#include "save.h"

/***********************************************************************************************/

void DeleteFile(char* fname){
	// If a file with fname exists and has permissions to be deleted, it will be deleted.
	char* systemCommand = malloc(sizeof(char)*(100 + sizeof(fname)));
	sprintf(systemCommand, "rm -f %s", fname);
	system(systemCommand);
	free(systemCommand);
}

void SaveEdges(struct tiling G){
	// Write the edge set to a file
	char* edgeFile = malloc(sizeof(char*) * 100);
	sprintf(edgeFile, "saved_tilings/%s_edges.txt", G.type);

	DeleteFile(edgeFile);
	
	FILE *edgeFid = fopen(edgeFile, "w");	
	int ei = 0;
	fprintf(edgeFid, "%d %d\n", G.v, G.e);
	for(ei = 0; ei < G.e; ei ++)
		fprintf(edgeFid, "%d %d\n", (G.E)[ei][0], (G.E)[ei][1]);

	fclose(edgeFid);

	free(edgeFile);
}

void SaveOpenEdges(struct tiling G){
	// Write the edge set to a file
	char* doEdgeFile = malloc(sizeof(char*) * 100);
	sprintf(doEdgeFile, "saved_tilings/%s_doEdges.txt", G.type);

	DeleteFile(doEdgeFile);
	
	FILE *doEdgeFid = fopen(doEdgeFile, "w");
		
	int ei = 0;
	int nOE = 0;
	for(ei = 0; ei < G.e; ei ++)
		nOE = nOE + G.doE[ei];

	fprintf(doEdgeFid, "%d\n", nOE);
	for(ei = 0; ei < G.e; ei ++)
		if((G.doE)[ei] == 1)
			fprintf(doEdgeFid, "%d\n", ei);

	fclose(doEdgeFid);

	free(doEdgeFile);
}

void SaveFaces(struct tiling G){
	// Write the edge set to a file
	char* faceFile = malloc(sizeof(char*) * 100);
	sprintf(faceFile, "saved_tilings/%s_faces.txt", G.type);

	DeleteFile(faceFile);
	FILE *faceFid = fopen(faceFile, "w");

	int ei = 0, fi = 0;
	fprintf(faceFid, "%d\n", G.f);
	for(fi = 0; fi < G.f; fi ++){
		for(ei = 0; ei <= (G.F)[fi][0]; ei ++)
			fprintf(faceFid, "%d ", (G.F)[fi][ei]);
	
		fprintf(faceFid, "\n");
	}

	fclose(faceFid);

	free(faceFile);
}

void SaveCoord(struct tiling G){
	// Write the vertex coordinates to a file
	char* coordFile = malloc(sizeof(char*) * 100);
	sprintf(coordFile, "saved_tilings/%s_coordinates.txt", G.type);

	DeleteFile(coordFile);
	FILE *coordFid = fopen(coordFile, "w");

	int vi = 0;
	for(vi = 0; vi < G.v; vi ++)
		fprintf(coordFid, "%f %f\n", G.Coord[vi][0], G.Coord[vi][1]);

	fclose(coordFid);

	free(coordFile);
}

void Save(struct tiling G, char* fname)
{
	// Change the tiling type to the filename
	G.type = fname;
	
	// Save a tiling to a file: save the edges, face and the open edges
	SaveEdges(G);
	SaveFaces(G);
	SaveOpenEdges(G);

	// Save the vertex coordinates to a file
	if(G.Coord != NULL)
		SaveCoord(G);
}

/***********************************************************************************************/
