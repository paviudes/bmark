#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <dirent.h>
#include <unistd.h>
#include <string.h>
#include "tiling.h"
#include "load.h"

/******************************************************************************************************/

#define DEBUGGING 0

void loadVertexCoords(struct tiling* pG){
  // Load Euclidean cooridinates of the vertices. This is useful for drawing the Tiling.
  // The coordinates are loaded only if the corresponding file exists. If not, nothing is loaded.
  char* coordFile = malloc(sizeof(char*) * 100);
  sprintf(coordFile, "tilings/%s_coordinates.txt", pG->type);

  if(access(coordFile, R_OK) == -1);
  else{
    (pG->Coord) = malloc(sizeof(float*) * (pG->v));

    FILE* coordFid = fopen(coordFile, "r");
	
    int vi = 0;
    for(vi = 0; vi < pG->v; vi ++){
      (pG->Coord)[vi] = malloc(sizeof(float) * 2);
      fscanf(coordFid, "%f %f", &((pG->Coord)[vi][0]), &((pG->Coord)[vi][1]));
    }
	
    fclose(coordFid);
  }

  free(coordFile);
}

void loadEdges(struct tiling* pG){
  char* edgeFile = malloc(sizeof(char*) * 100);
  sprintf(edgeFile, "tilings/%s_edges.txt", pG->type);

  FILE *edgeFid = fopen(edgeFile, "r");
  if (edgeFid == NULL)
    printf("Something went wrong when loading %s: %s\n", edgeFile, strerror(errno));
  else {
    int lNo = 0;
    fscanf(edgeFid, "%d %d", &(pG->v), &(pG->e));
    pG->E = malloc(sizeof(int*) * (pG->e));
		
    for(lNo = 0; lNo < pG->e; lNo ++){
      (pG->E)[lNo] = malloc(sizeof(int) * 2);
      fscanf(edgeFid, "%d %d", &((pG->E)[lNo][0]), &((pG->E)[lNo][1]));
    }
  }

  fclose(edgeFid);

  free(edgeFile);
}

void loadFaces(struct tiling* pG){
  char* faceFile = malloc(sizeof(char*) * 100);
  sprintf(faceFile, "tilings/%s_faces.txt", pG->type);

  FILE *faceFid = fopen(faceFile, "r");
  if (faceFid == NULL)
    printf("Something went wrong when loading %s: %s\n", faceFile, strerror(errno));
  else {

    int lNo = 0, cidx = 0, Fsize = 0;
    fscanf(faceFid, "%d", &(pG->f));
    pG->F = malloc(sizeof(int*) * (pG->f));

    for(lNo = 0; lNo < pG->f; lNo ++){
      fscanf(faceFid, "%d", &Fsize);
      (pG->F)[lNo] = malloc(sizeof(int) * (Fsize + 1));
      (pG->F)[lNo][0] = Fsize;
      for(cidx = 0; cidx < Fsize; cidx ++)
	fscanf(faceFid, "%d", &((pG->F)[lNo][1 + cidx]));
    }
  }
  fclose(faceFid);

  free(faceFile);
}


void loadOpenEdges(struct tiling* pG){
  if(DEBUGGING > 0){
    fprintf(stderr, "Going to read open edges.\n");
  }
  char* doEdgeFile = malloc(sizeof(char*) * 100);
  sprintf(doEdgeFile, "tilings/%s_doEdges.txt", pG->type);

  int ei = 0;
  pG->doE = malloc(sizeof(int) * (pG->e));
  for(ei = 0; ei < pG->e; ei ++)
    (pG->doE)[ei] = 0;

  if(access(doEdgeFile, R_OK) == -1);
  else{
      FILE *doEdgeFid = fopen(doEdgeFile, "r");
      int lNo = 0, nOE = 0;
      fscanf(doEdgeFid, "%d", &nOE);
				
      for(lNo = 0; lNo < nOE; lNo ++){
	       fscanf(doEdgeFid, "%d", &ei);
	       (pG->doE)[ei] = 1;
      }
      fclose(doEdgeFid);
    }

  free(doEdgeFile);
}

void ListTilings(char* dirName){
  // List all files in the directory named dirName
  DIR *dir = opendir(dirName);
  struct dirent *ent;
	
  if(dir != NULL){
    // Print the list of all edge files in the directory
    char* tileName;
    int ci = 0;
    while((ent = readdir(dir)) != NULL){
      if(strlen(ent->d_name) > 10){
      	tileName = ((ent->d_name) + (strlen(ent->d_name) - 10));
      	if(strncmp(tileName, "_edges.txt", 10) == 0){
      	  for(ci = 0; ci < strlen(ent->d_name) - 10; ci ++)
      	    printf("%c", (ent->d_name)[ci]);
      	  printf("\n");
      	}
      }
    }
    closedir(dir);
  }
}

struct tiling Load(char* tileType){
  // Construct a tiling by loading the necessary components from files
	struct tiling G;
  InitTiling(&G);
  // Load a tiling by reading its components from files
  G.length = 0;
  G.type = tileType;
  // Load edges, faces and the open edges
  loadEdges(&G);
  loadFaces(&G);
  loadOpenEdges(&G);
  // Load vertex Coordinates
  loadVertexCoords(&G);
  // Compute the rest of the parameters to define a tiling
  UpdateTiling(&G);
  return G;
}
/******************************************************************************************************/
