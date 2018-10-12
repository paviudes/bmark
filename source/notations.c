#include <stdio.h>
#include <stdlib.h>
#include "memory.h"
#include "notations.h"

void NoiseParameters(struct representation *repr){
	// Assign names for the different noise parameters
	int v, degfreedom[6] = {1, 1, 2, 2, 2, 4}, derived[6] = {0, 0, 4, 4, 4, 7};
	repr->nvars = degfreedom[repr->model];
	repr->nprops = derived[repr->model];
	repr->varnames = malloc(sizeof(char *) * repr->nvars);
	repr->propnames = malloc(sizeof(char *) * repr->nprops);
	repr->mname = malloc(sizeof(char *) * 100);
	for (v = 0; v < repr->nvars; v ++)
		(repr->varnames)[v] = malloc(sizeof(char) * 100);
	for (v = 0; v < repr->nprops; v ++)
		(repr->propnames)[v] = malloc(sizeof(char) * 100);

	if (repr->model == 1){
		sprintf(repr->mname, "uncorr");
		sprintf((repr->varnames)[0], "p");
	}
	else if (repr->model == 2){
		sprintf(repr->mname, "ball");
		sprintf((repr->varnames)[0], "q");
		sprintf((repr->varnames)[1], "r");
		sprintf((repr->propnames)[0], "<p>");
		sprintf((repr->propnames)[1], "<nc>");
		sprintf((repr->propnames)[2], "<cmax>");
		sprintf((repr->propnames)[3], "<c>");
	}
	else if (repr->model == 3){
		sprintf(repr->mname, "walk");
		sprintf((repr->varnames)[0], "q");
		sprintf((repr->varnames)[1], "d");
		sprintf((repr->propnames)[0], "<p>");
		sprintf((repr->propnames)[1], "<nc>");
		sprintf((repr->propnames)[2], "<cmax>");
		sprintf((repr->propnames)[3], "<c>");
	}
	else if (repr->model == 4){
		sprintf(repr->mname, "spread");
		sprintf((repr->varnames)[0], "q");
		sprintf((repr->varnames)[1], "a");
		sprintf((repr->propnames)[0], "<p>");
		sprintf((repr->propnames)[1], "<nc>");
		sprintf((repr->propnames)[2], "<cmax>");
		sprintf((repr->propnames)[3], "<c>");
	}
	else if (repr->model == 5){
		sprintf(repr->mname, "ising");
		sprintf((repr->varnames)[0], "J");
		sprintf((repr->varnames)[1], "B");
		sprintf((repr->varnames)[2], "T");
		sprintf((repr->varnames)[3], "l");
		sprintf((repr->propnames)[0], "<p>");
		sprintf((repr->propnames)[1], "<nc>");
		sprintf((repr->propnames)[2], "<cmax>");
		sprintf((repr->propnames)[3], "<c>");
		sprintf((repr->propnames)[4], "<m>");
		sprintf((repr->propnames)[5], "<u>");
		sprintf((repr->propnames)[6], "<sisj>");
	}
	else{
		sprintf(repr->mname, "unknown");
		sprintf((repr->varnames)[0], "X");
	}
}


void FamilyNames(struct representation *repr){
	// Give names to code families
	repr->fname = malloc(sizeof(char) * 100);
	if (repr->family == 1)
		sprintf(repr->fname, "Toric");
	else if(repr->family == 2)
		sprintf(repr->fname, "Hyperbolic_(6,6)_trivalent");
	else if(repr->family == 3)
		sprintf(repr->fname, "Hyperbolic_(6,6)_selfdual");
	else if(repr->family == 4)
		sprintf(repr->fname, "Bravyi_surface");
	else
		sprintf(repr->fname, "Unknown");
}

void ErrorTypes(struct representation *repr){
	// Assign names to error types
	repr->enames = malloc(sizeof(char *) * 3);
	int et;
	for (et = 0; et < 3; et ++)
		repr->enames[et] = malloc(sizeof(char) * 5);
	sprintf((repr->enames)[0], "X,Z");
	sprintf((repr->enames)[1], "X");
	sprintf((repr->enames)[2], "Z");
}

void FixRepresentation(struct representation *repr, int model, int family){
	// Initialize the representation data structure
	repr->model = model;
	repr->family = family;
	repr->nvars = 0;
	repr->varnames = NULL;
	repr->mname = NULL;
	ErrorTypes(repr);
	NoiseParameters(repr);
	// fprintf(stderr, "Noise parameter %s\n", (repr->varnames)[0]);
	FamilyNames(repr);
	// fprintf(stderr, "Family %d: %s\n", repr->family, repr->fname);
}

void FreeRep(struct representation *repr){
	// Free a representation structure
	free(repr->mname);
	free(repr->fname);
	FreeCharArray(repr->enames, 3);
	FreeCharArray(repr->varnames, repr->nvars);
	FreeCharArray(repr->propnames, repr->nprops);
}


void TilingFile(char* tilingFile, int family, int length){
	if(family == 1)
		sprintf(tilingFile, "%d_toric_square", length);
	else if(family == 2)
		sprintf(tilingFile, "6_6_%d_trivalent_hyp", length);
	else if(family == 3)
		sprintf(tilingFile, "6_%d_selfdual", length);
	else if (family == 4)
		sprintf(tilingFile, "bravyi_%d", length);
	else
		sprintf(tilingFile, "family%d_length%d", family, length);
}
