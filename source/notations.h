#ifndef NOTATIONS_H
#define NOTATIONS_H

struct representation
{
	/* Data required to verbally express a noise model. See noise.h for details. */
	// Index for the noise model.
	int model;
	// number of derived properties of the noise model
	int nprops;
	// Index for the family
	int family;
	// Name of the noise model.
	char *mname;
	// Number of parameters to describe a noise process.
	int nvars;
	// Names/symbols of the parameters.
	char **varnames;
	// Names/symbols for derived parameters of the noise models
	char **propnames;
	// Names to code families
	char *fname;
	// Names to types of errors -- X,Z and X and Z.
	char **enames;
};

// Initialize the representation data structure
extern void FixRepresentation(struct representation *repr, int model, int family);

// Free a representation structure
extern void FreeRep(struct representation *repr);

/*
	Assign a name to a tiling, that is derived from an index for a family and the length of the tiling.
	For a few standard families of tilings, we have the given names.
	The edges, faces, open edges and vertex-coordinates are stored in the files
			./../saved_tilings/tiling_folder/<tilingName>_edges.txt,
			./../saved_tilings/tiling_folder/<tilingName>_faces.txt,
			./../saved_tilings/tiling_folder/<tilingName>_doEdges.txt,
			./../saved_tilings/tiling_folder/<tilingName>_coordinates.txt,
	respectively.
*/
extern void TilingFile(char* tilingFile, int family, int length);

#endif /* NOTATIONS_H */
