#ifndef LOAD_H
#define LOAD_H

extern struct tiling Load(char* tileType);
/*
	Loading a tiling from a file.
	
	Input: name of a tiling
	Output: Tiling
	Algorithm: Read the list of edges, faces and open-edges from suitable files (derived from the name of the tiling).
			   From this compute the various other parameters that are used to define a Tiling.
*/

extern void ListTilings(char* dirName);
/*
	print the list of all files in a directory.
*/

#endif	/* LOAD_H */
