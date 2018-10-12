#ifndef GRID_H
#define GRID_H

//create a square tiling of size lx x ly
extern struct tiling Grid(int lx, int ly);

//Add a new face to the tiling given as an ordered list of vertices.
//length contains the length of the list.
extern void AddFace(struct tiling *G, int *vertexList, int length);

//
extern void RemoveIsolatedVertices(struct tiling *G);

#endif	/* GRID_H */
