#ifndef CELL_H
#define CELL_H


//return a tiling of a single face from its length 
//and the cordinates of its vertices
extern struct tiling ConstructCell(float **vertexCoord, int length);

//
extern struct tiling PeriodicTiling(struct tiling C, float *t, int n);


#endif	/* CELL_H */
