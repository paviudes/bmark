#ifndef HOLES_H
#define HOLES_H

//create a hole by replacing a faces by holes
//the second input is the indicator vector of the set of faces
extern void CreateHole(struct tiling *, int *);

//create a rectangle hole
//input: the tiling and the indices of 2 vertices i and j
//output: create a rectangle hole in the rectangle having
//vi and vj as opposite corners.
extern void CreateRectangleHole(struct tiling *, int, int);

//create a rectangle hole with open boundaries
//input: the tiling and the indices of 2 vertices i and j
//output: create a rectangle hole in the rectangle having
//vi and vj as opposite corners.
//Open the edges of the rectangle
extern void CreateOpenRectangleHole(struct tiling *, int, int);

//Open a subset of edges
//These egdes are given by their indicator vector
extern void OpenBoundary(struct tiling *, int *);

//Close a subset of edges
//These egdes are given by their indicator vector
extern void CloseBoundary(struct tiling *, int *);

#endif	/* HOLES_H */
