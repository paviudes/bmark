#ifndef DUAL_H
#define DUAL_H

//return the dual tiling
extern struct tiling DualTiling(struct tiling);

//return the array EToDE of indices of the dual edges.
//EToDE[k] is the index of the dual edge of the k-th
//edge ek of the tiling.
//EToDE[k] = -1, if ek is not an open edge. This indicates
//that ek has no dual edge.
extern int *EToDualE(struct tiling G);

#endif	/* DUAL_H */
