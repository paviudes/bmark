#ifndef HOMOLOGICAL_H
#define HOMOLOGICAL_H

// Explore a connected component of the induced erasure pattern on the dual graph having only non-open vertices
int ExploreNOVEdgeOrdering(int vertex, int* visited, int* erasure, struct tiling* pG, int* edgeOrdering);

// Explore a connected component of the erasure pattern on the dual graph having only non-open vertices
int ExploreNOV(int vertex, int* visited, int* erasure, struct tiling* pG);

// Explore a connected component of the erasure pattern having only non-closed edges
int ExploreNCE(int vertex, int* visited, struct tiling* pG);

/*
	Compute the connected components for a
		1. Graph, such that no connected component has a closed edge --- restriction: 1
		2. Subgraph induced by the erasure, such that no connected component has an open vertex --- restriction: 2

	Use erasure = NULL, edgeOrdering = NULL wherever applicable.
*/
extern int ConnectedComponents(struct tiling* pG, int* edgeOrdering, int* erasure, int restriction);

/*
	Determine the type of error that resulted in a decoding failure.
		
		Output: 0 -- No decoding failure
				1 -- Decoding failed for an error that was of Y type
				2 -- Decoding failed for an error that was of X-type
				3 -- Decoding failed for an error that was of Z-type

		Procedure:
			1. Remove any pending edges. 
			2. If anything remains in the erasure pattern, it is the cycle that is covered by it. Hence every non-open vertex has an even degree.
			3. Now that the erasure pattern is a cycle itself, we must determine if the cycle is trivial or not. This is done by compute the homological dimension.
*/
extern int LogicalType(struct tiling **pGs, struct simulation *pSim);

/*
	Walk on the erasure subgraph to explore all vertices with degree 1.
	If the vertex is not open and it has degree one, walk on it, to its neighbours.
*/
int DegreeOneWalk(int vertex, int* erasureDegrees, struct tiling* pT, int* erasure, int ewt);

/*
	Remove pending edges from an erasure subgraph.
		An edge e0 is called pending if there is some stabilizer generator G0 such that the only erased edge in the support of G0 is e0.
		The state of the qubit corresponding to the edge e0 can be recovered exactly as the eigenstate of the corresponding stabilizer generator.
		So, we can ignore these portions of the erasure pattern. Hence we say that we remove pending edges.
	Note: 
		VToN must be allocated before calling this function.
*/
int RemovePendingEdges(struct tiling *pT, int *erasure, int ewt);

/*
	Compute the Homological dimension of the erasure pattern.
	The homological dimension of the erasure pattern, denoted by h, is given by

		h = |E| - |V| + k1 - k2 + k3

	where |E| = number of edges in the erasure subgraph
		  |V| = number of non-open vertices
		  k1 = number of connected components in the erasure subgraph that have no open vertices
		  k2 = number of connected components in the dual of erasure subgraph that have no open vertices
		  k3 = number of connected components in the primal graph that have no closed edges

	isDual = 0 for testing on the primal graph and 1 for testing on the dual graph.
*/
int HomologicalDimension(struct tiling **pGs, struct simulation *pSim, int isDual);

#endif /* HOMOLOGICAL_H */
