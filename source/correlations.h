#ifndef CORRELATIONS_H
#define CORRELATIONS_H

/* 
	Design a correlated erasure pattern on n-bits in the following way.
	At every edge of the lattice, erase the edge with probability p and if the edge is erased, erase all the edges in a radius r around the erased edge.
	Algorithm:
		Function BubbleGrow(...) -- Erasing edges inside a ball of fixed radius around an edge
		Start from a vertex v0
		Recursion:
			1. For every neighbour of v0 that is not visited, erase the corresponding outgoing edge
			2. Call BubbleGrow on the endpoint of the erased edge, to explore erased edges that lie in a ball of a given radius.

*/
extern void BallErasure(struct simulation *pSim, struct tiling *pG, double center, double *cumulative);

// Infer one of the physical parameters of the ball erasures model, when provided with all other parameters and the empirical noise rate.
extern double BallInfer(double *available, int toinfer);

/*
	Design a correlated model of erasure where a random path originating from an erased edge, is also erased.
	The length of the random walk as well as the probability of choosing a root edge for the walk, are provided in "specs".
	Algorithm:
		Function Walk(...) -- Starting from a given edge, walk in any random direction, on unerased edges, for a fixed number of steps.
*/
extern void WalkErasure(struct simulation *pSim, struct tiling **pGs, double center, double *cumulative);

// Infer one of the Physical parameters of the Walk erasure model, when provided with the other and an empirical noise rate.
extern double WalkInfer(double *available, int toinfer);


/*
	Design a correlated erasure model where the erasure pattern is a subgraph that grows in an arbitrary direction
*/
extern void SpreadErasure(struct simulation *pSim, struct tiling **pGs, double center, double *cumulative);

/*
	Design an erasure pattern that corresponds to a low energy configuration of an Ising model.
	The Ising model is intialized with the lattice parameters of the code and few parameters of the noise model.
	Here, we will do several rounds of Markov chain Montecarlo (or Metropolis) to derive a configuration of lower energy than what is present already.
	The new configuration will be accepted as the new erasure pattern
*/
extern void RBIMErasure(struct simulation *pSim, struct tiling *pG, struct RandomBondIsingModel *prbim, int ncool);

// measure the cluster size of erasure patterns
// a cluster is a connected component where two edges are connected if they share a vertex either on the primal or on the dual lattice
extern void Cluster(struct tiling **pGs, struct simulation *pSim, double *clusterinfo);

#endif /* CORRELATIONS_H */
