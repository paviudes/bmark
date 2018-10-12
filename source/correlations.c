#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tiling.h"
#include "arrays.h"
#include "rgens.h"
#include "rbim.h"
#include "simulation.h"
#include "correlations.h"

void BubbleGrow(int root, int inEdge, int dist, int target, struct simulation *pSim, struct tiling *pG){
	// Explore erased edges that lie in a ball of a given radius around a given edge.
	if (dist <= target){
		(pSim->erasure)[1][inEdge] = 1;
		int vi;
		for (vi = 1; vi <= (pG->VToN)[root][0]; vi ++){
			if ((pG->doE)[(pG->VToE)[root][vi]] == 0)
				BubbleGrow((pG->VToN)[root][vi], (pG->VToE)[root][vi], dist + 1, target, pSim, pG);
		}
	}
}

void BallErasure(struct simulation *pSim, struct tiling *pG, double center, double *cumulative){
	// Erasing edges inside a ball of fixed radius around an edge
	int ei, radius;
	for (ei = 0; ei < pG->e; ei ++){
		if ((pG->doE)[ei] == 0){
			radius = (int) Random("cumulative", (int) cumulative[0], cumulative + 1);
			if (Random("uniform", 0, NULL) < center){
				BubbleGrow((pG->E)[ei][0], ei, 1, radius, pSim, pG);
				BubbleGrow((pG->E)[ei][1], ei, 1, radius, pSim, pG);
			}
		}
	}
}

double BallInfer(double *available, int toinfer){
	// Infer one of the physical parameters of the ball erasures model, when provided with all other parameters and the empirical noise rate.
	double infered = 0;
	if (toinfer == 0){
		// we need to infer q, the probability of creating the center of a ball. The data which is available in the place of q is the required empirical noise rate, p.
		// we have A q = p where A is the area of the ball. A is given by A = r^2 + 4 r + 2 where r is the radius of the ball, the second parameter of the noise mdoel.
		infered = available[0]/(double) (pow(available[1], 2) + 4 * available[1] + 2);
	}
	if (toinfer == 1){
		// we need to infer the radius of the ball which gives the desired p.
		// we have p = q (r^2 + 4r + 2). Hence we must solve for the equation r^2 + 4r + 2 - p/q = 0. If no solution is available, we will simply take r = 1.
		// Solution existence check: in a quadratic equation, ax^2 + bx + c = 0, if we have b^2 - 4 ac < 0, there is no real solution. This condition siplifies to 8 + 4p/q.
		infered = Max(sqrt(8 + 4 * available[1]/available[0])/(double) 2 - 2, 1);
	}
	return infered;
}




int Step(int inEdge, struct tiling **pGs, struct simulation *pSim, int *visited){
	// select the next edge to progress the random walk.
	// the next edge is one out of the set of all edges that either share a vertex or a face with the incoming edge.
	int dual, endpoint, ei, pi, incoming, nout = 0, next, pass;
	for (dual = 0; dual < 2; dual ++){
		if (dual == 0)
			incoming = inEdge;
		else
			incoming = (pSim->edgeOrderings)[0][inEdge];
		for (pi = 0; pi < 2; pi ++){
			endpoint = (pGs[dual]->E)[incoming][pi];
			for (ei = 1; ei <= (pGs[dual]->VToE)[endpoint][0]; ei ++){
				next = (pGs[dual]->VToE)[endpoint][ei];
				if (next != incoming){
					if ((pGs[dual]->doE)[next] == 0){
						pass = 0;
						if (visited != NULL){
							if (visited[next] == 0)
								pass = 1;
						}
						else
							pass = 1;
						if (pass == 1)
							nout ++;
					}
				}
			}
		}
	}
	if (nout == 0)
		return -1;

	int activated, nextstep = 0;
	activated = (int) (Random("uniform", 0, NULL) * nout);
	nout = 0;
	for (dual = 0; dual < 2; dual ++){
		if (dual == 0)
			incoming = inEdge;
		else
			incoming = (pSim->edgeOrderings)[0][inEdge];
		for (pi = 0; pi < 2; pi ++){
			endpoint = (pGs[dual]->E)[incoming][pi];
			for (ei = 1; ei <= (pGs[dual]->VToE)[endpoint][0]; ei ++){
				next = (pGs[dual]->VToE)[endpoint][ei];
				if (next != incoming){
					if ((pGs[dual]->doE)[next] == 0){
						pass = 0;
						if (visited != NULL){
							if (visited[next] == 0)
								pass = 1;
						}
						else
							pass = 1;

						if (pass == 1){
							if (nout == activated){
								if (dual == 0)
									nextstep = next;
								else
									nextstep = (pSim->edgeOrderings)[1][next];
							}
							nout ++;
						}
					}
				}
			}
		}
	}
	return nextstep;
}

void WalkErasure(struct simulation *pSim, struct tiling **pGs, double center, double *cumulative){
	// Perform a random walk starting at some edge. The length of the random walk is a random variate with given cumulative distribution.
	int length, nextstep = 0, ei, w;
	for (ei = 0; ei < pGs[0]->e; ei ++){
		if ((pGs[0]->doE)[ei] == 0){
			if (Random("uniform", 0, NULL) < center){
				length = (int) Random("cumulative", (int) cumulative[0], cumulative + 1);
				nextstep = ei;
				for (w = 0; w < length; w ++){
					if ((pGs[0]->doE)[nextstep] == 0){
						(pSim->erasure)[1][nextstep] = 1;
					}
					nextstep = Step(nextstep, pGs, pSim, NULL);
				}
			}
		}
	}
}

void SpreadErasure(struct simulation *pSim, struct tiling **pGs, double center, double *cumulative){
	// Grow a cluster of erased edges starting from a random edge. The size (area) of the cluster is a random variate with given cumulative distribution.
	int *active = malloc(sizeof(int) * ((int) cumulative[0] + 1));
	int *visited = malloc(sizeof(int) * pGs[0]->e);
	int activated, area, nextstep = 0, ei, w, randidx;
	for (ei = 0; ei < pGs[0]->e; ei ++){
		if ((pGs[0]->doE)[ei] == 0){
			if (Random("uniform", 0, NULL) < center){
				area = (int) Random("cumulative", (int) cumulative[0], cumulative + 1);
				FillIntArrayWithConstants(active, ((int) cumulative[0] + 1), 0);
				FillIntArrayWithConstants(visited, pGs[0]->e, 0);
				nextstep = ei;
				for (w = 0; w < area; w ++){
					if ((pGs[0]->doE)[nextstep] == 0){
						(pSim->erasure)[1][nextstep] = 1;
						visited[nextstep] = 1;
						active[++ active[0]] = nextstep;
					}
					do{
						randidx = (int) (Random("uniform", 0, NULL) * active[0]) + 1;
						activated = active[randidx];
						nextstep = Step(activated, pGs, pSim, visited);
						if (nextstep == -1){
							Swap(active, randidx, active[0]);
							active[0] --;
							if (active[0] == 0)
								break;
						}
					}while(nextstep == -1);
				}
			}
		}
	}
	free(active);
	free(visited);
}


double WalkInfer(double *available, int toinfer){
	// Infer one of the Physical parameters of the Walk erasure model, when provided with the other and an empirical noise rate.
	double infered = 0;
	if (toinfer == 0){
		// we need to infer q, the probability of choosing a starting point for the Random walk
		// p = q L where L is the length of the walk
		infered = available[0]/(double) available[1];
	}
	if (toinfer == 1){
		// we need to infer L where p = q L
		infered = available[1]/(double) available[0];
	}
	return infered;
}


void RBIMErasure(struct simulation *pSim, struct tiling *pG, struct RandomBondIsingModel *prbim, int ncool){
	// Design an erasure pattern that corresponds to a low energy configuration of an Ising model.
	Cool(prbim, ncool);
	// The alignment that is adopted by the minority of spins corresponds to erasures.
	int ei;
	for (ei = 0; ei < pG->e; ei ++)
		if ((pG->doE)[ei] == 0)
			(pSim->erasure)[1][ei] = (prbim->config)[ei];
}

int ExploreErasure(struct tiling **pGs, struct simulation *pSim, int edge, int *visited){
	// Given an edge, explore all the edges that either share a vertex or a face with the given edge.
	// from both the endpoints of the given edge, explore all outgoing unvisited edges.
	// from both the endpoints of the dual-edge, explore all the outgoing unvisited edges.
	int dual, incoming, ep, endpoint, ei, next;
	visited[edge] = 1;
	for (dual = 0; dual < 2; dual ++){
		if (dual == 1)
			incoming = (pSim->edgeOrderings)[0][edge];
		else
			incoming = edge;
		for (ep = 0; ep < 2; ep ++){
			endpoint = (pGs[dual]->E)[incoming][ep];
			for (ei = 1; ei <= (pGs[dual]->VToE)[endpoint][0]; ei ++){
				if (dual == 0)
					next = (pGs[dual]->VToE)[endpoint][ei];
				else
					next = (pSim->edgeOrderings)[1][(pGs[dual]->VToE)[endpoint][ei]];
				if ((pSim->erasure)[1][next] == 1)
					if (visited[next] == 0)
						return (1 + ExploreErasure(pGs, pSim, next, visited));
			}
		}
	}
	return 1;
}


void Cluster(struct tiling **pGs, struct simulation *pSim, double *clusterinfo){
	// measure the cluster size of erasure patterns
	// a cluster is a connected component where two edges are connected if they share a vertex either on the primal or on the dual lattice
	int *visited = malloc(sizeof(int) * pGs[0]->e);
	FillIntArrayWithConstants(visited, pGs[0]->e, 0);
	FillDoubleArrayWithConstants(clusterinfo, 3, 0);
	int ei, csize;
	for (ei = 0; ei < pGs[0]->e; ei ++){
		if ((pSim->erasure)[1][ei] == 1){
			if (visited[ei] == 0){
				csize = (double) ExploreErasure(pGs, pSim, ei, visited);
				clusterinfo[0] ++;
				if (clusterinfo[1] < csize)
					clusterinfo[1] = csize;
				clusterinfo[2] += csize;
			}
		}
	}
	free(visited);
	if (clusterinfo[0] > 0)
		clusterinfo[2] = clusterinfo[2]/clusterinfo[0];
}
