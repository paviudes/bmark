#include <stdio.h>
#include <stdlib.h>
#include "tiling.h"
#include "arrays.h"
#include "simulation.h"
#include "homological.h"

/*******************************************************************************************************/

int ExploreNOVEdgeOrdering(int vertex, int* visited, int* erasure, struct tiling* pG, int* edgeOrdering){
	// Count a finite contribution from the connected components having only non-open vertices
	int ei = 0, outEdge = 0, endPoint = 0, contribution = 1;
	visited[vertex] = 1;
	if(pG->doV[vertex])
		contribution = 0;
	for(ei = 1; ei <= pG->VToE[vertex][0]; ei ++){
		outEdge = pG->VToE[vertex][ei];
		endPoint = pG->VToN[vertex][ei];
		if(pG->doE[outEdge] == 0 && erasure[edgeOrdering[outEdge]] == 0){
			if(!visited[endPoint]){
				contribution = contribution * ExploreNOVEdgeOrdering(endPoint, visited, erasure, pG, edgeOrdering);
			}
		}
	}
	return contribution;
}


int ExploreNOV(int vertex, int* visited, int* erasure, struct tiling* pG){
	// Count a finite contribution from the connected components having only non-open vertices
	int ei = 0, outEdge = 0, endPoint = 0, contribution = 1;
	visited[vertex] = 1;
	if(pG->doV[vertex])
		contribution = 0;
	for(ei = 1; ei <= pG->VToE[vertex][0]; ei ++){
		outEdge = pG->VToE[vertex][ei];
		endPoint = pG->VToN[vertex][ei];
		if(pG->doE[outEdge] == 0 && erasure[outEdge] == 1){
			if(!visited[endPoint]){
				contribution = contribution * ExploreNOV(endPoint, visited, erasure, pG);
			}
		}
	}
	return contribution;
}


int ExploreNCE(int vertex, int* visited, struct tiling* pG){
	// Count a finite contribution from the connected components having only non-closed edges
	int ei = 0, outEdge = 0, endPoint = 0, contribution = 1;
	visited[vertex] = 1;
	for(ei = 1; ei <= pG->VToE[vertex][0]; ei ++){
		outEdge = pG->VToE[vertex][ei];
		endPoint = pG->VToN[vertex][ei];
		if(!visited[endPoint]){
			if(pG->dE[outEdge] && !pG->doE[outEdge]){
				contribution = 0;
			}
			contribution = contribution * ExploreNCE(endPoint, visited, pG);
		}
	}
	return contribution;
}


int ConnectedComponents(struct tiling* pG, int* edgeOrdering, int* erasure, int restriction){
	// Compute the connected components in a graph.
	int nComps = 0, vi = 0;
	int* visited = malloc(sizeof(int) * (pG->v));
	FillIntArrayWithConstants(visited, pG->v, 0);
	for(vi = 0; vi < pG->v; vi ++){
		if(visited[vi] == 0){
			if(restriction == 1)
				// Part 1
				nComps += ExploreNCE(vi, visited, pG);
			else if(restriction == 2)
				// Part 2
				nComps += ExploreNOV(vi, visited, erasure, pG);
			else if(restriction == 3)
				nComps += ExploreNOVEdgeOrdering(vi, visited, erasure, pG, edgeOrdering);
			else;
		}
	}
	free(visited);
	return nComps;
}


int DegreeOneWalk(int vertex, int* erasureDegrees, struct tiling* pT, int* erasure, int ewt){
	// Walk on the erasure subgraph to explore all vertices with degree 1.
	if((pT->doV[vertex] == 0) && (erasureDegrees[vertex] == 1)){
		erasureDegrees[vertex] --;
		int ei = 0, endPoint = 0, outEdge = 0, outErasedEdge = 0, outErasedEdgeIndex = 0;
		for(ei = 1; ei <= pT->VToE[vertex][0]; ei ++){
			outEdge = (pT->VToE)[vertex][ei];
			if(pT->doE[outEdge] == 0 && erasure[outEdge] == 1){
				outErasedEdge = outEdge;
				outErasedEdgeIndex = ei;
			}
		}
		erasure[outErasedEdge] = 0;
		ewt --;
		endPoint = (pT->VToN)[vertex][outErasedEdgeIndex];
		erasureDegrees[endPoint] --;
		return DegreeOneWalk(endPoint, erasureDegrees, pT, erasure, ewt);
	}
	return ewt;
}


int RemovePendingEdges(struct tiling *pT, int *erasure, int ewt){
	// Remove pending edges from an erasure subgraph.
	int vi = 0, ei = 0;
	int* erasureDegrees = malloc(sizeof(int) * pT->v);
	int* startWalk = malloc(sizeof(int) * pT->v);
	for(vi = 0; vi < pT->v; vi ++){
		erasureDegrees[vi] = 0;
		startWalk[vi] = 0;
	}
	for(vi = 0; vi < pT->v; vi ++){
		if(!pT->doV[vi]){
			for(ei = 1; ei <= pT->VToE[vi][0]; ei ++)
				erasureDegrees[vi] += erasure[(pT->VToE)[vi][ei]];
			startWalk[vi] = (erasureDegrees[vi] == 1);
		}
	}
	for(vi = 0; vi < pT->v; vi ++)
		if(startWalk[vi])
			ewt = DegreeOneWalk(vi, erasureDegrees, pT, erasure, ewt);
	
	free(startWalk);
	free(erasureDegrees);
	return ewt;
}


int HomologicalDimension(struct tiling **pGs, struct simulation *pSim, int isDual){
	// Compute the Homological dimension of the erasure pattern.
	int k1 = ConnectedComponents(pGs[isDual], NULL, (pSim->erasure)[1 + isDual], 2);
	int k2 = ConnectedComponents(pGs[1 + isDual], (pSim->edgeOrderings)[1 + isDual], (pSim->erasure)[1 + isDual], 3);
	int hdim = (pSim->erasure)[0][isDual] - (pSim->nOV)[isDual] + k1 - k2 + (pSim->comps)[isDual];
	return hdim;
}

int LogicalType(struct tiling **pGs, struct simulation *pSim){
	// Determine the type of logical error covered by the erasure pattern.
	int isLogZ = 0, isLogX = 0;
	if ((pSim->type == 0) || (pSim->type == 1)){
		(pSim->erasure)[0][0] = RemovePendingEdges(pGs[0], (pSim->erasure)[1], (pSim->erasure)[0][0]);
		if ((pSim->erasure)[0][0] > 0)
			isLogX = HomologicalDimension(pGs, pSim, 0);

	}
	if ((pSim->type == 0) || (pSim->type == 2)){
		(pSim->erasure)[0][1] = RemovePendingEdges(pGs[1], (pSim->erasure)[2], (pSim->erasure)[0][1]);
		if ((pSim->erasure)[0][1] > 0)
			isLogZ = HomologicalDimension(pGs, pSim, 1);
	}
	
	if (isLogX * isLogZ > 0)
		return 1;
	else if(isLogX)
		return 2;
	else if(isLogZ)
		return 3;
	else;
	return 0;
}
