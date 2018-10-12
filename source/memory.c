#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "memory.h"

void FreeIntArray(int **A, int length){
	// Free all the space allocated to an int array
	int ei;
	for (ei = 0; ei < length; ei ++)
		free(A[ei]);
	free(A);
}

void FreeFloatArray(float **A, int length){
	// Free all the space allocated to a float array
	int ei;
	for (ei = 0; ei < length; ei ++)
		free(A[ei]);
	free(A);
}

void FreeDoubleArray(double **A, int length){
	// Free all the space allocated to a double array
	int ei;
	for (ei = 0; ei < length; ei ++)
		free(A[ei]);
	free(A);
}

void FreeLongArray(long **A, int length){
	// Free all the space allocated to a long array
	int ei;
	for (ei = 0; ei < length; ei ++)
		free(A[ei]);
	free(A);
}

void FreeCharArray(char **A, int length){
	// Free all the space allocated to a char array
	int ei;
	for (ei = 0; ei < length; ei ++)
		free(A[ei]);
	free(A);
}
