#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "memory.h"
#include "arrays.h"

void Swap(int *arr, int from, int to){
	// swap two locations in an array
	int temp;
	temp = arr[from];
	arr[from] = arr[to];
	arr[to] = temp;
}

void Shift(int *arr, int start, int end, int direction, int nshifts){
	// shift the elements of an array (or a sub-array of an array) to the left (or right)
	// when direction = 1, the shift is done left-wards.
	int s, i;
	if (direction == 0)
		for (s = 0; s < nshifts; s ++)
			for (i = end; i > start; i --)
				arr[i] = arr[i - 1];
	else
		for (s = 0; s < nshifts; s ++)
			for (i = start; i < end; i ++)
				arr[i] = arr[i + 1];
}

double Min(double a, double b){
	// return the minimum of two numbers
	if (a < b)
		return a;
	return b;
}

double Max(double a, double b){
	// return the minimum of two numbers
	if (a < b)
		return b;
	return a;
}

void Reverse(double *arr, int start, int end){
	// Reverse an array or a subset of an array
	double temp;
	while (start < end - 1){
		// swap the array elements at locations start and end
		// inrease start by 1 and decrease end by 1
		temp = arr[start];
		arr[start] = arr[end];
		arr[end] = temp;
		start ++;
		end --;
	}
}

double* Slice(double **arr, int rows, int col){
	// Given a pointer-to-pointer, return a pointer to a list of number that are formed by taking a particular (col-th) column of the input (jagged) array.
	int ri;
	double *column = malloc(sizeof(double) * rows);
	for (ri = 0; ri < rows; ri ++)
		column[ri] = arr[ri][col];
	return column;
}

void Copy(int* original, int nelems, int *duplicate){
	// Copy the contents of an array into another
	int ei;
	for (ei = 0; ei < nelems; ei ++)
		duplicate[ei] = original[ei];
}


int Element(int elem, int *arr, int size){
	// Determine if a given element exists in an integer array. Return 1 if yes and 0 otherwise.
	int ei;
	for (ei = 0; ei < size; ei ++)
		if (arr[ei] == elem)
			return 1;
	return 0;
}

int Index(int elem, int *arr, int size, int *ignore){
	// Determine the index of an element in an array. If not found, return -1.
	int i;
	for (i = 0; i < size; i ++){
		if (arr[i] == elem){
			if (ignore != NULL){
				if (ignore[i] == 0)
					return i;
			}
			else
				return i;
		}
	}
	return -1;
}

int IsEqual(int *arr1, int *arr2, int size){
	// Determine if two arrays are equal up to permutation
	int i, index, *visited = malloc(sizeof(int) * size);
	FillIntArrayWithConstants(visited, size, 0);
	for (i = 0; i < size; i ++){
		index = Index(arr1[i], arr2, size, visited);
		if (index == -1)
			return 0;
		else
			visited[index] = 1;
	}
	free(visited);
	return 1;
}

void Eigenvalues(double **matrix, int dim, double *eigvals){
	// Find the eigenvalues of a matrix using Jacobi's iteration method
	if (dim == 2){
		// the eigenvalues of a 2x2 symmetric matrix {{a,b},{b,d}} are given by
		// l1 = 1/2 (a + d - sqrt(a^2 + 4 b^2 - 2 a d + d^2))
		// l2 = 1/2 (sqrt(a^2 - 2 a d + 4 b^2 + d^2) + a + d)
		eigvals[0] = (matrix[0][0] + matrix[1][1] - sqrt(pow(matrix[0][0], 2) + 4 * pow(matrix[0][1], 2) - 2 * matrix[0][0] * matrix[1][1] + pow(matrix[1][1], 2)))/(double) 2;
		eigvals[1] = (sqrt(pow(matrix[0][0], 2) - 2 * matrix[0][0] * matrix[1][1] + 4 * pow(matrix[0][1], 2) + pow(matrix[1][1], 2)) + matrix[0][0] + matrix[1][1])/(double) 2;
	}
	else{
		fprintf(stderr, "\033[2m[Positive Definite] Could not determine the eigenvalues of a %d x %d matrix.\033[0m\n", dim, dim);
	}
}	

int IsPositiveDefinite(double **matrix, int dim){
	// Test if a square matrix is positive definite
	double *eigvals = malloc(sizeof(double) * dim);
	Eigenvalues(matrix, dim, eigvals);
	int i, ispd = 1;
	for (i = 0; i < dim; i ++)
		if (eigvals[i] < 0)
			ispd = 0;
	free(eigvals);
	return ispd;
}

void Expand(double *range, double *expanded, int dir){
	// Given a range, compute all the numbers in that range.
	// if dir = 0, expand the range from largest to smallest. If dir = 1, expand the range from smallest to largest
	// There is a subtle issue with looping over floats, which is why the loop has to be over integers while quantity of interest is obtained by scaling the loop variable appropriately.
	// See http://stackoverflow.com/questions/40134984/looping-over-a-range-of-floats-in-c
	long scale, ivar;
	scale = (long) pow(10, ceil(fabs(log(range[2]))) + 1);
	int size = 0;
	if (dir == 1)
		for (ivar = (long) (range[0] * scale); ivar <= (long) (range[1] * scale); ivar += (long) (range[2] * scale))
			expanded[++ size] = (double) (ivar/(double) scale);
	else
		for (ivar = (long) (range[1] * scale); ivar >= (long) (range[0] * scale); ivar -= (long) (range[2] * scale))
			expanded[++ size] = (double) (ivar/(double) scale);
	// The negative part of the range is expanded according to the oppostive of "dir" while the positive part is expanded according to "dir"
	if (range[0] < 0){
		// there is a negative part
		int i, start = -1, end = - 1;
		for (i = 1; i < size; i ++){
			if (start == -1){
				if (expanded[i] < 0)
					start = i;
			}
			else{
				if (end == -1)
					if (expanded[i] > 0)
						end = i - 1;
			}
		}
		if (end == -1)
			end = size - 1;
		Reverse(expanded, start, end);
	}
}

int IntervalSize(double *range){
	// Compute the size of an interval defined by floating point numbers.
	long ivar, scale = (long) pow(10, ceil(fabs(log(range[2]))) + 1);
	int interval = 0;
	for (ivar = (long) (range[1] * scale); ivar >= (long) (range[0] * scale); ivar -= (long) (range[2] * scale))
		interval ++;
	return interval;
}

void FillIntArrayWithConstants(int* arr, int arrSize, int constant){
	// Fill an integer array with constants
	int i = 0;
	for(i = 0; i < arrSize; i ++){
		arr[i] = constant;
	}
}

void FillLongArrayWithConstants(long* arr, int arrSize, long constant){
	// Fill an integer array with constants
	int i = 0;
	for(i = 0; i < arrSize; i ++){
		arr[i] = constant;
	}
}

void FillFloatArrayWithConstants(float* arr, int arrSize, float constant){
	// Fill an integer array with constants
	int i = 0;
	for(i = 0; i < arrSize; i ++){
		arr[i] = constant;
	}
}

void FillDoubleArrayWithConstants(double* arr, int arrSize, double constant){
	// Fill an integer array with constants
	int i = 0;
	for(i = 0; i < arrSize; i ++){
		arr[i] = constant;
	}
}

void PrintIntArray(int* arr, int arrSize){
	// Print the contents of an integer array
	int i;
	for(i = 0; i < arrSize; i ++)
		printf(" %d", arr[i]);
	printf(".\n");
}

void PrintLongArray(long* arr, int arrSize){
	// Print the contents of an long array
	int i;
	for(i = 0; i < arrSize; i ++)
		printf(" %ld", arr[i]);
	printf(".\n");
}


void PrintDoubleArray(double* arr, int size){
	// Print the contents of an double array
	int i;
	for(i = 0; i < size; i ++)
		printf(" %g", arr[i]);
	printf(".\n");
}


int HammingWeight(unsigned int deci){
	// Compute the Hamming weight of the bit string that is used to represent a given decimal number
	// http://stackoverflow.com/questions/22081738/how-does-this-algorithm-to-count-the-number-of-set-bits-in-a-32-bit-integer-work
    deci = deci - ((deci >> 1) & 0x55555555);
    deci = (deci & 0x33333333) + ((deci >> 2) & 0x33333333);
    return (((deci + (deci >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

void IntToBinary(int deci, int nbits, int *bseq){
	// Convert an integer to binary
	int bi;
	for (bi = 0; bi < nbits; bi ++){
		bseq[nbits - bi - 1] = (deci % 2);
		deci = (int) (deci/2);
	}
}

int Sum(int *arr, int nelems){
	// Compute the sum of elements in an integer array
	int ei, sum = 0;
	for (ei = 0; ei < nelems; ei ++)
		sum += arr[ei];
	return sum;
}

double Prod(double *arr, int nelems){
	// Find the product of all the elements in an integer array
	int ei;
	double prod = 1;
	for (ei = 0 ; ei < nelems; ei ++)
		prod = prod * arr[ei];
	return prod;
}

long Combination(int total, int size){
	// Compute the number of subsets of a given size, in a superset of another given size.
	int ei;
	long nfact = 1, kfact = 1, diffFact = 1, nsubsets;
	// compute the binomial coefficient -- size of superset choose size of subset
	for (ei = 1; ei <= total; ei ++)
		nfact *= ei;
	for (ei = 1; ei <= size; ei ++)
		kfact *= ei;
	for (ei = 1; ei <= (total - size); ei ++)
		diffFact *= ei;
	nsubsets = (long) (nfact/(kfact * diffFact));
	return nsubsets;
}

void GenSubsets(int *superset, int length, int **collection){
	// Generate all possible collection of a set, all of which have a given number of elements.
	int ei, bi, idx;
	int *bseq = malloc(sizeof(int) * superset[0]);
	collection[0][0] = 1;
	collection[0][1] = length;
	for (bi = 0; bi < pow(2, superset[0]); bi ++){
		IntToBinary(bi, superset[0], bseq);
		if (Sum(bseq, superset[0]) == length){
			idx = 0;
			collection[collection[0][0]] = malloc(sizeof(int) * length);
			for (ei = 1; ei <= superset[0]; ei ++){
				if (bseq[ei - 1] == 1){
					collection[collection[0][0]][idx] = superset[ei];
					idx ++;
				}
			}
			collection[0][0] ++;
		}
	}
	collection[0][0] --;
	free(bseq);
}


void CartesianProduct(double **arr, double **cartProd){
	// Compute the cartersian product of a list of list of numbers.
	int si, id, prod = 1;
	cartProd[0][0] = Prod(Slice(arr + 1, (int) arr[0][0], 0), (int) arr[0][0]);
	cartProd[0][1] = arr[0][0];
	for (id = 0; id < (int) cartProd[0][0]; id ++){
		prod = 1;
		for (si = 1; si <= (int) arr[0][0]; si ++){
			cartProd[id + 1][si - 1] = arr[si][1 + ((id/prod) % (int) arr[si][0])];
			prod = prod * (int) arr[si][0];
		}
	}
}
