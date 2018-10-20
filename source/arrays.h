#ifndef ARRAYS_H
#define ARRAYS_H

// swap two locations in an array
extern void Swap(int *arr, int from, int to);

// shift the elements of an array (or a sub-array of an array) to the left (or right)
extern void Shift(int *arr, int start, int end, int direction, int nshifts);

// return the minimum of two real numbers
extern double Min(double a, double b);
// return the maximum of two real numbers
extern double Max(double a, double b);

/*
	Given a pointer-to-pointer, return a pointer to a list of number that are formed by taking a particular (col-th) column of the input (jagged) array.
	The col input indicates the column index which must be fetched from each row. This must be given keeping in mind any segmentation faults coming from accessing an unallocated memory location.
*/
extern double* Slice(double **arr, int rows, int col);

// Copy the contents of an array into another
extern void Copy(int* original, int nelems, int *duplicate);

// Determine if a given element esists in an integer array. Return 1 if yes and 0 otherwise.
extern int Element(int elem, int *arr, int size);

// Determine the index of an element in an array. If not found, return -1.
extern int Index(int elem, int *arr, int size, int *ignore);

// Determine if two arrays are equal up to permutation
extern int IsEqual(int *arr1, int *arr2, int size);

// Test if a square matrix is positive definite
extern int IsPositiveDefinite(double **matrix, int dim);

// Given a range, compute all the numbers in that range.
extern void Expand(double *range, double *expanded, int dir);

// Compute the size of an interval defined by floating point numbers.
extern int IntervalSize(double *range);

// Fill an integer array with constants
extern void FillIntArrayWithConstants(int* arr, int arrSize, int constant);

// Fill a long array with constants
extern void FillLongArrayWithConstants(long* arr, int arrSize, long constant);

// Fill a float array with constants
extern void FillFloatArrayWithConstants(float* arr, int arrSize, float constant);

// Fill a double array with constants
extern void FillDoubleArrayWithConstants(double* arr, int arrSize, double constant);

// Print the contents of an double array
extern void PrintDoubleArray(double* arr, int size);

// Print the contents of an integer array
extern void PrintIntArray(int* arr, int arrSize);

// Print the contents of an long array
extern void PrintLongArray(long* arr, int arrSize);

// Compute the Hamming weight of the bit string that is used to represent a given decimal number
// http://stackoverflow.com/questions/22081738/how-does-this-algorithm-to-count-the-number-of-set-bits-in-a-32-bit-integer-work
extern int HammingWeight(unsigned int deci);

// Compute the number of subsets of a given size, in a superset of another given size.
extern long Combination(int total, int size);

// Convert an integer to binary
extern void IntToBinary(int deci, int nbits, int *bseq);

// Compute the sum of elements in an integer array
extern int Sum(int *arr, int nelems);

// Find the product of all the elements in an integer array
extern double Prod(double *arr, int nelems);

// Compute the number of subsets of a given size, in a superset of another given size.
extern long Combination(int total, int size);

/* 
	Generate all possible collection of a set, all of which have a given number of elements.
	Every subset corresponds to a number from 0 to 2^(size of the superset)

*/
extern void GenSubsets(int *superset, int length, int **collection);

/*
	Compute the cartesian product of a list of list of numbers
	
	The input must be formatted as follows.
		arr[0][0] = number of lists described in arr
		arr[j][0] = (for j > 0) number of elements in the j-th list
		arr[j][1..] = elements of the set j.

	Output: A such that
		A[0][0] contains the number of multiples in the cartesian product
		A[0][1] contains the number of elements in each multiple of the cartesian product (equal to the number of lists in the input)
		A[j] contains the j-th element of the cartesian product.	

	Algorithm:
		Let there be M sets, arr[j][0] = aj and N = a1 * a2 ... aM.
		Every element of the cartesian product can be associated to a number between 1 and total number of elements in the cartesian product, in the following way.
			1. Take a number X between 0 and N. We will derive an element of the cartesian product, SX, from X.
			2. For k = 0 to M-1, do
				A. Divide X by ak, the remainder is the k-th element of SX.
				B. X -> X/ak
*/
extern void CartesianProduct(double **arr, double **cartProd);

#endif /* ARRAYS_H */
