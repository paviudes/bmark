#ifndef MEMORY_H
#define MEMORY_H

// Free all the space allocated to an int array
extern void FreeIntArray(int **A, int length);

// Free all the space allocated to a float array
extern void FreeFloatArray(float **A, int length);

// Free all the space allocated to a double array
extern void FreeDoubleArray(double **A, int length);

// Free all the space allocated to a long array
extern void FreeLongArray(long **A, int length);

// Free all the space allocated to a char array
extern void FreeCharArray(char **A, int length);

#endif	/* MEMORY_H */
