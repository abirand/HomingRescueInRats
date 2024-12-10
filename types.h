#ifndef TYPES_H_
#define TYPES_H_

// These should be reasonably portable 64-bit integer types
typedef unsigned long long int u64b;
typedef signed   long long int s64b;

// These should be reasonably portable 32-bit integer types
typedef unsigned int u32b;
typedef signed   int s32b;

typedef unsigned char byte;

/**
 * Gets the matrix entry (i,j) of a row-ordered, one-dimensional matrix of 'c' columns.
 * @param a The array of values. MUST be row-ordered and one-dimensional.
 * @param i The row index of the value to get.
 * @param j The column index of the value to get.
 * @param c The number of columns in the matrix.
 * @return The (i, j) entry of the matrix.
 */
#define VALUE(a,i,j,c)  ((a)[(i)*(c)+(j)])

#endif /*TYPES_H_*/
