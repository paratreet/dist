#ifndef DEFINES_H
#define DEFINES_H

#include "limits.h"
#include "charm++.h"

#include <stdint.h>
#include <cfloat>

#ifndef USE_DOUBLE_FP_REAL
typedef float Real;
#define REAL_MAX FLT_MAX
#else
typedef double Real;
#define REAL_MAX DBL_MAX
#endif
typedef uint64_t Key;

#define   PI         3.14159265358979323846

// Number of bits needed for a node or particle's Key (defaults to 64) 
#define TREE_KEY_BITS (sizeof(Key)*CHAR_BIT)

/*
 * When the floating point coordinates of a particle are integerized, 
 * there are BITS_PER_DIM bits available (out of a total of TREE_KEY_BITS)
 * to hold the integer equivalent of each floating point coordinate.
 */
#define BITS_PER_DIM 21
/*
 * Therefore, each integer coordinate of a particle's position can take one
 * of BOXES_PER_DIM values. 
 */
#define BOXES_PER_DIM (1<<(BITS_PER_DIM))

#define NDIMS 3

#define BRANCH_FACTOR 2
#define LOG_BRANCH_FACTOR 1

/*
 * When doing decomposition or tree building, we allow leaves to have up to
 * *_TOLERANCE times the number of particles specfied by user (through -ppc for 
 * decomposition and -b for tree building)
 */
#define DECOMP_TOLERANCE 1.2
#define BUCKET_TOLERANCE 1.2

#endif // DEFINES_H
