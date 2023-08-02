#ifndef QDLDL_ST_H
#define QDLDL_ST_H

#include "qdldl_types.h"
#include "types.h"

void QDLDL_solve_st(QDLDL_float x[16]);
void LDLSolve_st(QDLDL_float *x, QDLDL_float *b, c_float *bp);

#endif /* QDLDL_ST_H */
