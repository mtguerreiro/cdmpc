/**
 * @file mvops.h
 * @brief Provides a few vector and matrix operations.
 *
 *  Created on: 27.05.2021
 *      Author: mguerreiro
 */

#ifndef MVOPS_H_
#define MVOPS_H_

//=============================================================================
/*-------------------------------- Includes ---------------------------------*/
//=============================================================================
#include "stdint.h"

//=============================================================================

//=============================================================================
/*-------------------------------- Functions --------------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
/**
 * @brief Sums two vectors, element-wise.
 *
 * Both vectors and the output must have the same dimension.
 *
 * @param v1 Pointer to first vector.
 * @param v2 Pointer to second vector.
 * @param dim Number of elements to sum.
 * @param y Pointer to buffer where the sum is to be stored.
 */
void sumv(float *v1, float *v2, uint32_t dim, float *y);
//-----------------------------------------------------------------------------
/**
 * @brief Subtracts two vectors, element-wise.
 *
 * Both vectors and the output must have the same dimension.
 *
 * @param v1 Pointer to first vector.
 * @param v2 Pointer to second vector.
 * @param dim Number of elements to subtract.
 * @param y Pointer to buffer where the subtraction is to be stored.
 */
void subv(float *v1, float *v2, uint32_t dim, float *y);
//-----------------------------------------------------------------------------
/**
 * @brief Performs a matrix-vector multiplication.
 *
 * The number of columns of the matrix must be the same as the dimension of the
 * vector. The buffer should have be of appropriate size, to store the result.
 *
 * @param m Pointer to matrix.
 * @param mDim Number of lines of the matrix.
 * @param v Pointer to vector.
 * @param nDim Number of elements of the vector.
 * @param y Pointer to buffer where the result is to be stored.
 */
void mulmv(float *m, uint32_t mDim, float *v, uint32_t nDim, float* y);
//-----------------------------------------------------------------------------
//=============================================================================

#endif /* MVOPS_H_ */
