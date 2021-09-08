/**
 * @file qp.h
 * @brief Provides quadratic programming algorithms.
 *  Created on: 27.05.2021
 *      Author: mguerreiro
 */

#ifndef QP_H_
#define QP_H_

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
 * @brief Hildreth's quadratic programming algorithm.
 *
 * For the H matrix, we consider that the diagonal elements are already
 * inverted and with the sign flipped.
 *
 * @param H Pointer to H matrix.
 * @param K Pointer to K matrix.
 * @param nIter Maximum number of iterations.
 * @param lambda Buffer to hold the solution.
 * @param lambdaSize Size of the lambda vector.
 * @param tol Tolerance for stopping condition.
 * @return Number of iterations required to find solution.
 */
uint32_t qpHild(float *H, float *K, uint32_t nIter, float* lambda, uint32_t lambdaSize, float tol);
//-----------------------------------------------------------------------------
/**
 * @brief Hildreth's quadratic programming algorithm.
 *
 * For the H matrix, we consider that the diagonal elements are already
 * inverted and with the sign flipped.
 *
 * This function considers that there are only 4 variables to optimize.
 *
 * @param H Pointer to H matrix.
 * @param K Pointer to K matrix.
 * @param nIter Maximum number of iterations.
 * @param lambda Buffer to hold the solution.
 * @param tol Tolerance for stopping condition.
 * @return Number of iterations required to find solution.
 */
uint32_t qpHild4(float *H, float *K, uint32_t nIter, float* lambda, float tol);
//-----------------------------------------------------------------------------
//=============================================================================

#endif /* QP_H_ */
