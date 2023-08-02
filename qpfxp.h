/**
 * @file qpfxp.h
 * @brief Provides quadratic programming algorithms in fixed point.
 *  Created on: 12.07.2021
 *      Author: mguerreiro
 */

#ifndef QPFXP_H_
#define QPFXP_H_

//=============================================================================
/*-------------------------------- Includes ---------------------------------*/
//=============================================================================
#include "stdint.h"

#include "C:\Users\mguerreiro\Documents\Projects\C\psutils\fixedmath.h"
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
 * The algorithm runs for a fixed number of iterations.
 *
 * @param H Pointer to H matrix.
 * @param K Pointer to K matrix.
 * @param nIter Number of iterations.
 * @param lambda Buffer to hold the solution.
 * @param lambdaSize Size of the lambda vector.
 * @return Number of iterations required to find solution.
 */
uint32_t qpfxpHildFixedIter(fmint_t *H, fmint_t *K, uint32_t nIter, fmint_t* lambda, uint32_t lambdaSize);
//-----------------------------------------------------------------------------
uint32_t qpfxpHildFixedIterHLS(fmint_t *K_if, fmint_t *lambda_if);
//=============================================================================

#endif /* QPFXP_H_ */
