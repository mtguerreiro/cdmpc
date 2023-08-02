/**
 * @file qpfxp.h
 * @brief Provides quadratic programming algorithms in fixed point.
 *  Created on: 12.07.2021
 *      Author: mguerreiro
 */

#ifndef QPFXP_CPP_H_
#define QPFXP_CPP_H_

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
 * The algorithm runs for a fixed number of iterations.
 *
 * @param H Pointer to H matrix.
 * @param K Pointer to K matrix.
 * @param nIter Number of iterations.
 * @param lambda Buffer to hold the solution.
 * @param lambdaSize Size of the lambda vector.
 * @return Number of iterations required to find solution.
 */
//uint32_t qpfxpHildFixedIterCpp(fmint_t *H, fmint_t *K, uint32_t nIter, fmint_t* lambda, uint32_t lambdaSize);
//-----------------------------------------------------------------------------
//uint32_t qpfxpHildFixedIterHLSCpp(fmint_t *K, fmint_t* lambda);

#ifdef __cplusplus
extern "C"
#endif

uint32_t qpfxpHildFixedIterHLSCpp(float *K_float, float *lambda_float, uint32_t *t_exec);
//uint32_t qpfxpHildFixedIterHLSCpp(float *K_float, float *lambda_float);
//uint32_t qpfxpHildFixedIterHLSCpp(hls::stream<pkt > &A, hls::stream<pkt> &B)
//=============================================================================

#endif /* QPFXP_CPP_H_ */
