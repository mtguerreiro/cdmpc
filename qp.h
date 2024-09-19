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
 * The algorithm runs for a fixed number of iterations.
 *
 * @param H Pointer to H matrix.
 * @param K Pointer to K matrix.
 * @param nIter Number of iterations.
 * @param lambda Buffer to hold the solution.
 * @param lambdaSize Size of the lambda vector.
 * @return Number of iterations required to find solution.
 */
uint32_t qpHildFixedIter(float *H, float *K, uint32_t nIter, float* lambda, uint32_t lambdaSize);
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
/**
 * @brief Hildreth's quadratic programming algorithm.
 *
 * Runs Hildreth's QP procedure for a fixed number of iterations.
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
 */
void qpHild4FixedIter(float *H, float *K, uint32_t nIter, float* lambda);
//-----------------------------------------------------------------------------
/**
 * @brief Hildreth's quadratic programming algorithm.
 *
 * Runs Hildreth's QP procedure for a fixed number of iterations.
 *
 * For the H matrix, we consider that the diagonal elements are already
 * inverted and with the sign flipped.
 *
 * This function considers that there are only 8 variables to optimize.
 *
 * @param H Pointer to H matrix.
 * @param K Pointer to K matrix.
 * @param nIter Maximum number of iterations.
 * @param lambda Buffer to hold the solution.
 */
void qpHild8FixedIter(float *H, float *K, uint32_t nIter, float* lambda);
//-----------------------------------------------------------------------------
/**
 * @brief Initializes instance and configuration pointers for the IP core running the Hildreth algorithm and the AXI DMA.
 */
uint32_t qpHildInitFPGA(void);
//-----------------------------------------------------------------------------
/**
 * @brief Hildreth's quadratic programming algorithm.
 *
 * Runs Hildreth's QP procedure for the configured number of iterations on the HLS generated IP core.
 *
 * The H matrix is initialized on the FPGA.
 *
 * This function considers the configured number of variables in the IP core to optimize.
 *
 * @param K Pointer to K matrix.
 * @param lambda Buffer to hold the solution.
 */
uint32_t qpHildFixedIterFPGA(float *K, float* lambda);
//-----------------------------------------------------------------------------
//=============================================================================

#endif /* QP_H_ */
