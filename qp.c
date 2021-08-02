/*
 * qp.c
 *
 *  Created on: 27.05.2021
 *      Author: mguerreiro
 *
 *TODO: better initialize lambdap to zero.
 */

//=============================================================================
/*-------------------------------- Includes ---------------------------------*/
//=============================================================================
#include "qp.h"

//=============================================================================

//=============================================================================
/*-------------------------------- Functions --------------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
uint32_t qpHild(float *H, float *K, uint32_t nIter, float* lambda, uint32_t lambdaSize, float tol){

	uint32_t stopCond;
	float res;
	float lambdap[QP_CONFIG_MAX_LAMBDA];
	float *h;
	float *k;

	uint32_t n, i, j;

	for(i = 0; i < lambdaSize; i++){
		lambdap[i] = 0;
	}

	n = 0;
	k = K;
	while( n < nIter ){

		h = H;
		stopCond = 0;

		/* One iteration through the lambda vector */
		for(i = 0; i < lambdaSize; i++){

			res = lambdap[i];
			lambdap[i] = 0;
			for(j = 0; j < i; j++){
				lambdap[i] += h[j] * lambdap[j];
			}

			for(j = i + 1; j < lambdaSize; j++){
				lambdap[i] += h[j] * lambdap[j];
			}

			lambdap[i] = h[i] * (k[i] + lambdap[i]);
			if( lambdap[i] < 0 ) lambdap[i] = 0;

			if( stopCond == 0 ){
				if( res > lambdap[i] ) res = res - lambdap[i];
				else res = lambdap[i] - res;
				if( res > tol ) stopCond = 1;
			}

			h = h + lambdaSize;
		}

		n++;

		/*
		 * If stopCond is zero, all lambdas had a difference between the new
		 * and the previous value that is smaller than the tolerance; so we
		 * have reached the stopping condition.
		 */
		if( stopCond == 0 ) break;
	}

	for(j = 0; j < lambdaSize; j++){
		lambda[j] = lambdap[j];
	}

	return n;
}
//-----------------------------------------------------------------------------
//=============================================================================
