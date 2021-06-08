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
uint32_t qpHild(float *H, float *K, uint32_t nIter, float* lambda, uint32_t lambdaSize){

	float lambdap[lambdaSize];
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

		for(i = 0; i < lambdaSize; i++){

			lambdap[i] = 0;
			for(j = 0; j < i; j++){
				lambdap[i] += h[j] * lambdap[j];
			}

			for(j = i + 1; j < lambdaSize; j++){
				lambdap[i] += h[j] * lambdap[j];
			}

			lambdap[i] = h[i] * (k[i] + lambdap[i]);
			if( lambdap[i] < 0 ) lambdap[i] = 0;

			h = h + lambdaSize;
		}

		n++;

//		h = H;
//		k = K;
//
//		lambdap[0] = h[0] * (*k + h[1] * lambdap[1] + h[2] * lambdap[2] + h[3] * lambdap[3]);
//		if( lambdap[0] < 0 ) lambdap[0] = 0;
//
//		k++;
//		h = h + 4;
//		lambdap[1] = h[1] * (*k + h[0] * lambdap[0] + h[2] * lambdap[2] + h[3] * lambdap[3]);
//		if( lambdap[1] < 0 ) lambdap[1] = 0;
//
//		k++;
//		h = h + 4;
//		lambdap[2] = h[2] * (*k + h[0] * lambdap[0] + h[1] * lambdap[1] + h[3] * lambdap[3]);
//		if( lambdap[2] < 0 ) lambdap[2] = 0;
//
//		k++;
//		h = h + 4;
//		lambdap[3] = h[3] * (*k + h[0] * lambdap[0] + h[1] * lambdap[1] + h[2] * lambdap[2]);
//		if( lambdap[3] < 0 ) lambdap[3] = 0;
//
//		n++;
	}

	for(j = 0; j < lambdaSize; j++){
		lambda[j] = lambdap[j];
	}

	return n;
}
//-----------------------------------------------------------------------------
//=============================================================================
