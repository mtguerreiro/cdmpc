/*
 * qp.c
 *
 *  Created on: 27.05.2021
 *      Author: mguerreiro
 *
 *TODO: better initialize lambdap to zero.
 *TODO: implement stopping condition better. Is it possible do to it without having to save all previous values?
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

	float aux;

	uint32_t n, i, j;

	for(i = 0; i < lambdaSize; i++){
		lambdap[i] = 0;
		lambda[i] = 0;
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

		aux = 0;
		for(i = 0; i < lambdaSize; i++){
			if( lambdap[i] > lambda[i] ) aux += (lambdap[i] - lambda[i]);
			else aux += (lambda[i] - lambdap[i]);
		}
		if( aux < 1e-6 ){
			break;
		}
		else{
			for(j = 0; j < lambdaSize; j++){
				lambda[j] = lambdap[j];
			}
		}

	}

//	for(j = 0; j < lambdaSize; j++){
//		lambda[j] = lambdap[j];
//	}

	return n;
}
//-----------------------------------------------------------------------------
//=============================================================================
