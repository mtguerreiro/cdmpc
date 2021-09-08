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
	float *h;
	float *k;

	uint32_t n, i, j;

	for(i = 0; i < lambdaSize; i++){
		lambda[i] = 0;
	}

	n = 0;
	k = K;
	while( n < nIter ){

		h = H;
		stopCond = 0;

		/* One iteration through the lambda vector */
		for(i = 0; i < lambdaSize; i++){

			res = lambda[i];
			lambda[i] = 0;
			for(j = 0; j < i; j++){
				lambda[i] += h[j] * lambda[j];
			}

			for(j = i + 1; j < lambdaSize; j++){
				lambda[i] += h[j] * lambda[j];
			}

			lambda[i] = h[i] * (k[i] + lambda[i]);
			if( lambda[i] < 0 ) lambda[i] = 0;

			if( stopCond == 0 ){
				if( res > lambda[i] ) res = res - lambda[i];
				else res = lambda[i] - res;
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

	return n;
}
//-----------------------------------------------------------------------------
uint32_t qpHild4(float *H, float *K, uint32_t nIter, float* lambda, float tol){

    uint32_t stopCond;
    float res;
    float *h;
    float *k;

    uint32_t n;

    /* Initial value for lambda */
    lambda[0] = 0;
    lambda[1] = 0;
    lambda[2] = 0;
    lambda[3] = 0;

    n = 0;
    k = K;
    h = H;
    while( n < nIter ){

        stopCond = 0;

        res = lambda[0];
        lambda[0] = h[0] * (k[0] + h[1] * lambda[1] + h[2] * lambda[2] + h[3] * lambda[3]);
        if( lambda[0] < 0 ) lambda[0] = 0;

        if( stopCond == 0 ){
            if( res > lambda[0] ) res = res - lambda[0];
            else res = lambda[0] - res;
            if( res > tol ) stopCond = 1;
        }

        res = lambda[1];
        lambda[1] = h[5] * (k[1] + h[4] * lambda[0] + h[6] * lambda[2] + h[7] * lambda[3]);
        if( lambda[1] < 0 ) lambda[1] = 0;

        if( stopCond == 0 ){
            if( res > lambda[1] ) res = res - lambda[1];
            else res = lambda[1] - res;
            if( res > tol ) stopCond = 1;
        }

        res = lambda[2];
        lambda[2] = h[10] * (k[2] + h[8] * lambda[0] + h[9] * lambda[1] + h[11] * lambda[3]);
        if( lambda[2] < 0 ) lambda[2] = 0;

        if( stopCond == 0 ){
            if( res > lambda[2] ) res = res - lambda[2];
            else res = lambda[2] - res;
            if( res > tol ) stopCond = 1;
        }

        res = lambda[3];
        lambda[3] = h[15] * (k[3] + h[12] * lambda[0] + h[13] * lambda[1] + h[14] * lambda[2]);
        if( lambda[3] < 0 ) lambda[3] = 0;

        if( stopCond == 0 ){
            if( res > lambda[3] ) res = res - lambda[3];
            else res = lambda[3] - res;
            if( res > tol ) stopCond = 1;
        }

        n++;

        if( stopCond == 0 ) break;
    }

    return n;
}
//-----------------------------------------------------------------------------
void qpHild4FixedIter(float *H, float *K, uint32_t nIter, float* lambda){

    float *h;
    float *k;

    uint32_t n;

    /* Initial value for lambda */
    lambda[0] = 0;
    lambda[1] = 0;
    lambda[2] = 0;
    lambda[3] = 0;

    n = 0;
    k = K;
    h = H;
    while( n < nIter ){

        lambda[0] = h[0] * (k[0] + h[1] * lambda[1] + h[2] * lambda[2] + h[3] * lambda[3]);
        if( lambda[0] < 0 ) lambda[0] = 0;

        lambda[1] = h[5] * (k[1] + h[4] * lambda[0] + h[6] * lambda[2] + h[7] * lambda[3]);
        if( lambda[1] < 0 ) lambda[1] = 0;

        lambda[2] = h[10] * (k[2] + h[8] * lambda[0] + h[9] * lambda[1] + h[11] * lambda[3]);
        if( lambda[2] < 0 ) lambda[2] = 0;

        lambda[3] = h[15] * (k[3] + h[12] * lambda[0] + h[13] * lambda[1] + h[14] * lambda[2]);
        if( lambda[3] < 0 ) lambda[3] = 0;

        n++;
    }
}
//-----------------------------------------------------------------------------
//=============================================================================
