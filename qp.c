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
void qpHild8FixedIter(float *H, float *K, uint32_t nIter, float* lambda){

    float *h;
    float *k;

    uint32_t n;

    /* Initial value for lambda */
    lambda[0] = 0;
    lambda[1] = 0;
    lambda[2] = 0;
    lambda[3] = 0;
    lambda[4] = 0;
    lambda[5] = 0;
    lambda[6] = 0;
    lambda[7] = 0;

    n = 0;
    k = K;
    h = H;
    while( n < nIter ){

        lambda[0] = h[0] * (k[0] + h[1] * lambda[1] + h[2] * lambda[2] + h[3] * lambda[3] + h[4] * lambda[4] + h[5] * lambda[5] + h[6] * lambda[6] + h[7] * lambda[7]);
        if( lambda[0] < 0 ) lambda[0] = 0;

        lambda[1] = h[9] * (k[1] + h[8] * lambda[0] + h[10] * lambda[2] + h[11] * lambda[3] + h[12] * lambda[4] + h[13] * lambda[5] + h[14] * lambda[6] + h[15] * lambda[7]);
        if( lambda[1] < 0 ) lambda[1] = 0;

        lambda[2] = h[18] * (k[2] + h[16] * lambda[0] + h[17] * lambda[1] + h[19] * lambda[3] + h[20] * lambda[4] + h[21] * lambda[5] + h[22] * lambda[6] + h[23] * lambda[7]);
        if( lambda[2] < 0 ) lambda[2] = 0;

        lambda[3] = h[27] * (k[3] + h[24] * lambda[0] + h[25] * lambda[1] + h[26] * lambda[2] + h[28] * lambda[4] + h[29] * lambda[5] + h[30] * lambda[6] + h[31] * lambda[7]);
        if( lambda[3] < 0 ) lambda[3] = 0;

        lambda[4] = h[36] * (k[4] + h[32] * lambda[0] + h[33] * lambda[1] + h[34] * lambda[2] + h[35] * lambda[3] + h[37] * lambda[5] + h[38] * lambda[6] + h[39] * lambda[7]);
        if( lambda[4] < 0 ) lambda[4] = 0;

        lambda[5] = h[45] * (k[5] + h[40] * lambda[0] + h[41] * lambda[1] + h[42] * lambda[2] + h[43] * lambda[3] + h[44] * lambda[4] + h[46] * lambda[6] + h[47] * lambda[7]);
        if( lambda[5] < 0 ) lambda[5] = 0;

        lambda[6] = h[54] * (k[6] + h[48] * lambda[0] + h[49] * lambda[1] + h[50] * lambda[2] + h[51] * lambda[3] + h[52] * lambda[4] + h[53] * lambda[5] + h[55] * lambda[7]);
        if( lambda[6] < 0 ) lambda[6] = 0;

        lambda[7] = h[63] * (k[7] + h[56] * lambda[0] + h[57] * lambda[1] + h[58] * lambda[2] + h[59] * lambda[3] + h[60] * lambda[4] + h[61] * lambda[5] + h[62] * lambda[6]);
        if( lambda[7] < 0 ) lambda[7] = 0;

        n++;
    }
}
//-----------------------------------------------------------------------------
//=============================================================================
