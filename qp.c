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
uint32_t qpHild4(float *H, float *K, uint32_t nIter, float* lambda, float tol){

    uint32_t stopCond;
    float res;
    float lambdap[4] = {0};
    float *h;
    float *k;

    uint32_t n;

    n = 0;
    k = K;
    h = H;
    while( n < nIter ){

        stopCond = 0;

        res = lambdap[0];
        lambdap[0] = h[0] * (k[0] + h[1] * lambdap[1] + h[2] * lambdap[2] + h[3] * lambdap[3]);
        if( lambdap[0] < 0 ) lambdap[0] = 0;

        if( stopCond == 0 ){
            if( res > lambdap[0] ) res = res - lambdap[0];
            else res = lambdap[0] - res;
            if( res > tol ) stopCond = 1;
        }

        res = lambdap[1];
        lambdap[1] = h[5] * (k[1] + h[4] * lambdap[0] + h[6] * lambdap[2] + h[7] * lambdap[3]);
        if( lambdap[1] < 0 ) lambdap[1] = 0;

        if( stopCond == 0 ){
            if( res > lambdap[1] ) res = res - lambdap[1];
            else res = lambdap[1] - res;
            if( res > tol ) stopCond = 1;
        }

        res = lambdap[2];
        lambdap[2] = h[10] * (k[2] + h[8] * lambdap[0] + h[9] * lambdap[1] + h[11] * lambdap[3]);
        if( lambdap[2] < 0 ) lambdap[2] = 0;

        if( stopCond == 0 ){
            if( res > lambdap[2] ) res = res - lambdap[2];
            else res = lambdap[2] - res;
            if( res > tol ) stopCond = 1;
        }

        res = lambdap[3];
        lambdap[3] = h[15] * (k[3] + h[12] * lambdap[0] + h[13] * lambdap[1] + h[14] * lambdap[2]);
        if( lambdap[3] < 0 ) lambdap[3] = 0;

        if( stopCond == 0 ){
            if( res > lambdap[3] ) res = res - lambdap[3];
            else res = lambdap[3] - res;
            if( res > tol ) stopCond = 1;
        }

        n++;

        if( stopCond == 0 ) break;
    }

    lambda[0] = lambdap[0];
    lambda[1] = lambdap[1];
    lambda[2] = lambdap[2];
    lambda[3] = lambdap[3];

    return n;
}
//-----------------------------------------------------------------------------
//uint32_t qpHild4(float *H, float *K, uint32_t nIter, float* lambda, float tol){
//
//    uint32_t stopCond;
//    float res;
//    float lambdap[4] = {0};
//    float *h;
//    float *k;
//
//    uint32_t n, j;
//
//    n = 0;
//    k = K;
//    while( n < nIter ){
//
//        stopCond = 0;
//
//        h = H;
//        k = K;
//
//        res = lambdap[0];
//        lambdap[0] = h[0] * (*k + h[1] * lambdap[1] + h[2] * lambdap[2] + h[3] * lambdap[3]);
//        if( lambdap[0] < 0 ) lambdap[0] = 0;
//
//        if( stopCond == 0 ){
//            if( res > lambdap[0] ) res = res - lambdap[0];
//            else res = lambdap[0] - res;
//            if( res > tol ) stopCond = 1;
//        }
//
//        k++;
//        h = h + 4;
//        res = lambdap[1];
//        lambdap[1] = h[1] * (*k + h[0] * lambdap[0] + h[2] * lambdap[2] + h[3] * lambdap[3]);
//        if( lambdap[1] < 0 ) lambdap[1] = 0;
//
//        if( stopCond == 0 ){
//            if( res > lambdap[1] ) res = res - lambdap[1];
//            else res = lambdap[1] - res;
//            if( res > tol ) stopCond = 1;
//        }
//
//        k++;
//        h = h + 4;
//        res = lambdap[2];
//        lambdap[2] = h[2] * (*k + h[0] * lambdap[0] + h[1] * lambdap[1] + h[3] * lambdap[3]);
//        if( lambdap[2] < 0 ) lambdap[2] = 0;
//
//        if( stopCond == 0 ){
//            if( res > lambdap[2] ) res = res - lambdap[2];
//            else res = lambdap[2] - res;
//            if( res > tol ) stopCond = 1;
//        }
//
//        k++;
//        h = h + 4;
//        res = lambdap[3];
//        lambdap[3] = h[3] * (*k + h[0] * lambdap[0] + h[1] * lambdap[1] + h[2] * lambdap[2]);
//        if( lambdap[3] < 0 ) lambdap[3] = 0;
//
//        if( stopCond == 0 ){
//            if( res > lambdap[3] ) res = res - lambdap[3];
//            else res = lambdap[3] - res;
//            if( res > tol ) stopCond = 1;
//        }
//
//        n++;
//
//        if( stopCond == 0 ) break;
//    }
//
//    h = lambdap;
//    k = lambda;
//    j = 4;
//    while(j--){
//        *k++ = *h++;
//    }
//
//    return n;
//}
//-----------------------------------------------------------------------------
//=============================================================================
