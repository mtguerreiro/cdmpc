/*
 * dmpc_buck.c
 *
 *  Created on: 27.05.2021
 *      Author: mguerreiro
 */

//=============================================================================
/*-------------------------------- Includes ---------------------------------*/
//=============================================================================
#include "dmpc_buck.h"
#include "dmpc_buck_matrices.h"

#include "mvops.h"
#include "qp.h"
//=============================================================================


//=============================================================================
/*-------------------------------- Functions --------------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
float dmpcBuckOpt(float *x, float *x_1, float r, float u_1){

	uint32_t i;
	uint32_t n_iter;

	/* Augmented states */
	float xa[3];

	/* Auxiliary variables for intermediate computations */
	float auxm1[DMPC_BUCK_CONFIG_NC][DMPC_BUCK_CONFIG_NC];
	float auxm2[DMPC_BUCK_CONFIG_NC][DMPC_BUCK_CONFIG_NC];
	float aux1;

	/* Matrices and vectors */
	float Fj[DMPC_BUCK_CONFIG_NC];
	float Kj[DMPC_BUCK_CONFIG_NLAMBDA];
	float gam[DMPC_BUCK_CONFIG_NLAMBDA];
	float lambda[DMPC_BUCK_CONFIG_NLAMBDA];

	/* Optimal control increment */
	float DU;

	xa[0] = x[0] - x_1[0];
	xa[1] = x[1] - x_1[1];
	xa[2] = x[1];

	/* Computes Fj */
	mulmv(Fj_1, DMPC_BUCK_CONFIG_NC, &r, 1, (float *)auxm1);
	mulmv((float *)Fj_2, DMPC_BUCK_CONFIG_NC, xa, 3, (float *)auxm2);
	sumv((float *)auxm1, (float *)auxm2, DMPC_BUCK_CONFIG_NC, Fj);

	/* Computes gamma */
	mulmv((float *)Fx, DMPC_BUCK_CONFIG_NR, xa, 3, &aux1);
	for(i = 0; i < DMPC_BUCK_CONFIG_NR; i++){
		gam[i] = -DMPC_BUCK_CONFIG_U_MIN + u_1;
		gam[i + DMPC_BUCK_CONFIG_NR] = DMPC_BUCK_CONFIG_U_MAX - u_1;

		//mulmv((float *)Fx, DMPC_BUCK_CONFIG_NR, xa, 3, &aux1);
		mulmv((float *)Fx[i], DMPC_BUCK_CONFIG_NR, xa, 3, &aux1);
		gam[i + 2 * DMPC_BUCK_CONFIG_NR] = -(DMPC_BUCK_CONFIG_IL_MIN - x[0]) + aux1;
		gam[i + 3 * DMPC_BUCK_CONFIG_NR] = (DMPC_BUCK_CONFIG_IL_MAX - x[0]) - aux1;
	}

	/* Computes Kj */
	mulmv((float *)Kj_1, DMPC_BUCK_CONFIG_NLAMBDA, Fj, DMPC_BUCK_CONFIG_NC, (float *)auxm1);
	sumv(gam, (float *)auxm1, DMPC_BUCK_CONFIG_NLAMBDA, Kj);

	/* Opt */
	n_iter = qpHild((float *)Hj, Kj, 100, lambda, DMPC_BUCK_CONFIG_NLAMBDA);

	/* DU */
	mulmv(DU_1, 1, Fj, DMPC_BUCK_CONFIG_NC, &DU);
	mulmv(DU_2, 1, lambda, DMPC_BUCK_CONFIG_NLAMBDA, &aux1);
	DU = DU + aux1;


	return DU;
}
//-----------------------------------------------------------------------------
//=============================================================================
