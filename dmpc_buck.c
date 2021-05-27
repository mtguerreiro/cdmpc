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

	/* Augmented states */
	float xa[3];

	/* Auxiliary variables for intermediate computations */
	float auxm1[CONFIG_NC][CONFIG_NC];
	float auxm2[CONFIG_NC][CONFIG_NC];
	float aux1;

	/* Matrices and vectors */
	float Fj[CONFIG_NC];
	float Kj[4];
	float gam[4];
	float lambda[4];

	/* Optimal control increment */
	float DU;

	xa[0] = x[0] - x_1[0];
	xa[1] = x[1] - x_1[1];
	xa[2] = x[1];

	/* Computes Fj */
	mulmv(Fj_1, CONFIG_NC, &r, 1, (float *)auxm1);
	mulmv((float *)Fj_2, CONFIG_NC, xa, 3, (float *)auxm2);
	sumv((float *)auxm1, (float *)auxm2, CONFIG_NC, Fj);

	/* Computes gamma */
	mulmv((float *)Fx, 1, xa, 3, &aux1);
	gam[0] = -CONFIG_U_MIN + u_1;
	gam[1] = CONFIG_U_MAX - u_1;
	gam[2] = -(CONFIG_IL_MIN - x[0]) + aux1;
	gam[3] = (CONFIG_IL_MAX - x[0]) - aux1;

	/* Computes Kj */
	mulmv((float *)Kj_1, 4, Fj, CONFIG_NC, (float *)auxm1);
	sumv(gam, (float *)auxm1, 4, Kj);

	/* Opt */
	qpHild((float *)Hj, Kj, 4, lambda);

	/* DU */
	mulmv(DU_1, 1, Fj, CONFIG_NC, &DU);
	mulmv(DU_2, 1, lambda, 4, &aux1);
	DU = DU + aux1;


	return DU;
}
//-----------------------------------------------------------------------------
//=============================================================================
