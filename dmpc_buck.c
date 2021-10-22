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

/* Hildreth's QP */
#include "mvops.h"
#include "qp.h"

/* CVXGEN */
#include "cvxgen/solver.h"
Vars vars;
Params params = {.Ej = {0}, .Fj = {0}, .M = {0}, .y = {0}};
Workspace work;
Settings settings;
//=============================================================================

//=============================================================================
/*-------------------------------- Prototypes -------------------------------*/
//=============================================================================
static float dmpcBuckHildOpt(uint32_t* iters);

static void dmpcBuckCVXGENInitialize(void);
static float dmpcBuckCVXGENOpt(uint32_t* iters);
//=============================================================================

//=============================================================================
/*-------------------------------- Functions --------------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
void dmpcBuckInitialize(void){

	dmpcBuckCVXGENInitialize();
}
//-----------------------------------------------------------------------------
float dmpcBuckOpt(float *x, float *x_1, float r, float u_1, uint32_t* iters){

	uint32_t i;

	/* Optimal control increment */
	float du;

	/* Augmented states */
	float xa[3];

	/* Auxiliary variables for intermediate computations */
	float auxm1[DMPC_BUCK_CONFIG_NC];
	float auxm2[DMPC_BUCK_CONFIG_NC];


	xa[0] = x[0] - x_1[0];
	xa[1] = x[1] - x_1[1];
	xa[2] = x[1];

	/* Computes Fj and gam matrices */
	mulmv(Fj_1, DMPC_BUCK_CONFIG_NC, &r, 1, auxm1);
	mulmv((float *)Fj_2, DMPC_BUCK_CONFIG_NC, xa, 3, auxm2);
	sumv(auxm1, auxm2, DMPC_BUCK_CONFIG_NC, Fj);

	mulmv((float *)Fx, DMPC_BUCK_CONFIG_NR, xa, 3, (float *)auxm1);
	for(i = 0; i < DMPC_BUCK_CONFIG_NR; i++){
		gam[i] = -DMPC_BUCK_CONFIG_U_MIN + u_1;
		gam[i + DMPC_BUCK_CONFIG_NR] = DMPC_BUCK_CONFIG_U_MAX - u_1;

		gam[i + 2 * DMPC_BUCK_CONFIG_NR] = -(DMPC_BUCK_CONFIG_IL_MIN - x[0]) + auxm1[i];
		gam[i + 3 * DMPC_BUCK_CONFIG_NR] = (DMPC_BUCK_CONFIG_IL_MAX - x[0]) - auxm1[i];
	}

	du = dmpcBuckHildOpt(iters);
	//du = dmpcBuckCVXGENOpt(iters);

	return du;
}
//-----------------------------------------------------------------------------
//=============================================================================

//=============================================================================
/*---------------------------- Static functions -----------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
static float dmpcBuckHildOpt(uint32_t* iters){

	uint32_t n_iter;

	/* Auxiliary variables for intermediate computations */
	float auxm1[DMPC_BUCK_CONFIG_NLAMBDA];
	float aux1;

	/* Matrices and vectors */
	float Kj[DMPC_BUCK_CONFIG_NLAMBDA];
	float lambda[DMPC_BUCK_CONFIG_NLAMBDA];

	/* Optimal control increment */
	float DU;

	/* Computes Kj */
	mulmv((float *)Kj_1, DMPC_BUCK_CONFIG_NLAMBDA, Fj, DMPC_BUCK_CONFIG_NC, (float *)auxm1);
	sumv(gam, (float *)auxm1, DMPC_BUCK_CONFIG_NLAMBDA, Kj);

	/* Opt */
	n_iter = qpHild((float *)Hj, Kj, 10000, lambda, DMPC_BUCK_CONFIG_NLAMBDA, (float)1e-2);
	//n_iter = qpHild4((float *)Hj, Kj, 15, lambda, (float)1e-6);
	//qpHild4FixedIter((float *)Hj, Kj, 8, lambda);
	//n_iter = 15;

    if( iters != 0 ) *iters = n_iter;

	/* DU */
	mulmv(DU_1, 1, Fj, DMPC_BUCK_CONFIG_NC, &DU);
	mulmv(DU_2, 1, lambda, DMPC_BUCK_CONFIG_NLAMBDA, &aux1);
	DU = DU + aux1;

	return DU;
}
//-----------------------------------------------------------------------------
void dmpcBuckCVXGENInitialize(void){

	uint32_t i, j, k;
	float *src, *dst;

	set_defaults();  // Set basic algorithm parameters.
	setup_indexing();

	i = DMPC_BUCK_CONFIG_NP * DMPC_BUCK_CONFIG_NC;
	src = &Ej[0][0];
	dst = params.Ej;
	while(i--){
		*dst++ = *src++;
	}

	k = 0;
	for(j = 0; j < DMPC_BUCK_CONFIG_NC; j++){
		for(i = 0; i < DMPC_BUCK_CONFIG_NLAMBDA; i++){
			params.M[k++] = M[i][j];
		}
	}

}
//-----------------------------------------------------------------------------
static float dmpcBuckCVXGENOpt(uint32_t* iters){

	uint32_t i;
	uint32_t n_iter;

	float DU;

	/* We need to copy Fj and gam to CVXGEN settings */
	for(i = 0; i < DMPC_BUCK_CONFIG_NC; i++){
		params.Fj[i] = Fj[i];
	}
	for(i = 0; i < DMPC_BUCK_CONFIG_NLAMBDA; i++){
		params.y[i] = gam[i];
	}

	settings.eps = 1e-2;
	settings.resid_tol = 1e-2;
	settings.max_iters = 300;
	n_iter = solve();
	DU = (float)(*vars.du);
	if( iters != 0 ) *iters = n_iter;

	return DU;
}
//-----------------------------------------------------------------------------
//=============================================================================
