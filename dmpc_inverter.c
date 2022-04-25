/*
 * dmpc_inverter.c
 *
 *  Created on: 22.04.2022
 *      Author: mguerreiro
 */

//=============================================================================
/*-------------------------------- Includes ---------------------------------*/
//=============================================================================
#include "dmpc_inverter.h"
#include "dmpc_inverter_matrices.h"

/* Hildreth's QP */
#include "mvops.h"
#include "qp.h"

///* CVXGEN */
//#include "cvxgen/solver.h"
//Vars vars;
//Params params = {.Ej = {0}, .Fj = {0}, .M = {0}, .y = {0}};
//Workspace work;
//Settings settings;
//=============================================================================

#define DMPC_INVERTER_CONFIG_NC_x_NU	DMPC_INVERTER_CONFIG_NC * DMPC_INVERTER_CONFIG_NU

float DMPC_INVERTER_CONFIG_U_MIN[4] = {-375.27767497325675, -375.27767497325675, 162.5, 0.0};
float DMPC_INVERTER_CONFIG_U_MAX[4] = {375.27767497325675, 375.27767497325675, 162.5, 0.0};

#define DMPC_INVERTER_CONFIG_N_U_CONSTRAINTS		4

float DMPC_INVERTER_CONFIG_XM_MIN[2] = {-15.0, -15.0};
float DMPC_INVERTER_CONFIG_XM_MAX[2] = {15.0, 15.0};

#define DMPC_INVERTER_CONFIG_N_STATE_CONSTRAINTS	2

#define DMPC_INVERTER_CONFIG_N_STATES				6
#define DMPC_INVERTER_CONFIG_N_AUG_STATES			8

//=============================================================================
/*-------------------------------- Prototypes -------------------------------*/
//=============================================================================
static void dmpcInverterHildOpt(float *du, uint32_t* iters);

//static void dmpcInverterCVXGENInitialize(void);
//static float dmpcInverterCVXGENOpt(uint32_t* iters);
//=============================================================================

//=============================================================================
/*-------------------------------- Functions --------------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
//void dmpcInverterInitialize(void){
//
//	dmpcInverterCVXGENInitialize();
//}
//-----------------------------------------------------------------------------
void dmpcInverterOpt(float *x, float *x_1, float *r, float *u_1, float *du, uint32_t* iters){

	uint32_t i, j, k, w;

	/* Optimal control increment */
	//float du[DMPC_INVERTER_CONFIG_NC];

	/* Augmented states */
	float xa[8];

	/* Auxiliary variables for intermediate computations */
	float auxm1[DMPC_INVERTER_CONFIG_NC_x_NU];
	float auxm2[DMPC_INVERTER_CONFIG_NC_x_NU];

	xa[0] = x[0] - x_1[0];
	xa[1] = x[1] - x_1[1];
	xa[2] = x[2] - x_1[2];
	xa[3] = x[3] - x_1[3];
	xa[4] = x[4] - x_1[4];
	xa[5] = x[5] - x_1[5];
	xa[6] = x[0];
	xa[7] = x[1];

	//DMPC_INVERTER_CONFIG_U_MIN[2] = u_1[2];
	//DMPC_INVERTER_CONFIG_U_MAX[2] = u_1[2];

	//DMPC_INVERTER_CONFIG_U_MIN[3] = u_1[3];
	//DMPC_INVERTER_CONFIG_U_MAX[3] = u_1[3];
	/*
	 * Computes Fj matrix. This matrix is given by:
	 * Fj = Fj_1 * r + Fj_2 * xa,
	 *
	 * Fj_1 and Fj_2 are given by (both are computed off-line):
	 * Fj_1 = -Phi.T * R_s_bar,
	 * Fj_2 =  Phi.T * F
	 */
	mulmv((float *)Fj_1, DMPC_INVERTER_CONFIG_NC_x_NU, r, DMPC_INVERTER_CONFIG_NY, auxm1);
	mulmv((float *)Fj_2, DMPC_INVERTER_CONFIG_NC_x_NU, xa, DMPC_INVERTER_CONFIG_N_AUG_STATES, auxm2);
	sumv(auxm1, auxm2, DMPC_INVERTER_CONFIG_NC_x_NU, Fj);

	/*
	 * Computes the gam vector (or y vector). This vector holds the control
	 * and state inequalities.
	 */

	/* We start by assembling the control inequalities */
	j = 0;
	for(i = 0; i < DMPC_INVERTER_CONFIG_NR; i++){
		for(k = 0; k < DMPC_INVERTER_CONFIG_NU; k++){
			gam[j++] = -DMPC_INVERTER_CONFIG_U_MIN[k] + u_1[k];
		}
	}
	for(i = 0; i < DMPC_INVERTER_CONFIG_NR; i++){
		for(k = 0; k < DMPC_INVERTER_CONFIG_NU; k++){
			gam[j++] =  DMPC_INVERTER_CONFIG_U_MAX[k] - u_1[k];
		}
	}

	/* Now, the state inequalities */
	mulmv((float *)Fx, DMPC_INVERTER_CONFIG_NR * DMPC_INVERTER_CONFIG_N_STATE_CONSTRAINTS, xa, DMPC_INVERTER_CONFIG_N_STATES, auxm1);
	w = 0;
	for(i = 0; i < DMPC_INVERTER_CONFIG_NR; i++){
		for( k = 0; k < DMPC_INVERTER_CONFIG_N_STATE_CONSTRAINTS; k++){
			gam[j++] = -DMPC_INVERTER_CONFIG_XM_MIN[k] + x[k + 2] + auxm1[w++];
		}
	}
	w = 0;
	for(i = 0; i < DMPC_INVERTER_CONFIG_NR; i++){
		for( k = 0; k < DMPC_INVERTER_CONFIG_N_STATE_CONSTRAINTS; k++){
			gam[j++] =  DMPC_INVERTER_CONFIG_XM_MAX[k] - x[k + 2] - auxm1[w++];
		}
	}

	dmpcInverterHildOpt(du, iters);
	//du = dmpcBuckCVXGENOpt(iters);
}
//-----------------------------------------------------------------------------
//=============================================================================

//=============================================================================
/*---------------------------- Static functions -----------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
static void dmpcInverterHildOpt(float *du, uint32_t* iters){

	uint32_t n_iter;

	/* Auxiliary variables for intermediate computations */
	float auxm1[DMPC_INVERTER_CONFIG_NLAMBDA];
	float aux1;

	/* Matrices and vectors */
	float Kj[DMPC_INVERTER_CONFIG_NLAMBDA];
	float lambda[DMPC_INVERTER_CONFIG_NLAMBDA];
	//static float lambda[DMPC_INVERTER_CONFIG_NLAMBDA] = {0};

	/* Optimal control increment */
	//float DU;

	/* Computes Kj */
	mulmv((float *)Kj_1, DMPC_INVERTER_CONFIG_NLAMBDA, Fj, DMPC_INVERTER_CONFIG_NC_x_NU, (float *)auxm1);
	sumv(gam, (float *)auxm1, DMPC_INVERTER_CONFIG_NLAMBDA, Kj);

	/* Opt */
	n_iter = qpHild((float *)Hj, Kj, 20000, lambda, DMPC_INVERTER_CONFIG_NLAMBDA, (float)1e-4);
	//n_iter = qpHild4((float *)Hj, Kj, 15, lambda, (float)1e-6);
	//qpHild4FixedIter((float *)Hj, Kj, 8, lambda);
	//n_iter = 15;

    if( iters != 0 ) *iters = n_iter;

	/* DU */
	mulmv((float *)DU_1, DMPC_INVERTER_CONFIG_NU, Fj, DMPC_INVERTER_CONFIG_NC_x_NU, du);
	mulmv((float *)DU_2, DMPC_INVERTER_CONFIG_NU, lambda, DMPC_INVERTER_CONFIG_NLAMBDA, auxm1);
	sumv(du, auxm1, DMPC_INVERTER_CONFIG_NU, du);

}
//-----------------------------------------------------------------------------
//=============================================================================
