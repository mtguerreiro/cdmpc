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
	//float du[DIM_CONFIG_NC];

	/* Augmented states */
	float xa[8];

	/* Auxiliary variables for intermediate computations */
	float auxm1[DIM_CONFIG_NC_x_NU];
	float auxm2[DIM_CONFIG_NC_x_NU];

	xa[0] = x[0] - x_1[0];
	xa[1] = x[1] - x_1[1];
	xa[2] = x[2] - x_1[2];
	xa[3] = x[3] - x_1[3];
	xa[4] = x[4] - x_1[4];
	xa[5] = x[5] - x_1[5];
	xa[6] = x[0];
	xa[7] = x[1];

	//DIM_CONFIG_U_MIN[2] = u_1[2];
	//DIM_CONFIG_U_MAX[2] = u_1[2];

	//DIM_CONFIG_U_MIN[3] = u_1[3];
	//DIM_CONFIG_U_MAX[3] = u_1[3];
	/*
	 * Computes Fj matrix. This matrix is given by:
	 * Fj = Fj_1 * r + Fj_2 * xa,
	 *
	 * Fj_1 and Fj_2 are given by (both are computed off-line):
	 * Fj_1 = -Phi.T * R_s_bar,
	 * Fj_2 =  Phi.T * F
	 */
	mulmv((float *)DIM_Fj_1, DIM_CONFIG_NC_x_NU, r, DIM_CONFIG_NY, auxm1);
	mulmv((float *)DIM_Fj_2, DIM_CONFIG_NC_x_NU, xa, DIM_CONFIG_NXA, auxm2);
	sumv(auxm1, auxm2, DIM_CONFIG_NC_x_NU, DIM_Fj);

	/*
	 * Computes the gam vector (or y vector). This vector holds the control
	 * and state inequalities.
	 */

	/* We start by assembling the control inequalities */
	j = 0;
	for(i = 0; i < DIM_CONFIG_NR; i++){
		for(k = 0; k < DIM_CONFIG_NU; k++){
			DIM_gam[j++] = -DIM_CONFIG_U_MIN[k] + u_1[k];
		}
	}
	for(i = 0; i < DIM_CONFIG_NR; i++){
		for(k = 0; k < DIM_CONFIG_NU; k++){
			DIM_gam[j++] =  DIM_CONFIG_U_MAX[k] - u_1[k];
		}
	}

	/* Now, the state inequalities */
#if ( DIM_CONFIG_NXM_CTR != 0 )
    mulmv((float *)DIM_Fx, DIM_CONFIG_NR * DIM_CONFIG_NXM_CTR, xa, DIM_CONFIG_NXM, auxm1);
    w = 0;
    for(i = 0; i < DIM_CONFIG_NR; i++){
        for( k = 0; k < DIM_CONFIG_NXM_CTR; k++){
            DIM_gam[j++] = -DIM_CONFIG_XM_MIN[k] + x[DIM_CONFIG_XM_LIM_IDX[k]] + auxm1[w++];
        }
    }
    w = 0;
    for(i = 0; i < DIM_CONFIG_NR; i++){
        for( k = 0; k < DIM_CONFIG_NXM_CTR; k++){
            DIM_gam[j++] =  DIM_CONFIG_XM_MAX[k] - x[DIM_CONFIG_XM_LIM_IDX[k]] - auxm1[w++];
        }
    }   
#endif

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
	float auxm1[DIM_CONFIG_NLAMBDA];
	float aux1;

	/* Matrices and vectors */
	float Kj[DIM_CONFIG_NLAMBDA];
	float lambda[DIM_CONFIG_NLAMBDA];
	//static float lambda[DIM_CONFIG_NLAMBDA] = {0};

	/* Optimal control increment */
	//float DU;

	/* Computes Kj */
	mulmv((float *)DIM_Kj_1, DIM_CONFIG_NLAMBDA, DIM_Fj, DIM_CONFIG_NC_x_NU, (float *)auxm1);
	sumv(DIM_gam, (float *)auxm1, DIM_CONFIG_NLAMBDA, Kj);

	/* Opt */
	n_iter = qpHild((float *)DIM_Hj, Kj, 1000, lambda, DIM_CONFIG_NLAMBDA, (float)0.00005e-3);
	//n_iter = qpHildFixedIter((float *)DIM_Hj, Kj, 5, lambda, DIM_CONFIG_NLAMBDA);
	//n_iter = qpHild4((float *)Hj, Kj, 15, lambda, (float)1e-6);
	//qpHild4FixedIter((float *)Hj, Kj, 8, lambda);
	//n_iter = 15;

    if( iters != 0 ) *iters = n_iter;

	/* DU */
	mulmv((float *)DIM_DU_1, DIM_CONFIG_NU, DIM_Fj, DIM_CONFIG_NC_x_NU, du);
	mulmv((float *)DIM_DU_2, DIM_CONFIG_NU, lambda, DIM_CONFIG_NLAMBDA, auxm1);
	sumv(du, auxm1, DIM_CONFIG_NU, du);

}
//-----------------------------------------------------------------------------
//=============================================================================
