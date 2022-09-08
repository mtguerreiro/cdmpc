/*
 * dmpc.c
 *
 *  Created on: 08.09.2022
 *      Author: mguerreiro
 */

//=============================================================================
/*-------------------------------- Includes ---------------------------------*/
//=============================================================================
#include "dmpc.h"
#include "dmpc_matrices.h"

/* Hildreth's QP */
#include "mvops.h"
#include "qp.h"

//=============================================================================

//=============================================================================
/*-------------------------------- Prototypes -------------------------------*/
//=============================================================================
static uint32_t dmpcHildOpt(float *du);
//=============================================================================

//=============================================================================
/*-------------------------------- Functions --------------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
void dmpcInitialize(void){

}
//-----------------------------------------------------------------------------
uint32_t dmpcOpt(float *x, float *x_1, float *r, float *u_1, uint32_t *niters, float *du){

	uint32_t i, j, k, w;
    
    uint32_t iters;

	/* Augmented states */
	float xa[DMPC_CONFIG_NXA];

	/* Auxiliary variables for intermediate computations */
	float auxm1[DMPC_CONFIG_NC_x_NU];
	float auxm2[DMPC_CONFIG_NC_x_NU];

    /* Assembles augmented state vector */
    /*
    for(i = 0; i < DMPC_CONFIG_NXM; i++){
        xa[i] = x[i] - x_1[i];
    }
    for(i = DMPC_CONFIG_NXM; i < DMPC_CONFIG_NXA; i++){
        xa[i] = x[i - DMPC_CONFIG_NXM];
    }
    */
    xa[0] = x[0] - x_1[0];
    xa[1] = x[1] - x_1[1];
    xa[2] = x[1];
	/*
	 * Computes Fj matrix. This matrix is given by:
	 * Fj = Fj_1 * r + Fj_2 * xa,
	 *
	 * Fj_1 and Fj_2 are given by (both are computed off-line):
	 * Fj_1 = -Phi.T * R_s_bar,
	 * Fj_2 =  Phi.T * F
	 */
	mulmv((float *)DMPC_M_Fj_1, DMPC_CONFIG_NC_x_NU, r, DMPC_CONFIG_NY, auxm1);
	mulmv((float *)DMPC_M_Fj_2, DMPC_CONFIG_NC_x_NU, xa, DMPC_CONFIG_NXA, auxm2);
	sumv(auxm1, auxm2, DMPC_CONFIG_NC_x_NU, DMPC_M_Fj);
    
	/*
	 * Computes the gam vector (or y vector). This vector holds the control
	 * and state inequalities.
	 */

	/* We start by assembling the control inequalities */
	j = 0;
	for(i = 0; i < DMPC_CONFIG_NR; i++){
		for(k = 0; k < DMPC_CONFIG_NU; k++){
			DMPC_M_gam[j++] = -DMPC_CONFIG_U_MIN[k] + u_1[k];
		}
	}
	for(i = 0; i < DMPC_CONFIG_NR; i++){
		for(k = 0; k < DMPC_CONFIG_NU; k++){
			DMPC_M_gam[j++] =  DMPC_CONFIG_U_MAX[k] - u_1[k];
		}
	}

	/* Now, the state inequalities */
#if ( DMPC_CONFIG_NXM_CTR != 0 )
    mulmv((float *)DMPC_M_Fx, DMPC_CONFIG_NR * DMPC_CONFIG_NXM_CTR, xa, DMPC_CONFIG_NXM, auxm1);
    w = 0;
    for(i = 0; i < DMPC_CONFIG_NR; i++){
        for( k = 0; k < DMPC_CONFIG_NXM_CTR; k++){
            DMPC_M_gam[j++] = -DMPC_CONFIG_XM_MIN[k] + x[DMPC_CONFIG_XM_LIM_IDX[k]] + auxm1[w++];
        }
    }
    w = 0;
    for(i = 0; i < DMPC_CONFIG_NR; i++){
        for( k = 0; k < DMPC_CONFIG_NXM_CTR; k++){
            DMPC_M_gam[j++] =  DMPC_CONFIG_XM_MAX[k] - x[DMPC_CONFIG_XM_LIM_IDX[k]] - auxm1[w++];
        }
    }   
#endif

	iters = dmpcHildOpt(du);
    if( niters != 0 ) *niters = iters;

	return 0;
}
//-----------------------------------------------------------------------------
//=============================================================================

//=============================================================================
/*---------------------------- Static functions -----------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
static uint32_t dmpcHildOpt(float *du){

	uint32_t niter;

	/* Auxiliary variables for intermediate computations */
	float auxm1[DMPC_CONFIG_NLAMBDA];

	/* Matrices and vectors */
	float Kj[DMPC_CONFIG_NLAMBDA];
	float lambda[DMPC_CONFIG_NLAMBDA];
	//static float lambda[DIM_CONFIG_NLAMBDA] = {0};

	/* Computes Kj */
	mulmv((float *)DMPC_M_Kj_1, DMPC_CONFIG_NLAMBDA, DMPC_M_Fj, DMPC_CONFIG_NC_x_NU, (float *)auxm1);
	sumv(DMPC_M_gam, (float *)auxm1, DMPC_CONFIG_NLAMBDA, Kj);

	/* Opt */
	niter = qpHild((float *)DMPC_M_Hj, Kj, 1000, lambda, DMPC_CONFIG_NLAMBDA, (float)1e-6);
	//n_iter = qpHildFixedIter((float *)DIM_Hj, Kj, 5, lambda, DIM_CONFIG_NLAMBDA);
	//n_iter = qpHild4((float *)Hj, Kj, 15, lambda, (float)1e-6);
	//qpHild4FixedIter((float *)Hj, Kj, 8, lambda);
	//n_iter = 15;

	/* Optimal control increment */
	mulmv((float *)DMPC_M_DU_1, DMPC_CONFIG_NU, DMPC_M_Fj, DMPC_CONFIG_NC_x_NU, du);
	mulmv((float *)DMPC_M_DU_2, DMPC_CONFIG_NU, lambda, DMPC_CONFIG_NLAMBDA, auxm1);
	sumv(du, auxm1, DMPC_CONFIG_NU, du);
    
    return niter;
}
//-----------------------------------------------------------------------------
//=============================================================================
