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
#include "dmpc_defs.h"

#include "mvops.h"

#define DMPC_CONFIG_SOLVER_HILD

/* Hildreth's QP */
#ifdef DMPC_CONFIG_SOLVER_HILD
#include "qp.h"
#endif

/* OSQP */
#ifdef DMPC_CONFIG_SOLVER_OSQP
#include "osqp/workspace.h"
#include "osqp/osqp.h"
#endif
//=============================================================================

//=============================================================================
/*-------------------------------- Prototypes -------------------------------*/
//=============================================================================
#ifdef DMPC_CONFIG_SOLVER_HILD
static uint32_t dmpcHildOpt(float *du);
#endif

#ifdef DMPC_CONFIG_SOLVER_OSQP
static uint32_t dmpcOSQP(float *du);
#endif
//=============================================================================

//=============================================================================
/*-------------------------------- Functions --------------------------------*/
//=============================================================================
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
    for(i = 0; i < DMPC_CONFIG_NXM; i++){
        xa[i] = x[i] - x_1[i];
    }
    for(i = 0; i < (DMPC_CONFIG_NXA - DMPC_CONFIG_NXM); i++){
        xa[DMPC_CONFIG_NXM + i] = x[DMPC_CONFIG_Y_IDX[i]];
    }
    
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

#ifdef DMPC_CONFIG_SOLVER_HILD
	/*
	 * Computes the gam vector (or y vector). This vector holds the control
	 * and state inequalities.
	 */

	/* We start by assembling the control inequalities */
	j = 0;
#if ( DMPC_CONFIG_NU_CTR != 0 )
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
#endif

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
#endif

#ifdef DMPC_CONFIG_SOLVER_OSQP
	/*
	 * Computes the gam vector (or y vector). This vector holds the control
	 * and state inequalities.
	 */

	/* We start by assembling the control inequalities */
	j = 0;
#if ( DMPC_CONFIG_NU_CTR != 0 )
	for(i = 0; i < DMPC_CONFIG_NR; i++){
		for(k = 0; k < DMPC_CONFIG_NU; k++){
			ldata[j] = DMPC_CONFIG_U_MIN[k] - u_1[k];
			udata[j] = DMPC_CONFIG_U_MAX[k] - u_1[k];
			j++;
		}
	}
#endif

	/* Now, the state inequalities */
#if ( DMPC_CONFIG_NXM_CTR != 0 )
    mulmv((float *)DMPC_M_Fx, DMPC_CONFIG_NR * DMPC_CONFIG_NXM_CTR, xa, DMPC_CONFIG_NXM, auxm1);
    w = 0;
    for(i = 0; i < DMPC_CONFIG_NR; i++){
        for( k = 0; k < DMPC_CONFIG_NXM_CTR; k++){
            ldata[j] = DMPC_CONFIG_XM_MIN[k] - x[DMPC_CONFIG_XM_LIM_IDX[k]] - auxm1[w];
            udata[j] = DMPC_CONFIG_XM_MAX[k] - x[DMPC_CONFIG_XM_LIM_IDX[k]] - auxm1[w];
            j++;
            w++;
        }
    }
#endif

	iters = dmpcOSQP(du);
    if( niters != 0 ) *niters = iters;
#endif

	return 0;
}
//-----------------------------------------------------------------------------
void dmpcDelayComp(float *x_1, float *x, float *u){

    float aux1[DMPC_CONFIG_NXM], aux2[DMPC_CONFIG_NXM];

    mulmv((float *)DMPC_M_A, DMPC_CONFIG_NXM, x, DMPC_CONFIG_NXM, aux1);
    mulmv((float *)DMPC_M_B, DMPC_CONFIG_NXM, u, DMPC_CONFIG_NU+DMPC_CONFIG_ND, aux2);
    sumv(aux1, aux2, DMPC_CONFIG_NXM, x_1);
}
//-----------------------------------------------------------------------------
//=============================================================================

//=============================================================================
/*---------------------------- Static functions -----------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
#ifdef DMPC_CONFIG_SOLVER_HILD
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
	niter = qpHild((float *)DMPC_M_Hj, Kj, 200, lambda, DMPC_CONFIG_NLAMBDA, (float)1e-9);
	//niter = qpHildFixedIter((float *)DMPC_M_Hj, Kj, 30, lambda, DMPC_CONFIG_NLAMBDA);

	/* Optimal control increment */
	mulmv((float *)DMPC_M_DU_1, DMPC_CONFIG_NU, DMPC_M_Fj, DMPC_CONFIG_NC_x_NU, du);
	mulmv((float *)DMPC_M_DU_2, DMPC_CONFIG_NU, lambda, DMPC_CONFIG_NLAMBDA, auxm1);
	sumv(du, auxm1, DMPC_CONFIG_NU, du);
    
    return niter;
}
#endif
//-----------------------------------------------------------------------------
#ifdef DMPC_CONFIG_SOLVER_OSQP
static uint32_t dmpcOSQP(float *du){

    uint32_t i;

	osqp_update_bounds(&workspace, ldata, udata);
	osqp_update_lin_cost(&workspace, DMPC_M_Fj);

	osqp_solve(&workspace);
    
    for(i = 0; i < DMPC_CONFIG_NU; i++){
	    du[i] = workspace.solution->x[0];    
    }

	return workspace.info->iter;
}
#endif
//-----------------------------------------------------------------------------
//=============================================================================
