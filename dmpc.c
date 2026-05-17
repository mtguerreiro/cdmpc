
//=============================================================================
/*-------------------------------- Includes ---------------------------------*/
//=============================================================================
#include "dmpc.h"
#include "dmpc_data.h"

#include "mvops.h"

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
static uint32_t dmpcOptUnconstrained(float *x, float *x_1, float *r, float *u_1, uint32_t *niters, float *du);
uint32_t dmpcOptConstrained(float *x, float *x_1, float *r, float *u_1, uint32_t *niters, float *du);

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

    uint32_t status;

    if( (dmpc_data.l_u_cnt == 0) && (dmpc_data.l_x_cnt == 0) )
        status = dmpcOptUnconstrained(x, x_1, r, u_1, niters, du);
    else
        status = dmpcOptConstrained(x, x_1, r, u_1, niters, du);

    return status;
}
//-----------------------------------------------------------------------------
void dmpcDelayComp(float *x_1, float *x, float *u){

    mulmv((float *)dmpc_data.A, dmpc_data.n_xm, x, dmpc_data.n_xm, dmpc_data.auxm1);
    mulmv((float *)dmpc_data.B, dmpc_data.n_xm, u, dmpc_data.nu+dmpc_data.nd, dmpc_data.auxm2);
    sumv(dmpc_data.auxm1, dmpc_data.auxm2, dmpc_data.n_xm, x_1);
}
//-----------------------------------------------------------------------------
//=============================================================================

//=============================================================================
/*---------------------------- Static functions -----------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
static uint32_t dmpcOptUnconstrained(float *x, float *x_1, float *r, float *u_1, uint32_t *niters, float *du){

    uint32_t i;

    /* Assembles -dx state vector */
    for(i = 0; i < dmpc_data.n_xm; i++){
        dmpc_data.dx[i] = -(x[i] - x_1[i]);
    }

    /* Assembles -error vector */
    for(i = 0; i < dmpc_data.ny; i++){
        dmpc_data.e[i] = -( x[dmpc_data.y_idx[i]] - r[i] );
    }

    mulmv((float *)dmpc_data.Kx, dmpc_data.nu, dmpc_data.dx, dmpc_data.n_xm, dmpc_data.auxm1);
    mulmv((float *)dmpc_data.Ky, dmpc_data.nu, dmpc_data.e, dmpc_data.ny, dmpc_data.auxm2);
    sumv(dmpc_data.auxm1, dmpc_data.auxm2, dmpc_data.nu, du);

    if( niters != 0 ) *niters = 0;

    return 0;
}
//-----------------------------------------------------------------------------
uint32_t dmpcOptConstrained(float *x, float *x_1, float *r, float *u_1, uint32_t *niters, float *du){

    uint32_t i, j, k, w;

    uint32_t iters = 0;

    /* Assembles augmented state vector */
    for(i = 0; i < dmpc_data.n_xm; i++){
        dmpc_data.xa[i] = x[i] - x_1[i];
    }
    for(i = 0; i < (dmpc_data.n_xa - dmpc_data.n_xm); i++){
        dmpc_data.xa[dmpc_data.n_xm + i] = x[dmpc_data.y_idx[i]];
    }

    /*
     * Computes Fj matrix. This matrix is given by:
     * Fj = Fj_1 * r + Fj_2 * xa,
     *
     * Fj_1 and Fj_2 are given by (both are computed off-line):
     * Fj_1 = -Phi.T * R_s_bar,
     * Fj_2 =  Phi.T * F
     */
    mulmv((float *)dmpc_data.Fj_1, dmpc_data.u_size, r, dmpc_data.ny, dmpc_data.auxm1);
    mulmv((float *)dmpc_data.Fj_2, dmpc_data.u_size, dmpc_data.xa, dmpc_data.n_xa, dmpc_data.auxm2);
    sumv(dmpc_data.auxm1, dmpc_data.auxm2, dmpc_data.u_size, dmpc_data.Fj);

    /*
     * Computes the gam vector (or y vector). This vector holds the control
     * and state inequalities.
     */

    /* We start by assembling the control inequalities */
    j = 0;
    if(dmpc_data.l_u_cnt != 0 ){
        for(i = 0; i < dmpc_data.l_u_cnt; i++){
            for(k = 0; k < dmpc_data.nu; k++){
                dmpc_data.gam[j++] = -dmpc_data.u_min[k] + u_1[k];
            }
        }
        for(i = 0; i < dmpc_data.l_u_cnt; i++){
            for(k = 0; k < dmpc_data.nu; k++){
                dmpc_data.gam[j++] =  dmpc_data.u_max[k] - u_1[k];
            }
        }
    }

    /* Now, the state inequalities */
    if(dmpc_data.l_x_cnt != 0 ){
        mulmv((float *)dmpc_data.Fx, dmpc_data.l_x_cnt * dmpc_data.n_x_cnt, dmpc_data.xa, dmpc_data.n_xm, dmpc_data.auxm1);
        w = 0;
        for(i = 0; i < dmpc_data.l_x_cnt; i++){
            for( k = 0; k < dmpc_data.n_x_cnt; k++){
                dmpc_data.gam[j++] = -dmpc_data.x_min[k] + x[dmpc_data.x_cnt_idx[k]] + dmpc_data.auxm1[w++];
            }
        }
        w = 0;
        for(i = 0; i < dmpc_data.l_x_cnt; i++){
            for( k = 0; k < dmpc_data.n_x_cnt; k++){
                dmpc_data.gam[j++] =  dmpc_data.x_max[k] - x[dmpc_data.x_cnt_idx[k]] - dmpc_data.auxm1[w++];
            }
        }
    }

#ifdef DMPC_CONFIG_SOLVER_HILD
    iters = dmpcHildOpt(du);
#endif

    #ifdef DMPC_CONFIG_SOLVER_OSQP
    iters = dmpcOSQP(du);
#endif

    if( niters != 0 ) *niters = iters;

    return 0;
}
//-----------------------------------------------------------------------------
#ifdef DMPC_CONFIG_SOLVER_HILD
static uint32_t dmpcHildOpt(float *du){

	uint32_t niter;

	/* Computes Kj */
	mulmv((float *)dmpc_data.Kj_1, dmpc_data.n_lambda, dmpc_data.Fj, dmpc_data.u_size, dmpc_data.auxm1);
	sumv(dmpc_data.gam, dmpc_data.auxm1, dmpc_data.n_lambda, dmpc_data.Kj);

	/* Opt */
#if (DMPC_CONFIG_HILD_FIXED_ITER == 0)
	niter = qpHild((float *)dmpc_data.Hj, Kj, DMPC_CONFIG_HILD_N_ITER, lambda, dmpc_data.n_lambda, (float)DMPC_CONFIG_HILD_TOL);
#else
	niter = qpHildFixedIter((float *)dmpc_data.Hj, dmpc_data.Kj, DMPC_CONFIG_HILD_N_ITER, dmpc_data.lambda, dmpc_data.n_lambda);
#endif

	/* Optimal control increment */
	mulmv((float *)dmpc_data.DU_1, dmpc_data.nu, dmpc_data.Fj, dmpc_data.u_size, du);
	mulmv((float *)dmpc_data.DU_2, dmpc_data.nu, dmpc_data.lambda, dmpc_data.n_lambda, dmpc_data.auxm1);
	sumv(du, dmpc_data.auxm1, dmpc_data.nu, du);
    
    return niter;
}
#endif
//-----------------------------------------------------------------------------
#ifdef DMPC_CONFIG_SOLVER_OSQP
static uint32_t dmpcOSQP(float *du){

    uint32_t i, j, k, w;

    j = 0;
    if(dmpc_data.l_u_cnt != 0 ){
        for(i = 0; i < dmpc_data.l_u_cnt; i++){
            for(k = 0; k < dmpc_data.nu; k++){
                dmpc_data.ldata[j] = -dmpc_data.gam[k];
                dmpc_data.udata[j] = dmpc_data.gam[dmpc_data.nu*dmpc_data.l_u_cnt+k];
                j++;
            }
        }
    }

    if(dmpc_data.l_x_cnt != 0 ){
        w = 2*dmpc_data.nu*dmpc_data.l_u_cnt;
        for(i = 0; i < dmpc_data.l_x_cnt; i++){
            for( k = 0; k < dmpc_data.n_x_cnt; k++){
                dmpc_data.ldata[j] = -dmpc_data.gam[k + w];
                dmpc_data.udata[j] = dmpc_data.gam[dmpc_data.n_x_cnt*dmpc_data.l_x_cnt+k + w];
                j++;
            }
        }
    }

    osqp_update_data_vec(&solver, dmpc_data.Fj, dmpc_data.ldata, dmpc_data.udata);

    osqp_solve(&solver);

    for(i = 0; i < dmpc_data.nu; i++){
        du[i] = solver.solution->x[i];
    }

    return solver.info->iter;
}
#endif
//-----------------------------------------------------------------------------
//=============================================================================
