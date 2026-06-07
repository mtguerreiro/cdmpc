
//=============================================================================
/*-------------------------------- Includes ---------------------------------*/
//=============================================================================
#include "dmpc.h"

#include "mvops.h"
//=============================================================================

//=============================================================================
/*-------------------------------- Prototypes -------------------------------*/
//=============================================================================
static int32_t dmpcOptUnconstrained(dmpc_inst_t *inst);
static int32_t dmpcOptConstrained(dmpc_inst_t *inst);
//=============================================================================

//=============================================================================
/*-------------------------------- Functions --------------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
int32_t dmpcOpt(dmpc_inst_t *inst){

    int32_t status;

    if( (inst->prob_data->l_u_cnt == 0) && (inst->prob_data->l_x_cnt == 0) )
        status = dmpcOptUnconstrained(inst);
    else
        status = dmpcOptConstrained(inst);

    return status;
}
//-----------------------------------------------------------------------------
// void dmpcDelayComp(float *x_1, float *x, float *u){
//
//     mulmv((float *)inst->prob_data->A, inst->prob_data->n_xm, x, inst->prob_data->n_xm, inst->prob_data->auxm1);
//     mulmv((float *)inst->prob_data->B, inst->prob_data->n_xm, u, inst->prob_data->nu+inst->prob_data->nd, inst->prob_data->auxm2);
//     sumv(inst->prob_data->auxm1, inst->prob_data->auxm2, inst->prob_data->n_xm, x_1);
// }
//-----------------------------------------------------------------------------
//=============================================================================

//=============================================================================
/*---------------------------- Static functions -----------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
static int32_t dmpcOptUnconstrained(dmpc_inst_t *inst){

    uint32_t i;

    /* Assembles -dx state vector */
    for(i = 0; i < inst->prob_data->n_xm; i++){
        inst->prob_data->dx[i] = -(inst->prob_data->x[i] - inst->prob_data->x_1[i]);
    }

    /* Assembles -error vector */
    for(i = 0; i < inst->prob_data->ny; i++){
        inst->prob_data->e[i] = -( inst->prob_data->x[inst->prob_data->y_idx[i]] - inst->prob_data->r[i] );
    }

    mulmv((float *)inst->prob_data->Kx, inst->prob_data->nu, inst->prob_data->dx, inst->prob_data->n_xm, inst->prob_data->auxm1);
    mulmv((float *)inst->prob_data->Ky, inst->prob_data->nu, inst->prob_data->e, inst->prob_data->ny, inst->prob_data->auxm2);
    sumv(inst->prob_data->auxm1, inst->prob_data->auxm2, inst->prob_data->nu, inst->prob_data->du);

    inst->prob_data->n_iters = 0;

    return 0;
}
//-----------------------------------------------------------------------------
static int32_t dmpcOptConstrained(dmpc_inst_t *inst){

    int32_t status;
    uint32_t i, j, k, w;

    /* Assembles augmented state vector */
    for(i = 0; i < inst->prob_data->n_xm; i++){
        inst->prob_data->xa[i] = inst->prob_data->x[i] - inst->prob_data->x_1[i];
    }
    for(i = 0; i < (inst->prob_data->n_xa - inst->prob_data->n_xm); i++){
        inst->prob_data->xa[inst->prob_data->n_xm + i] = inst->prob_data->x[inst->prob_data->y_idx[i]];
    }

    /*
     * Computes Fj matrix. This matrix is given by:
     * Fj = Fj_1 * r + Fj_2 * xa,
     *
     * Fj_1 and Fj_2 are given by (both are computed off-line):
     * Fj_1 = -Phi.T * R_s_bar,
     * Fj_2 =  Phi.T * F
     */
    mulmv((float *)inst->prob_data->Fj_1, inst->prob_data->u_size, inst->prob_data->r, inst->prob_data->ny, inst->prob_data->auxm1);
    mulmv((float *)inst->prob_data->Fj_2, inst->prob_data->u_size, inst->prob_data->xa, inst->prob_data->n_xa, inst->prob_data->auxm2);
    sumv(inst->prob_data->auxm1, inst->prob_data->auxm2, inst->prob_data->u_size, inst->prob_data->Fj);

    /*
     * Computes the gam vector (or y vector). This vector holds the control
     * and state inequalities.
     */

    /* We start by assembling the control inequalities */
    j = 0;
    if(inst->prob_data->l_u_cnt != 0 ){
        for(i = 0; i < inst->prob_data->l_u_cnt; i++){
            for(k = 0; k < inst->prob_data->nu; k++){
                inst->prob_data->gam[j++] = -inst->prob_data->u_min[k] + inst->prob_data->u_1[k];
            }
        }
        for(i = 0; i < inst->prob_data->l_u_cnt; i++){
            for(k = 0; k < inst->prob_data->nu; k++){
                inst->prob_data->gam[j++] =  inst->prob_data->u_max[k] - inst->prob_data->u_1[k];
            }
        }
    }

    /* Now, the state inequalities */
    if(inst->prob_data->l_x_cnt != 0 ){
        mulmv((float *)inst->prob_data->Fx, inst->prob_data->l_x_cnt * inst->prob_data->n_x_cnt, inst->prob_data->xa, inst->prob_data->n_xm, inst->prob_data->auxm1);
        w = 0;
        for(i = 0; i < inst->prob_data->l_x_cnt; i++){
            for( k = 0; k < inst->prob_data->n_x_cnt; k++){
                inst->prob_data->gam[j++] = -inst->prob_data->x_min[k] + inst->prob_data->x[inst->prob_data->x_cnt_idx[k]] + inst->prob_data->auxm1[w++];
            }
        }
        w = 0;
        for(i = 0; i < inst->prob_data->l_x_cnt; i++){
            for( k = 0; k < inst->prob_data->n_x_cnt; k++){
                inst->prob_data->gam[j++] =  inst->prob_data->x_max[k] - inst->prob_data->x[inst->prob_data->x_cnt_idx[k]] - inst->prob_data->auxm1[w++];
            }
        }
    }

    status = inst->solve(inst->prob_data, inst->solver_data);

    return status;
}
//-----------------------------------------------------------------------------
//=============================================================================
