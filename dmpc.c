
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
    dmpc_data_t *pdata = inst->prob_data;

    /* Assembles -dx state vector */
    for(i = 0; i < pdata->n_xm; i++){
        pdata->dx[i] = -(pdata->x[i] - pdata->x_1[i]);
    }

    /* Assembles -error vector */
    for(i = 0; i < pdata->ny; i++){
        pdata->e[i] = -( pdata->x[pdata->y_idx[i]] - pdata->r[i] );
    }

    mulmv((float *)pdata->Kx, pdata->nu, pdata->dx, pdata->n_xm, pdata->auxm1);
    mulmv((float *)pdata->Ky, pdata->nu, pdata->e, pdata->ny, pdata->auxm2);
    sumv(pdata->auxm1, pdata->auxm2, pdata->nu, pdata->du);

    pdata->n_iters = 0;

    return 0;
}
//-----------------------------------------------------------------------------
static int32_t dmpcOptConstrained(dmpc_inst_t *inst){

    int32_t status;
    uint32_t i, j, k, w;
    dmpc_data_t *pdata = inst->prob_data;

    /* Assembles augmented state vector */
    for(i = 0; i < pdata->n_xm; i++){
        pdata->xa[i] = pdata->x[i] - pdata->x_1[i];
    }
    for(i = 0; i < (pdata->n_xa - pdata->n_xm); i++){
        pdata->xa[pdata->n_xm + i] = pdata->x[pdata->y_idx[i]];
    }

    /*
     * Computes Fj matrix. This matrix is given by:
     * Fj = Fj_1 * r + Fj_2 * xa,
     *
     * Fj_1 and Fj_2 are given by (both are computed off-line):
     * Fj_1 = -Phi.T * R_s_bar,
     * Fj_2 =  Phi.T * F
     */
    mulmv((float *)pdata->Fj_1, pdata->u_size, pdata->r, pdata->ny, pdata->auxm1);
    mulmv((float *)pdata->Fj_2, pdata->u_size, pdata->xa, pdata->n_xa, pdata->auxm2);
    sumv(pdata->auxm1, pdata->auxm2, pdata->u_size, pdata->Fj);

    /*
     * Computes the gam vector (or y vector). This vector holds the control
     * and state inequalities.
     */

    /* We start by assembling the control inequalities */
    j = 0;
    if(pdata->l_u_cnt != 0 ){
        for(i = 0; i < pdata->l_u_cnt; i++){
            for(k = 0; k < pdata->nu; k++){
                pdata->gam[j++] = -pdata->u_min[k] + pdata->u_1[k];
            }
        }
        for(i = 0; i < pdata->l_u_cnt; i++){
            for(k = 0; k < pdata->nu; k++){
                pdata->gam[j++] =  pdata->u_max[k] - pdata->u_1[k];
            }
        }
    }

    /* Now, the state inequalities */
    if(pdata->l_x_cnt != 0 ){
        mulmv((float *)pdata->Fx, pdata->l_x_cnt * pdata->n_x_cnt, pdata->xa, pdata->n_xm, pdata->auxm1);
        w = 0;
        for(i = 0; i < pdata->l_x_cnt; i++){
            for( k = 0; k < pdata->n_x_cnt; k++){
                pdata->gam[j++] = -pdata->x_min[k] + pdata->x[pdata->x_cnt_idx[k]] + pdata->auxm1[w++];
            }
        }
        w = 0;
        for(i = 0; i < pdata->l_x_cnt; i++){
            for( k = 0; k < pdata->n_x_cnt; k++){
                pdata->gam[j++] =  pdata->x_max[k] - pdata->x[pdata->x_cnt_idx[k]] - pdata->auxm1[w++];
            }
        }
    }

    status = inst->solve(inst->prob_data, inst->solver_data);

    return status;
}
//-----------------------------------------------------------------------------
//=============================================================================
