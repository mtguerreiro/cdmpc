
//=============================================================================
/*-------------------------------- Includes ---------------------------------*/
//=============================================================================
#include "dmpc_osqp.h"

#include "mvops.h"
#include "dmpc_data.h"

#include "osqp/workspace.h"
#include "osqp/osqp.h"
//=============================================================================

//=============================================================================
/*-------------------------------- Functions --------------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
int32_t dmpc_osqp_solve(float *du){

    int32_t niter;

    uint32_t i, j, k, w;

    j = 0;
    if(dmpc_data.l_u_cnt != 0 ){
        for(i = 0; i < dmpc_data.l_u_cnt; i++){
            for(k = 0; k < dmpc_data.nu; k++){
                dmpc_osqp_data.ldata[j] = -dmpc_data.gam[k];
                dmpc_osqp_data.udata[j] = dmpc_data.gam[dmpc_data.nu*dmpc_data.l_u_cnt+k];
                j++;
            }
        }
    }

    if(dmpc_data.l_x_cnt != 0 ){
        w = 2*dmpc_data.nu*dmpc_data.l_u_cnt;
        for(i = 0; i < dmpc_data.l_x_cnt; i++){
            for( k = 0; k < dmpc_data.n_x_cnt; k++){
                dmpc_osqp_data.ldata[j] = -dmpc_data.gam[k + w];
                dmpc_osqp_data.udata[j] = dmpc_data.gam[dmpc_data.n_x_cnt*dmpc_data.l_x_cnt+k + w];
                j++;
            }
        }
    }

    osqp_update_data_vec(&solver, dmpc_data.Fj, dmpc_osqp_data.ldata, dmpc_osqp_data.udata);

    osqp_solve(&solver);

    for(i = 0; i < dmpc_data.nu; i++){
        du[i] = solver.solution->x[i];
    }

    niter = (int32_t)solver.info->iter;

    return niter;
}
//-----------------------------------------------------------------------------
//=============================================================================
