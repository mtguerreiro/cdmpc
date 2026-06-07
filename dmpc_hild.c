
//=============================================================================
/*-------------------------------- Includes ---------------------------------*/
//=============================================================================
#include "dmpc_hild.h"

#include "mvops.h"
#include "dmpc_data.h"
#include "dmpc_data_hild.h"
//=============================================================================

//=============================================================================
/*-------------------------------- Prototypes -------------------------------*/
//=============================================================================
static uint32_t dmpc_hild_qp_solve(float *H, float *K, uint32_t n_iter, float* lambda, uint32_t lambda_size, float tol);
static void dmpc_hild_qp_solve_fixed_iter(float *H, float *K, uint32_t n_iter, float* lambda, uint32_t lambda_size);
//=============================================================================

//=============================================================================
/*-------------------------------- Functions --------------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
int32_t dmpc_hild_solve(dmpc_data_t *prob_data, void *hild_data){

    int32_t niter;
    dmpc_hild_data_t *hdata = (dmpc_hild_data_t *)hild_data;

    /* Computes Kj */
    mulmv((float *)hdata->Kj_1, hdata->n_lambda, prob_data->Fj, prob_data->u_size, hdata->aux);
    sumv(prob_data->gam, hdata->aux, hdata->n_lambda, hdata->Kj);

    /* Opt */
#if (DMPC_CONFIG_HILD_FIXED_ITER == 0)
    niter = (int32_t) dmpc_hild_qp_solve((float *)hdata->Hj, hdata->Kj, DMPC_CONFIG_HILD_N_ITER, hdata->lambda, hdata->n_lambda, (float)DMPC_CONFIG_HILD_TOL);
#else
    dmpc_hild_qp_solve_fixed_iter((float *)hdata->Hj, hdata->Kj, DMPC_CONFIG_HILD_N_ITER, hdata->lambda, hdata->n_lambda);
    niter = (int32_t)DMPC_CONFIG_HILD_N_ITER;
#endif

    /* Optimal control increment */
    mulmv((float *)hdata->DU_1, prob_data->nu, prob_data->Fj, prob_data->u_size, prob_data->du);
    mulmv((float *)hdata->DU_2, prob_data->nu, hdata->lambda, hdata->n_lambda, hdata->aux);
    sumv(prob_data->du, hdata->aux, prob_data->nu, prob_data->du);

    return niter;
}
//-----------------------------------------------------------------------------
//=============================================================================

//=============================================================================
/*---------------------------- Static functions -----------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
static uint32_t dmpc_hild_qp_solve(float *H, float *K, uint32_t n_iter, float* lambda, uint32_t lambda_size, float tol){

    uint32_t stopCond;
    float res;
    float *h;
    float *k;
    float acc;

    uint32_t n, i, j;

    for(i = 0; i < lambda_size; i++){
        lambda[i] = 0;
    }

    n = 0;
    k = K;
    while( n < n_iter ){

        h = H;
        stopCond = 0;

        /* One iteration through the lambda vector */
        for(i = 0; i < lambda_size; i++){

            res = lambda[i];
            acc = 0;
            lambda[i] = 0;
            for(j = 0; j < lambda_size; j++){
                acc += h[j] * lambda[j];
            }

            lambda[i] = h[i] * (k[i] + acc);
            if( lambda[i] < 0 ) lambda[i] = 0;

            if( stopCond == 0 ){
                if( res > lambda[i] ) res = res - lambda[i];
                else res = lambda[i] - res;
                if( res > tol ) stopCond = 1;
            }

            h = h + lambda_size;
        }

        n++;

        /*
         * If stopCond is zero, all lambdas had a difference between the new
         * and the previous value that is smaller than the tolerance; so we
         * have reached the stopping condition.
         */
        if( stopCond == 0 ) break;
    }

    return n;
}
//-----------------------------------------------------------------------------
static void dmpc_hild_qp_solve_fixed_iter(float *H, float *K, uint32_t n_iter, float* lambda, uint32_t lambda_size){

    float *h;
    float *k;
    float acc;

    uint32_t n, i, j;

    for(i = 0; i < lambda_size; i++){
        lambda[i] = 0;
    }

    n = 0;
    k = K;
    while( n < n_iter ){

        h = H;

        /* One iteration through the lambda vector */
        for(i = 0; i < lambda_size; i++){

            lambda[i] = 0;
            acc = 0;
            for(j = 0; j < lambda_size; j++){
                acc += h[j] * lambda[j];
            }

            lambda[i] = h[i] * (k[i] + acc);
            if( lambda[i] < 0 ) lambda[i] = 0;

            h = h + lambda_size;
        }

        n++;
    }
}
//-----------------------------------------------------------------------------
//=============================================================================
