
//=============================================================================
/*-------------------------------- Includes ---------------------------------*/
//=============================================================================
#include "dmpc_hild.h"
#include "mvops.h"
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

    uint32_t n_iters;
    dmpc_hild_data_t *hdata = (dmpc_hild_data_t *)hild_data;

    /* Computes Kj */
    mulmv(hdata->Kj_1, hdata->n_lambda, prob_data->Fj, prob_data->u_size, hdata->aux);
    sumv(prob_data->gam, hdata->aux, hdata->n_lambda, hdata->Kj);

    /* Opt */
    if( hdata->fixed_iter ){
        dmpc_hild_qp_solve_fixed_iter(hdata->Hj, hdata->Kj, hdata->n_iter, hdata->lambda, hdata->n_lambda);
        n_iters = hdata->n_iter;
    }
    else{
        n_iters = dmpc_hild_qp_solve(hdata->Hj, hdata->Kj, hdata->max_iter, hdata->lambda, hdata->n_lambda, hdata->tol);
    }

    prob_data->n_iters = n_iters;

    /* Optimal control increment */
    mulmv((float *)hdata->DU_1, prob_data->nu, prob_data->Fj, prob_data->u_size, prob_data->du);
    mulmv((float *)hdata->DU_2, prob_data->nu, hdata->lambda, hdata->n_lambda, hdata->aux);
    sumv(prob_data->du, hdata->aux, prob_data->nu, prob_data->du);

    return 0;
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
