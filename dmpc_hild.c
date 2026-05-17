
//=============================================================================
/*-------------------------------- Includes ---------------------------------*/
//=============================================================================
#include "dmpc_hild.h"

#include "mvops.h"
#include "dmpc_data.h"
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
int32_t dmpc_hild_solve(float *du){

    int32_t niter;

    /* Computes Kj */
    mulmv((float *)dmpc_hild_data.Kj_1, dmpc_hild_data.n_lambda, dmpc_data.Fj, dmpc_data.u_size, dmpc_hild_data.aux);
    sumv(dmpc_data.gam, dmpc_hild_data.aux, dmpc_hild_data.n_lambda, dmpc_hild_data.Kj);

    /* Opt */
#if (DMPC_CONFIG_HILD_FIXED_ITER == 0)
    niter = (int32_t) dmpc_hild_qp_solve((float *)dmpc_hild_data.Hj, dmpc_hild_data.Kj, DMPC_CONFIG_HILD_N_ITER, dmpc_hild_data.lambda, dmpc_hild_data.n_lambda, (float)DMPC_CONFIG_HILD_TOL);
#else
    dmpc_hild_qp_solve_fixed_iter((float *)dmpc_hild_data.Hj, dmpc_hild_data.Kj, DMPC_CONFIG_HILD_N_ITER, dmpc_hild_data.lambda, dmpc_hild_data.n_lambda);
    niter = (int32_t)DMPC_CONFIG_HILD_N_ITER;
#endif

    /* Optimal control increment */
    mulmv((float *)dmpc_hild_data.DU_1, dmpc_data.nu, dmpc_data.Fj, dmpc_data.u_size, du);
    mulmv((float *)dmpc_hild_data.DU_2, dmpc_data.nu, dmpc_hild_data.lambda, dmpc_hild_data.n_lambda, dmpc_hild_data.aux);
    sumv(du, dmpc_hild_data.aux, dmpc_data.nu, du);

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
