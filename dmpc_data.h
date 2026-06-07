
#ifndef DMPC_DATA_H_
#define DMPC_DATA_H_

#include "stdint.h"

typedef struct{
    float *x;
    float *x_1;
    float *r;
    float *u_1;
    float *du;
    uint32_t n_iters;
    uint32_t n_xm;
    uint32_t n_xa;
    uint32_t nu;
    uint32_t nd;
    uint32_t ny;
    uint32_t l_u_cnt;
    uint32_t n_x_cnt;
    uint32_t l_x_cnt;
    uint32_t u_size;
    float *u_min;
    float *u_max;
    uint32_t *u_cnt_idx;
    float *x_min;
    float *x_max;
    uint32_t *x_cnt_idx;
    uint32_t *y_idx;
    float *A;
    float *B;
    float *Kx;
    float *Ky;
    float *Ej;
    float *Fj;
    float *M;
    float *gam;
    float *Fj_1;
    float *Fj_2;
    float *Fx;
    float *xa;
    float *dx;
    float *e;
    float *auxm1;
    float *auxm2;
}dmpc_data_t;

typedef struct{
    uint32_t n_lambda;
    float *Kj_1;
    float *Hj;
    float *Kj;
    float *lambda;
    float *DU_1;
    float *DU_2;
    float *aux;
}dmpc_hild_data_t;

typedef struct{
    dmpc_data_t *prob_data;
    void *solver_data;
    int32_t (*solve)(dmpc_data_t *, void *);
}dmpc_inst_t;

#endif /* DMPC_DATA_H_ */
