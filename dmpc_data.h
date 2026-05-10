
#ifndef DMPC_DATA_H_
#define DMPC_DATA_H_

#include "stdint.h"

/* Number of model states and augmented states */
#define DMPC_CONFIG_NXM                6
#define DMPC_CONFIG_NXA                8

/* Length of prediction, control and constraint horizons */
#define DMPC_CONFIG_L_PRED             4
#define DMPC_CONFIG_L_CTL              4
#define DMPC_CONFIG_L_U_CNT            4
#define DMPC_CONFIG_L_X_CNT            4
#define DMPC_CONFIG_NLAMBDA            32

/* Number of inputs, outputs, and disturbances */
#define DMPC_CONFIG_NU                 2
#define DMPC_CONFIG_NY                 2
#define DMPC_CONFIG_ND                 0

/* Size of control vector */
#define DMPC_CONFIG_U_SIZE             (DMPC_CONFIG_L_CTL * DMPC_CONFIG_NU)

/* Input constraints */
#define DMPC_CONFIG_NU_CNT             2

/* State constraints */
#define DMPC_CONFIG_NXM_CNT            2

typedef struct{
    float u_min[DMPC_CONFIG_NU_CNT];
    float u_max[DMPC_CONFIG_NU_CNT];
    uint32_t u_cnt_idx[DMPC_CONFIG_NU_CNT];
    float x_min[DMPC_CONFIG_NXM_CNT];
    float x_max[DMPC_CONFIG_NXM_CNT];
    uint32_t x_cnt_idx[DMPC_CONFIG_NXM_CNT];
    uint32_t y_idx[DMPC_CONFIG_NXM_CNT];
    float A[DMPC_CONFIG_NXM][DMPC_CONFIG_NXM];
    float B[DMPC_CONFIG_NXM][DMPC_CONFIG_NU+DMPC_CONFIG_ND];
    float Kx[DMPC_CONFIG_NY][DMPC_CONFIG_NXM];
    float Ky[DMPC_CONFIG_NY][DMPC_CONFIG_NY];
    float Ej[DMPC_CONFIG_U_SIZE][DMPC_CONFIG_U_SIZE];
    float Fj[DMPC_CONFIG_U_SIZE];
    float M[2*((DMPC_CONFIG_NU_CNT*DMPC_CONFIG_L_U_CNT)+(DMPC_CONFIG_NXM_CNT*DMPC_CONFIG_L_X_CNT))][DMPC_CONFIG_U_SIZE];
    float gam[2*((DMPC_CONFIG_NU_CNT*DMPC_CONFIG_L_U_CNT)+(DMPC_CONFIG_NXM_CNT*DMPC_CONFIG_L_X_CNT))];
    float Fj_1[DMPC_CONFIG_U_SIZE][DMPC_CONFIG_NY];
    float Fj_2[DMPC_CONFIG_U_SIZE][DMPC_CONFIG_NXA];
    float Fx[DMPC_CONFIG_L_X_CNT * DMPC_CONFIG_NXM_CNT][DMPC_CONFIG_NXM];
    float Kj_1[DMPC_CONFIG_NLAMBDA][DMPC_CONFIG_U_SIZE];
    float Hj[DMPC_CONFIG_NLAMBDA][DMPC_CONFIG_NLAMBDA];
    float DU_1[DMPC_CONFIG_NU][DMPC_CONFIG_U_SIZE];
    float DU_2[DMPC_CONFIG_NU][DMPC_CONFIG_NLAMBDA];
}dmpc_data_t;

extern dmpc_data_t dmpc_data;

#endif /* DMPC_DATA_H_ */
