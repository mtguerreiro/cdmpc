#include "auxil_st.h"
#include "stdint.h"

void update_xzy_st(
    float x[8], float x_prev[8], float delta_x[8],
    float z[8], float z_prev[8],
    float y[8], float delta_y[8],
    float xz_tilde[16],
    float rho_vec[8], float rho_inv_vec[8],
    float l[8], float u[8]
    ){

    int32_t i;
    float alpha = (float)1.6;

    /* Updates x*/
    for(i = 0; i < 8; i++){
        x[i] = alpha * xz_tilde[i] + ( ((float)1.0f) - alpha) * x_prev[i];
        delta_x[i] = x[i] - x_prev[i];
    }

    /* Updates z*/
    for(i = 0; i < 8; i++){
        z[i] = alpha * xz_tilde[i + 8] + ( ((float)1.0f) - alpha) * z_prev[i] + rho_inv_vec[i] * y[i];
        if( z[i] > u[i] ) z[i] = u[i];
        else if( z[i] < l[i] ) z[i] = l[i];
    }

    /* Updates y*/
    for(i = 0; i < 8; i++){
        delta_y[i] = rho_vec[i] * (alpha * xz_tilde[i + 8] + ( ((float)1.0f) - alpha) * z_prev[i] - z[i]);
        y[i] += delta_y[i];
    }
}
