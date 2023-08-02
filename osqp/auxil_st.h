#ifndef AUXIL_ST_H
#define AUXIL_ST_H

void update_xzy_st(
    float x[8], float x_prev[8], float delta_x[8],
    float z[8], float z_prev[8],
    float y[8], float delta_y[8],
    float xz_tilde[16],
    float rho_vec[8], float rho_inv_vec[8],
    float l[8], float u[8]
    );

#endif /* AUXIL_ST_H */
