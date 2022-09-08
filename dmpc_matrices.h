/**
 * @file dmpc_matrices.h
 * @brief Header with data to run the DMPC algorithm.
 *
 * This file is generated automatically and should not be modified.
 *
 * The Hj matrix is already generated by flipping the sign and inverting its
 * diagonal elements, so that Hildreth's algorithm does not require any 
 * divisions.
 *
 *  Originally created on: 22.04.2022
 *      Author: mguerreiro
 */

#ifndef DMPC_MATRICES_H_
#define DMPC_MATRICES_H_

/* Scaling factor */
#define DMPC_CONFIG_SCALE			1.000000f

/* Number of model states and augmented states */
#define DMPC_CONFIG_NXM			2
#define DMPC_CONFIG_NXA			3

/* Prediction, control and constraint horizon */
#define DMPC_CONFIG_NP			5
#define DMPC_CONFIG_NC			5
#define DMPC_CONFIG_NR			5
#define DMPC_CONFIG_NLAMBDA		20

/* Number of inputs and outputs */
#define DMPC_CONFIG_NU			1
#define DMPC_CONFIG_NY			1

/* Number of external disturbances */
#define DMPC_CONFIG_ND			1

/* Size of control vector */
#define DMPC_CONFIG_NC_x_NU		DMPC_CONFIG_NC * DMPC_CONFIG_NU

/* Input constraints */
#define DMPC_CONFIG_NU_CTR		1
float DMPC_CONFIG_U_MIN[1] = {0.};
float DMPC_CONFIG_U_MAX[1] = {16.};
uint32_t DMPC_CONFIG_U_LIM_IDX[1] = {0};

/* State constraints */
#define DMPC_CONFIG_NXM_CTR		1
float DMPC_CONFIG_XM_MIN[1] = {-12.0};
float DMPC_CONFIG_XM_MAX[1] = {12.0};
uint32_t DMPC_CONFIG_XM_LIM_IDX[1] = {0};

/* Indexes of outputs */
uint32_t DMPC_CONFIG_Y_IDX[1] = {1};

/* A and B matrices for prediction */
float DMPC_M_A[2][2] = {{ 0.98622461,-0.41775318},
 { 0.03468545, 0.96742572}};
float DMPC_M_B[2][2] = {{ 0.41775318, 0.00750909},
 { 0.03257428,-0.03517406}};

/*
 * Matrices for QP solvers 
 *
 * The matrices were generated considering the following problem:
 *
 * min (1/2) * DU' * Ej * DU + DU' * Fj
 * DU
 *
 * s.t. M * DU <= gam
 *
 * The (1/2) term in from of DU' * Ej * DU needs to be considered in the QP
 * solver selected, or the solution will appear to be inconsistent.
 * Note that the Fj and gam matrices are usually updated online, while Ej
 * and M are static.
 */
float DMPC_M_Ej[5][5] = {{0.15995373,0.10042023,0.05973293,0.02914143,0.00929498},
 {0.10042023,0.07853072,0.04155739,0.02066878,0.00671959},
 {0.05973293,0.04155739,0.03597722,0.01331689,0.00445943},
 {0.02914143,0.02066878,0.01331689,0.0172355 ,0.0025596 },
 {0.00929498,0.00671959,0.00445943,0.0025596 ,0.01106108}};

float DMPC_M_Fj[5];

float DMPC_M_M[20][5] = {{-1.        ,-0.        ,-0.        ,-0.        ,-0.        },
 {-1.        ,-1.        ,-0.        ,-0.        ,-0.        },
 {-1.        ,-1.        ,-1.        ,-0.        ,-0.        },
 {-1.        ,-1.        ,-1.        ,-1.        ,-0.        },
 {-1.        ,-1.        ,-1.        ,-1.        ,-1.        },
 { 1.        , 0.        , 0.        , 0.        , 0.        },
 { 1.        , 1.        , 0.        , 0.        , 0.        },
 { 1.        , 1.        , 1.        , 0.        , 0.        },
 { 1.        , 1.        , 1.        , 1.        , 0.        },
 { 1.        , 1.        , 1.        , 1.        , 1.        },
 {-0.41775318,-0.        ,-0.        ,-0.        ,-0.        },
 {-0.39839046,-0.41775318,-0.        ,-0.        ,-0.        },
 {-0.37368451,-0.39839046,-0.41775318,-0.        ,-0.        },
 {-0.34417225,-0.37368451,-0.39839046,-0.41775318,-0.        },
 {-0.31044552,-0.34417225,-0.37368451,-0.39839046,-0.41775318},
 { 0.41775318, 0.        , 0.        , 0.        , 0.        },
 { 0.39839046, 0.41775318, 0.        , 0.        , 0.        },
 { 0.37368451, 0.39839046, 0.41775318, 0.        , 0.        },
 { 0.34417225, 0.37368451, 0.39839046, 0.41775318, 0.        },
 { 0.31044552, 0.34417225, 0.37368451, 0.39839046, 0.41775318}};

float DMPC_M_gam[20];

 /* Matrices for Hildreth's QP procedure */
float DMPC_M_Fj_1[5] = {-0.73968431,-0.45433713,-0.24805214,-0.11115172,-0.03257428};

float DMPC_M_Fj_2[5][3] = {{0.24198063,2.50376298,0.73968431},
 {0.16351302,1.64201081,0.45433713},
 {0.09832731,0.95149082,0.24805214},
 {0.04855841,0.45026209,0.11115172},
 {0.01569253,0.13877672,0.03257428}};

float DMPC_M_Fx[5][2] = {{ 0.98622461,-0.41775318},
 { 0.95814903,-0.81614363},
 { 0.91664185,-1.18982814},
 { 0.86274502,-1.53400039},
 { 0.79765288,-1.84444591}};

float DMPC_M_Kj_1[20][5] = {{-35.44815314, 35.6037436 , 15.83090353,  4.91243182,  0.63978617},
 {  0.15559046,-35.71325096, 34.96119   , 14.95632753,  4.00889132},
 { 15.98649398,-16.58296449,-46.44354615, 27.40079469,  9.02381895},
 { 20.8989258 , -6.53906877,-33.99907898,-61.37101301, 14.31928488},
 { 21.53871197, -3.16996362,-28.98415136,-56.07554707,-81.91936914},
 { 35.44815314,-35.6037436 ,-15.83090353, -4.91243182, -0.63978617},
 { -0.15559046, 35.71325096,-34.96119   ,-14.95632753, -4.00889132},
 {-15.98649398, 16.58296449, 46.44354615,-27.40079469, -9.02381895},
 {-20.8989258 ,  6.53906877, 33.99907898, 61.37101301,-14.31928488},
 {-21.53871197,  3.16996362, 28.98415136, 56.07554707, 81.91936914},
 {-14.80857862, 14.87357703,  6.61341026,  2.052184  ,  0.26727271},
 {  0.75137112,-15.60870944, 14.29861885,  6.15293531,  1.66233909},
 {  7.55117619, -7.11570458,-20.47000022, 11.03580758,  3.67629947},
 {  9.4633618 , -2.57904277,-14.63490491,-26.68297686,  5.68927608},
 {  9.38920235, -0.93476947,-11.8681833 ,-23.52146606,-34.86216844},
 { 14.80857862,-14.87357703, -6.61341026, -2.052184  , -0.26727271},
 { -0.75137112, 15.60870944,-14.29861885, -6.15293531, -1.66233909},
 { -7.55117619,  7.11570458, 20.47000022,-11.03580758, -3.67629947},
 { -9.4633618 ,  2.57904277, 14.63490491, 26.68297686, -5.68927608},
 { -9.38920235,  0.93476947, 11.8681833 , 23.52146606, 34.86216844}};

float DMPC_M_Hj[20][20] = {{-2.82102144e-02,-1.55590457e-01,-1.59864940e+01,-2.08989258e+01,
  -2.15387120e+01,-3.54481531e+01, 1.55590457e-01, 1.59864940e+01,
   2.08989258e+01, 2.15387120e+01, 1.48085786e+01,-7.51371117e-01,
  -7.55117619e+00,-9.46336180e+00,-9.38920235e+00,-1.48085786e+01,
   7.51371117e-01, 7.55117619e+00, 9.46336180e+00, 9.38920235e+00},
 {-1.55590457e-01,-2.81233351e-02, 5.96470505e-01,-1.43598570e+01,
  -1.83687483e+01, 1.55590457e-01,-3.55576605e+01,-5.96470505e-01,
   1.43598570e+01, 1.83687483e+01,-6.49984077e-02, 1.48573383e+01,
  -4.35471611e-01,-6.88431903e+00,-8.45443288e+00, 6.49984077e-02,
  -1.48573383e+01, 4.35471611e-01, 6.88431903e+00, 8.45443288e+00},
 {-1.59864940e+01, 5.96470505e-01,-2.12584959e-02, 1.96392220e+01,
   1.06154030e+01, 1.59864940e+01,-5.96470505e-01,-4.70400167e+01,
  -1.96392220e+01,-1.06154030e+01,-6.67840866e+00, 5.58719476e-01,
   2.00345286e+01, 7.75058588e+00, 3.41375041e+00, 6.67840866e+00,
  -5.58719476e-01,-2.00345286e+01,-7.75058588e+00,-3.41375041e+00},
 {-2.08989258e+01,-1.43598570e+01, 1.96392220e+01,-1.23441192e-02,
   6.66909501e+01, 2.08989258e+01, 1.43598570e+01,-1.96392220e+01,
  -8.10102350e+01,-6.66909501e+01,-8.73059266e+00,-5.59421583e+00,
   8.99872104e+00, 3.44335627e+01, 2.69352165e+01, 8.73059266e+00,
   5.59421583e+00,-8.99872104e+00,-3.44335627e+01,-2.69352165e+01},
 {-2.15387120e+01,-1.83687483e+01, 1.06154030e+01, 6.66909501e+01,
  -6.72900782e-03, 2.15387120e+01, 1.83687483e+01,-1.06154030e+01,
  -6.66909501e+01,-1.48610319e+02,-8.99786537e+00,-7.25655492e+00,
   5.32242156e+00, 2.87442867e+01, 6.17973849e+01, 8.99786537e+00,
   7.25655492e+00,-5.32242156e+00,-2.87442867e+01,-6.17973849e+01},
 {-3.54481531e+01, 1.55590457e-01, 1.59864940e+01, 2.08989258e+01,
   2.15387120e+01,-2.82102144e-02,-1.55590457e-01,-1.59864940e+01,
  -2.08989258e+01,-2.15387120e+01,-1.48085786e+01, 7.51371117e-01,
   7.55117619e+00, 9.46336180e+00, 9.38920235e+00, 1.48085786e+01,
  -7.51371117e-01,-7.55117619e+00,-9.46336180e+00,-9.38920235e+00},
 { 1.55590457e-01,-3.55576605e+01,-5.96470505e-01, 1.43598570e+01,
   1.83687483e+01,-1.55590457e-01,-2.81233351e-02, 5.96470505e-01,
  -1.43598570e+01,-1.83687483e+01, 6.49984077e-02,-1.48573383e+01,
   4.35471611e-01, 6.88431903e+00, 8.45443288e+00,-6.49984077e-02,
   1.48573383e+01,-4.35471611e-01,-6.88431903e+00,-8.45443288e+00},
 { 1.59864940e+01,-5.96470505e-01,-4.70400167e+01,-1.96392220e+01,
  -1.06154030e+01,-1.59864940e+01, 5.96470505e-01,-2.12584959e-02,
   1.96392220e+01, 1.06154030e+01, 6.67840866e+00,-5.58719476e-01,
  -2.00345286e+01,-7.75058588e+00,-3.41375041e+00,-6.67840866e+00,
   5.58719476e-01, 2.00345286e+01, 7.75058588e+00, 3.41375041e+00},
 { 2.08989258e+01, 1.43598570e+01,-1.96392220e+01,-8.10102350e+01,
  -6.66909501e+01,-2.08989258e+01,-1.43598570e+01, 1.96392220e+01,
  -1.23441192e-02, 6.66909501e+01, 8.73059266e+00, 5.59421583e+00,
  -8.99872104e+00,-3.44335627e+01,-2.69352165e+01,-8.73059266e+00,
  -5.59421583e+00, 8.99872104e+00, 3.44335627e+01, 2.69352165e+01},
 { 2.15387120e+01, 1.83687483e+01,-1.06154030e+01,-6.66909501e+01,
  -1.48610319e+02,-2.15387120e+01,-1.83687483e+01, 1.06154030e+01,
   6.66909501e+01,-6.72900782e-03, 8.99786537e+00, 7.25655492e+00,
  -5.32242156e+00,-2.87442867e+01,-6.17973849e+01,-8.99786537e+00,
  -7.25655492e+00, 5.32242156e+00, 2.87442867e+01, 6.17973849e+01},
 { 1.48085786e+01,-6.49984077e-02,-6.67840866e+00,-8.73059266e+00,
  -8.99786537e+00,-1.48085786e+01, 6.49984077e-02, 6.67840866e+00,
   8.73059266e+00, 8.99786537e+00,-1.61646707e-01,-3.13887672e-01,
  -3.15452785e+00,-3.95334947e+00,-3.92236912e+00,-6.18633078e+00,
   3.13887672e-01, 3.15452785e+00, 3.95334947e+00, 3.92236912e+00},
 {-7.51371117e-01, 1.48573383e+01, 5.58719476e-01,-5.59421583e+00,
  -7.25655492e+00, 7.51371117e-01,-1.48573383e+01,-5.58719476e-01,
   5.59421583e+00, 7.25655492e+00,-3.13887672e-01,-1.60739430e-01,
  -3.57083295e-02,-2.69270972e+00,-3.35006569e+00, 3.13887672e-01,
  -6.22124889e+00, 3.57083295e-02, 2.69270972e+00, 3.35006569e+00},
 {-7.55117619e+00,-4.35471611e-01, 2.00345286e+01, 8.99872104e+00,
   5.32242156e+00, 7.55117619e+00, 4.35471611e-01,-2.00345286e+01,
  -8.99872104e+00,-5.32242156e+00,-3.15452785e+00,-3.57083295e-02,
  -1.16761337e-01, 3.60493234e+00, 1.82177504e+00, 3.15452785e+00,
   3.57083295e-02,-8.56447886e+00,-3.60493234e+00,-1.82177504e+00},
 {-9.46336180e+00,-6.88431903e+00, 7.75058588e+00, 3.44335627e+01,
   2.87442867e+01, 9.46336180e+00, 6.88431903e+00,-7.75058588e+00,
  -3.44335627e+01,-2.87442867e+01,-3.95334947e+00,-2.69270972e+00,
   3.60493234e+00,-6.81012112e-02, 1.16721441e+01, 3.95334947e+00,
   2.69270972e+00,-3.60493234e+00,-1.46840267e+01,-1.16721441e+01},
 {-9.38920235e+00,-8.45443288e+00, 3.41375041e+00, 2.69352165e+01,
   6.17973849e+01, 9.38920235e+00, 8.45443288e+00,-3.41375041e+00,
  -2.69352165e+01,-6.17973849e+01,-3.92236912e+00,-3.35006569e+00,
   1.82177504e+00, 1.16721441e+01,-3.87952501e-02, 3.92236912e+00,
   3.35006569e+00,-1.82177504e+00,-1.16721441e+01,-2.57763514e+01},
 {-1.48085786e+01, 6.49984077e-02, 6.67840866e+00, 8.73059266e+00,
   8.99786537e+00, 1.48085786e+01,-6.49984077e-02,-6.67840866e+00,
  -8.73059266e+00,-8.99786537e+00,-6.18633078e+00, 3.13887672e-01,
   3.15452785e+00, 3.95334947e+00, 3.92236912e+00,-1.61646707e-01,
  -3.13887672e-01,-3.15452785e+00,-3.95334947e+00,-3.92236912e+00},
 { 7.51371117e-01,-1.48573383e+01,-5.58719476e-01, 5.59421583e+00,
   7.25655492e+00,-7.51371117e-01, 1.48573383e+01, 5.58719476e-01,
  -5.59421583e+00,-7.25655492e+00, 3.13887672e-01,-6.22124889e+00,
   3.57083295e-02, 2.69270972e+00, 3.35006569e+00,-3.13887672e-01,
  -1.60739430e-01,-3.57083295e-02,-2.69270972e+00,-3.35006569e+00},
 { 7.55117619e+00, 4.35471611e-01,-2.00345286e+01,-8.99872104e+00,
  -5.32242156e+00,-7.55117619e+00,-4.35471611e-01, 2.00345286e+01,
   8.99872104e+00, 5.32242156e+00, 3.15452785e+00, 3.57083295e-02,
  -8.56447886e+00,-3.60493234e+00,-1.82177504e+00,-3.15452785e+00,
  -3.57083295e-02,-1.16761337e-01, 3.60493234e+00, 1.82177504e+00},
 { 9.46336180e+00, 6.88431903e+00,-7.75058588e+00,-3.44335627e+01,
  -2.87442867e+01,-9.46336180e+00,-6.88431903e+00, 7.75058588e+00,
   3.44335627e+01, 2.87442867e+01, 3.95334947e+00, 2.69270972e+00,
  -3.60493234e+00,-1.46840267e+01,-1.16721441e+01,-3.95334947e+00,
  -2.69270972e+00, 3.60493234e+00,-6.81012112e-02, 1.16721441e+01},
 { 9.38920235e+00, 8.45443288e+00,-3.41375041e+00,-2.69352165e+01,
  -6.17973849e+01,-9.38920235e+00,-8.45443288e+00, 3.41375041e+00,
   2.69352165e+01, 6.17973849e+01, 3.92236912e+00, 3.35006569e+00,
  -1.82177504e+00,-1.16721441e+01,-2.57763514e+01,-3.92236912e+00,
  -3.35006569e+00, 1.82177504e+00, 1.16721441e+01,-3.87952501e-02}};

float DMPC_M_DU_1[5] = {-35.44815314, 35.6037436 , 15.83090353,  4.91243182,  0.63978617};

float DMPC_M_DU_2[20] = { 35.44815314, -0.15559046,-15.98649398,-20.8989258 ,-21.53871197,
 -35.44815314,  0.15559046, 15.98649398, 20.8989258 , 21.53871197,
  14.80857862, -0.75137112, -7.55117619, -9.4633618 , -9.38920235,
 -14.80857862,  0.75137112,  7.55117619,  9.4633618 ,  9.38920235};


#endif /* DMPC_MATRICES_H_ */
