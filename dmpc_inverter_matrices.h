/**
 * @file dmpc_inverter_matrices.h
 * @brief Header with data to run the DMPC algorithm for a two-level
 * inverter with LCL filter.
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

#ifndef DMPC_INVERTER_MATRICES_H_
#define DMPC_INVERTER_MATRICES_H_

/* Scaling factor */
#define DIM_CONFIG_SCALE			200.000000f

/* Number of model states and augmented states */
#define DIM_CONFIG_NXM			6
#define DIM_CONFIG_NXA			8

/* Prediction, control and constraint horizon */
#define DIM_CONFIG_NP			5
#define DIM_CONFIG_NC			5
#define DIM_CONFIG_NR			1
#define DIM_CONFIG_NLAMBDA		8

/* Number of inputs and outputs */
#define DIM_CONFIG_NU			2
#define DIM_CONFIG_NY			2

/* Size of control vector */
#define DIM_CONFIG_NC_x_NU		DIM_CONFIG_NC * DIM_CONFIG_NU

/* Input constraints */
#define DIM_CONFIG_NU_CTR		2
float DIM_CONFIG_U_MIN[2] = {-1.32680694,-1.32680694};
float DIM_CONFIG_U_MAX[2] = {1.32680694,1.32680694};
uint32_t DIM_CONFIG_U_LIM_IDX[2] = {0,1};

/* State constraints */
#define DIM_CONFIG_NXM_CTR		2
float DIM_CONFIG_XM_MIN[2] = {-0.075,-0.075};
float DIM_CONFIG_XM_MAX[2] = {0.075,0.075};
uint32_t DIM_CONFIG_XM_LIM_IDX[2] = {2,3};

/* A and B matrices for prediction */
float DIM_A[6][6] = {{ 5.11905185e-01, 3.22063444e-02, 4.75924325e-01, 2.99426206e-02,
  -8.11648780e-02,-5.10646129e-03},
 {-3.22063444e-02, 5.11905185e-01,-2.99426206e-02, 4.75924325e-01,
   5.10646129e-03,-8.11648780e-02},
 { 2.51985628e-01, 1.58535919e-02, 7.39747286e-01, 4.65409543e-02,
   4.31013897e-02, 2.71170959e-03},
 {-1.58535919e-02, 2.51985628e-01,-4.65409543e-02, 7.39747286e-01,
  -2.71170959e-03, 4.31013897e-02},
 { 7.30532601e+00, 4.59612155e-01,-7.32697764e+00,-4.60974360e-01,
   2.60613665e-01, 1.63964220e-02},
 {-4.59612155e-01, 7.30532601e+00, 4.60974360e-01,-7.32697764e+00,
  -1.63964220e-02, 2.60613665e-01}};
float DIM_B[6][4] = {{-1.00939177e-02,-4.70165160e-04, 9.13806134e-02, 2.57685402e-03},
 { 4.70165160e-04,-1.00939177e-02,-2.57685402e-03, 9.13806134e-02},
 {-5.32612596e-02,-1.59088687e-03, 1.00952411e-02, 4.70226939e-04},
 { 1.59088687e-03,-5.32612596e-02,-4.70226939e-04, 1.00952411e-02},
 { 2.55813379e-01, 1.03532214e-02, 4.82369470e-01, 1.95168774e-02},
 {-1.03532214e-02, 2.55813379e-01,-1.95168774e-02, 4.82369470e-01}};

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
float DIM_Ej[10][10] = {{ 8.63663382e-02, 2.10161705e-18, 6.44276499e-02,-9.42703517e-04,
   3.64856536e-02,-1.33330978e-03, 1.29331570e-02,-8.59721934e-04,
   1.82729724e-03,-1.91598092e-04},
 { 2.10161705e-18, 8.63663382e-02, 9.42703517e-04, 6.44276499e-02,
   1.33330978e-03, 3.64856536e-02, 8.59721934e-04, 1.29331570e-02,
   1.91598092e-04, 1.82729724e-03},
 { 6.44276499e-02, 9.42703517e-04, 5.33060742e-02, 1.77057870e-18,
   3.24071099e-02,-8.22571362e-04, 1.21619766e-02,-7.26228473e-04,
   1.77052614e-03,-1.78932547e-04},
 {-9.42703517e-04, 6.44276499e-02, 1.77057870e-18, 5.33060742e-02,
   8.22571362e-04, 3.24071099e-02, 7.26228473e-04, 1.21619766e-02,
   1.78932547e-04, 1.77052614e-03},
 { 3.64856536e-02, 1.33330978e-03, 3.24071099e-02, 8.22571362e-04,
   2.22921230e-02, 6.42620236e-19, 8.84619086e-03,-3.22968107e-04,
   1.34792945e-03,-1.07411487e-04},
 {-1.33330978e-03, 3.64856536e-02,-8.22571362e-04, 3.24071099e-02,
   6.42620236e-19, 2.22921230e-02, 3.22968107e-04, 8.84619086e-03,
   1.07411487e-04, 1.34792945e-03},
 { 1.29331570e-02, 8.59721934e-04, 1.21619766e-02, 7.26228473e-04,
   8.84619086e-03, 3.22968107e-04, 4.38513313e-03,-1.46255706e-20,
   6.20923128e-04,-2.70637346e-05},
 {-8.59721934e-04, 1.29331570e-02,-7.26228473e-04, 1.21619766e-02,
  -3.22968107e-04, 8.84619086e-03,-1.46255706e-20, 4.38513313e-03,
   2.70637346e-05, 6.20923128e-04},
 { 1.82729724e-03, 1.91598092e-04, 1.77052614e-03, 1.78932547e-04,
   1.34792945e-03, 1.07411487e-04, 6.20923128e-04, 2.70637346e-05,
   6.02108230e-04,-5.74096699e-21},
 {-1.91598092e-04, 1.82729724e-03,-1.78932547e-04, 1.77052614e-03,
  -1.07411487e-04, 1.34792945e-03,-2.70637346e-05, 6.20923128e-04,
  -5.74096699e-21, 6.02108230e-04}};

float DIM_Fj[10];

float DIM_M[8][10] = {{-1.        ,-0.        ,-0.        ,-0.        ,-0.        ,-0.        ,
  -0.        ,-0.        ,-0.        ,-0.        },
 {-0.        ,-1.        ,-0.        ,-0.        ,-0.        ,-0.        ,
  -0.        ,-0.        ,-0.        ,-0.        },
 { 1.        , 0.        , 0.        , 0.        , 0.        , 0.        ,
   0.        , 0.        , 0.        , 0.        },
 { 0.        , 1.        , 0.        , 0.        , 0.        , 0.        ,
   0.        , 0.        , 0.        , 0.        },
 { 0.05326126, 0.00159089,-0.        ,-0.        ,-0.        ,-0.        ,
  -0.        ,-0.        ,-0.        ,-0.        },
 {-0.00159089, 0.05326126,-0.        ,-0.        ,-0.        ,-0.        ,
  -0.        ,-0.        ,-0.        ,-0.        },
 {-0.05326126,-0.00159089, 0.        , 0.        , 0.        , 0.        ,
   0.        , 0.        , 0.        , 0.        },
 { 0.00159089,-0.05326126, 0.        , 0.        , 0.        , 0.        ,
   0.        , 0.        , 0.        , 0.        }};

float DIM_gam[8];

 /* Matrices for Hildreth's QP procedure */
float DIM_Fj_1[10][2] = {{ 5.58062802e-01,-7.60247478e-02},
 { 7.60247478e-02, 5.58062802e-01},
 { 3.78307407e-01,-4.86703722e-02},
 { 4.86703722e-02, 3.78307407e-01},
 { 2.04105805e-01,-2.28294570e-02},
 { 2.28294570e-02, 2.04105805e-01},
 { 7.13507090e-02,-6.00464101e-03},
 { 6.00464101e-03, 7.13507090e-02},
 { 1.00939177e-02,-4.70165160e-04},
 { 4.70165160e-04, 1.00939177e-02}};

float DIM_Fj_2[10][8] = {{-5.09636863e-01,-2.88229835e-02,-1.58964746e+00,-1.04072505e-02,
   1.85109579e-02,-6.94959973e-03,-5.58062802e-01, 7.60247478e-02},
 { 2.88229835e-02,-5.09636863e-01, 1.04072505e-02,-1.58964746e+00,
   6.94959973e-03, 1.85109579e-02,-7.60247478e-02,-5.58062802e-01},
 {-4.31003051e-01,-3.47686958e-02,-1.14993072e+00,-2.31602905e-02,
   6.82444847e-03,-5.27178262e-03,-3.78307407e-01, 4.86703722e-02},
 { 3.47686958e-02,-4.31003051e-01, 2.31602905e-02,-1.14993072e+00,
   5.27178262e-03, 6.82444847e-03,-4.86703722e-02,-3.78307407e-01},
 {-2.84886863e-01,-3.29049870e-02,-6.34974187e-01,-2.55046711e-02,
   2.11811280e-03,-3.04045187e-03,-2.04105805e-01, 2.28294570e-02},
 { 3.29049870e-02,-2.84886863e-01, 2.55046711e-02,-6.34974187e-01,
   3.04045187e-03, 2.11811280e-03,-2.28294570e-02,-2.04105805e-01},
 {-1.14670980e-01,-1.81713276e-02,-2.22524028e-01,-1.53153335e-02,
   7.92541403e-04,-1.02801512e-03,-7.13507090e-02, 6.00464101e-03},
 { 1.81713276e-02,-1.14670980e-01, 1.53153335e-02,-2.22524028e-01,
   1.02801512e-03, 7.92541403e-04,-6.00464101e-03,-7.13507090e-02},
 {-1.73519303e-02,-3.56979839e-03,-3.13322615e-02,-3.34445675e-03,
   1.42883590e-04,-1.32703803e-04,-1.00939177e-02, 4.70165160e-04},
 { 3.56979839e-03,-1.73519303e-02, 3.34445675e-03,-3.13322615e-02,
   1.32703803e-04, 1.42883590e-04,-4.70165160e-04,-1.00939177e-02}};

float DIM_Fx[2][6] = {{ 0.25198563, 0.01585359, 0.73974729, 0.04654095, 0.04310139, 0.00271171},
 {-0.01585359, 0.25198563,-0.04654095, 0.73974729,-0.00271171, 0.04310139}};

float DIM_Kj_1[8][10] = {{-1.78413556e+02, 6.73300491e-14, 3.21277926e+02,-1.00674173e+01,
  -1.50526494e+02, 1.39463180e+01,-5.99430574e+01, 6.89571888e+00,
  -4.28572468e+00, 4.26973138e-01},
 { 6.55595403e-14,-1.78413556e+02, 1.00674173e+01, 3.21277926e+02,
  -1.39463180e+01,-1.50526494e+02,-6.89571888e+00,-5.99430574e+01,
  -4.26973138e-01,-4.28572468e+00},
 { 1.78413556e+02,-6.73300491e-14,-3.21277926e+02, 1.00674173e+01,
   1.50526494e+02,-1.39463180e+01, 5.99430574e+01,-6.89571888e+00,
   4.28572468e+00,-4.26973138e-01},
 {-6.55595403e-14, 1.78413556e+02,-1.00674173e+01,-3.21277926e+02,
   1.39463180e+01, 1.50526494e+02, 6.89571888e+00, 5.99430574e+01,
   4.26973138e-01, 4.28572468e+00},
 { 9.50253074e+00, 2.83835784e-01,-1.71276831e+01, 2.50864907e-02,
   8.03941770e+00,-5.03327840e-01, 3.20361305e+00,-2.71912050e-01,
   2.28942360e-01,-1.59230240e-02},
 {-2.83835784e-01, 9.50253074e+00,-2.50864907e-02,-1.71276831e+01,
   5.03327840e-01, 8.03941770e+00, 2.71912050e-01, 3.20361305e+00,
   1.59230240e-02, 2.28942360e-01},
 {-9.50253074e+00,-2.83835784e-01, 1.71276831e+01,-2.50864907e-02,
  -8.03941770e+00, 5.03327840e-01,-3.20361305e+00, 2.71912050e-01,
  -2.28942360e-01, 1.59230240e-02},
 { 2.83835784e-01,-9.50253074e+00, 2.50864907e-02, 1.71276831e+01,
  -5.03327840e-01,-8.03941770e+00,-2.71912050e-01,-3.20361305e+00,
  -1.59230240e-02,-2.28942360e-01}};

float DIM_Hj[8][8] = {{-5.60495525e-03,-6.73300491e-14,-1.78413556e+02, 6.73300491e-14,
  -9.50253074e+00, 2.83835784e-01, 9.50253074e+00,-2.83835784e-01},
 {-6.55595403e-14,-5.60495525e-03, 6.55595403e-14,-1.78413556e+02,
  -2.83835784e-01,-9.50253074e+00, 2.83835784e-01, 9.50253074e+00},
 {-1.78413556e+02, 6.73300491e-14,-5.60495525e-03,-6.73300491e-14,
   9.50253074e+00,-2.83835784e-01,-9.50253074e+00, 2.83835784e-01},
 { 6.55595403e-14,-1.78413556e+02,-6.55595403e-14,-5.60495525e-03,
   2.83835784e-01, 9.50253074e+00,-2.83835784e-01,-9.50253074e+00},
 {-9.50253074e+00,-2.83835784e-01, 9.50253074e+00, 2.83835784e-01,
  -1.97406744e+00,-1.80094680e-16,-5.06568307e-01, 1.80094680e-16},
 { 2.83835784e-01,-9.50253074e+00,-2.83835784e-01, 9.50253074e+00,
  -1.74531766e-16,-1.97406744e+00, 1.74531766e-16,-5.06568307e-01},
 { 9.50253074e+00, 2.83835784e-01,-9.50253074e+00,-2.83835784e-01,
  -5.06568307e-01, 1.80094680e-16,-1.97406744e+00,-1.80094680e-16},
 {-2.83835784e-01, 9.50253074e+00, 2.83835784e-01,-9.50253074e+00,
   1.74531766e-16,-5.06568307e-01,-1.74531766e-16,-1.97406744e+00}};

float DIM_DU_1[2][10] = {{-1.78413556e+02, 6.73300491e-14, 3.21277926e+02,-1.00674173e+01,
  -1.50526494e+02, 1.39463180e+01,-5.99430574e+01, 6.89571888e+00,
  -4.28572468e+00, 4.26973138e-01},
 { 6.55595403e-14,-1.78413556e+02, 1.00674173e+01, 3.21277926e+02,
  -1.39463180e+01,-1.50526494e+02,-6.89571888e+00,-5.99430574e+01,
  -4.26973138e-01,-4.28572468e+00}};

float DIM_DU_2[2][8] = {{ 1.78413556e+02,-6.73300491e-14,-1.78413556e+02, 6.73300491e-14,
  -9.50253074e+00, 2.83835784e-01, 9.50253074e+00,-2.83835784e-01},
 {-6.55595403e-14, 1.78413556e+02, 6.55595403e-14,-1.78413556e+02,
  -2.83835784e-01,-9.50253074e+00, 2.83835784e-01, 9.50253074e+00}};


#endif /* DMPC_INVERTER_MATRICES_H_ */
