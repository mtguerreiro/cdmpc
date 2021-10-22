/**
 * @file dmpc_buck_matrices.h
 * @brief Header with data to run the DMPC algorithm for a buck converter.
 *
 * This file is generated automatically and should not be modified.
 *
 * The Hj matrix is already generated by flipping the sign and inverting its
 * diagonal elements, so that Hildreth's algorithm does not require any 
 * divisions.
 *
 *  Originally created on: 27.05.2021
 *      Author: mguerreiro
 */

#ifndef DMPC_BUCK_MATRICES_H_
#define DMPC_BUCK_MATRICES_H_

/* Prediction, control and restriction horizon */
#define DMPC_BUCK_CONFIG_NP			6
#define DMPC_BUCK_CONFIG_NC			6
#define DMPC_BUCK_CONFIG_NR			4
#define DMPC_BUCK_CONFIG_NLAMBDA	16

/* Constraints */
#define DMPC_BUCK_CONFIG_IL_MIN		-10
#define DMPC_BUCK_CONFIG_IL_MAX		10
#define DMPC_BUCK_CONFIG_U_MIN		0
#define DMPC_BUCK_CONFIG_U_MAX		1

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
float Ej[6][6] = {{74.5318609 ,46.44163231,30.84199385,18.18834553, 8.79660161, 2.78156124},
 {46.44163231,43.87740395,22.811542  ,13.64075874, 6.68887115, 2.14417575},
 {30.84199385,22.811542  ,25.66207009, 9.55191259, 4.77635059, 1.5608265 },
 {18.18834553,13.64075874, 9.55191259,16.00988509, 3.09917513, 1.04345323},
 { 8.79660161, 6.68887115, 4.77635059, 3.09917513,11.69606178, 0.60363587},
 { 2.78156124, 2.14417575, 1.5608265 , 1.04345323, 0.60363587,10.25239667}};

float Fj[6];

float M[16][6] = {{ -1.        , -0.        , -0.        , -0.        , -0.        ,
   -0.        },
 { -1.        , -1.        , -0.        , -0.        , -0.        ,
   -0.        },
 { -1.        , -1.        , -1.        , -0.        , -0.        ,
   -0.        },
 { -1.        , -1.        , -1.        , -1.        , -0.        ,
   -0.        },
 {  1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
    0.        },
 {  1.        ,  1.        ,  0.        ,  0.        ,  0.        ,
    0.        },
 {  1.        ,  1.        ,  1.        ,  0.        ,  0.        ,
    0.        },
 {  1.        ,  1.        ,  1.        ,  1.        ,  0.        ,
    0.        },
 { -6.66607032,  0.        ,  0.        ,  0.        ,  0.        ,
    0.        },
 {-12.99278635, -6.66607032,  0.        ,  0.        ,  0.        ,
    0.        },
 {-18.90553232,-12.99278635, -6.66607032,  0.        ,  0.        ,
    0.        },
 {-24.3396114 ,-18.90553232,-12.99278635, -6.66607032,  0.        ,
    0.        },
 {  6.66607032,  0.        ,  0.        ,  0.        ,  0.        ,
    0.        },
 { 12.99278635,  6.66607032,  0.        ,  0.        ,  0.        ,
    0.        },
 { 18.90553232, 12.99278635,  6.66607032,  0.        ,  0.        ,
    0.        },
 { 24.3396114 , 18.90553232, 12.99278635,  6.66607032,  0.        ,
    0.        }};

float gam[16];

 /* Matrices for Hildreth's QP procedure */
float Fj_1[6] = {-16.69227755,-11.15563083, -6.88768825, -3.78089167, -1.70391711,
  -0.50239096};

float Fj_2[6][3] = {{ 6.74369522,62.18772113,16.69227755},
 { 4.88955505,43.87904879,11.15563083},
 { 3.27570949,28.4520598 , 6.88768825},
 { 1.95049214,16.32768834, 3.78089167},
 { 0.95306163, 7.66223847, 1.70391711},
 { 0.30461254, 2.3446102 , 0.50239096}};

float Fx[4][3] = {{ 0.9802431 ,-0.4133303 , 0.        },
 { 1.92746113,-1.21246104, 0.        },
 { 2.8295578 ,-2.36566796, 0.        },
 { 3.67572411,-3.83771326, 0.        }};

float Kj_1[16][6] = {{-0.04658414, 0.03488136, 0.02037572, 0.01020682, 0.00401225, 0.00096657},
 {-0.01170278,-0.0386159 , 0.0386021 , 0.02117757, 0.0093725 , 0.00266715},
 { 0.00867293,-0.02038952,-0.04626708, 0.03183897, 0.01534702, 0.00481086},
 { 0.01887975,-0.00941876,-0.03560568,-0.05917504, 0.02104345, 0.00705185},
 { 0.04658414,-0.03488136,-0.02037572,-0.01020682,-0.00401225,-0.00096657},
 { 0.01170278, 0.0386159 ,-0.0386021 ,-0.02117757,-0.0093725 ,-0.00266715},
 {-0.00867293, 0.02038952, 0.04626708,-0.03183897,-0.01534702,-0.00481086},
 {-0.01887975, 0.00941876, 0.03560568, 0.05917504,-0.02104345,-0.00705185},
 {-0.31053316, 0.23252158, 0.13582596, 0.06803937, 0.02674592, 0.00644322},
 {-0.37273622,-0.03673187, 0.38623569, 0.20574683, 0.0878621 , 0.02389465},
 {-0.29166601,-0.1739852 , 0.05628139, 0.40657562, 0.1853249 , 0.05465891},
 {-0.14161257,-0.23056279,-0.19110094,-0.01234662, 0.31459349, 0.09846769},
 { 0.31053316,-0.23252158,-0.13582596,-0.06803937,-0.02674592,-0.00644322},
 { 0.37273622, 0.03673187,-0.38623569,-0.20574683,-0.0878621 ,-0.02389465},
 { 0.29166601, 0.1739852 ,-0.05628139,-0.40657562,-0.1853249 ,-0.05465891},
 { 0.14161257, 0.23056279, 0.19110094, 0.01234662,-0.31459349,-0.09846769}};

float Hj[16][16] = {{-2.14665328e+01, 1.17027843e-02,-8.67293276e-03,-1.88797512e-02,
  -4.65841415e-02,-1.17027843e-02, 8.67293276e-03, 1.88797512e-02,
   3.10533163e-01, 3.72736218e-01, 2.91666008e-01, 1.41612568e-01,
  -3.10533163e-01,-3.72736218e-01,-2.91666008e-01,-1.41612568e-01},
 { 1.17027843e-02,-1.98733331e+01, 1.17165850e-02,-9.46098816e-03,
  -1.17027843e-02,-5.03186856e-02,-1.17165850e-02, 9.46098816e-03,
   7.80115829e-02, 4.09468090e-01, 4.65651205e-01, 3.72175355e-01,
  -7.80115829e-02,-4.09468090e-01,-4.65651205e-01,-3.72175355e-01},
 {-8.67293276e-03, 1.17165850e-02,-1.72462379e+01, 2.61446914e-02,
   8.67293276e-03,-1.17165850e-02,-5.79836602e-02,-2.61446914e-02,
  -5.78143797e-02, 2.32323971e-02, 4.09369815e-01, 5.63276292e-01,
   5.78143797e-02,-2.32323971e-02,-4.09369815e-01,-5.63276292e-01},
 {-1.88797512e-02,-9.46098816e-03, 2.61446914e-02,-1.17206177e+01,
   1.88797512e-02, 9.46098816e-03,-2.61446914e-02,-8.53197351e-02,
  -1.25853749e-01,-1.82514437e-01, 2.79419336e-03, 5.75622912e-01,
   1.25853749e-01, 1.82514437e-01,-2.79419336e-03,-5.75622912e-01},
 {-4.65841415e-02,-1.17027843e-02, 8.67293276e-03, 1.88797512e-02,
  -2.14665328e+01, 1.17027843e-02,-8.67293276e-03,-1.88797512e-02,
  -3.10533163e-01,-3.72736218e-01,-2.91666008e-01,-1.41612568e-01,
   3.10533163e-01, 3.72736218e-01, 2.91666008e-01, 1.41612568e-01},
 {-1.17027843e-02,-5.03186856e-02,-1.17165850e-02, 9.46098816e-03,
   1.17027843e-02,-1.98733331e+01, 1.17165850e-02,-9.46098816e-03,
  -7.80115829e-02,-4.09468090e-01,-4.65651205e-01,-3.72175355e-01,
   7.80115829e-02, 4.09468090e-01, 4.65651205e-01, 3.72175355e-01},
 { 8.67293276e-03,-1.17165850e-02,-5.79836602e-02,-2.61446914e-02,
  -8.67293276e-03, 1.17165850e-02,-1.72462379e+01, 2.61446914e-02,
   5.78143797e-02,-2.32323971e-02,-4.09369815e-01,-5.63276292e-01,
  -5.78143797e-02, 2.32323971e-02, 4.09369815e-01, 5.63276292e-01},
 { 1.88797512e-02, 9.46098816e-03,-2.61446914e-02,-8.53197351e-02,
  -1.88797512e-02,-9.46098816e-03, 2.61446914e-02,-1.17206177e+01,
   1.25853749e-01, 1.82514437e-01,-2.79419336e-03,-5.75622912e-01,
  -1.25853749e-01,-1.82514437e-01, 2.79419336e-03, 5.75622912e-01},
 { 3.10533163e-01, 7.80115829e-02,-5.78143797e-02,-1.25853749e-01,
  -3.10533163e-01,-7.80115829e-02, 5.78143797e-02, 1.25853749e-01,
  -4.83083408e-01, 2.48468584e+00, 1.94426612e+00, 9.43999337e-01,
  -2.07003590e+00,-2.48468584e+00,-1.94426612e+00,-9.43999337e-01},
 { 3.72736218e-01, 4.09468090e-01, 2.32323971e-02,-1.82514437e-01,
  -3.72736218e-01,-4.09468090e-01,-2.32323971e-02, 1.82514437e-01,
   2.48468584e+00,-1.96550952e-01, 4.94935169e+00, 3.37688959e+00,
  -2.48468584e+00,-5.08773929e+00,-4.94935169e+00,-3.37688959e+00},
 { 2.91666008e-01, 4.65651205e-01, 4.09369815e-01, 2.79419336e-03,
  -2.91666008e-01,-4.65651205e-01,-4.09369815e-01,-2.79419336e-03,
   1.94426612e+00, 4.94935169e+00,-1.35144670e-01, 6.94680630e+00,
  -1.94426612e+00,-4.94935169e+00,-7.39947793e+00,-6.94680630e+00},
 { 1.41612568e-01, 3.72175355e-01, 5.63276292e-01, 5.75622912e-01,
  -1.41612568e-01,-3.72175355e-01,-5.63276292e-01,-5.75622912e-01,
   9.43999337e-01, 3.37688959e+00, 6.94680630e+00,-9.64232362e-02,
  -9.43999337e-01,-3.37688959e+00,-6.94680630e+00,-1.03709442e+01},
 {-3.10533163e-01,-7.80115829e-02, 5.78143797e-02, 1.25853749e-01,
   3.10533163e-01, 7.80115829e-02,-5.78143797e-02,-1.25853749e-01,
  -2.07003590e+00,-2.48468584e+00,-1.94426612e+00,-9.43999337e-01,
  -4.83083408e-01, 2.48468584e+00, 1.94426612e+00, 9.43999337e-01},
 {-3.72736218e-01,-4.09468090e-01,-2.32323971e-02, 1.82514437e-01,
   3.72736218e-01, 4.09468090e-01, 2.32323971e-02,-1.82514437e-01,
  -2.48468584e+00,-5.08773929e+00,-4.94935169e+00,-3.37688959e+00,
   2.48468584e+00,-1.96550952e-01, 4.94935169e+00, 3.37688959e+00},
 {-2.91666008e-01,-4.65651205e-01,-4.09369815e-01,-2.79419336e-03,
   2.91666008e-01, 4.65651205e-01, 4.09369815e-01, 2.79419336e-03,
  -1.94426612e+00,-4.94935169e+00,-7.39947793e+00,-6.94680630e+00,
   1.94426612e+00, 4.94935169e+00,-1.35144670e-01, 6.94680630e+00},
 {-1.41612568e-01,-3.72175355e-01,-5.63276292e-01,-5.75622912e-01,
   1.41612568e-01, 3.72175355e-01, 5.63276292e-01, 5.75622912e-01,
  -9.43999337e-01,-3.37688959e+00,-6.94680630e+00,-1.03709442e+01,
   9.43999337e-01, 3.37688959e+00, 6.94680630e+00,-9.64232362e-02}};

float DU_1[6] = {-0.04658414, 0.03488136, 0.02037572, 0.01020682, 0.00401225, 0.00096657};

float DU_2[16] = { 0.04658414, 0.01170278,-0.00867293,-0.01887975,-0.04658414,-0.01170278,
  0.00867293, 0.01887975, 0.31053316, 0.37273622, 0.29166601, 0.14161257,
 -0.31053316,-0.37273622,-0.29166601,-0.14161257};


#endif /* DMPC_BUCK_MATRICES_H_ */
