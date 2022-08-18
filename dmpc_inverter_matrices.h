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

/* Number of model states and augmented states */
#define DIM_CONFIG_NXM			6
#define DIM_CONFIG_NXA			8

/* Prediction, control and constraint horizon */
#define DIM_CONFIG_NP			5
#define DIM_CONFIG_NC			5
#define DIM_CONFIG_NR			1
#define DIM_CONFIG_NLAMBDA		12

/* Number of inputs and outputs */
#define DIM_CONFIG_NU			4
#define DIM_CONFIG_NY			2

/* Size of control vector */
#define DIM_CONFIG_NC_x_NU		DIM_CONFIG_NC * DIM_CONFIG_NU

/* Input constraints */
#define DIM_CONFIG_NU_CTR		4
float DIM_CONFIG_U_MIN[4] = {-375.27767497,-375.27767497, 162.5       ,   0.        };
float DIM_CONFIG_U_MAX[4] = {375.27767497,375.27767497,162.5       ,  0.        };
uint32_t DIM_CONFIG_U_LIM_IDX[4] = {0,1,2,3};

/* State constraints */
#define DIM_CONFIG_NXM_CTR		2
float DIM_CONFIG_XM_MIN[2] = {-15,-15};
float DIM_CONFIG_XM_MAX[2] = {15,15};
uint32_t DIM_CONFIG_XM_LIM_IDX[2] = {2,3};

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
float DIM_Ej[20][20] = {{ 8.93665937e-02, 1.21213170e-17,-7.41456087e-02, 3.04220386e-03,
   6.66138672e-02,-9.72205252e-04,-5.28342800e-02, 5.80545157e-03,
   3.75783493e-02,-1.37567420e-03,-4.52790031e-02, 5.74317201e-03,
   1.32219792e-02,-8.81535933e-04,-3.57875684e-02, 4.29163515e-03,
   1.84214519e-03,-1.93781969e-04,-1.69429504e-02, 2.09776080e-03},
 { 1.21213170e-17, 8.93665937e-02,-3.04220386e-03,-7.41456087e-02,
   9.72205252e-04, 6.66138672e-02,-5.80545157e-03,-5.28342800e-02,
   1.37567420e-03, 3.75783493e-02,-5.74317201e-03,-4.52790031e-02,
   8.81535933e-04, 1.32219792e-02,-4.29163515e-03,-3.57875684e-02,
   1.93781969e-04, 1.84214519e-03,-2.09776080e-03,-1.69429504e-02},
 {-7.41456087e-02,-3.04220386e-03, 1.00778742e+01, 7.44268487e-18,
  -5.64786351e-02,-8.52349656e-04, 4.69854092e-02,-2.61001238e-03,
  -3.48546848e-02, 8.28728239e-04, 3.34986479e-02,-2.74402660e-03,
  -1.36034468e-02, 9.68342359e-04, 3.08808840e-02,-3.11657305e-03,
  -2.02136306e-03, 2.43817039e-04, 1.85860075e-02,-2.58920039e-03},
 { 3.04220386e-03,-7.41456087e-02, 7.44268487e-18, 1.00778742e+01,
   8.52349656e-04,-5.64786351e-02, 2.61001238e-03, 4.69854092e-02,
  -8.28728239e-04,-3.48546848e-02, 2.74402660e-03, 3.34986479e-02,
  -9.68342359e-04,-1.36034468e-02, 3.11657305e-03, 3.08808840e-02,
  -2.43817039e-04,-2.02136306e-03, 2.58920039e-03, 1.85860075e-02},
 { 6.66138672e-02, 9.72205252e-04,-5.64786351e-02, 8.52349656e-04,
   5.50800720e-02, 5.86199074e-18,-3.64630503e-02, 3.61622910e-03,
   3.33821698e-02,-8.51968431e-04,-3.35529420e-02, 3.91409598e-03,
   1.24406140e-02,-7.46818761e-04,-3.12315572e-02, 3.40588757e-03,
   1.78615119e-03,-1.81360190e-04,-1.64290573e-02, 1.97380727e-03},
 {-9.72205252e-04, 6.66138672e-02,-8.52349656e-04,-5.64786351e-02,
   5.86199074e-18, 5.50800720e-02,-3.61622910e-03,-3.64630503e-02,
   8.51968431e-04, 3.33821698e-02,-3.91409598e-03,-3.35529420e-02,
   7.46818761e-04, 1.24406140e-02,-3.40588757e-03,-3.12315572e-02,
   1.81360190e-04, 1.78615119e-03,-1.97380727e-03,-1.64290573e-02},
 {-5.28342800e-02,-5.80545157e-03, 4.69854092e-02, 2.61001238e-03,
  -3.64630503e-02,-3.61622910e-03, 1.00364496e+01, 2.11416734e-18,
  -1.99573930e-02,-1.54086114e-03, 2.58259441e-02,-2.08512077e-04,
  -7.23761678e-03,-2.83270720e-04, 1.80989510e-02, 5.99453787e-05,
  -1.04663678e-03, 7.35636067e-06, 9.64374187e-03,-2.45066884e-04},
 { 5.80545157e-03,-5.28342800e-02,-2.61001238e-03, 4.69854092e-02,
   3.61622910e-03,-3.64630503e-02, 2.11416734e-18, 1.00364496e+01,
   1.54086114e-03,-1.99573930e-02, 2.08512077e-04, 2.58259441e-02,
   2.83270720e-04,-7.23761678e-03,-5.99453787e-05, 1.80989510e-02,
  -7.35636067e-06,-1.04663678e-03, 2.45066884e-04, 9.64374187e-03},
 { 3.75783493e-02, 1.37567420e-03,-3.48546848e-02,-8.28728239e-04,
   3.33821698e-02, 8.51968431e-04,-1.99573930e-02, 1.54086114e-03,
   2.28703256e-02, 1.80522336e-18,-1.77682697e-02, 1.85067730e-03,
   9.01559121e-03,-3.30613348e-04,-1.99294691e-02, 1.69797998e-03,
   1.35415367e-03,-1.08287357e-04,-1.24604821e-02, 1.22725453e-03},
 {-1.37567420e-03, 3.75783493e-02, 8.28728239e-04,-3.48546848e-02,
  -8.51968431e-04, 3.33821698e-02,-1.54086114e-03,-1.99573930e-02,
   1.80522336e-18, 2.28703256e-02,-1.85067730e-03,-1.77682697e-02,
   3.30613348e-04, 9.01559121e-03,-1.69797998e-03,-1.99294691e-02,
   1.08287357e-04, 1.35415367e-03,-1.22725453e-03,-1.24604821e-02},
 {-4.52790031e-02,-5.74317201e-03, 3.34986479e-02, 2.74402660e-03,
  -3.35529420e-02,-3.91409598e-03, 2.58259441e-02, 2.08512077e-04,
  -1.77682697e-02,-1.85067730e-03, 1.00255022e+01, 1.03499330e-18,
  -5.78626868e-03,-5.07824831e-04, 1.77973026e-02, 3.30978516e-04,
  -7.67950707e-04,-4.61836167e-05, 7.08465691e-03, 2.95518835e-04},
 { 5.74317201e-03,-4.52790031e-02,-2.74402660e-03, 3.34986479e-02,
   3.91409598e-03,-3.35529420e-02,-2.08512077e-04, 2.58259441e-02,
   1.85067730e-03,-1.77682697e-02, 1.03499330e-18, 1.00255022e+01,
   5.07824831e-04,-5.78626868e-03,-3.30978516e-04, 1.77973026e-02,
   4.61836167e-05,-7.67950707e-04,-2.95518835e-04, 7.08465691e-03},
 { 1.32219792e-02, 8.81535933e-04,-1.36034468e-02,-9.68342359e-04,
   1.24406140e-02, 7.46818761e-04,-7.23761678e-03, 2.83270720e-04,
   9.01559121e-03, 3.30613348e-04,-5.78626868e-03, 5.07824831e-04,
   4.42864169e-03, 2.39642515e-19,-7.42628184e-03, 3.94705035e-04,
   6.18381101e-04,-2.70144605e-05,-5.69394205e-03, 3.53683845e-04},
 {-8.81535933e-04, 1.32219792e-02, 9.68342359e-04,-1.36034468e-02,
  -7.46818761e-04, 1.24406140e-02,-2.83270720e-04,-7.23761678e-03,
  -3.30613348e-04, 9.01559121e-03,-5.07824831e-04,-5.78626868e-03,
   2.39642515e-19, 4.42864169e-03,-3.94705035e-04,-7.42628184e-03,
   2.70144605e-05, 6.18381101e-04,-3.53683845e-04,-5.69394205e-03},
 {-3.57875684e-02,-4.29163515e-03, 3.08808840e-02, 3.11657305e-03,
  -3.12315572e-02,-3.40588757e-03, 1.80989510e-02,-5.99453787e-05,
  -1.99294691e-02,-1.69797998e-03, 1.77973026e-02,-3.30978516e-04,
  -7.42628184e-03,-3.94705035e-04, 1.00195875e+01, 5.29144694e-19,
  -1.05319349e-03,-1.51209590e-05, 9.70797054e-03,-3.90440034e-05},
 { 4.29163515e-03,-3.57875684e-02,-3.11657305e-03, 3.08808840e-02,
   3.40588757e-03,-3.12315572e-02, 5.99453787e-05, 1.80989510e-02,
   1.69797998e-03,-1.99294691e-02, 3.30978516e-04, 1.77973026e-02,
   3.94705035e-04,-7.42628184e-03, 5.29144694e-19, 1.00195875e+01,
   1.51209590e-05,-1.05319349e-03, 3.90440034e-05, 9.70797054e-03},
 { 1.84214519e-03, 1.93781969e-04,-2.02136306e-03,-2.43817039e-04,
   1.78615119e-03, 1.81360190e-04,-1.04663678e-03,-7.35636067e-06,
   1.35415367e-03, 1.08287357e-04,-7.67950707e-04, 4.61836167e-05,
   6.18381101e-04, 2.70144605e-05,-1.05319349e-03, 1.51209590e-05,
   6.00069946e-04, 6.06130173e-21,-9.22166534e-04, 1.69495670e-05},
 {-1.93781969e-04, 1.84214519e-03, 2.43817039e-04,-2.02136306e-03,
  -1.81360190e-04, 1.78615119e-03, 7.35636067e-06,-1.04663678e-03,
  -1.08287357e-04, 1.35415367e-03,-4.61836167e-05,-7.67950707e-04,
  -2.70144605e-05, 6.18381101e-04,-1.51209590e-05,-1.05319349e-03,
   6.06130173e-21, 6.00069946e-04,-1.69495670e-05,-9.22166534e-04},
 {-1.69429504e-02,-2.09776080e-03, 1.85860075e-02, 2.58920039e-03,
  -1.64290573e-02,-1.97380727e-03, 9.64374187e-03, 2.45066884e-04,
  -1.24604821e-02,-1.22725453e-03, 7.08465691e-03,-2.95518835e-04,
  -5.69394205e-03,-3.53683845e-04, 9.70797054e-03, 3.90440034e-05,
  -9.22166534e-04,-1.69495670e-05, 1.00085008e+01, 8.24199638e-20},
 { 2.09776080e-03,-1.69429504e-02,-2.58920039e-03, 1.85860075e-02,
   1.97380727e-03,-1.64290573e-02,-2.45066884e-04, 9.64374187e-03,
   1.22725453e-03,-1.24604821e-02, 2.95518835e-04, 7.08465691e-03,
   3.53683845e-04,-5.69394205e-03,-3.90440034e-05, 9.70797054e-03,
   1.69495670e-05,-9.22166534e-04, 8.24199638e-20, 1.00085008e+01}};

float DIM_Fj[20];

float DIM_M[12][20] = {{-1.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00},
 {-0.00000000e+00,-1.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00},
 {-0.00000000e+00,-0.00000000e+00,-1.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00},
 {-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-1.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00},
 { 1.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
 { 0.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
 { 0.00000000e+00, 0.00000000e+00, 1.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
 { 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
 { 5.34946300e-02, 1.60032361e-03,-9.99260520e-03,-4.66676790e-04,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00},
 {-1.60032361e-03, 5.34946300e-02, 4.66676790e-04,-9.99260520e-03,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00,
  -0.00000000e+00,-0.00000000e+00,-0.00000000e+00,-0.00000000e+00},
 {-5.34946300e-02,-1.60032361e-03, 9.99260520e-03, 4.66676790e-04,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
 { 1.60032361e-03,-5.34946300e-02,-4.66676790e-04, 9.99260520e-03,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00}};

float DIM_gam[12];

 /* Matrices for Hildreth's QP procedure */
float DIM_Fj_1[20][2] = {{ 5.66889900e-01,-7.75571854e-02},
 { 7.75571854e-02, 5.66889900e-01},
 {-5.79489019e-01, 4.43681426e-02},
 {-4.43681426e-02,-5.79489019e-01},
 { 3.83843972e-01,-4.96159980e-02},
 { 4.96159980e-02, 3.83843972e-01},
 {-3.78780415e-01, 1.05948593e-02},
 {-1.05948593e-02,-3.78780415e-01},
 { 2.06331464e-01,-2.31763307e-02},
 { 2.31763307e-02, 2.06331464e-01},
 {-2.74301542e-01, 4.97928413e-03},
 {-4.97928413e-03,-2.74301542e-01},
 { 7.16158141e-02,-6.04806229e-03},
 { 6.04806229e-03, 7.16158141e-02},
 {-1.97401520e-01, 6.00966220e-03},
 {-6.00966220e-03,-1.97401520e-01},
 { 9.99260520e-03,-4.66676790e-04},
 { 4.66676790e-04, 9.99260520e-03},
 {-9.21630964e-02, 2.60800967e-03},
 {-2.60800967e-03,-9.21630964e-02}};

float DIM_Fj_2[20][8] = {{-5.31658102e-01,-3.08652621e-02,-1.64652804e+00,-1.01891367e-02,
   1.81992378e-02,-7.23455948e-03,-5.66889900e-01, 7.75571854e-02},
 { 3.08652621e-02,-5.31658102e-01, 1.01891367e-02,-1.64652804e+00,
   7.23455948e-03, 1.81992378e-02,-7.75571854e-02,-5.66889900e-01},
 { 5.47883118e-01, 4.17638324e-02, 1.41118960e+00, 7.15540182e-02,
  -2.74872286e-02, 4.03427195e-03, 5.79489019e-01,-4.43681426e-02},
 {-4.17638324e-02, 5.47883118e-01,-7.15540182e-02, 1.41118960e+00,
  -4.03427195e-03,-2.74872286e-02, 4.43681426e-02, 5.79489019e-01},
 {-4.50823619e-01,-3.69629723e-02,-1.18910179e+00,-2.34359322e-02,
   6.33800738e-03,-5.48850292e-03,-3.83843972e-01, 4.96159980e-02},
 { 3.69629723e-02,-4.50823619e-01, 2.34359322e-02,-1.18910179e+00,
   5.48850292e-03, 6.33800738e-03,-4.96159980e-02,-3.83843972e-01},
 { 2.95942510e-01, 4.65122985e-02, 1.02083581e+00, 1.17132346e-01,
  -1.95175277e-02, 2.53255181e-03, 3.78780415e-01,-1.05948593e-02},
 {-4.65122985e-02, 2.95942510e-01,-1.17132346e-01, 1.02083581e+00,
  -2.53255181e-03,-1.95175277e-02, 1.05948593e-02, 3.78780415e-01},
 {-2.97452123e-01,-3.47492769e-02,-6.53789565e-01,-2.60100176e-02,
   1.80225234e-03,-3.15690083e-03,-2.06331464e-01, 2.31763307e-02},
 { 3.47492769e-02,-2.97452123e-01, 2.60100176e-02,-6.53789565e-01,
   3.15690083e-03, 1.80225234e-03,-2.31763307e-02,-2.06331464e-01},
 { 2.38171857e-01, 4.69450564e-02, 8.28863466e-01, 1.08998570e-01,
  -6.56177284e-03, 3.19406320e-03, 2.74301542e-01,-4.97928413e-03},
 {-4.69450564e-02, 2.38171857e-01,-1.08998570e-01, 8.28863466e-01,
  -3.19406320e-03,-6.56177284e-03, 4.97928413e-03, 2.74301542e-01},
 {-1.18969843e-01,-1.90180744e-02,-2.27375996e-01,-1.55830762e-02,
   6.87690949e-04,-1.06148415e-03,-7.16158141e-02, 6.04806229e-03},
 { 1.90180744e-02,-1.18969843e-01, 1.55830762e-02,-2.27375996e-01,
   1.06148415e-03, 6.87690949e-04,-6.04806229e-03,-7.16158141e-02},
 { 2.62374020e-01, 5.27612432e-02, 6.20477767e-01, 7.60569374e-02,
  -8.07142824e-04, 3.15433995e-03, 1.97401520e-01,-6.00966220e-03},
 {-5.27612432e-02, 2.62374020e-01,-7.60569374e-02, 6.20477767e-01,
  -3.15433995e-03,-8.07142824e-04, 6.00966220e-03, 1.97401520e-01},
 {-1.77521862e-02,-3.67811795e-03,-3.15678299e-02,-3.36408587e-03,
   1.27532932e-04,-1.35752435e-04,-9.99260520e-03, 4.66676790e-04},
 { 3.67811795e-03,-1.77521862e-02, 3.36408587e-03,-3.15678299e-02,
   1.35752435e-04, 1.27532932e-04,-4.66676790e-04,-9.99260520e-03},
 { 1.62967306e-01, 3.69014804e-02, 2.90334688e-01, 3.63476609e-02,
  -1.19823735e-03, 1.22938734e-03, 9.21630964e-02,-2.60800967e-03},
 {-3.69014804e-02, 1.62967306e-01,-3.63476609e-02, 2.90334688e-01,
  -1.22938734e-03,-1.19823735e-03, 2.60800967e-03, 9.21630964e-02}};

float DIM_Fx[2][6] = {{ 0.25425333, 0.01599626, 0.74377339, 0.04679426, 0.04343702, 0.00273283},
 {-0.01599626, 0.25425333,-0.04679426, 0.74377339,-0.00273283, 0.04343702}};

float DIM_Kj_1[12][20] = {{-1.74790567e+02,-9.51597663e-14,-1.17200226e-01, 9.62488105e-03,
   3.16729036e+02,-9.98490369e+00,-1.14180542e-01, 8.41429746e-03,
  -1.51475211e+02, 1.40412469e+01,-3.24627802e-02, 2.77939846e-03,
  -5.88621375e+01, 6.80007231e+00, 1.74730554e-02,-2.99687302e-03,
  -4.11382805e+00, 4.14199484e-01, 1.89198366e-03,-2.25686395e-04},
 {-8.48499093e-14,-1.74790567e+02,-9.62488105e-03,-1.17200226e-01,
   9.98490369e+00, 3.16729036e+02,-8.41429746e-03,-1.14180542e-01,
  -1.40412469e+01,-1.51475211e+02,-2.77939846e-03,-3.24627802e-02,
  -6.80007231e+00,-5.88621375e+01, 2.99687302e-03, 1.74730554e-02,
  -4.14199484e-01,-4.11382805e+00, 2.25686395e-04, 1.89198366e-03},
 {-1.17200226e-01,-9.62488105e-03,-9.98859666e-02,-6.48074347e-21,
   1.17366657e-01,-3.58445961e-04, 2.17006798e-05, 4.86734188e-07,
  -8.85322603e-02, 1.24447851e-02,-2.29700054e-05, 2.34158161e-06,
  -1.02445284e-01, 1.56675951e-02,-6.07127707e-07,-4.06492526e-07,
  -2.03839237e-02, 4.00409916e-03, 9.35820657e-06,-2.01756126e-06},
 { 9.62488105e-03,-1.17200226e-01, 5.01067163e-21,-9.98859666e-02,
   3.58445961e-04, 1.17366657e-01,-4.86734188e-07, 2.17006798e-05,
  -1.24447851e-02,-8.85322603e-02,-2.34158161e-06,-2.29700054e-05,
  -1.56675951e-02,-1.02445284e-01, 4.06492526e-07,-6.07127707e-07,
  -4.00409916e-03,-2.03839237e-02, 2.01756126e-06, 9.35820657e-06},
 { 1.74790567e+02, 9.51597663e-14, 1.17200226e-01,-9.62488105e-03,
  -3.16729036e+02, 9.98490369e+00, 1.14180542e-01,-8.41429746e-03,
   1.51475211e+02,-1.40412469e+01, 3.24627802e-02,-2.77939846e-03,
   5.88621375e+01,-6.80007231e+00,-1.74730554e-02, 2.99687302e-03,
   4.11382805e+00,-4.14199484e-01,-1.89198366e-03, 2.25686395e-04},
 { 8.48499093e-14, 1.74790567e+02, 9.62488105e-03, 1.17200226e-01,
  -9.98490369e+00,-3.16729036e+02, 8.41429746e-03, 1.14180542e-01,
   1.40412469e+01, 1.51475211e+02, 2.77939846e-03, 3.24627802e-02,
   6.80007231e+00, 5.88621375e+01,-2.99687302e-03,-1.74730554e-02,
   4.14199484e-01, 4.11382805e+00,-2.25686395e-04,-1.89198366e-03},
 { 1.17200226e-01, 9.62488105e-03, 9.98859666e-02, 6.48074347e-21,
  -1.17366657e-01, 3.58445961e-04,-2.17006798e-05,-4.86734188e-07,
   8.85322603e-02,-1.24447851e-02, 2.29700054e-05,-2.34158161e-06,
   1.02445284e-01,-1.56675951e-02, 6.07127707e-07, 4.06492526e-07,
   2.03839237e-02,-4.00409916e-03,-9.35820657e-06, 2.01756126e-06},
 {-9.62488105e-03, 1.17200226e-01,-5.01067163e-21, 9.98859666e-02,
  -3.58445961e-04,-1.17366657e-01, 4.86734188e-07,-2.17006798e-05,
   1.24447851e-02, 8.85322603e-02, 2.34158161e-06, 2.29700054e-05,
   1.56675951e-02, 1.02445284e-01,-4.06492526e-07, 6.07127707e-07,
   4.00409916e-03, 2.03839237e-02,-2.01756126e-06,-9.35820657e-06},
 { 9.34919005e+00, 2.79570599e-01, 5.28686460e-03,-3.73935624e-04,
  -1.69581087e+01, 2.73209641e-02, 6.12172806e-03,-2.67378921e-04,
   8.12469041e+00,-5.08638912e-01, 1.74080173e-03,-9.67192598e-05,
   3.15865957e+00,-2.69460133e-01,-9.39516476e-04, 1.32349725e-04,
   2.20525005e-01,-1.55434933e-02,-1.01477683e-04, 9.02943062e-06},
 {-2.79570599e-01, 9.34919005e+00, 3.73935624e-04, 5.28686460e-03,
  -2.73209641e-02,-1.69581087e+01, 2.67378921e-04, 6.12172806e-03,
   5.08638912e-01, 8.12469041e+00, 9.67192598e-05, 1.74080173e-03,
   2.69460133e-01, 3.15865957e+00,-1.32349725e-04,-9.39516476e-04,
   1.55434933e-02, 2.20525005e-01,-9.02943062e-06,-1.01477683e-04},
 {-9.34919005e+00,-2.79570599e-01,-5.28686460e-03, 3.73935624e-04,
   1.69581087e+01,-2.73209641e-02,-6.12172806e-03, 2.67378921e-04,
  -8.12469041e+00, 5.08638912e-01,-1.74080173e-03, 9.67192598e-05,
  -3.15865957e+00, 2.69460133e-01, 9.39516476e-04,-1.32349725e-04,
  -2.20525005e-01, 1.55434933e-02, 1.01477683e-04,-9.02943062e-06},
 { 2.79570599e-01,-9.34919005e+00,-3.73935624e-04,-5.28686460e-03,
   2.73209641e-02, 1.69581087e+01,-2.67378921e-04,-6.12172806e-03,
  -5.08638912e-01,-8.12469041e+00,-9.67192598e-05,-1.74080173e-03,
  -2.69460133e-01,-3.15865957e+00, 1.32349725e-04, 9.39516476e-04,
  -1.55434933e-02,-2.20525005e-01, 9.02943062e-06, 1.01477683e-04}};

float DIM_Hj[12][12] = {{-5.72113254e-03, 9.51597663e-14, 1.17200226e-01,-9.62488105e-03,
  -1.74790567e+02,-9.51597663e-14,-1.17200226e-01, 9.62488105e-03,
  -9.34919005e+00, 2.79570599e-01, 9.34919005e+00,-2.79570599e-01},
 { 8.48499093e-14,-5.72113254e-03, 9.62488105e-03, 1.17200226e-01,
  -8.48499093e-14,-1.74790567e+02,-9.62488105e-03,-1.17200226e-01,
  -2.79570599e-01,-9.34919005e+00, 2.79570599e-01, 9.34919005e+00},
 { 1.17200226e-01, 9.62488105e-03,-1.00114164e+01, 6.48074347e-21,
  -1.17200226e-01,-9.62488105e-03,-9.98859666e-02,-6.48074347e-21,
  -5.28686460e-03,-3.73935624e-04, 5.28686460e-03, 3.73935624e-04},
 {-9.62488105e-03, 1.17200226e-01,-5.01067163e-21,-1.00114164e+01,
   9.62488105e-03,-1.17200226e-01, 5.01067163e-21,-9.98859666e-02,
   3.73935624e-04,-5.28686460e-03,-3.73935624e-04, 5.28686460e-03},
 {-1.74790567e+02,-9.51597663e-14,-1.17200226e-01, 9.62488105e-03,
  -5.72113254e-03, 9.51597663e-14, 1.17200226e-01,-9.62488105e-03,
   9.34919005e+00,-2.79570599e-01,-9.34919005e+00, 2.79570599e-01},
 {-8.48499093e-14,-1.74790567e+02,-9.62488105e-03,-1.17200226e-01,
   8.48499093e-14,-5.72113254e-03, 9.62488105e-03, 1.17200226e-01,
   2.79570599e-01, 9.34919005e+00,-2.79570599e-01,-9.34919005e+00},
 {-1.17200226e-01,-9.62488105e-03,-9.98859666e-02,-6.48074347e-21,
   1.17200226e-01, 9.62488105e-03,-1.00114164e+01, 6.48074347e-21,
   5.28686460e-03, 3.73935624e-04,-5.28686460e-03,-3.73935624e-04},
 { 9.62488105e-03,-1.17200226e-01, 5.01067163e-21,-9.98859666e-02,
  -9.62488105e-03, 1.17200226e-01,-5.01067163e-21,-1.00114164e+01,
  -3.73935624e-04, 5.28686460e-03, 3.73935624e-04,-5.28686460e-03},
 {-9.34919005e+00,-2.79570599e-01,-5.28686460e-03, 3.73935624e-04,
   9.34919005e+00, 2.79570599e-01, 5.28686460e-03,-3.73935624e-04,
  -1.99789737e+00, 1.66013672e-17,-5.00526211e-01,-1.66013672e-17},
 { 2.79570599e-01,-9.34919005e+00,-3.73935624e-04,-5.28686460e-03,
  -2.79570599e-01, 9.34919005e+00, 3.73935624e-04, 5.28686460e-03,
  -1.70048741e-17,-1.99789737e+00, 1.70048741e-17,-5.00526211e-01},
 { 9.34919005e+00, 2.79570599e-01, 5.28686460e-03,-3.73935624e-04,
  -9.34919005e+00,-2.79570599e-01,-5.28686460e-03, 3.73935624e-04,
  -5.00526211e-01,-1.66013672e-17,-1.99789737e+00, 1.66013672e-17},
 {-2.79570599e-01, 9.34919005e+00, 3.73935624e-04, 5.28686460e-03,
   2.79570599e-01,-9.34919005e+00,-3.73935624e-04,-5.28686460e-03,
   1.70048741e-17,-5.00526211e-01,-1.70048741e-17,-1.99789737e+00}};

float DIM_DU_1[4][20] = {{-1.74790567e+02,-9.51597663e-14,-1.17200226e-01, 9.62488105e-03,
   3.16729036e+02,-9.98490369e+00,-1.14180542e-01, 8.41429746e-03,
  -1.51475211e+02, 1.40412469e+01,-3.24627802e-02, 2.77939846e-03,
  -5.88621375e+01, 6.80007231e+00, 1.74730554e-02,-2.99687302e-03,
  -4.11382805e+00, 4.14199484e-01, 1.89198366e-03,-2.25686395e-04},
 {-8.48499093e-14,-1.74790567e+02,-9.62488105e-03,-1.17200226e-01,
   9.98490369e+00, 3.16729036e+02,-8.41429746e-03,-1.14180542e-01,
  -1.40412469e+01,-1.51475211e+02,-2.77939846e-03,-3.24627802e-02,
  -6.80007231e+00,-5.88621375e+01, 2.99687302e-03, 1.74730554e-02,
  -4.14199484e-01,-4.11382805e+00, 2.25686395e-04, 1.89198366e-03},
 {-1.17200226e-01,-9.62488105e-03,-9.98859666e-02,-6.48074347e-21,
   1.17366657e-01,-3.58445961e-04, 2.17006798e-05, 4.86734188e-07,
  -8.85322603e-02, 1.24447851e-02,-2.29700054e-05, 2.34158161e-06,
  -1.02445284e-01, 1.56675951e-02,-6.07127707e-07,-4.06492526e-07,
  -2.03839237e-02, 4.00409916e-03, 9.35820657e-06,-2.01756126e-06},
 { 9.62488105e-03,-1.17200226e-01, 5.01067163e-21,-9.98859666e-02,
   3.58445961e-04, 1.17366657e-01,-4.86734188e-07, 2.17006798e-05,
  -1.24447851e-02,-8.85322603e-02,-2.34158161e-06,-2.29700054e-05,
  -1.56675951e-02,-1.02445284e-01, 4.06492526e-07,-6.07127707e-07,
  -4.00409916e-03,-2.03839237e-02, 2.01756126e-06, 9.35820657e-06}};

float DIM_DU_2[4][12] = {{ 1.74790567e+02, 9.51597663e-14, 1.17200226e-01,-9.62488105e-03,
  -1.74790567e+02,-9.51597663e-14,-1.17200226e-01, 9.62488105e-03,
  -9.34919005e+00, 2.79570599e-01, 9.34919005e+00,-2.79570599e-01},
 { 8.48499093e-14, 1.74790567e+02, 9.62488105e-03, 1.17200226e-01,
  -8.48499093e-14,-1.74790567e+02,-9.62488105e-03,-1.17200226e-01,
  -2.79570599e-01,-9.34919005e+00, 2.79570599e-01, 9.34919005e+00},
 { 1.17200226e-01, 9.62488105e-03, 9.98859666e-02, 6.48074347e-21,
  -1.17200226e-01,-9.62488105e-03,-9.98859666e-02,-6.48074347e-21,
  -5.28686460e-03,-3.73935624e-04, 5.28686460e-03, 3.73935624e-04},
 {-9.62488105e-03, 1.17200226e-01,-5.01067163e-21, 9.98859666e-02,
   9.62488105e-03,-1.17200226e-01, 5.01067163e-21,-9.98859666e-02,
   3.73935624e-04,-5.28686460e-03,-3.73935624e-04, 5.28686460e-03}};


#endif /* DMPC_INVERTER_MATRICES_H_ */