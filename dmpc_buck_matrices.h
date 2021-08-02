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
#define DMPC_BUCK_CONFIG_NP			20
#define DMPC_BUCK_CONFIG_NC			2
#define DMPC_BUCK_CONFIG_NR			1
#define DMPC_BUCK_CONFIG_NLAMBDA	4

/* Constraints */
#define DMPC_BUCK_CONFIG_IL_MIN		-8
#define DMPC_BUCK_CONFIG_IL_MAX		8
#define DMPC_BUCK_CONFIG_U_MIN		0
#define DMPC_BUCK_CONFIG_U_MAX		1

 /* Matrices */
float Fj_1[2] = {-68.57386378,-60.44625265};

float Fj_2[2][3] = {{116.95950028,707.83900161, 68.57386378},
 {105.79037349,633.55081889, 60.44625265}};

float Fx[3] = { 0.98236791,-0.20888864, 0.        };

float Kj_1[4][2] = {{-0.0230602 , 0.02386833},
 { 0.0230602 ,-0.02386833},
 {-0.06797821, 0.07036046},
 { 0.06797821,-0.07036046}};

float Hj[4][4] = {{-4.33647650e+01,-2.30601965e-02, 6.79782069e-02,-6.79782069e-02},
 {-2.30601965e-02,-4.33647650e+01,-6.79782069e-02, 6.79782069e-02},
 { 6.79782069e-02,-6.79782069e-02,-4.99026482e+00,-2.00390167e-01},
 {-6.79782069e-02, 6.79782069e-02,-2.00390167e-01,-4.99026482e+00}};

float DU_1[2] = {-0.0230602 , 0.02386833};

float DU_2[4] = { 0.0230602 ,-0.0230602 , 0.06797821,-0.06797821};


#endif /* DMPC_BUCK_MATRICES_H_ */
