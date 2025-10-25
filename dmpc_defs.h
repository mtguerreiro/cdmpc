/**
 * @file dmpc_defs.h
 * @brief Header with definitions to aid the DMPC algorithm.
 *
 * This file is generated automatically and should not be modified.
 *
 */

#ifndef DMPC_DEFS_H_
#define DMPC_DEFS_H_

/* Scaling factor */
#define DMPC_CONFIG_SCALE              1.000000f

/* Number of model states and augmented states */
#define DMPC_CONFIG_NXM                2
#define DMPC_CONFIG_NXA                3

/* Length of prediction, control and constraint horizons */
#define DMPC_CONFIG_L_PRED             10
#define DMPC_CONFIG_L_CTL              4
#define DMPC_CONFIG_L_U_CNT            4
#define DMPC_CONFIG_L_X_CNT            4
#define DMPC_CONFIG_NLAMBDA            16

/* Number of inputs, outputs, and disturbances */
#define DMPC_CONFIG_NU                 1
#define DMPC_CONFIG_NY                 1
#define DMPC_CONFIG_ND                 0

/* Size of control vector */
#define DMPC_CONFIG_U_SIZE             (DMPC_CONFIG_L_CTL * DMPC_CONFIG_NU)

/* Input constraints */
#define DMPC_CONFIG_NU_CNT             1

/* State constraints */
#define DMPC_CONFIG_NXM_CNT            1

/* Solver settings */
#define DMPC_CONFIG_HILD_TOL           1e-06
#define DMPC_CONFIG_HILD_N_ITER        200
#define DMPC_CONFIG_HILD_FIXED_ITER    1

#if !defined(DMPC_CONFIG_SOLVER_HILD) && !defined(DMPC_CONFIG_SOLVER_OSQP)
#define DMPC_CONFIG_SOLVER_HILD
#endif

#endif /* DMPC_DEFS_H_ */
