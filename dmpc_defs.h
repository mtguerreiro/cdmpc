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
#define DMPC_CONFIG_SCALE      200.000000f

/* Number of model states and augmented states */
#define DMPC_CONFIG_NXM        6
#define DMPC_CONFIG_NXA        8

/* Prediction, control and constraint horizon */
#define DMPC_CONFIG_NP         4
#define DMPC_CONFIG_NC         4
#define DMPC_CONFIG_NR         4
#define DMPC_CONFIG_NLAMBDA    32

/* Number of inputs, outputs, and disturbances */
#define DMPC_CONFIG_NU         2
#define DMPC_CONFIG_NY         2
#define DMPC_CONFIG_ND         2

/* Size of control vector */
#define DMPC_CONFIG_NC_x_NU    (DMPC_CONFIG_NC * DMPC_CONFIG_NU)

/* Input constraints */
#define DMPC_CONFIG_NU_CTR     2

/* Input constraints */
#define DMPC_CONFIG_NXM_CTR    2

#endif /* DMPC_DEFS_H_ */
