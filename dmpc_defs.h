/**
 * @file dmpc_defs.h
 * @brief Header with definitions to aid the DMPC algorithm.
 *
 * This file is generated automatically and should not be modified.
 *
 *  Originally created on: 21.10.2022
 *      Author: mguerreiro
 */

#ifndef DMPC_DEFS_H_
#define DMPC_DEFS_H_

/* Scaling factor */
#define DMPC_CONFIG_SCALE			1.000000f

/* Number of model states and augmented states */
#define DMPC_CONFIG_NXM			2
#define DMPC_CONFIG_NXA			3

/* Prediction, control and constraint horizon */
#define DMPC_CONFIG_NP			5
#define DMPC_CONFIG_NC			5
#define DMPC_CONFIG_NR			1
#define DMPC_CONFIG_NLAMBDA		4

/* Number of inputs and outputs */
#define DMPC_CONFIG_NU			1
#define DMPC_CONFIG_NY			1

/* Number of external disturbances */
#define DMPC_CONFIG_ND			0

/* Size of control vector */
#define DMPC_CONFIG_NC_x_NU		DMPC_CONFIG_NC * DMPC_CONFIG_NU

/* Input constraints */
#define DMPC_CONFIG_NU_CTR		1

/* State constraints */
#define DMPC_CONFIG_NXM_CTR		1

#endif /* DMPC_DEFS_H_ */
