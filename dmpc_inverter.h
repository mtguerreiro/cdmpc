/*
 * dmpc_inverter.h
 *
 *  Created on: 22.04.2022
 *      Author: mguerreiro
 */

#ifndef DMPC_INVERTER_H_
#define DMPC_INVERTER_H_

//=============================================================================
/*-------------------------------- Includes ---------------------------------*/
//=============================================================================
#include "stdint.h"

//=============================================================================

//=============================================================================
/*-------------------------------- Functions --------------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
void dmpcInverterInitialize(void);
//-----------------------------------------------------------------------------
void dmpcInverterOpt(float *x, float *x_1, float *r, float *u_1, float *du, uint32_t* iters);
//-----------------------------------------------------------------------------
//=============================================================================

#endif /* DMPC_INVERTER_H_ */