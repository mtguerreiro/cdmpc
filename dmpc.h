/*
 * dmpc.h
 *
 *  Created on: 08.09.2022
 *      Author: mguerreiro
 */

#ifndef DMPC_H_
#define DMPC_H_

//=============================================================================
/*-------------------------------- Includes ---------------------------------*/
//=============================================================================
#include "stdint.h"

//=============================================================================

//=============================================================================
/*-------------------------------- Functions --------------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------
uint32_t dmpcOpt(float *x, float *x_1, float *r, float *u_1, uint32_t *niters, float *du);
//-----------------------------------------------------------------------------
void dmpcDelayComp(float *x_1, float *x, float *u);
//-----------------------------------------------------------------------------
//=============================================================================

#endif /* DMPC_H_ */
