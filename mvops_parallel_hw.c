/*
 * mvops_parallel_hw.c
 *
 *  Created on: 03.11.2021
 *      Author: mguerreiro
 */

//=============================================================================
/*-------------------------------- Includes ---------------------------------*/
//=============================================================================
#include "mvops_parallel_hw.h"

/* Libs */
#include "mvops.h"

/* Only required for automatic definition of flags */
#include "platdefs.h"

//=============================================================================

typedef void (*mvopsParallelOp_t)(float *x, float *y, float *z, uint32_t dim);



#define MVOPS_PAR_IPC_FLAG_MAIN    PLAT_IPC_FLAG_MVOPS_PAR_MAIN
#define MVOPS_PAR_IPC_FLAG_SEC     PLAT_IPC_FLAG_MVOPS_PAR_SEC

//=============================================================================
/*-------------------------------- Functions --------------------------------*/
//=============================================================================
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void mvopsParallelHWSecondarySignal(uint32_t op, float *x, float *y, uint32_t dim){

    uint32_t data;

    data = dim & 0x00FFFFFF;
    data |= ((op & 0xFF) << 24);

    HWREG(IPC_BASE + IPC_O_SENDCOM) = data;
    HWREG(IPC_BASE + IPC_O_SENDDATA) = (uint32_t)(x);
    HWREG(IPC_BASE + IPC_O_SENDADDR) = (uint32_t)(y);

    HWREG(IPC_BASE + IPC_O_SET) = 1UL << MVOPS_PAR_IPC_FLAG_SEC;
}
//-----------------------------------------------------------------------------
//=============================================================================
