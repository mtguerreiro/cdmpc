/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xgerc.c
 *
 * MATLAB Coder version            : 5.0
 * C/C++ source code generated on  : 13-Sep-2021 17:45:18
 */

/* Include Files */
#include "xgerc.h"
#include "qpm2.h"

/* Function Definitions */

/*
 * Arguments    : int m
 *                int n
 *                float alpha1
 *                int ix0
 *                const float y[10]
 *                float A[100]
 *                int ia0
 * Return Type  : void
 */
void xgerc(int m, int n, float alpha1, int ix0, const float y[10], float A[100],
           int ia0)
{
  int jA;
  int jy;
  int j;
  float temp;
  int ix;
  int i;
  int ijA;
  if (alpha1 != 0.0F) {
    jA = ia0;
    jy = 0;
    for (j = 0; j < n; j++) {
      if (y[jy] != 0.0F) {
        temp = y[jy] * alpha1;
        ix = ix0;
        i = m + jA;
        for (ijA = jA; ijA < i; ijA++) {
          A[ijA - 1] += A[ix - 1] * temp;
          ix++;
        }
      }

      jy++;
      jA += 10;
    }
  }
}

/*
 * File trailer for xgerc.c
 *
 * [EOF]
 */
