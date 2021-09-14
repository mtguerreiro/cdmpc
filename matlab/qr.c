/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: qr.c
 *
 * MATLAB Coder version            : 5.0
 * C/C++ source code generated on  : 13-Sep-2021 17:45:18
 */

/* Include Files */
#include "qr.h"
#include "qpm2.h"
#include "xgerc.h"
#include "xnrm2.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static float rt_hypotf(float u0, float u1);

/* Function Definitions */

/*
 * Arguments    : float u0
 *                float u1
 * Return Type  : float
 */
static float rt_hypotf(float u0, float u1)
{
  float y;
  float a;
  float b;
  a = fabsf(u0);
  b = fabsf(u1);
  if (a < b) {
    a /= b;
    y = b * sqrtf(a * a + 1.0F);
  } else if (a > b) {
    b /= a;
    y = a * sqrtf(b * b + 1.0F);
  } else {
    y = a * 1.41421354F;
  }

  return y;
}

/*
 * Arguments    : const float A[100]
 *                float Q[100]
 *                float R[100]
 * Return Type  : void
 */
void qr(const float A[100], float Q[100], float R[100])
{
  float b_A[100];
  int i;
  float tau[10];
  float work[10];
  int itau;
  int ii;
  int ix0;
  float atmp;
  int b_i;
  int knt;
  float c;
  float beta1;
  int lastv;
  int lastc;
  bool exitg2;
  int k;
  int ia;
  int ix;
  long exitg1;
  int i1;
  memcpy(&b_A[0], &A[0], 100U * sizeof(float));
  for (i = 0; i < 10; i++) {
    tau[i] = 0.0F;
    work[i] = 0.0F;
  }

  for (i = 0; i < 10; i++) {
    ii = i * 10 + i;
    if (i + 1 < 10) {
      atmp = b_A[ii];
      ix0 = ii + 2;
      tau[i] = 0.0F;
      c = xnrm2(9 - i, b_A, ii + 2);
      if (c != 0.0F) {
        beta1 = rt_hypotf(b_A[ii], c);
        if (b_A[ii] >= 0.0F) {
          beta1 = -beta1;
        }

        if (fabsf(beta1) < 9.86076132E-32F) {
          knt = -1;
          b_i = (ii - i) + 10;
          do {
            knt++;
            for (k = ix0; k <= b_i; k++) {
              b_A[k - 1] *= 1.01412048E+31F;
            }

            beta1 *= 1.01412048E+31F;
            atmp *= 1.01412048E+31F;
          } while (!(fabsf(beta1) >= 9.86076132E-32F));

          beta1 = rt_hypotf(atmp, xnrm2(9 - i, b_A, ii + 2));
          if (atmp >= 0.0F) {
            beta1 = -beta1;
          }

          tau[i] = (beta1 - atmp) / beta1;
          c = 1.0F / (atmp - beta1);
          for (k = ix0; k <= b_i; k++) {
            b_A[k - 1] *= c;
          }

          for (k = 0; k <= knt; k++) {
            beta1 *= 9.86076132E-32F;
          }

          atmp = beta1;
        } else {
          tau[i] = (beta1 - b_A[ii]) / beta1;
          c = 1.0F / (b_A[ii] - beta1);
          b_i = (ii - i) + 10;
          for (k = ix0; k <= b_i; k++) {
            b_A[k - 1] *= c;
          }

          atmp = beta1;
        }
      }

      b_A[ii] = 1.0F;
      if (tau[i] != 0.0F) {
        lastv = 10 - i;
        knt = (ii - i) + 9;
        while ((lastv > 0) && (b_A[knt] == 0.0F)) {
          lastv--;
          knt--;
        }

        lastc = 9 - i;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          knt = (ii + (lastc - 1) * 10) + 10;
          ia = knt;
          do {
            exitg1 = 0L;
            if (ia + 1 <= knt + lastv) {
              if (b_A[ia] != 0.0F) {
                exitg1 = 1L;
              } else {
                ia++;
              }
            } else {
              lastc--;
              exitg1 = 2L;
            }
          } while (exitg1 == 0L);

          if (exitg1 == 1L) {
            exitg2 = true;
          }
        }
      } else {
        lastv = 0;
        lastc = 0;
      }

      if (lastv > 0) {
        knt = ii + 11;
        if (lastc != 0) {
          if (0 <= lastc - 1) {
            memset(&work[0], 0, lastc * sizeof(float));
          }

          ix0 = 0;
          b_i = (ii + 10 * (lastc - 1)) + 11;
          for (k = knt; k <= b_i; k += 10) {
            ix = ii;
            c = 0.0F;
            i1 = (k + lastv) - 1;
            for (ia = k; ia <= i1; ia++) {
              c += b_A[ia - 1] * b_A[ix];
              ix++;
            }

            work[ix0] += c;
            ix0++;
          }
        }

        xgerc(lastv, lastc, -tau[i], ii + 1, work, b_A, ii + 11);
      }

      b_A[ii] = atmp;
    } else {
      tau[9] = 0.0F;
    }
  }

  itau = 9;
  for (ix0 = 0; ix0 < 10; ix0++) {
    for (i = 0; i <= ix0; i++) {
      knt = i + 10 * ix0;
      R[knt] = b_A[knt];
    }

    b_i = ix0 + 2;
    if (b_i <= 10) {
      memset(&R[(ix0 * 10 + b_i) + -1], 0, (11 - b_i) * sizeof(float));
    }

    work[ix0] = 0.0F;
  }

  for (i = 9; i >= 0; i--) {
    ii = (i + i * 10) + 11;
    if (i + 1 < 10) {
      b_A[ii - 11] = 1.0F;
      if (tau[itau] != 0.0F) {
        lastv = 10 - i;
        knt = ii - i;
        while ((lastv > 0) && (b_A[knt - 2] == 0.0F)) {
          lastv--;
          knt--;
        }

        lastc = 9 - i;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          knt = ii + (lastc - 1) * 10;
          ia = knt;
          do {
            exitg1 = 0L;
            if (ia <= (knt + lastv) - 1) {
              if (b_A[ia - 1] != 0.0F) {
                exitg1 = 1L;
              } else {
                ia++;
              }
            } else {
              lastc--;
              exitg1 = 2L;
            }
          } while (exitg1 == 0L);

          if (exitg1 == 1L) {
            exitg2 = true;
          }
        }
      } else {
        lastv = 0;
        lastc = 0;
      }

      if (lastv > 0) {
        if (lastc != 0) {
          if (0 <= lastc - 1) {
            memset(&work[0], 0, lastc * sizeof(float));
          }

          ix0 = 0;
          b_i = ii + 10 * (lastc - 1);
          for (k = ii; k <= b_i; k += 10) {
            ix = ii;
            c = 0.0F;
            i1 = (k + lastv) - 1;
            for (ia = k; ia <= i1; ia++) {
              c += b_A[ia - 1] * b_A[ix - 11];
              ix++;
            }

            work[ix0] += c;
            ix0++;
          }
        }

        xgerc(lastv, lastc, -tau[itau], ii - 10, work, b_A, ii);
      }

      ix0 = ii - 9;
      b_i = (ii - i) - 1;
      for (k = ix0; k <= b_i; k++) {
        b_A[k - 1] *= -tau[itau];
      }
    }

    b_A[ii - 11] = 1.0F - tau[itau];
    for (ix0 = 0; ix0 < i; ix0++) {
      b_A[(ii - ix0) - 12] = 0.0F;
    }

    itau--;
  }

  for (ix0 = 0; ix0 < 10; ix0++) {
    for (i = 0; i < 10; i++) {
      knt = i + 10 * ix0;
      Q[knt] = b_A[knt];
    }
  }
}

/*
 * File trailer for qr.c
 *
 * [EOF]
 */
