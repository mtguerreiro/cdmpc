/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: qpm2.c
 *
 * MATLAB Coder version            : 5.0
 * C/C++ source code generated on  : 13-Sep-2021 17:45:18
 */

/* Include Files */
#include "qpm2.h"
#include "qr.h"
#include <math.h>
#include <string.h>

/* Function Definitions */

/*
 * Arguments    : const float f[10]
 *                const float b[4]
 *                float x[10]
 *                float *status
 * Return Type  : void
 */
void qpm2(const float f[10], const float b[4], float x[10], float *status)
{
  float b_b[4];
  int iA1[4];
  float lam[4];
  float rMin;
  bool cTolComputed;
  float cTol[4];
  int iC[4];
  int nA;
  float Xnorm0;
  float scale;
  int i;
  float r[10];
  bool exitg2;
  float absxk;
  int RLinv_tmp;
  static const float Hinv[100] = { 0.0172031429F, -0.0259538591F, 0.00240378059F,
    0.00716733327F, 0.00268285442F, -0.00115655281F, -0.00189037377F,
    -0.000928911497F, 5.37549495E-5F, 0.000420423399F, -0.0259538591F,
    0.0563482977F, -0.0295797382F, -0.00838870183F, 0.00316519383F,
    0.00445540575F, 0.00162811589F, -0.000660337857F, -0.00106465735F,
    4.54776637E-5F, 0.00240378059F, -0.0295797382F, 0.0566798709F, -0.028573161F,
    -0.00800702348F, 0.00300231576F, 0.00418721745F, 0.00147864874F,
    -0.000659043668F, -0.000937960285F, 0.00716733327F, -0.00838870183F,
    -0.028573161F, 0.0596176162F, -0.0275616366F, -0.008544554F, 0.0023613465F,
    0.00418354059F, 0.00162381551F, -0.00188117812F, 0.00268285442F,
    0.00316519383F, -0.00800702348F, -0.0275616366F, 0.0598375127F,
    -0.027858723F, -0.00854671653F, 0.00298701413F, 0.00443615578F,
    -0.0011118853F, -0.00115655281F, 0.00445540575F, 0.00300231576F,
    -0.008544554F, -0.027858723F, 0.0598409846F, -0.0275577828F, -0.00802704599F,
    0.00312671857F, 0.00275294809F, -0.00189037377F, 0.00162811589F,
    0.00418721745F, 0.0023613465F, -0.00854671653F, -0.0275577828F,
    0.0596211962F, -0.0285849012F, -0.00839225575F, 0.0071820938F,
    -0.000928911497F, -0.000660337857F, 0.00147864874F, 0.00418354059F,
    0.00298701413F, -0.00802704599F, -0.0285849012F, 0.0567295365F, -0.02951831F,
    0.00228004134F, 5.37549495E-5F, -0.00106465735F, -0.000659043668F,
    0.00162381551F, 0.00443615578F, 0.00312671857F, -0.00839225575F,
    -0.02951831F, 0.0564342476F, -0.0261241514F, 0.000420423399F, 4.54776637E-5F,
    -0.000937960285F, -0.00188117812F, -0.0011118853F, 0.00275294809F,
    0.0071820938F, 0.00228004134F, -0.0261241514F, 0.0175386332F };

  int kNext;
  float t;
  int k;
  float Xnorm;
  float z[10];
  static const float Ac[40] = { 1.0F, -1.0F, 24.4311924F, -24.4311924F, 0.0F,
    -0.0F, -0.0F, -0.0F, 0.0F, -0.0F, -0.0F, -0.0F, 0.0F, -0.0F, -0.0F, -0.0F,
    0.0F, -0.0F, -0.0F, -0.0F, 0.0F, -0.0F, -0.0F, -0.0F, 0.0F, -0.0F, -0.0F,
    -0.0F, 0.0F, -0.0F, -0.0F, -0.0F, 0.0F, -0.0F, -0.0F, -0.0F, 0.0F, -0.0F,
    -0.0F, -0.0F };

  long exitg1;
  bool guard1 = false;
  int Stat;
  float RLinv[100];
  int iSave;
  float QQ[100];
  float RR[100];
  int kDrop;
  long exitg3;
  bool isT1Inf;
  bool tempOK;
  static const float Linv[100] = { 0.0082094F, -0.0878326F, 0.0803394F,
    0.017533F, -0.0333682F, -0.0354325F, -0.0159355F, 0.0002282F, 0.005137F,
    0.0031746F, 0.0F, 0.0889459F, -0.1994441F, 0.0518489F, 0.070822F, 0.0225759F,
    -0.0111123F, -0.0162763F, -0.0075313F, 0.0003434F, 0.0F, 0.0F, 0.120322F,
    -0.189681098F, 0.0392104F, 0.0621935F, 0.0200985F, -0.0110763F, -0.0155334F,
    -0.0070825F, 0.0F, 0.0F, 0.0F, 0.121633403F, -0.201697797F, 0.0206892F,
    0.0549891F, 0.0202011F, -0.0089011F, -0.0142047F, 0.0F, 0.0F, 0.0F, 0.0F,
    0.126233697F, -0.200090706F, 0.0181014F, 0.0550187F, 0.0210016F, -0.0083958F,
    0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.131219104F, -0.1971872F, 0.0180706F,
    0.0545992F, 0.0207874F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.132221103F,
    -0.197216406F, 0.017418F, 0.0542317F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
    0.0F, 0.132247F, -0.197342F, 0.0172165F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
    0.0F, 0.0F, 0.132369906F, -0.197262406F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
    0.0F, 0.0F, 0.0F, 0.132433504F };

  int j;
  bool exitg4;
  float TL[100];
  float H[100];
  int qY;
  float D[100];
  *status = 1.0F;
  b_b[0] = -b[0];
  iA1[0] = 0;
  lam[0] = 0.0F;
  b_b[1] = -b[1];
  iA1[1] = 0;
  lam[1] = 0.0F;
  b_b[2] = -b[2];
  iA1[2] = 0;
  lam[2] = 0.0F;
  b_b[3] = -b[3];
  iA1[3] = 0;
  lam[3] = 0.0F;
  rMin = 0.0F;
  cTolComputed = false;
  cTol[0] = 1.0F;
  iC[0] = 0;
  cTol[1] = 1.0F;
  iC[1] = 0;
  cTol[2] = 1.0F;
  iC[2] = 0;
  cTol[3] = 1.0F;
  iC[3] = 0;
  nA = 0;
  Xnorm0 = 0.0F;
  scale = 1.29246971E-26F;
  for (i = 0; i < 10; i++) {
    r[i] = 0.0F;
    absxk = 0.0F;
    for (RLinv_tmp = 0; RLinv_tmp < 10; RLinv_tmp++) {
      absxk += -Hinv[i + 10 * RLinv_tmp] * f[RLinv_tmp];
    }

    x[i] = absxk;
    absxk = fabsf(absxk);
    if (absxk > scale) {
      t = scale / absxk;
      Xnorm0 = Xnorm0 * t * t + 1.0F;
      scale = absxk;
    } else {
      t = absxk / scale;
      Xnorm0 += t * t;
    }
  }

  Xnorm0 = scale * sqrtf(Xnorm0);
  exitg2 = false;
  while ((!exitg2) && (*status <= 5.0F)) {
    scale = -1.0E-6F;
    kNext = -1;
    for (i = 0; i < 4; i++) {
      if (!cTolComputed) {
        for (k = 0; k < 10; k++) {
          z[k] = fabsf(Ac[i + (k << 2U)] * x[k]);
        }

        absxk = z[0];
        for (k = 0; k < 9; k++) {
          Xnorm = z[k + 1];
          if (absxk < Xnorm) {
            absxk = Xnorm;
          }
        }

        cTol[i] = fmaxf(cTol[i], absxk);
      }

      if (iA1[i] == 0) {
        Xnorm = 0.0F;
        for (RLinv_tmp = 0; RLinv_tmp < 10; RLinv_tmp++) {
          Xnorm += Ac[i + (RLinv_tmp << 2U)] * x[RLinv_tmp];
        }

        absxk = (Xnorm - b_b[i]) / cTol[i];
        if (absxk < scale) {
          scale = absxk;
          kNext = i;
        }
      }
    }

    cTolComputed = true;
    if (kNext + 1 <= 0) {
      exitg2 = true;
    } else if (*status == 5.0F) {
      *status = 0.0F;
      exitg2 = true;
    } else {
      do {
        exitg1 = 0L;
        if ((kNext + 1 > 0) && (*status <= 5.0F)) {
          guard1 = false;
          if (nA == 0) {
            for (RLinv_tmp = 0; RLinv_tmp < 10; RLinv_tmp++) {
              Xnorm = 0.0F;
              for (iSave = 0; iSave < 10; iSave++) {
                Xnorm += Hinv[RLinv_tmp + 10 * iSave] * Ac[kNext + (iSave << 2U)];
              }

              z[RLinv_tmp] = Xnorm;
            }

            guard1 = true;
          } else {
            Stat = 1;
            memset(&RLinv[0], 0, 100U * sizeof(float));
            for (i = 0; i < nA; i++) {
              for (RLinv_tmp = 0; RLinv_tmp < 10; RLinv_tmp++) {
                Xnorm = 0.0F;
                for (iSave = 0; iSave < 10; iSave++) {
                  Xnorm += Linv[RLinv_tmp + 10 * iSave] * Ac[(iC[i] + (iSave <<
                    2U)) - 1];
                }

                RLinv[RLinv_tmp + 10 * i] = Xnorm;
              }
            }

            qr(RLinv, QQ, RR);
            i = 0;
            do {
              exitg3 = 0L;
              if (i <= nA - 1) {
                if (fabsf(RR[i + 10 * i]) < 1.0E-12F) {
                  Stat = -2;
                  exitg3 = 1L;
                } else {
                  i++;
                }
              } else {
                for (i = 0; i < 10; i++) {
                  for (j = 0; j < 10; j++) {
                    absxk = 0.0F;
                    for (RLinv_tmp = 0; RLinv_tmp < 10; RLinv_tmp++) {
                      absxk += Linv[RLinv_tmp + 10 * i] * QQ[RLinv_tmp + 10 * j];
                    }

                    TL[i + 10 * j] = absxk;
                  }
                }

                memset(&RLinv[0], 0, 100U * sizeof(float));
                for (j = nA; j >= 1; j--) {
                  iSave = 10 * (j - 1);
                  kDrop = (j + iSave) - 1;
                  RLinv[kDrop] = 1.0F;
                  for (k = j; k <= nA; k++) {
                    RLinv_tmp = (j + 10 * (k - 1)) - 1;
                    RLinv[RLinv_tmp] /= RR[kDrop];
                  }

                  if (j > 1) {
                    for (i = 0; i <= j - 2; i++) {
                      for (k = j; k <= nA; k++) {
                        kDrop = 10 * (k - 1);
                        RLinv_tmp = i + kDrop;
                        RLinv[RLinv_tmp] -= RR[i + iSave] * RLinv[(j + kDrop) -
                          1];
                      }
                    }
                  }
                }

                for (i = 0; i < 10; i++) {
                  if (i + 1 <= 10) {
                    if (nA > 32766) {
                      qY = MAX_int16_T;
                    } else {
                      qY = nA + 1;
                    }
                  }

                  for (j = i + 1; j < 11; j++) {
                    iSave = i + 10 * (j - 1);
                    H[iSave] = 0.0F;
                    for (k = qY; k < 11; k++) {
                      kDrop = 10 * (k - 1);
                      H[iSave] -= TL[i + kDrop] * TL[(j + kDrop) - 1];
                    }

                    H[(j + 10 * i) - 1] = H[iSave];
                  }
                }

                for (j = 0; j < nA; j++) {
                  for (i = 0; i < 10; i++) {
                    iSave = i + 10 * j;
                    D[iSave] = 0.0F;
                    for (k = j + 1; k <= nA; k++) {
                      kDrop = 10 * (k - 1);
                      D[iSave] += TL[i + kDrop] * RLinv[j + kDrop];
                    }
                  }
                }

                exitg3 = 1L;
              }
            } while (exitg3 == 0L);

            if (Stat <= 0) {
              *status = -2.0F;
              exitg1 = 1L;
            } else {
              for (RLinv_tmp = 0; RLinv_tmp < 100; RLinv_tmp++) {
                TL[RLinv_tmp] = -H[RLinv_tmp];
              }

              for (RLinv_tmp = 0; RLinv_tmp < 10; RLinv_tmp++) {
                Xnorm = 0.0F;
                for (iSave = 0; iSave < 10; iSave++) {
                  Xnorm += TL[RLinv_tmp + 10 * iSave] * Ac[kNext + (iSave << 2U)];
                }

                z[RLinv_tmp] = Xnorm;
              }

              for (i = 0; i < nA; i++) {
                Xnorm = 0.0F;
                for (RLinv_tmp = 0; RLinv_tmp < 10; RLinv_tmp++) {
                  Xnorm += Ac[kNext + (RLinv_tmp << 2U)] * D[RLinv_tmp + 10 * i];
                }

                r[i] = Xnorm;
              }

              guard1 = true;
            }
          }

          if (guard1) {
            kDrop = 0;
            scale = 0.0F;
            isT1Inf = true;
            tempOK = true;
            if (nA > 0) {
              iSave = 0;
              exitg4 = false;
              while ((!exitg4) && (iSave <= nA - 1)) {
                if (r[iSave] >= 1.0E-12F) {
                  tempOK = false;
                  exitg4 = true;
                } else {
                  iSave++;
                }
              }
            }

            if ((nA != 0) && (!tempOK)) {
              for (i = 0; i < nA; i++) {
                if (r[i] > 1.0E-12F) {
                  absxk = lam[iC[i] - 1] / r[i];
                  if ((kDrop == 0) || (absxk < rMin)) {
                    rMin = absxk;
                    kDrop = i + 1;
                  }
                }
              }

              if (kDrop > 0) {
                scale = rMin;
                isT1Inf = false;
              }
            }

            absxk = 0.0F;
            for (k = 0; k < 10; k++) {
              absxk += z[k] * Ac[kNext + (k << 2U)];
            }

            if (absxk <= 0.0F) {
              absxk = 0.0F;
              tempOK = true;
            } else {
              Xnorm = 0.0F;
              for (RLinv_tmp = 0; RLinv_tmp < 10; RLinv_tmp++) {
                Xnorm += Ac[kNext + (RLinv_tmp << 2U)] * x[RLinv_tmp];
              }

              absxk = (b_b[kNext] - Xnorm) / absxk;
              tempOK = false;
            }

            if (isT1Inf && tempOK) {
              *status = -1.0F;
              exitg1 = 1L;
            } else {
              if (tempOK) {
                t = scale;
              } else if (isT1Inf) {
                t = absxk;
              } else {
                t = fminf(scale, absxk);
              }

              for (i = 0; i < nA; i++) {
                lam[iC[i] - 1] -= t * r[i];
                if (lam[iC[i] - 1] < 0.0F) {
                  lam[iC[i] - 1] = 0.0F;
                }
              }

              lam[kNext] += t;
              if (t == scale) {
                iA1[iC[kDrop - 1] - 1] = 0;
                if (kDrop < nA) {
                  RLinv_tmp = nA - 1;
                  for (i = kDrop; i <= RLinv_tmp; i++) {
                    iC[i - 1] = iC[i];
                  }
                }

                iC[nA - 1] = 0;
                nA--;
              }

              if (!tempOK) {
                for (RLinv_tmp = 0; RLinv_tmp < 10; RLinv_tmp++) {
                  x[RLinv_tmp] += t * z[RLinv_tmp];
                }

                if (t == absxk) {
                  if (nA == 10) {
                    *status = -1.0F;
                    exitg1 = 1L;
                  } else {
                    if (nA > 32766) {
                      nA = MAX_int16_T;
                    } else {
                      nA++;
                    }

                    iC[nA - 1] = kNext + 1;
                    i = nA - 1;
                    exitg4 = false;
                    while ((!exitg4) && (i + 1 > 1)) {
                      RLinv_tmp = iC[i - 1];
                      if (iC[i] > RLinv_tmp) {
                        exitg4 = true;
                      } else {
                        iSave = iC[i];
                        iC[i] = RLinv_tmp;
                        iC[i - 1] = iSave;
                        i--;
                      }
                    }

                    iA1[kNext] = 1;
                    kNext = -1;
                    (*status)++;
                  }
                } else {
                  (*status)++;
                }
              } else {
                (*status)++;
              }
            }
          }
        } else {
          Xnorm = 0.0F;
          scale = 1.29246971E-26F;
          for (k = 0; k < 10; k++) {
            absxk = fabsf(x[k]);
            if (absxk > scale) {
              t = scale / absxk;
              Xnorm = Xnorm * t * t + 1.0F;
              scale = absxk;
            } else {
              t = absxk / scale;
              Xnorm += t * t;
            }
          }

          Xnorm = scale * sqrtf(Xnorm);
          if (fabsf(Xnorm - Xnorm0) > 0.001F) {
            Xnorm0 = Xnorm;
            cTol[0] = fmaxf(fabsf(b_b[0]), 1.0F);
            cTol[1] = fmaxf(fabsf(b_b[1]), 1.0F);
            cTol[2] = fmaxf(fabsf(b_b[2]), 1.0F);
            cTol[3] = fmaxf(fabsf(b_b[3]), 1.0F);
            cTolComputed = false;
          }

          exitg1 = 2L;
        }
      } while (exitg1 == 0L);

      if (exitg1 == 1L) {
        exitg2 = true;
      }
    }
  }
}

/*
 * File trailer for qpm2.c
 *
 * [EOF]
 */
