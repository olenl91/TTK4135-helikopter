/*
 * helikopter04_3.c
 *
 * Real-Time Workshop code generation for Simulink model "helikopter04_3.mdl".
 *
 * Model version              : 1.73
 * Real-Time Workshop version : 7.5  (R2010a)  25-Jan-2010
 * C source code generated on : Fri Apr 17 11:28:31 2015
 *
 * Target selection: quarc_windows.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "helikopter04_3.h"
#include "helikopter04_3_private.h"
#include <stdio.h>
#include "helikopter04_3_dt.h"

/* Block signals (auto storage) */
BlockIO_helikopter04_3 helikopter04_3_B;

/* Continuous states */
ContinuousStates_helikopter04_3 helikopter04_3_X;

/* Block states (auto storage) */
D_Work_helikopter04_3 helikopter04_3_DWork;

/* Real-time model */
RT_MODEL_helikopter04_3 helikopter04_3_M_;
RT_MODEL_helikopter04_3 *helikopter04_3_M = &helikopter04_3_M_;

/*
 * Writes out MAT-file header.  Returns success or failure.
 * Returns:
 *      0 - success
 *      1 - failure
 */
int_T rt_WriteMat4FileHeader(FILE *fp, int32_T m, int32_T n, const char *name)
{
  typedef enum { ELITTLE_ENDIAN, EBIG_ENDIAN } ByteOrder;

  int16_T one = 1;
  ByteOrder byteOrder = (*((int8_T *)&one)==1) ? ELITTLE_ENDIAN : EBIG_ENDIAN;
  int32_T type = (byteOrder == ELITTLE_ENDIAN) ? 0: 1000;
  int32_T imagf = 0;
  int32_T name_len = strlen(name) + 1;
  if ((fwrite(&type, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&m, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&n, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&imagf, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&name_len, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(name, sizeof(char), name_len, fp) == 0)) {
    return(1);
  } else {
    return(0);
  }
}                                      /* end rt_WriteMat4FileHeader */

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 5;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  helikopter04_3_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helikopter04_3_output(int_T tid)
{
  /* local block i/o variables */
  real_T rtb_HILReadEncoder_o1;
  real_T rtb_HILReadEncoder_o2;
  real_T rtb_HILReadEncoder_o3;
  real_T rtb_FromWorkspace[2];
  real_T rtb_Saturation_k;
  real_T rtb_Gain[6];
  real_T rtb_Gain_c[2];
  real_T rtb_Sum;
  real_T rtb_Gain1;
  real_T tmp[6];
  int32_T tmp_0;
  real_T rtb_Gain_0[6];
  int32_T tmp_1;
  if (rtmIsMajorTimeStep(helikopter04_3_M)) {
    /* set solver stop time */
    if (!(helikopter04_3_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helikopter04_3_M->solverInfo,
                            ((helikopter04_3_M->Timing.clockTickH0 + 1) *
        helikopter04_3_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helikopter04_3_M->solverInfo,
                            ((helikopter04_3_M->Timing.clockTick0 + 1) *
        helikopter04_3_M->Timing.stepSize0 +
        helikopter04_3_M->Timing.clockTickH0 *
        helikopter04_3_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helikopter04_3_M)) {
    helikopter04_3_M->Timing.t[0] = rtsiGetT(&helikopter04_3_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helikopter04_3_M)) {
    /* S-Function (hil_read_encoder_block): '<S2>/HIL Read Encoder' */

    /* S-Function Block: helikopter04_3/Heli 3D/HIL Read Encoder (hil_read_encoder_block) */
    {
      t_error result = hil_read_encoder(helikopter04_3_DWork.HILInitialize_Card,
        helikopter04_3_P.HILReadEncoder_Channels, 3,
        &helikopter04_3_DWork.HILReadEncoder_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter04_3_M, _rt_error_message);
      } else {
        rtb_HILReadEncoder_o1 = helikopter04_3_DWork.HILReadEncoder_Buffer[0];
        rtb_HILReadEncoder_o2 = helikopter04_3_DWork.HILReadEncoder_Buffer[1];
        rtb_HILReadEncoder_o3 = helikopter04_3_DWork.HILReadEncoder_Buffer[2];
      }
    }

    /* Gain: '<S2>/Kalibrer-Elev' */
    helikopter04_3_B.KalibrerElev = helikopter04_3_P.KalibrerElev_Gain *
      rtb_HILReadEncoder_o3;

    /* Sum: '<Root>/Add' incorporates:
     *  Constant: '<Root>/Constant'
     */
    helikopter04_3_B.Add = helikopter04_3_B.KalibrerElev +
      helikopter04_3_P.Constant_Value;
  }

  /* FromWorkspace: '<Root>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *)
      helikopter04_3_DWork.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopter04_3_DWork.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex = helikopter04_3_DWork.FromWorkspace1_IWORK.PrevIndex;
    real_T t = helikopter04_3_M->Timing.t[0];

    /* get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[79]) {
      currTimeIndex = 78;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helikopter04_3_DWork.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T i1;
            real_T *y0 = helikopter04_3_B.FromWorkspace1;
            for (i1=0; i1 < 6; i1++) {
              y0[i1] = pDataValues[currTimeIndex];
              pDataValues += 80;
            }
          }
        } else {
          {
            int_T i1;
            real_T *y0 = helikopter04_3_B.FromWorkspace1;
            for (i1=0; i1 < 6; i1++) {
              y0[i1] = pDataValues[currTimeIndex + 1];
              pDataValues += 80;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T i1;
          real_T *y0 = helikopter04_3_B.FromWorkspace1;
          for (i1=0; i1 < 6; i1++) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            y0[i1] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 80;
          }
        }
      }
    }
  }

  if (rtmIsMajorTimeStep(helikopter04_3_M)) {
    /* Gain: '<S2>/Kalibrer-Pitch' */
    helikopter04_3_B.KalibrerPitch = helikopter04_3_P.KalibrerPitch_Gain *
      rtb_HILReadEncoder_o2;
  }

  /* TransferFcn: '<S2>/Vandring Lavpass' */
  helikopter04_3_B.VandringLavpass = helikopter04_3_P.VandringLavpass_C*
    helikopter04_3_X.VandringLavpass_CSTATE;

  /* Sum: '<Root>/Add1' incorporates:
   *  Constant: '<Root>/Constant1'
   */
  helikopter04_3_B.Add1 = helikopter04_3_P.Constant1_Value +
    helikopter04_3_B.VandringLavpass;
  if (rtmIsMajorTimeStep(helikopter04_3_M)) {
    /* ToFile: '<Root>/To File' */
    if (rtmIsMajorTimeStep(helikopter04_3_M)) {
      if (!(++helikopter04_3_DWork.ToFile_IWORK.Decimation % 1) &&
          (helikopter04_3_DWork.ToFile_IWORK.Count*2)+1 < 100000000 ) {
        FILE *fp = (FILE *) helikopter04_3_DWork.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          helikopter04_3_DWork.ToFile_IWORK.Decimation = 0;
          u[0] = helikopter04_3_M->Timing.t[1];
          u[1] = helikopter04_3_B.Add1;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(helikopter04_3_M,
                              "Error writing to MAT-file travel.mat");
            return;
          }

          if (((++helikopter04_3_DWork.ToFile_IWORK.Count)*2)+1 >= 100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file travel.mat.\n");
          }
        }
      }
    }

    /* Gain: '<S2>/Kalibrer -Vandring' */
    helikopter04_3_B.KalibrerVandring = helikopter04_3_P.KalibrerVandring_Gain *
      rtb_HILReadEncoder_o1;
  }

  /* TransferFcn: '<S2>/Vandring Deriv' */
  helikopter04_3_B.VandringDeriv = helikopter04_3_P.VandringDeriv_D*
    helikopter04_3_B.KalibrerVandring;
  helikopter04_3_B.VandringDeriv += helikopter04_3_P.VandringDeriv_C*
    helikopter04_3_X.VandringDeriv_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn4' */
  helikopter04_3_B.TransferFcn4 = helikopter04_3_P.TransferFcn4_D*
    helikopter04_3_B.KalibrerPitch;
  helikopter04_3_B.TransferFcn4 += helikopter04_3_P.TransferFcn4_C*
    helikopter04_3_X.TransferFcn4_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn5' */
  helikopter04_3_B.TransferFcn5 = helikopter04_3_P.TransferFcn5_D*
    helikopter04_3_B.KalibrerElev;
  helikopter04_3_B.TransferFcn5 += helikopter04_3_P.TransferFcn5_C*
    helikopter04_3_X.TransferFcn5_CSTATE;

  /* Gain: '<Root>/Gain' incorporates:
   *  SignalConversion: '<Root>/TmpSignal ConversionAtGainInport1'
   */
  tmp[0] = helikopter04_3_B.Add1;
  tmp[1] = helikopter04_3_B.VandringDeriv;
  tmp[2] = helikopter04_3_B.KalibrerPitch;
  tmp[3] = helikopter04_3_B.TransferFcn4;
  tmp[4] = helikopter04_3_B.Add;
  tmp[5] = helikopter04_3_B.TransferFcn5;

  /* Gain: '<S3>/Gain' incorporates:
   *  Sum: '<S3>/Sum1'
   */
  for (tmp_0 = 0; tmp_0 < 6; tmp_0++) {
    rtb_Gain[tmp_0] = 0.0;
    for (tmp_1 = 0; tmp_1 < 6; tmp_1++) {
      rtb_Gain[tmp_0] += helikopter04_3_P.Gain_Gain[6 * tmp_1 + tmp_0] *
        tmp[tmp_1];
    }

    rtb_Gain_0[tmp_0] = rtb_Gain[tmp_0] - helikopter04_3_B.FromWorkspace1[tmp_0];
  }

  for (tmp_0 = 0; tmp_0 < 2; tmp_0++) {
    rtb_Gain_c[tmp_0] = 0.0;
    for (tmp_1 = 0; tmp_1 < 6; tmp_1++) {
      rtb_Gain_c[tmp_0] += helikopter04_3_P.Gain_Gain_h[(tmp_1 << 1) + tmp_0] *
        rtb_Gain_0[tmp_1];
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helikopter04_3_DWork.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopter04_3_DWork.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helikopter04_3_DWork.FromWorkspace_IWORK.PrevIndex;
    real_T t = helikopter04_3_M->Timing.t[0];

    /* get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[79]) {
      currTimeIndex = 78;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helikopter04_3_DWork.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_FromWorkspace[0] = pDataValues[currTimeIndex];
          pDataValues += 80;
          rtb_FromWorkspace[1] = pDataValues[currTimeIndex];
          pDataValues += 80;
        } else {
          rtb_FromWorkspace[0] = pDataValues[currTimeIndex + 1];
          pDataValues += 80;
          rtb_FromWorkspace[1] = pDataValues[currTimeIndex + 1];
          pDataValues += 80;
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_FromWorkspace[0] = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 80;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_FromWorkspace[1] = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 80;
      }
    }
  }

  /* Sum: '<S3>/Sum2' */
  helikopter04_3_B.Sum2[0] = rtb_FromWorkspace[0] - rtb_Gain_c[0];
  helikopter04_3_B.Sum2[1] = rtb_FromWorkspace[1] - rtb_Gain_c[1];
  if (rtmIsMajorTimeStep(helikopter04_3_M)) {
  }

  /* Integrator: '<S1>/Integrator'
   *
   * Regarding '<S1>/Integrator':
   *  Limited Integrator
   */
  if (helikopter04_3_X.Integrator_CSTATE >= helikopter04_3_P.Integrator_UpperSat
      ) {
    helikopter04_3_X.Integrator_CSTATE = helikopter04_3_P.Integrator_UpperSat;
  } else if (helikopter04_3_X.Integrator_CSTATE <=
             helikopter04_3_P.Integrator_LowerSat ) {
    helikopter04_3_X.Integrator_CSTATE = helikopter04_3_P.Integrator_LowerSat;
  }

  rtb_Saturation_k = helikopter04_3_X.Integrator_CSTATE;

  /* Sum: '<S1>/Sum' */
  rtb_Sum = helikopter04_3_B.Sum2[1] - rtb_Gain[4];

  /* Gain: '<S1>/K_ei' */
  helikopter04_3_B.K_ei = helikopter04_3_P.K_ei_Gain * rtb_Sum;

  /* Sum: '<S1>/Sum1' incorporates:
   *  Gain: '<S1>/K_ed'
   *  Gain: '<S1>/K_ep'
   */
  rtb_Saturation_k = (helikopter04_3_P.K_ep_Gain * rtb_Sum + rtb_Saturation_k) +
    helikopter04_3_P.K_ed_Gain * rtb_Gain[5];

  /* Saturate: '<S1>/Saturation' */
  rtb_Saturation_k = rt_SATURATE(rtb_Saturation_k,
    helikopter04_3_P.Saturation_LowerSat, helikopter04_3_P.Saturation_UpperSat);
  if (rtmIsMajorTimeStep(helikopter04_3_M)) {
  }

  /* Sum: '<S4>/Sum' incorporates:
   *  Gain: '<S4>/K_pd'
   *  Gain: '<S4>/K_pp'
   *  Saturate: '<S4>/Saturation'
   *  Sum: '<S4>/Sum1'
   */
  rtb_Sum = (rt_SATURATE(helikopter04_3_B.Sum2[0],
              helikopter04_3_P.Saturation_LowerSat_i,
              helikopter04_3_P.Saturation_UpperSat_k) - rtb_Gain[2]) *
    helikopter04_3_P.K_pp_Gain - helikopter04_3_P.K_pd_Gain * rtb_Gain[3];

  /* Gain: '<S5>/Gain2' incorporates:
   *  Sum: '<S5>/Sum4'
   */
  rtb_Gain1 = (rtb_Saturation_k - rtb_Sum) * helikopter04_3_P.Gain2_Gain;

  /* Saturate: '<S2>/Sat B' */
  helikopter04_3_B.SatB = rt_SATURATE(rtb_Gain1, helikopter04_3_P.SatB_LowerSat,
    helikopter04_3_P.SatB_UpperSat);
  if (rtmIsMajorTimeStep(helikopter04_3_M)) {
  }

  /* Gain: '<S5>/Gain1' incorporates:
   *  Sum: '<S5>/Sum3'
   */
  rtb_Gain1 = (rtb_Sum + rtb_Saturation_k) * helikopter04_3_P.Gain1_Gain;

  /* Saturate: '<S2>/Sat' */
  helikopter04_3_B.Sat = rt_SATURATE(rtb_Gain1, helikopter04_3_P.Sat_LowerSat,
    helikopter04_3_P.Sat_UpperSat);
  if (rtmIsMajorTimeStep(helikopter04_3_M)) {
    /* S-Function (hil_write_analog_block): '<S2>/HIL Write Analog' */

    /* S-Function Block: helikopter04_3/Heli 3D/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helikopter04_3_DWork.HILWriteAnalog_Buffer[0] = helikopter04_3_B.SatB;
      helikopter04_3_DWork.HILWriteAnalog_Buffer[1] = helikopter04_3_B.Sat;
      result = hil_write_analog(helikopter04_3_DWork.HILInitialize_Card,
        helikopter04_3_P.HILWriteAnalog_Channels, 2,
        &helikopter04_3_DWork.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter04_3_M, _rt_error_message);
      }
    }
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Model update function */
void helikopter04_3_update(int_T tid)
{
  if (rtmIsMajorTimeStep(helikopter04_3_M)) {
    rt_ertODEUpdateContinuousStates(&helikopter04_3_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++helikopter04_3_M->Timing.clockTick0)) {
    ++helikopter04_3_M->Timing.clockTickH0;
  }

  helikopter04_3_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helikopter04_3_M->solverInfo);
  if (rtmIsMajorTimeStep(helikopter04_3_M)) {
    /* Update absolute timer for sample time: [0.001s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helikopter04_3_M->Timing.clockTick1)) {
      ++helikopter04_3_M->Timing.clockTickH1;
    }

    helikopter04_3_M->Timing.t[1] = helikopter04_3_M->Timing.clockTick1 *
      helikopter04_3_M->Timing.stepSize1 + helikopter04_3_M->Timing.clockTickH1 *
      helikopter04_3_M->Timing.stepSize1 * 4294967296.0;
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Derivatives for root system: '<Root>' */
void helikopter04_3_derivatives(void)
{
  /* Derivatives for TransferFcn: '<S2>/Vandring Lavpass' */
  {
    ((StateDerivatives_helikopter04_3 *) helikopter04_3_M->ModelData.derivs)
      ->VandringLavpass_CSTATE = helikopter04_3_B.KalibrerVandring;
    ((StateDerivatives_helikopter04_3 *) helikopter04_3_M->ModelData.derivs)
      ->VandringLavpass_CSTATE += (helikopter04_3_P.VandringLavpass_A)*
      helikopter04_3_X.VandringLavpass_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Vandring Deriv' */
  {
    ((StateDerivatives_helikopter04_3 *) helikopter04_3_M->ModelData.derivs)
      ->VandringDeriv_CSTATE = helikopter04_3_B.KalibrerVandring;
    ((StateDerivatives_helikopter04_3 *) helikopter04_3_M->ModelData.derivs)
      ->VandringDeriv_CSTATE += (helikopter04_3_P.VandringDeriv_A)*
      helikopter04_3_X.VandringDeriv_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn4' */
  {
    ((StateDerivatives_helikopter04_3 *) helikopter04_3_M->ModelData.derivs)
      ->TransferFcn4_CSTATE = helikopter04_3_B.KalibrerPitch;
    ((StateDerivatives_helikopter04_3 *) helikopter04_3_M->ModelData.derivs)
      ->TransferFcn4_CSTATE += (helikopter04_3_P.TransferFcn4_A)*
      helikopter04_3_X.TransferFcn4_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn5' */
  {
    ((StateDerivatives_helikopter04_3 *) helikopter04_3_M->ModelData.derivs)
      ->TransferFcn5_CSTATE = helikopter04_3_B.KalibrerElev;
    ((StateDerivatives_helikopter04_3 *) helikopter04_3_M->ModelData.derivs)
      ->TransferFcn5_CSTATE += (helikopter04_3_P.TransferFcn5_A)*
      helikopter04_3_X.TransferFcn5_CSTATE;
  }

  /* Derivatives for Integrator: '<S1>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helikopter04_3_X.Integrator_CSTATE <=
            helikopter04_3_P.Integrator_LowerSat );
    usat = ( helikopter04_3_X.Integrator_CSTATE >=
            helikopter04_3_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helikopter04_3_B.K_ei > 0)) ||
        (usat && (helikopter04_3_B.K_ei < 0)) ) {
      ((StateDerivatives_helikopter04_3 *) helikopter04_3_M->ModelData.derivs)
        ->Integrator_CSTATE = helikopter04_3_B.K_ei;
    } else {
      /* in saturation */
      ((StateDerivatives_helikopter04_3 *) helikopter04_3_M->ModelData.derivs)
        ->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helikopter04_3_initialize(boolean_T firstTime)
{
  (void)firstTime;

  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helikopter04_3_P.Integrator_UpperSat = rtInf;
  helikopter04_3_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helikopter04_3_M, 0,
                sizeof(RT_MODEL_helikopter04_3));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helikopter04_3_M->solverInfo,
                          &helikopter04_3_M->Timing.simTimeStep);
    rtsiSetTPtr(&helikopter04_3_M->solverInfo, &rtmGetTPtr(helikopter04_3_M));
    rtsiSetStepSizePtr(&helikopter04_3_M->solverInfo,
                       &helikopter04_3_M->Timing.stepSize0);
    rtsiSetdXPtr(&helikopter04_3_M->solverInfo,
                 &helikopter04_3_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helikopter04_3_M->solverInfo,
                         &helikopter04_3_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helikopter04_3_M->solverInfo,
      &helikopter04_3_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helikopter04_3_M->solverInfo, (&rtmGetErrorStatus
      (helikopter04_3_M)));
    rtsiSetRTModelPtr(&helikopter04_3_M->solverInfo, helikopter04_3_M);
  }

  rtsiSetSimTimeStep(&helikopter04_3_M->solverInfo, MAJOR_TIME_STEP);
  helikopter04_3_M->ModelData.intgData.f[0] = helikopter04_3_M->ModelData.odeF[0];
  helikopter04_3_M->ModelData.contStates = ((real_T *) &helikopter04_3_X);
  rtsiSetSolverData(&helikopter04_3_M->solverInfo, (void *)
                    &helikopter04_3_M->ModelData.intgData);
  rtsiSetSolverName(&helikopter04_3_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helikopter04_3_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helikopter04_3_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helikopter04_3_M->Timing.sampleTimes =
      (&helikopter04_3_M->Timing.sampleTimesArray[0]);
    helikopter04_3_M->Timing.offsetTimes =
      (&helikopter04_3_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helikopter04_3_M->Timing.sampleTimes[0] = (0.0);
    helikopter04_3_M->Timing.sampleTimes[1] = (0.001);

    /* task offsets */
    helikopter04_3_M->Timing.offsetTimes[0] = (0.0);
    helikopter04_3_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helikopter04_3_M, &helikopter04_3_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helikopter04_3_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helikopter04_3_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helikopter04_3_M, -1);
  helikopter04_3_M->Timing.stepSize0 = 0.001;
  helikopter04_3_M->Timing.stepSize1 = 0.001;

  /* external mode info */
  helikopter04_3_M->Sizes.checksums[0] = (2614143017U);
  helikopter04_3_M->Sizes.checksums[1] = (996152973U);
  helikopter04_3_M->Sizes.checksums[2] = (2204923576U);
  helikopter04_3_M->Sizes.checksums[3] = (3660220029U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helikopter04_3_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helikopter04_3_M->extModeInfo,
      &helikopter04_3_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helikopter04_3_M->extModeInfo,
                        helikopter04_3_M->Sizes.checksums);
    rteiSetTPtr(helikopter04_3_M->extModeInfo, rtmGetTPtr(helikopter04_3_M));
  }

  helikopter04_3_M->solverInfoPtr = (&helikopter04_3_M->solverInfo);
  helikopter04_3_M->Timing.stepSize = (0.001);
  rtsiSetFixedStepSize(&helikopter04_3_M->solverInfo, 0.001);
  rtsiSetSolverMode(&helikopter04_3_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helikopter04_3_M->ModelData.blockIO = ((void *) &helikopter04_3_B);

  {
    int_T i;
    for (i = 0; i < 6; i++) {
      helikopter04_3_B.FromWorkspace1[i] = 0.0;
    }

    helikopter04_3_B.KalibrerElev = 0.0;
    helikopter04_3_B.Add = 0.0;
    helikopter04_3_B.KalibrerPitch = 0.0;
    helikopter04_3_B.VandringLavpass = 0.0;
    helikopter04_3_B.Add1 = 0.0;
    helikopter04_3_B.KalibrerVandring = 0.0;
    helikopter04_3_B.VandringDeriv = 0.0;
    helikopter04_3_B.TransferFcn4 = 0.0;
    helikopter04_3_B.TransferFcn5 = 0.0;
    helikopter04_3_B.Sum2[0] = 0.0;
    helikopter04_3_B.Sum2[1] = 0.0;
    helikopter04_3_B.K_ei = 0.0;
    helikopter04_3_B.SatB = 0.0;
    helikopter04_3_B.Sat = 0.0;
  }

  /* parameters */
  helikopter04_3_M->ModelData.defaultParam = ((real_T *)&helikopter04_3_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helikopter04_3_X;
    helikopter04_3_M->ModelData.contStates = (x);
    (void) memset((void *)&helikopter04_3_X, 0,
                  sizeof(ContinuousStates_helikopter04_3));
  }

  /* states (dwork) */
  helikopter04_3_M->Work.dwork = ((void *) &helikopter04_3_DWork);
  (void) memset((void *)&helikopter04_3_DWork, 0,
                sizeof(D_Work_helikopter04_3));
  helikopter04_3_DWork.HILInitialize_AIMinimums[0] = 0.0;
  helikopter04_3_DWork.HILInitialize_AIMinimums[1] = 0.0;
  helikopter04_3_DWork.HILInitialize_AIMinimums[2] = 0.0;
  helikopter04_3_DWork.HILInitialize_AIMinimums[3] = 0.0;
  helikopter04_3_DWork.HILInitialize_AIMaximums[0] = 0.0;
  helikopter04_3_DWork.HILInitialize_AIMaximums[1] = 0.0;
  helikopter04_3_DWork.HILInitialize_AIMaximums[2] = 0.0;
  helikopter04_3_DWork.HILInitialize_AIMaximums[3] = 0.0;
  helikopter04_3_DWork.HILInitialize_AOMinimums[0] = 0.0;
  helikopter04_3_DWork.HILInitialize_AOMinimums[1] = 0.0;
  helikopter04_3_DWork.HILInitialize_AOMinimums[2] = 0.0;
  helikopter04_3_DWork.HILInitialize_AOMinimums[3] = 0.0;
  helikopter04_3_DWork.HILInitialize_AOMaximums[0] = 0.0;
  helikopter04_3_DWork.HILInitialize_AOMaximums[1] = 0.0;
  helikopter04_3_DWork.HILInitialize_AOMaximums[2] = 0.0;
  helikopter04_3_DWork.HILInitialize_AOMaximums[3] = 0.0;
  helikopter04_3_DWork.HILInitialize_AOVoltages[0] = 0.0;
  helikopter04_3_DWork.HILInitialize_AOVoltages[1] = 0.0;
  helikopter04_3_DWork.HILInitialize_AOVoltages[2] = 0.0;
  helikopter04_3_DWork.HILInitialize_AOVoltages[3] = 0.0;
  helikopter04_3_DWork.HILInitialize_FilterFrequency[0] = 0.0;
  helikopter04_3_DWork.HILInitialize_FilterFrequency[1] = 0.0;
  helikopter04_3_DWork.HILInitialize_FilterFrequency[2] = 0.0;
  helikopter04_3_DWork.HILInitialize_FilterFrequency[3] = 0.0;
  helikopter04_3_DWork.HILWriteAnalog_Buffer[0] = 0.0;
  helikopter04_3_DWork.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helikopter04_3_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 15;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }
}

/* Model terminate function */
void helikopter04_3_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopter04_3/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    hil_task_stop_all(helikopter04_3_DWork.HILInitialize_Card);
    hil_task_delete_all(helikopter04_3_DWork.HILInitialize_Card);
    hil_monitor_stop_all(helikopter04_3_DWork.HILInitialize_Card);
    hil_monitor_delete_all(helikopter04_3_DWork.HILInitialize_Card);
    is_switching = false;
    if ((helikopter04_3_P.HILInitialize_AOTerminate && !is_switching) ||
        (helikopter04_3_P.HILInitialize_AOExit && is_switching)) {
      helikopter04_3_DWork.HILInitialize_AOVoltages[0] =
        helikopter04_3_P.HILInitialize_AOFinal;
      helikopter04_3_DWork.HILInitialize_AOVoltages[1] =
        helikopter04_3_P.HILInitialize_AOFinal;
      helikopter04_3_DWork.HILInitialize_AOVoltages[2] =
        helikopter04_3_P.HILInitialize_AOFinal;
      helikopter04_3_DWork.HILInitialize_AOVoltages[3] =
        helikopter04_3_P.HILInitialize_AOFinal;
      result = hil_write_analog(helikopter04_3_DWork.HILInitialize_Card,
        &helikopter04_3_P.HILInitialize_AOChannels[0], 4U,
        &helikopter04_3_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter04_3_M, _rt_error_message);
      }
    }

    hil_close(helikopter04_3_DWork.HILInitialize_Card);
    helikopter04_3_DWork.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helikopter04_3_DWork.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      const char *fileName = "travel.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopter04_3_M, "Error closing MAT-file travel.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helikopter04_3_M,
                          "Error reopening MAT-file travel.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2, helikopter04_3_DWork.ToFile_IWORK.Count,
           "travel")) {
        rtmSetErrorStatus(helikopter04_3_M,
                          "Error writing header for travel to MAT-file travel.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopter04_3_M, "Error closing MAT-file travel.mat");
        return;
      }

      helikopter04_3_DWork.ToFile_PWORK.FilePtr = (NULL);
    }
  }
}

/*========================================================================*
 * Start of GRT compatible call interface                                 *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  helikopter04_3_output(tid);
}

void MdlUpdate(int_T tid)
{
  helikopter04_3_update(tid);
}

void MdlInitializeSizes(void)
{
  helikopter04_3_M->Sizes.numContStates = (5);/* Number of continuous states */
  helikopter04_3_M->Sizes.numY = (0);  /* Number of model outputs */
  helikopter04_3_M->Sizes.numU = (0);  /* Number of model inputs */
  helikopter04_3_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helikopter04_3_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helikopter04_3_M->Sizes.numBlocks = (52);/* Number of blocks */
  helikopter04_3_M->Sizes.numBlockIO = (14);/* Number of block outputs */
  helikopter04_3_M->Sizes.numBlockPrms = (158);/* Sum of parameter "widths" */
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
  /* InitializeConditions for TransferFcn: '<S2>/Vandring Lavpass' */
  helikopter04_3_X.VandringLavpass_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Vandring Deriv' */
  helikopter04_3_X.VandringDeriv_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn4' */
  helikopter04_3_X.TransferFcn4_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn5' */
  helikopter04_3_X.TransferFcn5_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S1>/Integrator' */
  helikopter04_3_X.Integrator_CSTATE = helikopter04_3_P.Integrator_IC;
}

void MdlStart(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopter04_3/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q4", "0", &helikopter04_3_DWork.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter04_3_M, _rt_error_message);
      return;
    }

    is_switching = false;
    if ((helikopter04_3_P.HILInitialize_CKPStart && !is_switching) ||
        (helikopter04_3_P.HILInitialize_CKPEnter && is_switching)) {
      result = hil_set_clock_mode(helikopter04_3_DWork.HILInitialize_Card,
        (t_clock *) &helikopter04_3_P.HILInitialize_CKChannels[0], 2U,
        (t_clock_mode *) &helikopter04_3_P.HILInitialize_CKModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter04_3_M, _rt_error_message);
        return;
      }
    }

    result = hil_watchdog_clear(helikopter04_3_DWork.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter04_3_M, _rt_error_message);
      return;
    }

    if ((helikopter04_3_P.HILInitialize_AIPStart && !is_switching) ||
        (helikopter04_3_P.HILInitialize_AIPEnter && is_switching)) {
      helikopter04_3_DWork.HILInitialize_AIMinimums[0] =
        helikopter04_3_P.HILInitialize_AILow;
      helikopter04_3_DWork.HILInitialize_AIMinimums[1] =
        helikopter04_3_P.HILInitialize_AILow;
      helikopter04_3_DWork.HILInitialize_AIMinimums[2] =
        helikopter04_3_P.HILInitialize_AILow;
      helikopter04_3_DWork.HILInitialize_AIMinimums[3] =
        helikopter04_3_P.HILInitialize_AILow;
      helikopter04_3_DWork.HILInitialize_AIMaximums[0] =
        helikopter04_3_P.HILInitialize_AIHigh;
      helikopter04_3_DWork.HILInitialize_AIMaximums[1] =
        helikopter04_3_P.HILInitialize_AIHigh;
      helikopter04_3_DWork.HILInitialize_AIMaximums[2] =
        helikopter04_3_P.HILInitialize_AIHigh;
      helikopter04_3_DWork.HILInitialize_AIMaximums[3] =
        helikopter04_3_P.HILInitialize_AIHigh;
      result = hil_set_analog_input_ranges
        (helikopter04_3_DWork.HILInitialize_Card,
         &helikopter04_3_P.HILInitialize_AIChannels[0], 4U,
         &helikopter04_3_DWork.HILInitialize_AIMinimums[0],
         &helikopter04_3_DWork.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter04_3_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter04_3_P.HILInitialize_AOPStart && !is_switching) ||
        (helikopter04_3_P.HILInitialize_AOPEnter && is_switching)) {
      helikopter04_3_DWork.HILInitialize_AOMinimums[0] =
        helikopter04_3_P.HILInitialize_AOLow;
      helikopter04_3_DWork.HILInitialize_AOMinimums[1] =
        helikopter04_3_P.HILInitialize_AOLow;
      helikopter04_3_DWork.HILInitialize_AOMinimums[2] =
        helikopter04_3_P.HILInitialize_AOLow;
      helikopter04_3_DWork.HILInitialize_AOMinimums[3] =
        helikopter04_3_P.HILInitialize_AOLow;
      helikopter04_3_DWork.HILInitialize_AOMaximums[0] =
        helikopter04_3_P.HILInitialize_AOHigh;
      helikopter04_3_DWork.HILInitialize_AOMaximums[1] =
        helikopter04_3_P.HILInitialize_AOHigh;
      helikopter04_3_DWork.HILInitialize_AOMaximums[2] =
        helikopter04_3_P.HILInitialize_AOHigh;
      helikopter04_3_DWork.HILInitialize_AOMaximums[3] =
        helikopter04_3_P.HILInitialize_AOHigh;
      result = hil_set_analog_output_ranges
        (helikopter04_3_DWork.HILInitialize_Card,
         &helikopter04_3_P.HILInitialize_AOChannels[0], 4U,
         &helikopter04_3_DWork.HILInitialize_AOMinimums[0],
         &helikopter04_3_DWork.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter04_3_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter04_3_P.HILInitialize_AOStart && !is_switching) ||
        (helikopter04_3_P.HILInitialize_AOEnter && is_switching)) {
      helikopter04_3_DWork.HILInitialize_AOVoltages[0] =
        helikopter04_3_P.HILInitialize_AOInitial;
      helikopter04_3_DWork.HILInitialize_AOVoltages[1] =
        helikopter04_3_P.HILInitialize_AOInitial;
      helikopter04_3_DWork.HILInitialize_AOVoltages[2] =
        helikopter04_3_P.HILInitialize_AOInitial;
      helikopter04_3_DWork.HILInitialize_AOVoltages[3] =
        helikopter04_3_P.HILInitialize_AOInitial;
      result = hil_write_analog(helikopter04_3_DWork.HILInitialize_Card,
        &helikopter04_3_P.HILInitialize_AOChannels[0], 4U,
        &helikopter04_3_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter04_3_M, _rt_error_message);
        return;
      }
    }

    if (helikopter04_3_P.HILInitialize_AOReset) {
      helikopter04_3_DWork.HILInitialize_AOVoltages[0] =
        helikopter04_3_P.HILInitialize_AOWatchdog;
      helikopter04_3_DWork.HILInitialize_AOVoltages[1] =
        helikopter04_3_P.HILInitialize_AOWatchdog;
      helikopter04_3_DWork.HILInitialize_AOVoltages[2] =
        helikopter04_3_P.HILInitialize_AOWatchdog;
      helikopter04_3_DWork.HILInitialize_AOVoltages[3] =
        helikopter04_3_P.HILInitialize_AOWatchdog;
      result = hil_watchdog_set_analog_expiration_state
        (helikopter04_3_DWork.HILInitialize_Card,
         &helikopter04_3_P.HILInitialize_AOChannels[0], 4U,
         &helikopter04_3_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter04_3_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter04_3_P.HILInitialize_EIPStart && !is_switching) ||
        (helikopter04_3_P.HILInitialize_EIPEnter && is_switching)) {
      helikopter04_3_DWork.HILInitialize_QuadratureModes[0] =
        helikopter04_3_P.HILInitialize_EIQuadrature;
      helikopter04_3_DWork.HILInitialize_QuadratureModes[1] =
        helikopter04_3_P.HILInitialize_EIQuadrature;
      helikopter04_3_DWork.HILInitialize_QuadratureModes[2] =
        helikopter04_3_P.HILInitialize_EIQuadrature;
      helikopter04_3_DWork.HILInitialize_QuadratureModes[3] =
        helikopter04_3_P.HILInitialize_EIQuadrature;
      result = hil_set_encoder_quadrature_mode
        (helikopter04_3_DWork.HILInitialize_Card,
         &helikopter04_3_P.HILInitialize_EIChannels[0], 4U,
         (t_encoder_quadrature_mode *)
         &helikopter04_3_DWork.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter04_3_M, _rt_error_message);
        return;
      }

      helikopter04_3_DWork.HILInitialize_FilterFrequency[0] =
        helikopter04_3_P.HILInitialize_EIFrequency;
      helikopter04_3_DWork.HILInitialize_FilterFrequency[1] =
        helikopter04_3_P.HILInitialize_EIFrequency;
      helikopter04_3_DWork.HILInitialize_FilterFrequency[2] =
        helikopter04_3_P.HILInitialize_EIFrequency;
      helikopter04_3_DWork.HILInitialize_FilterFrequency[3] =
        helikopter04_3_P.HILInitialize_EIFrequency;
      result = hil_set_encoder_filter_frequency
        (helikopter04_3_DWork.HILInitialize_Card,
         &helikopter04_3_P.HILInitialize_EIChannels[0], 4U,
         &helikopter04_3_DWork.HILInitialize_FilterFrequency[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter04_3_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter04_3_P.HILInitialize_EIStart && !is_switching) ||
        (helikopter04_3_P.HILInitialize_EIEnter && is_switching)) {
      helikopter04_3_DWork.HILInitialize_InitialEICounts[0] =
        helikopter04_3_P.HILInitialize_EIInitial;
      helikopter04_3_DWork.HILInitialize_InitialEICounts[1] =
        helikopter04_3_P.HILInitialize_EIInitial;
      helikopter04_3_DWork.HILInitialize_InitialEICounts[2] =
        helikopter04_3_P.HILInitialize_EIInitial;
      helikopter04_3_DWork.HILInitialize_InitialEICounts[3] =
        helikopter04_3_P.HILInitialize_EIInitial;
      result = hil_set_encoder_counts(helikopter04_3_DWork.HILInitialize_Card,
        &helikopter04_3_P.HILInitialize_EIChannels[0], 4U,
        &helikopter04_3_DWork.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter04_3_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for FromWorkspace: '<Root>/From Workspace1' */
  {
    static real_T pTimeValues[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75 } ;

    static real_T pDataValues[] = { 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1394217997731140E+000, 3.1350768885809868E+000,
      3.1285586271465977E+000, 3.1198722531193908E+000, 3.1090220041114693E+000,
      3.0960047105800568E+000, 3.0808112100938687E+000, 3.0634368226320010E+000,
      3.0438898258081375E+000, 3.0221851157849957E+000, 2.9983258824117573E+000,
      2.9722942491822626E+000, 2.9440678284623765E+000, 2.9136495138154226E+000,
      2.8810738597005487E+000, 2.8463713284604064E+000, 2.8095235150528790E+000,
      2.7704689341777318E+000, 2.7291721291949060E+000, 2.6856845968098235E+000,
      2.6401079769891509E+000, 2.5924710823603316E+000, 2.5426589849314691E+000,
      2.4905160095021688E+000, 2.4360490135968731E+000, 2.3794850127197460E+000,
      2.3210413243708832E+000, 2.2608044050614380E+000, 2.1991055460646520E+000,
      2.1362563179091500E+000, 2.0725216520400829E+000, 2.0081303810128501E+000,
      1.9432774855701971E+000, 1.8781224445675542E+000, 1.8127892318653918E+000,
      1.7473682946919253E+000, 1.6819187630156478E+000, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.8367099231598242E-040,
      -8.6834152667160051E-003, -1.7379644768507573E-002,
      -2.6073045737556168E-002, -3.4745496108826490E-002,
      -4.3400996031685270E-002, -5.2069174125650361E-002,
      -6.0774001944753679E-002, -6.9497549847470586E-002,
      -7.8187987295453343E-002, -8.6818840092566879E-002,
      -9.5436933492952555E-002, -1.0412653291797888E-001,
      -1.1290568287954611E-001, -1.2167325858781605E-001,
      -1.3030261645949476E-001, -1.3881012496056960E-001,
      -1.4739125363010980E-001, -1.5621832350058990E-001,
      -1.6518721993130209E-001, -1.7395012954033079E-001,
      -1.8230647928269100E-001, -1.9054757851527776E-001,
      -1.9924838971545011E-001, -2.0857190171720022E-001,
      -2.1786798362118320E-001, -2.2625600350850725E-001,
      -2.3377475339545190E-001, -2.4094767723778118E-001,
      -2.4679543598714376E-001, -2.5139691262200825E-001,
      -2.5493866347626715E-001, -2.5756508410893109E-001,
      -2.5941158177061258E-001, -2.6062016401057225E-001,
      -2.6133285080864976E-001, -2.6168374869386535E-001,
      -2.6179812670511021E-001, -2.6179843716457196E-001, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.2371940531501560E-001,
      5.2449226522867942E-001, 5.2432166904727329E-001, 5.2305808389414943E-001,
      5.2203575817445902E-001, 5.2280041171552982E-001, 5.2501085215454946E-001,
      5.2613990918527032E-001, 5.2414293136831158E-001, 5.2054922578581475E-001,
      5.1977967331581920E-001, 5.2409238801976654E-001, 5.2949341426373930E-001,
      5.2879534088274383E-001, 5.2045906293678912E-001, 5.1311000983376964E-001,
      5.1755023406163125E-001, 5.3238358885838299E-001, 5.4093751833211512E-001,
      5.2851391627678912E-001, 5.0399323114812800E-001, 4.9704217253963540E-001,
      5.2476859939877540E-001, 5.6232530761489607E-001, 5.6067092693067522E-001,
      5.0590334012912763E-001, 4.5347540092853317E-001, 4.3261773089142064E-001,
      3.5269357051031563E-001, 2.7752704814422291E-001, 2.1361222447534892E-001,
      1.5840627329182805E-001, 1.1136708628897693E-001, 7.2892744679867838E-002,
      4.2983998185102096E-002, 2.1163565961865898E-002, 6.8984359483390363E-003,
      1.8724619253156656E-005, -4.6212980953686808E-005, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 2.0948776212600624E+000, 3.0914396546553109E-003,
      -6.8238472562417562E-004, -5.0543406124955031E-003,
      -4.0893028787616541E-003, 3.0586141642830215E-003, 8.8417617560786598E-003,
      4.5162281228834531E-003, -7.9879112678346437E-003,
      -1.4374822329987451E-002, -3.0782098799823749E-003,
      1.7250858815789472E-002, 2.1604104975891086E-002, -2.7922935239816800E-003,
      -3.3345111783818833E-002, -2.9396212412077342E-002,
      1.7760896911446694E-002, 5.9333419187006735E-002, 3.4215717894928444E-002,
      -4.9694408221304290E-002, -9.8082740514644523E-002,
      -2.7804234433970662E-002, 1.1090570743656032E-001, 1.5022683286448246E-001,
      -6.6175227368833684E-003, -2.1907034720619042E-001,
      -2.0971175680237789E-001, -8.3430680148449968E-002,
      -3.1969664152442007E-001, -3.0066608946437096E-001,
      -2.5565929467549603E-001, -2.2082380473408339E-001,
      -1.8815674801140453E-001, -1.5389736643643639E-001,
      -1.1963498597906293E-001, -8.7281728892944804E-002,
      -5.7060520054107451E-002, -2.7518845316343517E-002,
      -2.5975040082737389E-004, -3.3858189980697474E-005, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.6230631511411558E-007,
      2.0519883929218080E-006, 3.2549673722427490E-006, 4.6084741248338291E-006,
      5.3543671063659680E-006, 6.4038457086196035E-006, 6.7930522320859838E-006,
      3.4381087788040977E-006, 7.0715399215680290E-010, 1.8275940999382850E-007,
      2.9535957566465125E-009, 6.6721473224711066E-009, 6.3834036055766117E-006,
      1.1974842754824259E-005, 1.9292479518598220E-005, 2.3072678322406391E-005,
      2.1026547133725250E-005, 1.1677960122280505E-005, 7.2275471144823889E-006,
      2.1498694547023179E-005, 6.3258858328507306E-005, 1.8375573001642236E-004,
      5.1802412550542144E-004, 1.4003521352033868E-003, 3.5948677337676945E-003,
      8.6718479019526964E-003, 1.9370852072355813E-002, 3.9361659253524653E-002,
      7.1589410701868417E-002, 1.1494758550612363E-001, 1.6061847404069390E-001,
      1.9311196289573584E-001, 1.9809533605803428E-001, 1.7234276672327364E-001,
      1.2667017794412069E-001, 7.8479709666989278E-002, 4.0947061067476552E-002,
      1.7989086300724356E-002, 6.6564262842517029E-003, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 3.4492252604564835E-006, 4.7587283112308755E-006,
      4.8119159172836524E-006, 5.4140270103644068E-006, 2.9835719261285005E-006,
      4.1979144090145047E-006, 1.5568260938656153E-006, -1.3419773813127593E-005,
      -1.3749606499247820E-005, 7.2820902400668816E-007,
      -7.1922325694889534E-007, 1.4874206263330140E-008, 2.5506925833016493E-005,
      2.2365756596990622E-005, 2.9270547055095901E-005, 1.5120795215232712E-005,
      -8.1845247547246807E-006, -3.7394348045778946E-005,
      -1.7801652031192440E-005, 5.7084589730163219E-005, 1.6704065512593646E-004,
      4.8198748675166010E-004, 1.3370735819559962E-003, 3.5293120387918614E-003,
      8.7780623942572307E-003, 2.0307920672740010E-002, 4.2796016681612467E-002,
      7.9963228724675331E-002, 1.2891100579337506E-001, 1.7343269921702087E-001,
      1.8268355413828100E-001, 1.2997395542016782E-001, 1.9933492649193713E-002,
      -1.0301027733904246E-001, -1.8269035511661189E-001,
      -1.9276187310852563E-001, -1.5013059439805090E-001,
      -9.1831899067008796E-002, -4.5330640065890618E-002,
      -2.5069889165425543E-002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helikopter04_3_DWork.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues;
    helikopter04_3_DWork.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues;
    helikopter04_3_DWork.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    const char *fileName = "travel.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helikopter04_3_M, "Error creating .mat file travel.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"travel")) {
      rtmSetErrorStatus(helikopter04_3_M,
                        "Error writing mat file header to file travel.mat");
      return;
    }

    helikopter04_3_DWork.ToFile_IWORK.Count = 0;
    helikopter04_3_DWork.ToFile_IWORK.Decimation = -1;
    helikopter04_3_DWork.ToFile_PWORK.FilePtr = fp;
  }

  /* Start for FromWorkspace: '<Root>/From Workspace' */
  {
    static real_T pTimeValues[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75 } ;

    static real_T pDataValues[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      3.7032554293723352E-001, 5.4649448585890109E-004, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      4.9115474020613392E-001, 3.9696048514241261E-001, 3.7946697320814621E-001,
      3.0749896527571235E-001, 2.3849054570050182E-001, 1.8035049636413078E-001,
      1.3120080543833929E-001, 9.0218408825738686E-002, 5.7463369086396071E-002,
      3.2897028018029503E-002, 1.6298875428463359E-002, 6.8525181308787958E-003,
      1.2739280319643262E-005, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      4.3060921107147092E-006, 2.1840154069412119E-006, 9.2021778851098495E-007,
      1.7465739680675887E-006, -1.8093996979121004E-006, 2.5047191598777486E-006,
      -2.0320008368614268E-006, -1.7735491457205686E-005,
      -1.7913945637730224E-006, 1.6268363590669265E-005,
      -1.6909800888396110E-006, 8.2231593881368330E-007, 3.1827560837382800E-005,
      1.4056717500084493E-007, 1.2892755622963500E-005, -1.1669599682003458E-005,
      -2.4537007017071110E-005, -3.5197775590171502E-005,
      2.0849431316705462E-005, 9.1956837343895270E-005, 1.4716651640303519E-004,
      4.2218004082900719E-004, 1.1513040880867034E-003, 2.9702168738407063E-003,
      7.1723516404172561E-003, 1.5947903481565686E-002, 3.1708851754331004E-002,
      5.4181123231240809E-002, 7.5998669516456285E-002, 8.0494858267672872E-002,
      4.7142900562051247E-002, -2.3904406810149180E-002,
      -9.8779891086025112E-002, -1.2878819577573633E-001,
      -9.3797053772550570E-002, -2.2453470911493616E-002,
      3.6647736687743925E-002, 5.7624072924248815E-002, 4.7995185649841884E-002,
      2.0081246002859061E-002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helikopter04_3_DWork.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues;
    helikopter04_3_DWork.FromWorkspace_PWORK.DataPtr = (void *) pDataValues;
    helikopter04_3_DWork.FromWorkspace_IWORK.PrevIndex = 0;
  }

  MdlInitialize();
}

void MdlTerminate(void)
{
  helikopter04_3_terminate();
}

RT_MODEL_helikopter04_3 *helikopter04_3(void)
{
  helikopter04_3_initialize(1);
  return helikopter04_3_M;
}

/*========================================================================*
 * End of GRT compatible call interface                                   *
 *========================================================================*/
