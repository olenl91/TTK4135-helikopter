/*
 * helikopter03_3.c
 *
 * Real-Time Workshop code generation for Simulink model "helikopter03_3.mdl".
 *
 * Model version              : 1.64
 * Real-Time Workshop version : 7.5  (R2010a)  25-Jan-2010
 * C source code generated on : Wed Mar 04 16:41:38 2015
 *
 * Target selection: quarc_windows.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "helikopter03_3.h"
#include "helikopter03_3_private.h"
#include <stdio.h>
#include "helikopter03_3_dt.h"

/* Block signals (auto storage) */
BlockIO_helikopter03_3 helikopter03_3_B;

/* Continuous states */
ContinuousStates_helikopter03_3 helikopter03_3_X;

/* Block states (auto storage) */
D_Work_helikopter03_3 helikopter03_3_DWork;

/* Real-time model */
RT_MODEL_helikopter03_3 helikopter03_3_M_;
RT_MODEL_helikopter03_3 *helikopter03_3_M = &helikopter03_3_M_;

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
  helikopter03_3_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helikopter03_3_output(int_T tid)
{
  /* local block i/o variables */
  real_T rtb_HILReadEncoder_o1;
  real_T rtb_HILReadEncoder_o2;
  real_T rtb_HILReadEncoder_o3;
  real_T rtb_Saturation_k;
  real_T rtb_Gain2;
  real_T rtb_Gain[6];
  real_T rtb_Gain1;
  real_T tmp[6];
  int32_T tmp_0;
  int32_T tmp_1;
  if (rtmIsMajorTimeStep(helikopter03_3_M)) {
    /* set solver stop time */
    if (!(helikopter03_3_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helikopter03_3_M->solverInfo,
                            ((helikopter03_3_M->Timing.clockTickH0 + 1) *
        helikopter03_3_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helikopter03_3_M->solverInfo,
                            ((helikopter03_3_M->Timing.clockTick0 + 1) *
        helikopter03_3_M->Timing.stepSize0 +
        helikopter03_3_M->Timing.clockTickH0 *
        helikopter03_3_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helikopter03_3_M)) {
    helikopter03_3_M->Timing.t[0] = rtsiGetT(&helikopter03_3_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helikopter03_3_M)) {
    /* S-Function (hil_read_encoder_block): '<S2>/HIL Read Encoder' */

    /* S-Function Block: helikopter03_3/Heli 3D/HIL Read Encoder (hil_read_encoder_block) */
    {
      t_error result = hil_read_encoder(helikopter03_3_DWork.HILInitialize_Card,
        helikopter03_3_P.HILReadEncoder_Channels, 3,
        &helikopter03_3_DWork.HILReadEncoder_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_3_M, _rt_error_message);
      } else {
        rtb_HILReadEncoder_o1 = helikopter03_3_DWork.HILReadEncoder_Buffer[0];
        rtb_HILReadEncoder_o2 = helikopter03_3_DWork.HILReadEncoder_Buffer[1];
        rtb_HILReadEncoder_o3 = helikopter03_3_DWork.HILReadEncoder_Buffer[2];
      }
    }

    /* Gain: '<S2>/Kalibrer-Elev' */
    helikopter03_3_B.KalibrerElev = helikopter03_3_P.KalibrerElev_Gain *
      rtb_HILReadEncoder_o3;

    /* Sum: '<Root>/Add' incorporates:
     *  Constant: '<Root>/Constant'
     */
    helikopter03_3_B.Add = helikopter03_3_B.KalibrerElev +
      helikopter03_3_P.Constant_Value;
  }

  /* FromWorkspace: '<Root>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *)
      helikopter03_3_DWork.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopter03_3_DWork.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex = helikopter03_3_DWork.FromWorkspace1_IWORK.PrevIndex;
    real_T t = helikopter03_3_M->Timing.t[0];

    /* get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[139]) {
      currTimeIndex = 138;
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

    helikopter03_3_DWork.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          helikopter03_3_B.FromWorkspace1[0] = pDataValues[currTimeIndex];
          pDataValues += 140;
          helikopter03_3_B.FromWorkspace1[1] = pDataValues[currTimeIndex];
          pDataValues += 140;
          helikopter03_3_B.FromWorkspace1[2] = pDataValues[currTimeIndex];
          pDataValues += 140;
          helikopter03_3_B.FromWorkspace1[3] = pDataValues[currTimeIndex];
          pDataValues += 140;
        } else {
          helikopter03_3_B.FromWorkspace1[0] = pDataValues[currTimeIndex + 1];
          pDataValues += 140;
          helikopter03_3_B.FromWorkspace1[1] = pDataValues[currTimeIndex + 1];
          pDataValues += 140;
          helikopter03_3_B.FromWorkspace1[2] = pDataValues[currTimeIndex + 1];
          pDataValues += 140;
          helikopter03_3_B.FromWorkspace1[3] = pDataValues[currTimeIndex + 1];
          pDataValues += 140;
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        helikopter03_3_B.FromWorkspace1[0] = (real_T) rtInterpolate(d1, d2, f1,
          f2);
        pDataValues += 140;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        helikopter03_3_B.FromWorkspace1[1] = (real_T) rtInterpolate(d1, d2, f1,
          f2);
        pDataValues += 140;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        helikopter03_3_B.FromWorkspace1[2] = (real_T) rtInterpolate(d1, d2, f1,
          f2);
        pDataValues += 140;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        helikopter03_3_B.FromWorkspace1[3] = (real_T) rtInterpolate(d1, d2, f1,
          f2);
        pDataValues += 140;
      }
    }
  }

  if (rtmIsMajorTimeStep(helikopter03_3_M)) {
    /* Gain: '<S2>/Kalibrer-Pitch' */
    helikopter03_3_B.KalibrerPitch = helikopter03_3_P.KalibrerPitch_Gain *
      rtb_HILReadEncoder_o2;
  }

  /* TransferFcn: '<S2>/Vandring Lavpass' */
  helikopter03_3_B.VandringLavpass = helikopter03_3_P.VandringLavpass_C*
    helikopter03_3_X.VandringLavpass_CSTATE;

  /* Sum: '<Root>/Add1' incorporates:
   *  Constant: '<Root>/Constant1'
   */
  helikopter03_3_B.Add1 = helikopter03_3_P.Constant1_Value +
    helikopter03_3_B.VandringLavpass;
  if (rtmIsMajorTimeStep(helikopter03_3_M)) {
    /* ToFile: '<Root>/To File' */
    if (rtmIsMajorTimeStep(helikopter03_3_M)) {
      if (!(++helikopter03_3_DWork.ToFile_IWORK.Decimation % 1) &&
          (helikopter03_3_DWork.ToFile_IWORK.Count*2)+1 < 100000000 ) {
        FILE *fp = (FILE *) helikopter03_3_DWork.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          helikopter03_3_DWork.ToFile_IWORK.Decimation = 0;
          u[0] = helikopter03_3_M->Timing.t[1];
          u[1] = helikopter03_3_B.Add1;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(helikopter03_3_M,
                              "Error writing to MAT-file travel.mat");
            return;
          }

          if (((++helikopter03_3_DWork.ToFile_IWORK.Count)*2)+1 >= 100000000) {
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
    helikopter03_3_B.KalibrerVandring = helikopter03_3_P.KalibrerVandring_Gain *
      rtb_HILReadEncoder_o1;
  }

  /* TransferFcn: '<S2>/Vandring Deriv' */
  helikopter03_3_B.VandringDeriv = helikopter03_3_P.VandringDeriv_D*
    helikopter03_3_B.KalibrerVandring;
  helikopter03_3_B.VandringDeriv += helikopter03_3_P.VandringDeriv_C*
    helikopter03_3_X.VandringDeriv_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn4' */
  helikopter03_3_B.TransferFcn4 = helikopter03_3_P.TransferFcn4_D*
    helikopter03_3_B.KalibrerPitch;
  helikopter03_3_B.TransferFcn4 += helikopter03_3_P.TransferFcn4_C*
    helikopter03_3_X.TransferFcn4_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn5' */
  helikopter03_3_B.TransferFcn5 = helikopter03_3_P.TransferFcn5_D*
    helikopter03_3_B.KalibrerElev;
  helikopter03_3_B.TransferFcn5 += helikopter03_3_P.TransferFcn5_C*
    helikopter03_3_X.TransferFcn5_CSTATE;
  if (rtmIsMajorTimeStep(helikopter03_3_M)) {
  }

  /* Integrator: '<S1>/Integrator'
   *
   * Regarding '<S1>/Integrator':
   *  Limited Integrator
   */
  if (helikopter03_3_X.Integrator_CSTATE >= helikopter03_3_P.Integrator_UpperSat
      ) {
    helikopter03_3_X.Integrator_CSTATE = helikopter03_3_P.Integrator_UpperSat;
  } else if (helikopter03_3_X.Integrator_CSTATE <=
             helikopter03_3_P.Integrator_LowerSat ) {
    helikopter03_3_X.Integrator_CSTATE = helikopter03_3_P.Integrator_LowerSat;
  }

  rtb_Saturation_k = helikopter03_3_X.Integrator_CSTATE;

  /* Gain: '<Root>/Gain' incorporates:
   *  SignalConversion: '<Root>/TmpSignal ConversionAtGainInport1'
   */
  tmp[0] = helikopter03_3_B.Add1;
  tmp[1] = helikopter03_3_B.VandringDeriv;
  tmp[2] = helikopter03_3_B.KalibrerPitch;
  tmp[3] = helikopter03_3_B.TransferFcn4;
  tmp[4] = helikopter03_3_B.Add;
  tmp[5] = helikopter03_3_B.TransferFcn5;
  for (tmp_0 = 0; tmp_0 < 6; tmp_0++) {
    rtb_Gain[tmp_0] = 0.0;
    for (tmp_1 = 0; tmp_1 < 6; tmp_1++) {
      rtb_Gain[tmp_0] += helikopter03_3_P.Gain_Gain[6 * tmp_1 + tmp_0] *
        tmp[tmp_1];
    }
  }

  /* Sum: '<S1>/Sum' incorporates:
   *  Constant: '<Root>/elevation'
   */
  rtb_Gain2 = helikopter03_3_P.elevation_Value - rtb_Gain[4];

  /* Gain: '<S1>/K_ei' */
  helikopter03_3_B.K_ei = helikopter03_3_P.K_ei_Gain * rtb_Gain2;

  /* Sum: '<S1>/Sum1' incorporates:
   *  Gain: '<S1>/K_ed'
   *  Gain: '<S1>/K_ep'
   */
  rtb_Saturation_k = (helikopter03_3_P.K_ep_Gain * rtb_Gain2 + rtb_Saturation_k)
    + helikopter03_3_P.K_ed_Gain * rtb_Gain[5];

  /* Saturate: '<S1>/Saturation' */
  rtb_Saturation_k = rt_SATURATE(rtb_Saturation_k,
    helikopter03_3_P.Saturation_LowerSat, helikopter03_3_P.Saturation_UpperSat);

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helikopter03_3_DWork.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopter03_3_DWork.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helikopter03_3_DWork.FromWorkspace_IWORK.PrevIndex;
    real_T t = helikopter03_3_M->Timing.t[0];

    /* get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[139]) {
      currTimeIndex = 138;
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

    helikopter03_3_DWork.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Gain2 = pDataValues[currTimeIndex];
        } else {
          rtb_Gain2 = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Gain2 = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 140;
      }
    }
  }

  if (rtmIsMajorTimeStep(helikopter03_3_M)) {
  }

  /* Sum: '<S3>/Sum2' incorporates:
   *  Gain: '<S3>/Gain'
   *  Sum: '<S3>/Sum1'
   */
  rtb_Gain1 = rtb_Gain2 - ((((rtb_Gain[0] - helikopter03_3_B.FromWorkspace1[0]) *
    helikopter03_3_P.Gain_Gain_a[0] + (rtb_Gain[1] -
    helikopter03_3_B.FromWorkspace1[1]) * helikopter03_3_P.Gain_Gain_a[1]) +
    (rtb_Gain[2] - helikopter03_3_B.FromWorkspace1[2]) *
    helikopter03_3_P.Gain_Gain_a[2]) + (rtb_Gain[3] -
    helikopter03_3_B.FromWorkspace1[3]) * helikopter03_3_P.Gain_Gain_a[3]);

  /* Sum: '<S4>/Sum' incorporates:
   *  Gain: '<S4>/K_pd'
   *  Gain: '<S4>/K_pp'
   *  Saturate: '<S4>/Saturation'
   *  Sum: '<S4>/Sum1'
   */
  rtb_Gain1 = (rt_SATURATE(rtb_Gain1, helikopter03_3_P.Saturation_LowerSat_a,
    helikopter03_3_P.Saturation_UpperSat_h) - rtb_Gain[2]) *
    helikopter03_3_P.K_pp_Gain - helikopter03_3_P.K_pd_Gain * rtb_Gain[3];

  /* Gain: '<S5>/Gain2' incorporates:
   *  Sum: '<S5>/Sum4'
   */
  rtb_Gain2 = (rtb_Saturation_k - rtb_Gain1) * helikopter03_3_P.Gain2_Gain;

  /* Saturate: '<S2>/Sat B' */
  helikopter03_3_B.SatB = rt_SATURATE(rtb_Gain2, helikopter03_3_P.SatB_LowerSat,
    helikopter03_3_P.SatB_UpperSat);
  if (rtmIsMajorTimeStep(helikopter03_3_M)) {
  }

  /* Gain: '<S5>/Gain1' incorporates:
   *  Sum: '<S5>/Sum3'
   */
  rtb_Gain1 = (rtb_Gain1 + rtb_Saturation_k) * helikopter03_3_P.Gain1_Gain;

  /* Saturate: '<S2>/Sat' */
  helikopter03_3_B.Sat = rt_SATURATE(rtb_Gain1, helikopter03_3_P.Sat_LowerSat,
    helikopter03_3_P.Sat_UpperSat);
  if (rtmIsMajorTimeStep(helikopter03_3_M)) {
    /* S-Function (hil_write_analog_block): '<S2>/HIL Write Analog' */

    /* S-Function Block: helikopter03_3/Heli 3D/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helikopter03_3_DWork.HILWriteAnalog_Buffer[0] = helikopter03_3_B.SatB;
      helikopter03_3_DWork.HILWriteAnalog_Buffer[1] = helikopter03_3_B.Sat;
      result = hil_write_analog(helikopter03_3_DWork.HILInitialize_Card,
        helikopter03_3_P.HILWriteAnalog_Channels, 2,
        &helikopter03_3_DWork.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_3_M, _rt_error_message);
      }
    }
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Model update function */
void helikopter03_3_update(int_T tid)
{
  if (rtmIsMajorTimeStep(helikopter03_3_M)) {
    rt_ertODEUpdateContinuousStates(&helikopter03_3_M->solverInfo);
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
  if (!(++helikopter03_3_M->Timing.clockTick0)) {
    ++helikopter03_3_M->Timing.clockTickH0;
  }

  helikopter03_3_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helikopter03_3_M->solverInfo);
  if (rtmIsMajorTimeStep(helikopter03_3_M)) {
    /* Update absolute timer for sample time: [0.001s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helikopter03_3_M->Timing.clockTick1)) {
      ++helikopter03_3_M->Timing.clockTickH1;
    }

    helikopter03_3_M->Timing.t[1] = helikopter03_3_M->Timing.clockTick1 *
      helikopter03_3_M->Timing.stepSize1 + helikopter03_3_M->Timing.clockTickH1 *
      helikopter03_3_M->Timing.stepSize1 * 4294967296.0;
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Derivatives for root system: '<Root>' */
void helikopter03_3_derivatives(void)
{
  /* Derivatives for TransferFcn: '<S2>/Vandring Lavpass' */
  {
    ((StateDerivatives_helikopter03_3 *) helikopter03_3_M->ModelData.derivs)
      ->VandringLavpass_CSTATE = helikopter03_3_B.KalibrerVandring;
    ((StateDerivatives_helikopter03_3 *) helikopter03_3_M->ModelData.derivs)
      ->VandringLavpass_CSTATE += (helikopter03_3_P.VandringLavpass_A)*
      helikopter03_3_X.VandringLavpass_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Vandring Deriv' */
  {
    ((StateDerivatives_helikopter03_3 *) helikopter03_3_M->ModelData.derivs)
      ->VandringDeriv_CSTATE = helikopter03_3_B.KalibrerVandring;
    ((StateDerivatives_helikopter03_3 *) helikopter03_3_M->ModelData.derivs)
      ->VandringDeriv_CSTATE += (helikopter03_3_P.VandringDeriv_A)*
      helikopter03_3_X.VandringDeriv_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn4' */
  {
    ((StateDerivatives_helikopter03_3 *) helikopter03_3_M->ModelData.derivs)
      ->TransferFcn4_CSTATE = helikopter03_3_B.KalibrerPitch;
    ((StateDerivatives_helikopter03_3 *) helikopter03_3_M->ModelData.derivs)
      ->TransferFcn4_CSTATE += (helikopter03_3_P.TransferFcn4_A)*
      helikopter03_3_X.TransferFcn4_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn5' */
  {
    ((StateDerivatives_helikopter03_3 *) helikopter03_3_M->ModelData.derivs)
      ->TransferFcn5_CSTATE = helikopter03_3_B.KalibrerElev;
    ((StateDerivatives_helikopter03_3 *) helikopter03_3_M->ModelData.derivs)
      ->TransferFcn5_CSTATE += (helikopter03_3_P.TransferFcn5_A)*
      helikopter03_3_X.TransferFcn5_CSTATE;
  }

  /* Derivatives for Integrator: '<S1>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helikopter03_3_X.Integrator_CSTATE <=
            helikopter03_3_P.Integrator_LowerSat );
    usat = ( helikopter03_3_X.Integrator_CSTATE >=
            helikopter03_3_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helikopter03_3_B.K_ei > 0)) ||
        (usat && (helikopter03_3_B.K_ei < 0)) ) {
      ((StateDerivatives_helikopter03_3 *) helikopter03_3_M->ModelData.derivs)
        ->Integrator_CSTATE = helikopter03_3_B.K_ei;
    } else {
      /* in saturation */
      ((StateDerivatives_helikopter03_3 *) helikopter03_3_M->ModelData.derivs)
        ->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helikopter03_3_initialize(boolean_T firstTime)
{
  (void)firstTime;

  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helikopter03_3_P.Integrator_UpperSat = rtInf;
  helikopter03_3_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helikopter03_3_M, 0,
                sizeof(RT_MODEL_helikopter03_3));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helikopter03_3_M->solverInfo,
                          &helikopter03_3_M->Timing.simTimeStep);
    rtsiSetTPtr(&helikopter03_3_M->solverInfo, &rtmGetTPtr(helikopter03_3_M));
    rtsiSetStepSizePtr(&helikopter03_3_M->solverInfo,
                       &helikopter03_3_M->Timing.stepSize0);
    rtsiSetdXPtr(&helikopter03_3_M->solverInfo,
                 &helikopter03_3_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helikopter03_3_M->solverInfo,
                         &helikopter03_3_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helikopter03_3_M->solverInfo,
      &helikopter03_3_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helikopter03_3_M->solverInfo, (&rtmGetErrorStatus
      (helikopter03_3_M)));
    rtsiSetRTModelPtr(&helikopter03_3_M->solverInfo, helikopter03_3_M);
  }

  rtsiSetSimTimeStep(&helikopter03_3_M->solverInfo, MAJOR_TIME_STEP);
  helikopter03_3_M->ModelData.intgData.f[0] = helikopter03_3_M->ModelData.odeF[0];
  helikopter03_3_M->ModelData.contStates = ((real_T *) &helikopter03_3_X);
  rtsiSetSolverData(&helikopter03_3_M->solverInfo, (void *)
                    &helikopter03_3_M->ModelData.intgData);
  rtsiSetSolverName(&helikopter03_3_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helikopter03_3_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helikopter03_3_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helikopter03_3_M->Timing.sampleTimes =
      (&helikopter03_3_M->Timing.sampleTimesArray[0]);
    helikopter03_3_M->Timing.offsetTimes =
      (&helikopter03_3_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helikopter03_3_M->Timing.sampleTimes[0] = (0.0);
    helikopter03_3_M->Timing.sampleTimes[1] = (0.001);

    /* task offsets */
    helikopter03_3_M->Timing.offsetTimes[0] = (0.0);
    helikopter03_3_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helikopter03_3_M, &helikopter03_3_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helikopter03_3_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helikopter03_3_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helikopter03_3_M, -1);
  helikopter03_3_M->Timing.stepSize0 = 0.001;
  helikopter03_3_M->Timing.stepSize1 = 0.001;

  /* external mode info */
  helikopter03_3_M->Sizes.checksums[0] = (678985645U);
  helikopter03_3_M->Sizes.checksums[1] = (561742273U);
  helikopter03_3_M->Sizes.checksums[2] = (3765221545U);
  helikopter03_3_M->Sizes.checksums[3] = (1163615554U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helikopter03_3_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helikopter03_3_M->extModeInfo,
      &helikopter03_3_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helikopter03_3_M->extModeInfo,
                        helikopter03_3_M->Sizes.checksums);
    rteiSetTPtr(helikopter03_3_M->extModeInfo, rtmGetTPtr(helikopter03_3_M));
  }

  helikopter03_3_M->solverInfoPtr = (&helikopter03_3_M->solverInfo);
  helikopter03_3_M->Timing.stepSize = (0.001);
  rtsiSetFixedStepSize(&helikopter03_3_M->solverInfo, 0.001);
  rtsiSetSolverMode(&helikopter03_3_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helikopter03_3_M->ModelData.blockIO = ((void *) &helikopter03_3_B);

  {
    helikopter03_3_B.KalibrerElev = 0.0;
    helikopter03_3_B.Add = 0.0;
    helikopter03_3_B.FromWorkspace1[0] = 0.0;
    helikopter03_3_B.FromWorkspace1[1] = 0.0;
    helikopter03_3_B.FromWorkspace1[2] = 0.0;
    helikopter03_3_B.FromWorkspace1[3] = 0.0;
    helikopter03_3_B.KalibrerPitch = 0.0;
    helikopter03_3_B.VandringLavpass = 0.0;
    helikopter03_3_B.Add1 = 0.0;
    helikopter03_3_B.KalibrerVandring = 0.0;
    helikopter03_3_B.VandringDeriv = 0.0;
    helikopter03_3_B.TransferFcn4 = 0.0;
    helikopter03_3_B.TransferFcn5 = 0.0;
    helikopter03_3_B.K_ei = 0.0;
    helikopter03_3_B.SatB = 0.0;
    helikopter03_3_B.Sat = 0.0;
  }

  /* parameters */
  helikopter03_3_M->ModelData.defaultParam = ((real_T *)&helikopter03_3_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helikopter03_3_X;
    helikopter03_3_M->ModelData.contStates = (x);
    (void) memset((void *)&helikopter03_3_X, 0,
                  sizeof(ContinuousStates_helikopter03_3));
  }

  /* states (dwork) */
  helikopter03_3_M->Work.dwork = ((void *) &helikopter03_3_DWork);
  (void) memset((void *)&helikopter03_3_DWork, 0,
                sizeof(D_Work_helikopter03_3));
  helikopter03_3_DWork.HILInitialize_AIMinimums[0] = 0.0;
  helikopter03_3_DWork.HILInitialize_AIMinimums[1] = 0.0;
  helikopter03_3_DWork.HILInitialize_AIMinimums[2] = 0.0;
  helikopter03_3_DWork.HILInitialize_AIMinimums[3] = 0.0;
  helikopter03_3_DWork.HILInitialize_AIMaximums[0] = 0.0;
  helikopter03_3_DWork.HILInitialize_AIMaximums[1] = 0.0;
  helikopter03_3_DWork.HILInitialize_AIMaximums[2] = 0.0;
  helikopter03_3_DWork.HILInitialize_AIMaximums[3] = 0.0;
  helikopter03_3_DWork.HILInitialize_AOMinimums[0] = 0.0;
  helikopter03_3_DWork.HILInitialize_AOMinimums[1] = 0.0;
  helikopter03_3_DWork.HILInitialize_AOMinimums[2] = 0.0;
  helikopter03_3_DWork.HILInitialize_AOMinimums[3] = 0.0;
  helikopter03_3_DWork.HILInitialize_AOMaximums[0] = 0.0;
  helikopter03_3_DWork.HILInitialize_AOMaximums[1] = 0.0;
  helikopter03_3_DWork.HILInitialize_AOMaximums[2] = 0.0;
  helikopter03_3_DWork.HILInitialize_AOMaximums[3] = 0.0;
  helikopter03_3_DWork.HILInitialize_AOVoltages[0] = 0.0;
  helikopter03_3_DWork.HILInitialize_AOVoltages[1] = 0.0;
  helikopter03_3_DWork.HILInitialize_AOVoltages[2] = 0.0;
  helikopter03_3_DWork.HILInitialize_AOVoltages[3] = 0.0;
  helikopter03_3_DWork.HILInitialize_FilterFrequency[0] = 0.0;
  helikopter03_3_DWork.HILInitialize_FilterFrequency[1] = 0.0;
  helikopter03_3_DWork.HILInitialize_FilterFrequency[2] = 0.0;
  helikopter03_3_DWork.HILInitialize_FilterFrequency[3] = 0.0;
  helikopter03_3_DWork.HILWriteAnalog_Buffer[0] = 0.0;
  helikopter03_3_DWork.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helikopter03_3_M->SpecialInfo.mappingInfo = (&dtInfo);
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
void helikopter03_3_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopter03_3/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    hil_task_stop_all(helikopter03_3_DWork.HILInitialize_Card);
    hil_task_delete_all(helikopter03_3_DWork.HILInitialize_Card);
    hil_monitor_stop_all(helikopter03_3_DWork.HILInitialize_Card);
    hil_monitor_delete_all(helikopter03_3_DWork.HILInitialize_Card);
    is_switching = false;
    if ((helikopter03_3_P.HILInitialize_AOTerminate && !is_switching) ||
        (helikopter03_3_P.HILInitialize_AOExit && is_switching)) {
      helikopter03_3_DWork.HILInitialize_AOVoltages[0] =
        helikopter03_3_P.HILInitialize_AOFinal;
      helikopter03_3_DWork.HILInitialize_AOVoltages[1] =
        helikopter03_3_P.HILInitialize_AOFinal;
      helikopter03_3_DWork.HILInitialize_AOVoltages[2] =
        helikopter03_3_P.HILInitialize_AOFinal;
      helikopter03_3_DWork.HILInitialize_AOVoltages[3] =
        helikopter03_3_P.HILInitialize_AOFinal;
      result = hil_write_analog(helikopter03_3_DWork.HILInitialize_Card,
        &helikopter03_3_P.HILInitialize_AOChannels[0], 4U,
        &helikopter03_3_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_3_M, _rt_error_message);
      }
    }

    hil_close(helikopter03_3_DWork.HILInitialize_Card);
    helikopter03_3_DWork.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helikopter03_3_DWork.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      const char *fileName = "travel.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopter03_3_M, "Error closing MAT-file travel.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helikopter03_3_M,
                          "Error reopening MAT-file travel.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2, helikopter03_3_DWork.ToFile_IWORK.Count,
           "travel")) {
        rtmSetErrorStatus(helikopter03_3_M,
                          "Error writing header for travel to MAT-file travel.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopter03_3_M, "Error closing MAT-file travel.mat");
        return;
      }

      helikopter03_3_DWork.ToFile_PWORK.FilePtr = (NULL);
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
  helikopter03_3_output(tid);
}

void MdlUpdate(int_T tid)
{
  helikopter03_3_update(tid);
}

void MdlInitializeSizes(void)
{
  helikopter03_3_M->Sizes.numContStates = (5);/* Number of continuous states */
  helikopter03_3_M->Sizes.numY = (0);  /* Number of model outputs */
  helikopter03_3_M->Sizes.numU = (0);  /* Number of model inputs */
  helikopter03_3_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helikopter03_3_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helikopter03_3_M->Sizes.numBlocks = (52);/* Number of blocks */
  helikopter03_3_M->Sizes.numBlockIO = (13);/* Number of block outputs */
  helikopter03_3_M->Sizes.numBlockPrms = (151);/* Sum of parameter "widths" */
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
  /* InitializeConditions for TransferFcn: '<S2>/Vandring Lavpass' */
  helikopter03_3_X.VandringLavpass_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Vandring Deriv' */
  helikopter03_3_X.VandringDeriv_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn4' */
  helikopter03_3_X.TransferFcn4_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn5' */
  helikopter03_3_X.TransferFcn5_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S1>/Integrator' */
  helikopter03_3_X.Integrator_CSTATE = helikopter03_3_P.Integrator_IC;
}

void MdlStart(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopter03_3/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q4", "0", &helikopter03_3_DWork.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter03_3_M, _rt_error_message);
      return;
    }

    is_switching = false;
    if ((helikopter03_3_P.HILInitialize_CKPStart && !is_switching) ||
        (helikopter03_3_P.HILInitialize_CKPEnter && is_switching)) {
      result = hil_set_clock_mode(helikopter03_3_DWork.HILInitialize_Card,
        (t_clock *) &helikopter03_3_P.HILInitialize_CKChannels[0], 2U,
        (t_clock_mode *) &helikopter03_3_P.HILInitialize_CKModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_3_M, _rt_error_message);
        return;
      }
    }

    result = hil_watchdog_clear(helikopter03_3_DWork.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter03_3_M, _rt_error_message);
      return;
    }

    if ((helikopter03_3_P.HILInitialize_AIPStart && !is_switching) ||
        (helikopter03_3_P.HILInitialize_AIPEnter && is_switching)) {
      helikopter03_3_DWork.HILInitialize_AIMinimums[0] =
        helikopter03_3_P.HILInitialize_AILow;
      helikopter03_3_DWork.HILInitialize_AIMinimums[1] =
        helikopter03_3_P.HILInitialize_AILow;
      helikopter03_3_DWork.HILInitialize_AIMinimums[2] =
        helikopter03_3_P.HILInitialize_AILow;
      helikopter03_3_DWork.HILInitialize_AIMinimums[3] =
        helikopter03_3_P.HILInitialize_AILow;
      helikopter03_3_DWork.HILInitialize_AIMaximums[0] =
        helikopter03_3_P.HILInitialize_AIHigh;
      helikopter03_3_DWork.HILInitialize_AIMaximums[1] =
        helikopter03_3_P.HILInitialize_AIHigh;
      helikopter03_3_DWork.HILInitialize_AIMaximums[2] =
        helikopter03_3_P.HILInitialize_AIHigh;
      helikopter03_3_DWork.HILInitialize_AIMaximums[3] =
        helikopter03_3_P.HILInitialize_AIHigh;
      result = hil_set_analog_input_ranges
        (helikopter03_3_DWork.HILInitialize_Card,
         &helikopter03_3_P.HILInitialize_AIChannels[0], 4U,
         &helikopter03_3_DWork.HILInitialize_AIMinimums[0],
         &helikopter03_3_DWork.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_3_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter03_3_P.HILInitialize_AOPStart && !is_switching) ||
        (helikopter03_3_P.HILInitialize_AOPEnter && is_switching)) {
      helikopter03_3_DWork.HILInitialize_AOMinimums[0] =
        helikopter03_3_P.HILInitialize_AOLow;
      helikopter03_3_DWork.HILInitialize_AOMinimums[1] =
        helikopter03_3_P.HILInitialize_AOLow;
      helikopter03_3_DWork.HILInitialize_AOMinimums[2] =
        helikopter03_3_P.HILInitialize_AOLow;
      helikopter03_3_DWork.HILInitialize_AOMinimums[3] =
        helikopter03_3_P.HILInitialize_AOLow;
      helikopter03_3_DWork.HILInitialize_AOMaximums[0] =
        helikopter03_3_P.HILInitialize_AOHigh;
      helikopter03_3_DWork.HILInitialize_AOMaximums[1] =
        helikopter03_3_P.HILInitialize_AOHigh;
      helikopter03_3_DWork.HILInitialize_AOMaximums[2] =
        helikopter03_3_P.HILInitialize_AOHigh;
      helikopter03_3_DWork.HILInitialize_AOMaximums[3] =
        helikopter03_3_P.HILInitialize_AOHigh;
      result = hil_set_analog_output_ranges
        (helikopter03_3_DWork.HILInitialize_Card,
         &helikopter03_3_P.HILInitialize_AOChannels[0], 4U,
         &helikopter03_3_DWork.HILInitialize_AOMinimums[0],
         &helikopter03_3_DWork.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_3_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter03_3_P.HILInitialize_AOStart && !is_switching) ||
        (helikopter03_3_P.HILInitialize_AOEnter && is_switching)) {
      helikopter03_3_DWork.HILInitialize_AOVoltages[0] =
        helikopter03_3_P.HILInitialize_AOInitial;
      helikopter03_3_DWork.HILInitialize_AOVoltages[1] =
        helikopter03_3_P.HILInitialize_AOInitial;
      helikopter03_3_DWork.HILInitialize_AOVoltages[2] =
        helikopter03_3_P.HILInitialize_AOInitial;
      helikopter03_3_DWork.HILInitialize_AOVoltages[3] =
        helikopter03_3_P.HILInitialize_AOInitial;
      result = hil_write_analog(helikopter03_3_DWork.HILInitialize_Card,
        &helikopter03_3_P.HILInitialize_AOChannels[0], 4U,
        &helikopter03_3_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_3_M, _rt_error_message);
        return;
      }
    }

    if (helikopter03_3_P.HILInitialize_AOReset) {
      helikopter03_3_DWork.HILInitialize_AOVoltages[0] =
        helikopter03_3_P.HILInitialize_AOWatchdog;
      helikopter03_3_DWork.HILInitialize_AOVoltages[1] =
        helikopter03_3_P.HILInitialize_AOWatchdog;
      helikopter03_3_DWork.HILInitialize_AOVoltages[2] =
        helikopter03_3_P.HILInitialize_AOWatchdog;
      helikopter03_3_DWork.HILInitialize_AOVoltages[3] =
        helikopter03_3_P.HILInitialize_AOWatchdog;
      result = hil_watchdog_set_analog_expiration_state
        (helikopter03_3_DWork.HILInitialize_Card,
         &helikopter03_3_P.HILInitialize_AOChannels[0], 4U,
         &helikopter03_3_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_3_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter03_3_P.HILInitialize_EIPStart && !is_switching) ||
        (helikopter03_3_P.HILInitialize_EIPEnter && is_switching)) {
      helikopter03_3_DWork.HILInitialize_QuadratureModes[0] =
        helikopter03_3_P.HILInitialize_EIQuadrature;
      helikopter03_3_DWork.HILInitialize_QuadratureModes[1] =
        helikopter03_3_P.HILInitialize_EIQuadrature;
      helikopter03_3_DWork.HILInitialize_QuadratureModes[2] =
        helikopter03_3_P.HILInitialize_EIQuadrature;
      helikopter03_3_DWork.HILInitialize_QuadratureModes[3] =
        helikopter03_3_P.HILInitialize_EIQuadrature;
      result = hil_set_encoder_quadrature_mode
        (helikopter03_3_DWork.HILInitialize_Card,
         &helikopter03_3_P.HILInitialize_EIChannels[0], 4U,
         (t_encoder_quadrature_mode *)
         &helikopter03_3_DWork.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_3_M, _rt_error_message);
        return;
      }

      helikopter03_3_DWork.HILInitialize_FilterFrequency[0] =
        helikopter03_3_P.HILInitialize_EIFrequency;
      helikopter03_3_DWork.HILInitialize_FilterFrequency[1] =
        helikopter03_3_P.HILInitialize_EIFrequency;
      helikopter03_3_DWork.HILInitialize_FilterFrequency[2] =
        helikopter03_3_P.HILInitialize_EIFrequency;
      helikopter03_3_DWork.HILInitialize_FilterFrequency[3] =
        helikopter03_3_P.HILInitialize_EIFrequency;
      result = hil_set_encoder_filter_frequency
        (helikopter03_3_DWork.HILInitialize_Card,
         &helikopter03_3_P.HILInitialize_EIChannels[0], 4U,
         &helikopter03_3_DWork.HILInitialize_FilterFrequency[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_3_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter03_3_P.HILInitialize_EIStart && !is_switching) ||
        (helikopter03_3_P.HILInitialize_EIEnter && is_switching)) {
      helikopter03_3_DWork.HILInitialize_InitialEICounts[0] =
        helikopter03_3_P.HILInitialize_EIInitial;
      helikopter03_3_DWork.HILInitialize_InitialEICounts[1] =
        helikopter03_3_P.HILInitialize_EIInitial;
      helikopter03_3_DWork.HILInitialize_InitialEICounts[2] =
        helikopter03_3_P.HILInitialize_EIInitial;
      helikopter03_3_DWork.HILInitialize_InitialEICounts[3] =
        helikopter03_3_P.HILInitialize_EIInitial;
      result = hil_set_encoder_counts(helikopter03_3_DWork.HILInitialize_Card,
        &helikopter03_3_P.HILInitialize_EIChannels[0], 4U,
        &helikopter03_3_DWork.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_3_M, _rt_error_message);
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
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75 } ;

    static real_T pDataValues[] = { 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897913E+000,
      3.1415926535897927E+000, 3.1415926535897927E+000, 3.1394222997917800E+000,
      3.1350815921957644E+000, 3.1285705308017406E+000, 3.1198891156097086E+000,
      3.1090373466196675E+000, 3.0960152238316159E+000, 3.0808227472455596E+000,
      3.0634599168614938E+000, 3.0439267326794210E+000, 3.0222231946993410E+000,
      2.9983493029212491E+000, 2.9723050573451522E+000, 2.9440904579710474E+000,
      2.9137055047989331E+000, 2.8811501978288097E+000, 2.8464245370606784E+000,
      2.8095285224945408E+000, 2.7704621541303949E+000, 2.7292254319682399E+000,
      2.6858183560080784E+000, 2.6402409262499105E+000, 2.5924931426937294E+000,
      2.5425750053395406E+000, 2.4904865141873502E+000, 2.4362276692371458E+000,
      2.3797984704889354E+000, 2.3211989179427128E+000, 2.2604290115984877E+000,
      2.1977172483364606E+000, 2.1334432702748529E+000, 2.0679612132912251E+000,
      2.0015548081893528E+000, 1.9344719336352727E+000, 1.8669514417887767E+000,
      1.7992239523182436E+000, 1.7315048943112468E+000, 1.6639921306455814E+000,
      1.5968675936909547E+000, 1.5302991300241484E+000, 1.4644412579663528E+000,
      1.3994354427196429E+000, 1.3354104965427982E+000, 1.2724831605823368E+000,
      1.2107587127132009E+000, 1.1503315242354899E+000, 1.0912855801453745E+000,
      1.0336949902844907E+000, 9.7762449691323561E-001, 9.2312997231865457E-001,
      8.7025890187446098E-001, 8.1905085219992069E-001, 7.6953792523043596E-001,
      7.2174519829485173E-001, 6.7569114973921240E-001, 6.3138806968698635E-001,
      5.8884245575856187E-001, 5.4805539369627421E-001, 5.0902292286687312E-001,
      4.7173638663003337E-001, 4.3618276760029684E-001, 4.0234500787702743E-001,
      3.7020231435823975E-001, 3.3973044928828261E-001, 3.1090200621998593E-001,
      2.8368667160107425E-001, 2.5805147222201225E-001, 2.3396100878746184E-001,
      2.1137767589606027E-001, 1.9026186873397738E-001, 1.7057217680723505E-001,
      1.5226556505578373E-001, 1.3529754270642416E-001, 1.1962232023071502E-001,
      1.0519295478464642E-001, 9.1961484532090243E-002, 7.9879052288913038E-002,
      6.8896018917483148E-002, 5.8962066801908353E-002, 5.0026293653795136E-002,
      4.2037297157437162E-002, 3.4943251603656800E-002, 2.8691977633841244E-002,
      2.3231003988769743E-002, 1.8507617318872349E-002, 1.4468900038706647E-002,
      1.1061771185645356E-002, 8.2330515430311876E-003, 5.9295310459889236E-003,
      4.0979375865115462E-003, 2.6847508832659487E-003, 1.6361357655773150E-003,
      8.9855327876641616E-004, 4.1984491399408869E-004, 1.4859738111082838E-004,
      2.9593628253519906E-005, -3.6648624472349444E-016, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -9.8725394512401298E-016,
      -3.2521775393214009E-015, -8.6814151920363170E-003,
      -1.7362830384067166E-002, -2.6044245576099221E-002,
      -3.4725660768130384E-002, -4.3407075960162925E-002,
      -5.2088491152195007E-002, -6.0769906344228276E-002,
      -6.9451321536259977E-002, -7.8132736728293156E-002,
      -8.6814151920324989E-002, -9.5495567112358321E-002,
      -1.0417698230439051E-001, -1.1285839749642305E-001,
      -1.2153981268845668E-001, -1.3022122788048907E-001,
      -1.3890264307252054E-001, -1.4758405826455248E-001,
      -1.5626547345658520E-001, -1.6494688864861748E-001,
      -1.7362830384065031E-001, -1.8230971903268325E-001,
      -1.9099113422471528E-001, -1.9967254941674709E-001,
      -2.0835396460877875E-001, -2.1703537980081100E-001,
      -2.2571679499284419E-001, -2.3439821018487830E-001,
      -2.4307962537691086E-001, -2.5084705304811206E-001,
      -2.5709591224641487E-001, -2.6192822793450771E-001,
      -2.6562562040748805E-001, -2.6833149821632413E-001,
      -2.7008196738598633E-001, -2.7090995788212818E-001,
      -2.7087623202798311E-001, -2.7005105466267065E-001,
      -2.6849814781851228E-001, -2.6627385466722625E-001,
      -2.6343148823118057E-001, -2.6002326098683454E-001,
      -2.5609978470738642E-001, -2.5170934384185234E-001,
      -2.4689779147653637E-001, -2.4170875391084240E-001,
      -2.3618377636045082E-001, -2.3036235944353814E-001,
      -2.2428197348502971E-001, -2.1797809837831356E-001,
      -2.1148428177676665E-001, -2.0483219869815092E-001,
      -1.9805170787794080E-001, -1.9117090774233272E-001,
      -1.8421619422256183E-001, -1.7721232020890682E-001,
      -1.7018245571369678E-001, -1.6314824824915117E-001,
      -1.5612988331760169E-001, -1.4914614494735659E-001,
      -1.4221447611894431E-001, -1.3535103889307662E-001,
      -1.2857077407514456E-001, -1.2188746027982675E-001,
      -1.1531377227318869E-001, -1.0886133847565090E-001,
      -1.0254079751624659E-001, -9.6361853738200956E-002,
      -9.0333331565602984E-002, -8.4463228648330044E-002,
      -7.8758767706967811E-002, -7.3226447005802831E-002,
      -6.7872089397440202E-002, -6.2700889902834914E-002,
      -5.7717461784272880E-002, -5.2925881010223066E-002,
      -4.8329728972711629E-002, -4.3932133485720899E-002,
      -3.9735808462297559E-002, -3.5743092592449267E-002,
      -3.1955985985432436E-002, -2.8376182215122614E-002,
      -2.5005095879257738E-002, -2.1843894580289420E-002,
      -1.8893546679588570E-002, -1.6154869120664276E-002,
      -1.3628515412246647E-002, -1.1314878570455529E-002,
      -9.2140819881688669E-003, -7.3263738379079702E-003,
      -5.6527468129810386E-003, -4.1944604707541749E-003,
      -2.9503299472431370E-003, -1.9148334590913159E-003,
      -1.0849901315341563E-003, -4.7601501142772334E-004,
      -1.1837451301833624E-004, 2.1091971947117046E-016, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.5207646309492972E-017,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 4.6847380619717482E-001, 3.7688498392758313E-001,
      2.9144955305350478E-001, 2.2299937613104578E-001, 1.6319854266678158E-001,
      1.0557535766734777E-001, 4.9938264717986569E-002, -2.0340941592853616E-003,
      -4.9768597466189106E-002, -9.3659858932440165E-002,
      -1.3415291687168684E-001, -1.7143052748838297E-001,
      -2.0555906757441414E-001, -2.3663507971526906E-001,
      -2.6479893089776274E-001, -2.9019726294377529E-001,
      -3.1296437917459669E-001, -3.3322579517315770E-001,
      -3.5110482594364412E-001, -3.6672392376329049E-001,
      -3.8020313674556966E-001, -3.9165900332134146E-001,
      -4.0120446702485635E-001, -4.0894907257457797E-001,
      -4.1499898881087305E-001, -4.1945689764159083E-001,
      -4.2242189514907968E-001, -4.2398944882763667E-001,
      -4.2425138462680079E-001, -4.2329588017330699E-001,
      -4.2120746201784126E-001, -4.1806701224708759E-001,
      -4.1395178647347003E-001, -4.0893544178826940E-001,
      -4.0308807294221216E-001, -3.9647625592195934E-001,
      -3.8916309856034403E-001, -3.8120829775544990E-001,
      -3.7266820271937240E-001, -3.6359588366841639E-001,
      -3.5404120551003304E-001, -3.4405090624932305E-001,
      -3.3366867973388625E-001, -3.2293526180233850E-001,
      -3.1188851861787953E-001, -3.0056353756195420E-001,
      -2.8899272422487260E-001, -2.7720590780096982E-001,
      -2.6523044476416463E-001, -2.5309129855851598E-001,
      -2.4081109986342791E-001, -2.2841026936644959E-001,
      -2.1590729501529812E-001, -2.0331900258785626E-001,
      -1.9066028901355456E-001, -1.7794317104119103E-001,
      -1.6517678107685357E-001, -1.5237097629769486E-001,
      -1.3954146769211381E-001, -1.2670451694051832E-001,
      -1.1385259825731046E-001, -1.0094080765361217E-001,
      -8.7953049862493082E-002, -7.5036754306186504E-002,
      -6.2453491894992483E-002, -5.0050013809521440E-002,
      -3.6728876594930675E-002, -2.1570230536064366E-002,
      -7.1394753859541681E-003, 2.9895240250835913E-016, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0943951023932001E+000,
      -2.5701586533746604E-015, -3.2546616274269247E-015,
      3.8352143484655727E-015, 4.0240874122224388E-015, -1.6568083679919893E-016,
      3.2893122382051069E-015, 9.3961415130138102E-016, 2.0032965005122033E-015,
      -2.7932301093418985E-015, -2.0694358262954436E-015,
      -9.7228555379101986E-016, 8.1922713598803085E-017,
      -2.0927656540033175E-016, -2.3545635770191572E-015,
      -1.0006881253598346E-015, 2.9878186960571047E-015, 1.1837456474657640E-015,
      -1.6714393322179605E-015, 2.3844394659036018E-015,
      -1.1607239597141367E-015, -3.8229763492255061E-015,
      2.9466870618571244E-015, 3.5596451765774081E-015, 3.7789435723811561E-015,
      5.1141894736445666E-016, -3.8142349233177695E-015,
      -1.9927470011136630E-015, -2.2049987760449175E-001,
      -3.6635528907836157E-001, -3.4174172349631210E-001,
      -2.7380070768982373E-001, -2.3920333385705428E-001,
      -2.3049273999773329E-001, -2.2254837179744011E-001,
      -2.0788943550908590E-001, -1.9093801322761142E-001,
      -1.7556504586500291E-001, -1.6197223175698353E-001,
      -1.4911044246677455E-001, -1.3651416034412092E-001,
      -1.2430404856342320E-001, -1.1265540472997154E-001,
      -1.0159332818404462E-001, -9.1068464923279968E-002,
      -8.1045663994241693E-002, -7.1516123081945951E-002,
      -6.2476391278583294E-002, -5.3916851929118520E-002,
      -4.5823466303086993E-002, -3.8181854814058350E-002,
      -3.0978422198887731E-002, -2.4199664945182463E-002,
      -1.7831635322873497E-002, -1.1859990029954744E-002,
      -6.2702147142223597E-003, -1.0477431966515757E-003,
      3.8220178139750623E-003, 8.3536726218582727E-003, 1.2561799083010866E-002,
      1.6460903094468966E-002, 2.0065378740800430E-002, 2.3389475384230210E-002,
      2.6447268081017373E-002, 2.9252629446455607E-002, 3.1819203219575896E-002,
      3.4160380144311796E-002, 3.6289276203822898E-002, 3.8218712633535129E-002,
      3.9961197042840596E-002, 4.1528906061746182E-002, 4.2933671726190975E-002,
      4.4186972737834194E-002, 4.5299924223700418E-002, 4.6283253348326878E-002,
      4.7147265695612360E-002, 4.7901852147220168E-002, 4.8556584822594569E-002,
      4.9120794780353146E-002, 4.9603321987911821E-002, 5.0011897404603427E-002,
      5.0353169709766374E-002, 5.0634854297206731E-002, 5.0868471889453870E-002,
      5.1065559857348924E-002, 5.1223219116638019E-002, 5.1318034422326797E-002,
      5.1347803006383648E-002, 5.1407674732831976E-002, 5.1647162414796839E-002,
      5.1951031164477385E-002, 5.1665182225226688E-002, 5.0333049644776105E-002,
      4.9613912341886085E-002, 5.3284548858361327E-002, 6.0634584235465853E-002,
      5.7723020600442002E-002, 2.8557901543820052E-002, 6.5388600108723035E-017,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helikopter03_3_DWork.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues;
    helikopter03_3_DWork.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues;
    helikopter03_3_DWork.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    const char *fileName = "travel.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helikopter03_3_M, "Error creating .mat file travel.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"travel")) {
      rtmSetErrorStatus(helikopter03_3_M,
                        "Error writing mat file header to file travel.mat");
      return;
    }

    helikopter03_3_DWork.ToFile_IWORK.Count = 0;
    helikopter03_3_DWork.ToFile_IWORK.Decimation = -1;
    helikopter03_3_DWork.ToFile_PWORK.FilePtr = fp;
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
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75 } ;

    static real_T pDataValues[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      3.7024024484653367E-001, -2.7162992997560800E-015, 5.2359877559829671E-001,
      5.2359877559829870E-001, 5.2359877559829882E-001, 5.2359877559829937E-001,
      5.2359877559829793E-001, 5.2359877559829993E-001, 5.2359877559829981E-001,
      5.2359877559829937E-001, 5.2359877559830070E-001, 5.2359877559829904E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829959E-001,
      5.2359877559829793E-001, 5.2359877559829771E-001, 5.2359877559829737E-001,
      5.2359877559829671E-001, 5.2359877559829793E-001, 5.2359877559829826E-001,
      5.2359877559829826E-001, 5.2359877559829948E-001, 5.2359877559830059E-001,
      5.2359877559829504E-001, 5.2359877559829349E-001, 5.2359877559829593E-001,
      5.2359877559830048E-001, 4.8461953592206514E-001, 4.5883569829058313E-001,
      4.0806183367252036E-001, 3.2848339965229489E-001, 2.4916397819031194E-001,
      1.8225363126437860E-001, 1.2385717695678342E-001, 6.8825350270967131E-002,
      1.6184873733103913E-002, -3.3069902776905170E-002,
      -7.8401513326008462E-002, -1.2001911018593596E-001,
      -1.5828543899851308E-001, -1.9340458640541913E-001,
      -2.2547391772988584E-001, -2.5459441253583331E-001,
      -2.8089771317263612E-001, -3.0452424759280000E-001,
      -3.2560676307345021E-001, -3.4427016515744602E-001,
      -3.6063606884847754E-001, -3.7482444470339304E-001,
      -3.8695279885989470E-001, -3.9713526642316266E-001,
      -4.0548240382115214E-001, -4.1210129013868968E-001,
      -4.1709555865461934E-001, -4.2056532547757031E-001,
      -4.2260711172890358E-001, -4.2331380514911848E-001,
      -4.2277464998711872E-001, -4.2107524684443154E-001,
      -4.1829755796720264E-001, -4.1451992090341355E-001,
      -4.0981707231032344E-001, -4.0426018113728274E-001,
      -3.9791688977993450E-001, -3.9085136232983109E-001,
      -3.8312433944835478E-001, -3.7479319943343092E-001,
      -3.6591202500152320E-001, -3.5653167531508534E-001,
      -3.4669986273715397E-001, -3.3646123364461727E-001,
      -3.2585745271812905E-001, -3.1492729090088550E-001,
      -3.0370671804338456E-001, -2.9222899974001171E-001,
      -2.8052479310369827E-001, -2.6862223520114137E-001,
      -2.5654703299254944E-001, -2.4432258722175745E-001,
      -2.3197016191472775E-001, -2.1950900242744692E-001,
      -2.0695623280581102E-001, -1.9432664223244864E-001,
      -1.8163308809849968E-001, -1.6888809964329596E-001,
      -1.5610494854255536E-001, -1.4329388137148311E-001,
      -1.3045378883995895E-001, -1.1757450224738200E-001,
      -1.0466886665080181E-001, -9.1807607477937980E-002,
      -7.9055339682088069E-002, -6.6266170841651179E-002,
      -5.3034025437439264E-002, -3.9331232387690269E-002,
      -2.6524791770644388E-002, -1.6521859076540663E-002,
      -7.1394753859537457E-003, -7.1394753859537457E-003, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0 } ;

    helikopter03_3_DWork.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues;
    helikopter03_3_DWork.FromWorkspace_PWORK.DataPtr = (void *) pDataValues;
    helikopter03_3_DWork.FromWorkspace_IWORK.PrevIndex = 0;
  }

  MdlInitialize();
}

void MdlTerminate(void)
{
  helikopter03_3_terminate();
}

RT_MODEL_helikopter03_3 *helikopter03_3(void)
{
  helikopter03_3_initialize(1);
  return helikopter03_3_M;
}

/*========================================================================*
 * End of GRT compatible call interface                                   *
 *========================================================================*/
