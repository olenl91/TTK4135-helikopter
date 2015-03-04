/*
 * helikopter03.c
 *
 * Real-Time Workshop code generation for Simulink model "helikopter03.mdl".
 *
 * Model version              : 1.56
 * Real-Time Workshop version : 7.5  (R2010a)  25-Jan-2010
 * C source code generated on : Wed Mar 04 16:57:08 2015
 *
 * Target selection: quarc_windows.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "helikopter03.h"
#include "helikopter03_private.h"
#include <stdio.h>
#include "helikopter03_dt.h"

/* Block signals (auto storage) */
BlockIO_helikopter03 helikopter03_B;

/* Continuous states */
ContinuousStates_helikopter03 helikopter03_X;

/* Block states (auto storage) */
D_Work_helikopter03 helikopter03_DWork;

/* Real-time model */
RT_MODEL_helikopter03 helikopter03_M_;
RT_MODEL_helikopter03 *helikopter03_M = &helikopter03_M_;

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
  helikopter03_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helikopter03_output(int_T tid)
{
  /* local block i/o variables */
  real_T rtb_HILReadEncoder_o1;
  real_T rtb_HILReadEncoder_o2;
  real_T rtb_HILReadEncoder_o3;
  real_T rtb_VandringDeriv;
  real_T rtb_Gain2;
  real_T rtb_Saturation_b;
  real_T rtb_Gain1;
  real_T rtb_Gain[6];
  real_T tmp[6];
  int32_T tmp_0;
  int32_T tmp_1;
  if (rtmIsMajorTimeStep(helikopter03_M)) {
    /* set solver stop time */
    if (!(helikopter03_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helikopter03_M->solverInfo,
                            ((helikopter03_M->Timing.clockTickH0 + 1) *
        helikopter03_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helikopter03_M->solverInfo,
                            ((helikopter03_M->Timing.clockTick0 + 1) *
        helikopter03_M->Timing.stepSize0 + helikopter03_M->Timing.clockTickH0 *
        helikopter03_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helikopter03_M)) {
    helikopter03_M->Timing.t[0] = rtsiGetT(&helikopter03_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helikopter03_M)) {
    /* S-Function (hil_read_encoder_block): '<S2>/HIL Read Encoder' */

    /* S-Function Block: helikopter03/Heli 3D/HIL Read Encoder (hil_read_encoder_block) */
    {
      t_error result = hil_read_encoder(helikopter03_DWork.HILInitialize_Card,
        helikopter03_P.HILReadEncoder_Channels, 3,
        &helikopter03_DWork.HILReadEncoder_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_M, _rt_error_message);
      } else {
        rtb_HILReadEncoder_o1 = helikopter03_DWork.HILReadEncoder_Buffer[0];
        rtb_HILReadEncoder_o2 = helikopter03_DWork.HILReadEncoder_Buffer[1];
        rtb_HILReadEncoder_o3 = helikopter03_DWork.HILReadEncoder_Buffer[2];
      }
    }

    /* Gain: '<S2>/Kalibrer-Elev' */
    helikopter03_B.KalibrerElev = helikopter03_P.KalibrerElev_Gain *
      rtb_HILReadEncoder_o3;

    /* Sum: '<Root>/Add' incorporates:
     *  Constant: '<Root>/Constant'
     */
    helikopter03_B.Add = helikopter03_B.KalibrerElev +
      helikopter03_P.Constant_Value;

    /* Gain: '<S2>/Kalibrer-Pitch' */
    helikopter03_B.KalibrerPitch = helikopter03_P.KalibrerPitch_Gain *
      rtb_HILReadEncoder_o2;
  }

  /* TransferFcn: '<S2>/Vandring Lavpass' */
  helikopter03_B.VandringLavpass = helikopter03_P.VandringLavpass_C*
    helikopter03_X.VandringLavpass_CSTATE;

  /* Sum: '<Root>/Add1' incorporates:
   *  Constant: '<Root>/Constant1'
   */
  helikopter03_B.Add1 = helikopter03_P.Constant1_Value +
    helikopter03_B.VandringLavpass;
  if (rtmIsMajorTimeStep(helikopter03_M)) {
    /* ToFile: '<Root>/To File' */
    if (rtmIsMajorTimeStep(helikopter03_M)) {
      if (!(++helikopter03_DWork.ToFile_IWORK.Decimation % 1) &&
          (helikopter03_DWork.ToFile_IWORK.Count*2)+1 < 100000000 ) {
        FILE *fp = (FILE *) helikopter03_DWork.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          helikopter03_DWork.ToFile_IWORK.Decimation = 0;
          u[0] = helikopter03_M->Timing.t[1];
          u[1] = helikopter03_B.Add1;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(helikopter03_M,
                              "Error writing to MAT-file travel.mat");
            return;
          }

          if (((++helikopter03_DWork.ToFile_IWORK.Count)*2)+1 >= 100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file travel.mat.\n");
          }
        }
      }
    }

    /* ToFile: '<Root>/To File1' */
    if (rtmIsMajorTimeStep(helikopter03_M)) {
      if (!(++helikopter03_DWork.ToFile1_IWORK.Decimation % 1) &&
          (helikopter03_DWork.ToFile1_IWORK.Count*2)+1 < 100000000 ) {
        FILE *fp = (FILE *) helikopter03_DWork.ToFile1_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          helikopter03_DWork.ToFile1_IWORK.Decimation = 0;
          u[0] = helikopter03_M->Timing.t[1];
          u[1] = helikopter03_B.Add;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(helikopter03_M,
                              "Error writing to MAT-file elevation.mat");
            return;
          }

          if (((++helikopter03_DWork.ToFile1_IWORK.Count)*2)+1 >= 100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file elevation.mat.\n");
          }
        }
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helikopter03_DWork.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopter03_DWork.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helikopter03_DWork.FromWorkspace_IWORK.PrevIndex;
    real_T t = helikopter03_M->Timing.t[0];

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

    helikopter03_DWork.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          helikopter03_B.FromWorkspace = pDataValues[currTimeIndex];
        } else {
          helikopter03_B.FromWorkspace = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        helikopter03_B.FromWorkspace = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 140;
      }
    }
  }

  if (rtmIsMajorTimeStep(helikopter03_M)) {
  }

  /* Integrator: '<S1>/Integrator'
   *
   * Regarding '<S1>/Integrator':
   *  Limited Integrator
   */
  if (helikopter03_X.Integrator_CSTATE >= helikopter03_P.Integrator_UpperSat ) {
    helikopter03_X.Integrator_CSTATE = helikopter03_P.Integrator_UpperSat;
  } else if (helikopter03_X.Integrator_CSTATE <=
             helikopter03_P.Integrator_LowerSat ) {
    helikopter03_X.Integrator_CSTATE = helikopter03_P.Integrator_LowerSat;
  }

  rtb_Saturation_b = helikopter03_X.Integrator_CSTATE;
  if (rtmIsMajorTimeStep(helikopter03_M)) {
    /* Gain: '<S2>/Kalibrer -Vandring' */
    helikopter03_B.KalibrerVandring = helikopter03_P.KalibrerVandring_Gain *
      rtb_HILReadEncoder_o1;
  }

  /* TransferFcn: '<S2>/Vandring Deriv' */
  rtb_VandringDeriv = helikopter03_P.VandringDeriv_D*
    helikopter03_B.KalibrerVandring;
  rtb_VandringDeriv += helikopter03_P.VandringDeriv_C*
    helikopter03_X.VandringDeriv_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn4' */
  rtb_Gain2 = helikopter03_P.TransferFcn4_D*helikopter03_B.KalibrerPitch;
  rtb_Gain2 += helikopter03_P.TransferFcn4_C*helikopter03_X.TransferFcn4_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn5' */
  rtb_Gain1 = helikopter03_P.TransferFcn5_D*helikopter03_B.KalibrerElev;
  rtb_Gain1 += helikopter03_P.TransferFcn5_C*helikopter03_X.TransferFcn5_CSTATE;

  /* Gain: '<Root>/Gain' incorporates:
   *  SignalConversion: '<Root>/TmpSignal ConversionAtGainInport1'
   */
  tmp[0] = helikopter03_B.Add1;
  tmp[1] = rtb_VandringDeriv;
  tmp[2] = helikopter03_B.KalibrerPitch;
  tmp[3] = rtb_Gain2;
  tmp[4] = helikopter03_B.Add;
  tmp[5] = rtb_Gain1;
  for (tmp_0 = 0; tmp_0 < 6; tmp_0++) {
    rtb_Gain[tmp_0] = 0.0;
    for (tmp_1 = 0; tmp_1 < 6; tmp_1++) {
      rtb_Gain[tmp_0] += helikopter03_P.Gain_Gain[6 * tmp_1 + tmp_0] * tmp[tmp_1];
    }
  }

  /* Sum: '<S1>/Sum' incorporates:
   *  Constant: '<Root>/elevation'
   */
  rtb_Gain1 = helikopter03_P.elevation_Value - rtb_Gain[4];

  /* Gain: '<S1>/K_ei' */
  helikopter03_B.K_ei = helikopter03_P.K_ei_Gain * rtb_Gain1;

  /* Sum: '<S1>/Sum1' incorporates:
   *  Gain: '<S1>/K_ed'
   *  Gain: '<S1>/K_ep'
   */
  rtb_Saturation_b = (helikopter03_P.K_ep_Gain * rtb_Gain1 + rtb_Saturation_b) +
    helikopter03_P.K_ed_Gain * rtb_Gain[5];

  /* Saturate: '<S1>/Saturation' */
  rtb_Saturation_b = rt_SATURATE(rtb_Saturation_b,
    helikopter03_P.Saturation_LowerSat, helikopter03_P.Saturation_UpperSat);
  if (rtmIsMajorTimeStep(helikopter03_M)) {
  }

  /* Sum: '<S3>/Sum' incorporates:
   *  Gain: '<S3>/K_pd'
   *  Gain: '<S3>/K_pp'
   *  Saturate: '<S3>/Saturation'
   *  Sum: '<S3>/Sum1'
   */
  rtb_Gain1 = (rt_SATURATE(helikopter03_B.FromWorkspace,
    helikopter03_P.Saturation_LowerSat_p, helikopter03_P.Saturation_UpperSat_g)
               - rtb_Gain[2]) * helikopter03_P.K_pp_Gain -
    helikopter03_P.K_pd_Gain * rtb_Gain[3];

  /* Gain: '<S4>/Gain2' incorporates:
   *  Sum: '<S4>/Sum4'
   */
  rtb_Gain2 = (rtb_Saturation_b - rtb_Gain1) * helikopter03_P.Gain2_Gain;

  /* Saturate: '<S2>/Sat B' */
  helikopter03_B.SatB = rt_SATURATE(rtb_Gain2, helikopter03_P.SatB_LowerSat,
    helikopter03_P.SatB_UpperSat);
  if (rtmIsMajorTimeStep(helikopter03_M)) {
  }

  /* Gain: '<S4>/Gain1' incorporates:
   *  Sum: '<S4>/Sum3'
   */
  rtb_Gain1 = (rtb_Gain1 + rtb_Saturation_b) * helikopter03_P.Gain1_Gain;

  /* Saturate: '<S2>/Sat' */
  helikopter03_B.Sat = rt_SATURATE(rtb_Gain1, helikopter03_P.Sat_LowerSat,
    helikopter03_P.Sat_UpperSat);
  if (rtmIsMajorTimeStep(helikopter03_M)) {
    /* S-Function (hil_write_analog_block): '<S2>/HIL Write Analog' */

    /* S-Function Block: helikopter03/Heli 3D/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helikopter03_DWork.HILWriteAnalog_Buffer[0] = helikopter03_B.SatB;
      helikopter03_DWork.HILWriteAnalog_Buffer[1] = helikopter03_B.Sat;
      result = hil_write_analog(helikopter03_DWork.HILInitialize_Card,
        helikopter03_P.HILWriteAnalog_Channels, 2,
        &helikopter03_DWork.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_M, _rt_error_message);
      }
    }
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Model update function */
void helikopter03_update(int_T tid)
{
  if (rtmIsMajorTimeStep(helikopter03_M)) {
    rt_ertODEUpdateContinuousStates(&helikopter03_M->solverInfo);
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
  if (!(++helikopter03_M->Timing.clockTick0)) {
    ++helikopter03_M->Timing.clockTickH0;
  }

  helikopter03_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helikopter03_M->solverInfo);
  if (rtmIsMajorTimeStep(helikopter03_M)) {
    /* Update absolute timer for sample time: [0.001s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helikopter03_M->Timing.clockTick1)) {
      ++helikopter03_M->Timing.clockTickH1;
    }

    helikopter03_M->Timing.t[1] = helikopter03_M->Timing.clockTick1 *
      helikopter03_M->Timing.stepSize1 + helikopter03_M->Timing.clockTickH1 *
      helikopter03_M->Timing.stepSize1 * 4294967296.0;
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Derivatives for root system: '<Root>' */
void helikopter03_derivatives(void)
{
  /* Derivatives for TransferFcn: '<S2>/Vandring Lavpass' */
  {
    ((StateDerivatives_helikopter03 *) helikopter03_M->ModelData.derivs)
      ->VandringLavpass_CSTATE = helikopter03_B.KalibrerVandring;
    ((StateDerivatives_helikopter03 *) helikopter03_M->ModelData.derivs)
      ->VandringLavpass_CSTATE += (helikopter03_P.VandringLavpass_A)*
      helikopter03_X.VandringLavpass_CSTATE;
  }

  /* Derivatives for Integrator: '<S1>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helikopter03_X.Integrator_CSTATE <=
            helikopter03_P.Integrator_LowerSat );
    usat = ( helikopter03_X.Integrator_CSTATE >=
            helikopter03_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helikopter03_B.K_ei > 0)) ||
        (usat && (helikopter03_B.K_ei < 0)) ) {
      ((StateDerivatives_helikopter03 *) helikopter03_M->ModelData.derivs)
        ->Integrator_CSTATE = helikopter03_B.K_ei;
    } else {
      /* in saturation */
      ((StateDerivatives_helikopter03 *) helikopter03_M->ModelData.derivs)
        ->Integrator_CSTATE = 0.0;
    }
  }

  /* Derivatives for TransferFcn: '<S2>/Vandring Deriv' */
  {
    ((StateDerivatives_helikopter03 *) helikopter03_M->ModelData.derivs)
      ->VandringDeriv_CSTATE = helikopter03_B.KalibrerVandring;
    ((StateDerivatives_helikopter03 *) helikopter03_M->ModelData.derivs)
      ->VandringDeriv_CSTATE += (helikopter03_P.VandringDeriv_A)*
      helikopter03_X.VandringDeriv_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn4' */
  {
    ((StateDerivatives_helikopter03 *) helikopter03_M->ModelData.derivs)
      ->TransferFcn4_CSTATE = helikopter03_B.KalibrerPitch;
    ((StateDerivatives_helikopter03 *) helikopter03_M->ModelData.derivs)
      ->TransferFcn4_CSTATE += (helikopter03_P.TransferFcn4_A)*
      helikopter03_X.TransferFcn4_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn5' */
  {
    ((StateDerivatives_helikopter03 *) helikopter03_M->ModelData.derivs)
      ->TransferFcn5_CSTATE = helikopter03_B.KalibrerElev;
    ((StateDerivatives_helikopter03 *) helikopter03_M->ModelData.derivs)
      ->TransferFcn5_CSTATE += (helikopter03_P.TransferFcn5_A)*
      helikopter03_X.TransferFcn5_CSTATE;
  }
}

/* Model initialize function */
void helikopter03_initialize(boolean_T firstTime)
{
  (void)firstTime;

  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helikopter03_P.Integrator_UpperSat = rtInf;
  helikopter03_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helikopter03_M, 0,
                sizeof(RT_MODEL_helikopter03));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helikopter03_M->solverInfo,
                          &helikopter03_M->Timing.simTimeStep);
    rtsiSetTPtr(&helikopter03_M->solverInfo, &rtmGetTPtr(helikopter03_M));
    rtsiSetStepSizePtr(&helikopter03_M->solverInfo,
                       &helikopter03_M->Timing.stepSize0);
    rtsiSetdXPtr(&helikopter03_M->solverInfo, &helikopter03_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helikopter03_M->solverInfo,
                         &helikopter03_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helikopter03_M->solverInfo,
      &helikopter03_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helikopter03_M->solverInfo, (&rtmGetErrorStatus
      (helikopter03_M)));
    rtsiSetRTModelPtr(&helikopter03_M->solverInfo, helikopter03_M);
  }

  rtsiSetSimTimeStep(&helikopter03_M->solverInfo, MAJOR_TIME_STEP);
  helikopter03_M->ModelData.intgData.f[0] = helikopter03_M->ModelData.odeF[0];
  helikopter03_M->ModelData.contStates = ((real_T *) &helikopter03_X);
  rtsiSetSolverData(&helikopter03_M->solverInfo, (void *)
                    &helikopter03_M->ModelData.intgData);
  rtsiSetSolverName(&helikopter03_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helikopter03_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helikopter03_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helikopter03_M->Timing.sampleTimes =
      (&helikopter03_M->Timing.sampleTimesArray[0]);
    helikopter03_M->Timing.offsetTimes =
      (&helikopter03_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helikopter03_M->Timing.sampleTimes[0] = (0.0);
    helikopter03_M->Timing.sampleTimes[1] = (0.001);

    /* task offsets */
    helikopter03_M->Timing.offsetTimes[0] = (0.0);
    helikopter03_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helikopter03_M, &helikopter03_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helikopter03_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helikopter03_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helikopter03_M, -1);
  helikopter03_M->Timing.stepSize0 = 0.001;
  helikopter03_M->Timing.stepSize1 = 0.001;

  /* external mode info */
  helikopter03_M->Sizes.checksums[0] = (3583586733U);
  helikopter03_M->Sizes.checksums[1] = (1339983735U);
  helikopter03_M->Sizes.checksums[2] = (945622544U);
  helikopter03_M->Sizes.checksums[3] = (3654633784U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helikopter03_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helikopter03_M->extModeInfo,
      &helikopter03_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helikopter03_M->extModeInfo,
                        helikopter03_M->Sizes.checksums);
    rteiSetTPtr(helikopter03_M->extModeInfo, rtmGetTPtr(helikopter03_M));
  }

  helikopter03_M->solverInfoPtr = (&helikopter03_M->solverInfo);
  helikopter03_M->Timing.stepSize = (0.001);
  rtsiSetFixedStepSize(&helikopter03_M->solverInfo, 0.001);
  rtsiSetSolverMode(&helikopter03_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helikopter03_M->ModelData.blockIO = ((void *) &helikopter03_B);

  {
    helikopter03_B.KalibrerElev = 0.0;
    helikopter03_B.Add = 0.0;
    helikopter03_B.KalibrerPitch = 0.0;
    helikopter03_B.VandringLavpass = 0.0;
    helikopter03_B.Add1 = 0.0;
    helikopter03_B.FromWorkspace = 0.0;
    helikopter03_B.KalibrerVandring = 0.0;
    helikopter03_B.K_ei = 0.0;
    helikopter03_B.SatB = 0.0;
    helikopter03_B.Sat = 0.0;
  }

  /* parameters */
  helikopter03_M->ModelData.defaultParam = ((real_T *)&helikopter03_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helikopter03_X;
    helikopter03_M->ModelData.contStates = (x);
    (void) memset((void *)&helikopter03_X, 0,
                  sizeof(ContinuousStates_helikopter03));
  }

  /* states (dwork) */
  helikopter03_M->Work.dwork = ((void *) &helikopter03_DWork);
  (void) memset((void *)&helikopter03_DWork, 0,
                sizeof(D_Work_helikopter03));
  helikopter03_DWork.HILInitialize_AIMinimums[0] = 0.0;
  helikopter03_DWork.HILInitialize_AIMinimums[1] = 0.0;
  helikopter03_DWork.HILInitialize_AIMinimums[2] = 0.0;
  helikopter03_DWork.HILInitialize_AIMinimums[3] = 0.0;
  helikopter03_DWork.HILInitialize_AIMaximums[0] = 0.0;
  helikopter03_DWork.HILInitialize_AIMaximums[1] = 0.0;
  helikopter03_DWork.HILInitialize_AIMaximums[2] = 0.0;
  helikopter03_DWork.HILInitialize_AIMaximums[3] = 0.0;
  helikopter03_DWork.HILInitialize_AOMinimums[0] = 0.0;
  helikopter03_DWork.HILInitialize_AOMinimums[1] = 0.0;
  helikopter03_DWork.HILInitialize_AOMinimums[2] = 0.0;
  helikopter03_DWork.HILInitialize_AOMinimums[3] = 0.0;
  helikopter03_DWork.HILInitialize_AOMaximums[0] = 0.0;
  helikopter03_DWork.HILInitialize_AOMaximums[1] = 0.0;
  helikopter03_DWork.HILInitialize_AOMaximums[2] = 0.0;
  helikopter03_DWork.HILInitialize_AOMaximums[3] = 0.0;
  helikopter03_DWork.HILInitialize_AOVoltages[0] = 0.0;
  helikopter03_DWork.HILInitialize_AOVoltages[1] = 0.0;
  helikopter03_DWork.HILInitialize_AOVoltages[2] = 0.0;
  helikopter03_DWork.HILInitialize_AOVoltages[3] = 0.0;
  helikopter03_DWork.HILInitialize_FilterFrequency[0] = 0.0;
  helikopter03_DWork.HILInitialize_FilterFrequency[1] = 0.0;
  helikopter03_DWork.HILInitialize_FilterFrequency[2] = 0.0;
  helikopter03_DWork.HILInitialize_FilterFrequency[3] = 0.0;
  helikopter03_DWork.HILWriteAnalog_Buffer[0] = 0.0;
  helikopter03_DWork.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helikopter03_M->SpecialInfo.mappingInfo = (&dtInfo);
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
void helikopter03_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopter03/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    hil_task_stop_all(helikopter03_DWork.HILInitialize_Card);
    hil_task_delete_all(helikopter03_DWork.HILInitialize_Card);
    hil_monitor_stop_all(helikopter03_DWork.HILInitialize_Card);
    hil_monitor_delete_all(helikopter03_DWork.HILInitialize_Card);
    is_switching = false;
    if ((helikopter03_P.HILInitialize_AOTerminate && !is_switching) ||
        (helikopter03_P.HILInitialize_AOExit && is_switching)) {
      helikopter03_DWork.HILInitialize_AOVoltages[0] =
        helikopter03_P.HILInitialize_AOFinal;
      helikopter03_DWork.HILInitialize_AOVoltages[1] =
        helikopter03_P.HILInitialize_AOFinal;
      helikopter03_DWork.HILInitialize_AOVoltages[2] =
        helikopter03_P.HILInitialize_AOFinal;
      helikopter03_DWork.HILInitialize_AOVoltages[3] =
        helikopter03_P.HILInitialize_AOFinal;
      result = hil_write_analog(helikopter03_DWork.HILInitialize_Card,
        &helikopter03_P.HILInitialize_AOChannels[0], 4U,
        &helikopter03_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_M, _rt_error_message);
      }
    }

    hil_close(helikopter03_DWork.HILInitialize_Card);
    helikopter03_DWork.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helikopter03_DWork.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      const char *fileName = "travel.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopter03_M, "Error closing MAT-file travel.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helikopter03_M, "Error reopening MAT-file travel.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2, helikopter03_DWork.ToFile_IWORK.Count,
           "travel")) {
        rtmSetErrorStatus(helikopter03_M,
                          "Error writing header for travel to MAT-file travel.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopter03_M, "Error closing MAT-file travel.mat");
        return;
      }

      helikopter03_DWork.ToFile_PWORK.FilePtr = (NULL);
    }
  }

  /* Terminate for ToFile: '<Root>/To File1' */
  {
    FILE *fp = (FILE *) helikopter03_DWork.ToFile1_PWORK.FilePtr;
    if (fp != (NULL)) {
      const char *fileName = "elevation.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopter03_M, "Error closing MAT-file elevation.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helikopter03_M,
                          "Error reopening MAT-file elevation.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2, helikopter03_DWork.ToFile1_IWORK.Count,
           "elevation")) {
        rtmSetErrorStatus(helikopter03_M,
                          "Error writing header for elevation to MAT-file elevation.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopter03_M, "Error closing MAT-file elevation.mat");
        return;
      }

      helikopter03_DWork.ToFile1_PWORK.FilePtr = (NULL);
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
  helikopter03_output(tid);
}

void MdlUpdate(int_T tid)
{
  helikopter03_update(tid);
}

void MdlInitializeSizes(void)
{
  helikopter03_M->Sizes.numContStates = (5);/* Number of continuous states */
  helikopter03_M->Sizes.numY = (0);    /* Number of model outputs */
  helikopter03_M->Sizes.numU = (0);    /* Number of model inputs */
  helikopter03_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helikopter03_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helikopter03_M->Sizes.numBlocks = (47);/* Number of blocks */
  helikopter03_M->Sizes.numBlockIO = (10);/* Number of block outputs */
  helikopter03_M->Sizes.numBlockPrms = (147);/* Sum of parameter "widths" */
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
  /* InitializeConditions for TransferFcn: '<S2>/Vandring Lavpass' */
  helikopter03_X.VandringLavpass_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S1>/Integrator' */
  helikopter03_X.Integrator_CSTATE = helikopter03_P.Integrator_IC;

  /* InitializeConditions for TransferFcn: '<S2>/Vandring Deriv' */
  helikopter03_X.VandringDeriv_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn4' */
  helikopter03_X.TransferFcn4_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn5' */
  helikopter03_X.TransferFcn5_CSTATE = 0.0;
}

void MdlStart(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopter03/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q4", "0", &helikopter03_DWork.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter03_M, _rt_error_message);
      return;
    }

    is_switching = false;
    if ((helikopter03_P.HILInitialize_CKPStart && !is_switching) ||
        (helikopter03_P.HILInitialize_CKPEnter && is_switching)) {
      result = hil_set_clock_mode(helikopter03_DWork.HILInitialize_Card,
        (t_clock *) &helikopter03_P.HILInitialize_CKChannels[0], 2U,
        (t_clock_mode *) &helikopter03_P.HILInitialize_CKModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_M, _rt_error_message);
        return;
      }
    }

    result = hil_watchdog_clear(helikopter03_DWork.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter03_M, _rt_error_message);
      return;
    }

    if ((helikopter03_P.HILInitialize_AIPStart && !is_switching) ||
        (helikopter03_P.HILInitialize_AIPEnter && is_switching)) {
      helikopter03_DWork.HILInitialize_AIMinimums[0] =
        helikopter03_P.HILInitialize_AILow;
      helikopter03_DWork.HILInitialize_AIMinimums[1] =
        helikopter03_P.HILInitialize_AILow;
      helikopter03_DWork.HILInitialize_AIMinimums[2] =
        helikopter03_P.HILInitialize_AILow;
      helikopter03_DWork.HILInitialize_AIMinimums[3] =
        helikopter03_P.HILInitialize_AILow;
      helikopter03_DWork.HILInitialize_AIMaximums[0] =
        helikopter03_P.HILInitialize_AIHigh;
      helikopter03_DWork.HILInitialize_AIMaximums[1] =
        helikopter03_P.HILInitialize_AIHigh;
      helikopter03_DWork.HILInitialize_AIMaximums[2] =
        helikopter03_P.HILInitialize_AIHigh;
      helikopter03_DWork.HILInitialize_AIMaximums[3] =
        helikopter03_P.HILInitialize_AIHigh;
      result = hil_set_analog_input_ranges(helikopter03_DWork.HILInitialize_Card,
        &helikopter03_P.HILInitialize_AIChannels[0], 4U,
        &helikopter03_DWork.HILInitialize_AIMinimums[0],
        &helikopter03_DWork.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter03_P.HILInitialize_AOPStart && !is_switching) ||
        (helikopter03_P.HILInitialize_AOPEnter && is_switching)) {
      helikopter03_DWork.HILInitialize_AOMinimums[0] =
        helikopter03_P.HILInitialize_AOLow;
      helikopter03_DWork.HILInitialize_AOMinimums[1] =
        helikopter03_P.HILInitialize_AOLow;
      helikopter03_DWork.HILInitialize_AOMinimums[2] =
        helikopter03_P.HILInitialize_AOLow;
      helikopter03_DWork.HILInitialize_AOMinimums[3] =
        helikopter03_P.HILInitialize_AOLow;
      helikopter03_DWork.HILInitialize_AOMaximums[0] =
        helikopter03_P.HILInitialize_AOHigh;
      helikopter03_DWork.HILInitialize_AOMaximums[1] =
        helikopter03_P.HILInitialize_AOHigh;
      helikopter03_DWork.HILInitialize_AOMaximums[2] =
        helikopter03_P.HILInitialize_AOHigh;
      helikopter03_DWork.HILInitialize_AOMaximums[3] =
        helikopter03_P.HILInitialize_AOHigh;
      result = hil_set_analog_output_ranges
        (helikopter03_DWork.HILInitialize_Card,
         &helikopter03_P.HILInitialize_AOChannels[0], 4U,
         &helikopter03_DWork.HILInitialize_AOMinimums[0],
         &helikopter03_DWork.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter03_P.HILInitialize_AOStart && !is_switching) ||
        (helikopter03_P.HILInitialize_AOEnter && is_switching)) {
      helikopter03_DWork.HILInitialize_AOVoltages[0] =
        helikopter03_P.HILInitialize_AOInitial;
      helikopter03_DWork.HILInitialize_AOVoltages[1] =
        helikopter03_P.HILInitialize_AOInitial;
      helikopter03_DWork.HILInitialize_AOVoltages[2] =
        helikopter03_P.HILInitialize_AOInitial;
      helikopter03_DWork.HILInitialize_AOVoltages[3] =
        helikopter03_P.HILInitialize_AOInitial;
      result = hil_write_analog(helikopter03_DWork.HILInitialize_Card,
        &helikopter03_P.HILInitialize_AOChannels[0], 4U,
        &helikopter03_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_M, _rt_error_message);
        return;
      }
    }

    if (helikopter03_P.HILInitialize_AOReset) {
      helikopter03_DWork.HILInitialize_AOVoltages[0] =
        helikopter03_P.HILInitialize_AOWatchdog;
      helikopter03_DWork.HILInitialize_AOVoltages[1] =
        helikopter03_P.HILInitialize_AOWatchdog;
      helikopter03_DWork.HILInitialize_AOVoltages[2] =
        helikopter03_P.HILInitialize_AOWatchdog;
      helikopter03_DWork.HILInitialize_AOVoltages[3] =
        helikopter03_P.HILInitialize_AOWatchdog;
      result = hil_watchdog_set_analog_expiration_state
        (helikopter03_DWork.HILInitialize_Card,
         &helikopter03_P.HILInitialize_AOChannels[0], 4U,
         &helikopter03_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter03_P.HILInitialize_EIPStart && !is_switching) ||
        (helikopter03_P.HILInitialize_EIPEnter && is_switching)) {
      helikopter03_DWork.HILInitialize_QuadratureModes[0] =
        helikopter03_P.HILInitialize_EIQuadrature;
      helikopter03_DWork.HILInitialize_QuadratureModes[1] =
        helikopter03_P.HILInitialize_EIQuadrature;
      helikopter03_DWork.HILInitialize_QuadratureModes[2] =
        helikopter03_P.HILInitialize_EIQuadrature;
      helikopter03_DWork.HILInitialize_QuadratureModes[3] =
        helikopter03_P.HILInitialize_EIQuadrature;
      result = hil_set_encoder_quadrature_mode
        (helikopter03_DWork.HILInitialize_Card,
         &helikopter03_P.HILInitialize_EIChannels[0], 4U,
         (t_encoder_quadrature_mode *)
         &helikopter03_DWork.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_M, _rt_error_message);
        return;
      }

      helikopter03_DWork.HILInitialize_FilterFrequency[0] =
        helikopter03_P.HILInitialize_EIFrequency;
      helikopter03_DWork.HILInitialize_FilterFrequency[1] =
        helikopter03_P.HILInitialize_EIFrequency;
      helikopter03_DWork.HILInitialize_FilterFrequency[2] =
        helikopter03_P.HILInitialize_EIFrequency;
      helikopter03_DWork.HILInitialize_FilterFrequency[3] =
        helikopter03_P.HILInitialize_EIFrequency;
      result = hil_set_encoder_filter_frequency
        (helikopter03_DWork.HILInitialize_Card,
         &helikopter03_P.HILInitialize_EIChannels[0], 4U,
         &helikopter03_DWork.HILInitialize_FilterFrequency[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter03_P.HILInitialize_EIStart && !is_switching) ||
        (helikopter03_P.HILInitialize_EIEnter && is_switching)) {
      helikopter03_DWork.HILInitialize_InitialEICounts[0] =
        helikopter03_P.HILInitialize_EIInitial;
      helikopter03_DWork.HILInitialize_InitialEICounts[1] =
        helikopter03_P.HILInitialize_EIInitial;
      helikopter03_DWork.HILInitialize_InitialEICounts[2] =
        helikopter03_P.HILInitialize_EIInitial;
      helikopter03_DWork.HILInitialize_InitialEICounts[3] =
        helikopter03_P.HILInitialize_EIInitial;
      result = hil_set_encoder_counts(helikopter03_DWork.HILInitialize_Card,
        &helikopter03_P.HILInitialize_EIChannels[0], 4U,
        &helikopter03_DWork.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter03_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    const char *fileName = "travel.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helikopter03_M, "Error creating .mat file travel.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"travel")) {
      rtmSetErrorStatus(helikopter03_M,
                        "Error writing mat file header to file travel.mat");
      return;
    }

    helikopter03_DWork.ToFile_IWORK.Count = 0;
    helikopter03_DWork.ToFile_IWORK.Decimation = -1;
    helikopter03_DWork.ToFile_PWORK.FilePtr = fp;
  }

  /* Start for ToFile: '<Root>/To File1' */
  {
    const char *fileName = "elevation.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helikopter03_M, "Error creating .mat file elevation.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"elevation")) {
      rtmSetErrorStatus(helikopter03_M,
                        "Error writing mat file header to file elevation.mat");
      return;
    }

    helikopter03_DWork.ToFile1_IWORK.Count = 0;
    helikopter03_DWork.ToFile1_IWORK.Decimation = -1;
    helikopter03_DWork.ToFile1_PWORK.FilePtr = fp;
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

    helikopter03_DWork.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues;
    helikopter03_DWork.FromWorkspace_PWORK.DataPtr = (void *) pDataValues;
    helikopter03_DWork.FromWorkspace_IWORK.PrevIndex = 0;
  }

  MdlInitialize();
}

void MdlTerminate(void)
{
  helikopter03_terminate();
}

RT_MODEL_helikopter03 *helikopter03(void)
{
  helikopter03_initialize(1);
  return helikopter03_M;
}

/*========================================================================*
 * End of GRT compatible call interface                                   *
 *========================================================================*/
