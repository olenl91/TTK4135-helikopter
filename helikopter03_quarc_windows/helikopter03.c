/*
 * helikopter03.c
 *
 * Real-Time Workshop code generation for Simulink model "helikopter03.mdl".
 *
 * Model version              : 1.50
 * Real-Time Workshop version : 7.5  (R2010a)  25-Jan-2010
 * C source code generated on : Wed Feb 25 17:50:09 2015
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
  if (rtmIsMajorTimeStep(helikopter03_M)) {
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
    } else if (t >= pTimeValues[300]) {
      currTimeIndex = 299;
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
        pDataValues += 301;
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
  tmp[0] = helikopter03_B.VandringLavpass;
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
  helikopter03_M->Sizes.checksums[0] = (1031523503U);
  helikopter03_M->Sizes.checksums[1] = (4033324845U);
  helikopter03_M->Sizes.checksums[2] = (2788055736U);
  helikopter03_M->Sizes.checksums[3] = (597412017U);

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
  helikopter03_M->Sizes.numBlocks = (45);/* Number of blocks */
  helikopter03_M->Sizes.numBlockIO = (9);/* Number of block outputs */
  helikopter03_M->Sizes.numBlockPrms = (146);/* Sum of parameter "widths" */
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

  /* Start for FromWorkspace: '<Root>/From Workspace' */
  {
    static real_T pTimeValues[] = { 0.0, 0.1, 0.2, 3.0000000000000004E-001, 0.4,
      0.5, 6.0000000000000009E-001, 7.0000000000000007E-001, 0.8, 0.9, 1.0, 1.1,
      1.2000000000000002E+000, 1.3, 1.4000000000000001E+000, 1.5, 1.6,
      1.7000000000000002E+000, 1.8, 1.9000000000000001E+000, 2.0, 2.1, 2.2,
      2.3000000000000003E+000, 2.4000000000000004E+000, 2.5, 2.6, 2.7,
      2.8000000000000003E+000, 2.9000000000000004E+000, 3.0, 3.1, 3.2,
      3.3000000000000003E+000, 3.4000000000000004E+000, 3.5, 3.6, 3.7,
      3.8000000000000003E+000, 3.9000000000000004E+000, 4.0,
      4.1000000000000005E+000, 4.2, 4.3, 4.4, 4.5, 4.6000000000000005E+000, 4.7,
      4.8000000000000007E+000, 4.9, 5.0, 5.1000000000000005E+000, 5.2,
      5.3000000000000007E+000, 5.4, 5.5, 5.6000000000000005E+000, 5.7,
      5.8000000000000007E+000, 5.9, 6.0, 6.1000000000000005E+000, 6.2,
      6.3000000000000007E+000, 6.4, 6.5, 6.6000000000000005E+000, 6.7,
      6.8000000000000007E+000, 6.9, 7.0, 7.1000000000000005E+000, 7.2,
      7.3000000000000007E+000, 7.4, 7.5, 7.6000000000000005E+000, 7.7,
      7.8000000000000007E+000, 7.9, 8.0, 8.1, 8.2000000000000011E+000, 8.3, 8.4,
      8.5, 8.6, 8.7000000000000011E+000, 8.8, 8.9, 9.0, 9.1,
      9.2000000000000011E+000, 9.3, 9.4, 9.5, 9.6000000000000014E+000,
      9.7000000000000011E+000, 9.8, 9.9, 10.0, 1.0100000000000001E+001,
      1.0200000000000001E+001, 10.3, 10.4, 10.5, 1.0600000000000001E+001,
      1.0700000000000001E+001, 10.8, 10.9, 11.0, 1.1100000000000001E+001,
      1.1200000000000001E+001, 11.3, 11.4, 11.5, 1.1600000000000001E+001,
      1.1700000000000001E+001, 11.8, 11.9, 12.0, 1.2100000000000001E+001,
      1.2200000000000001E+001, 12.3, 12.4, 12.5, 1.2600000000000001E+001,
      1.2700000000000001E+001, 12.8, 12.9, 13.0, 1.3100000000000001E+001,
      1.3200000000000001E+001, 13.3, 13.4, 13.5, 1.3600000000000001E+001,
      1.3700000000000001E+001, 13.8, 13.9, 14.0, 1.4100000000000001E+001,
      1.4200000000000001E+001, 14.3, 14.4, 14.5, 1.4600000000000001E+001,
      1.4700000000000001E+001, 14.8, 14.9, 15.0, 15.1, 15.2,
      1.5299999999999999E+001, 1.5399999999999999E+001, 15.5, 15.6, 15.7,
      1.5799999999999999E+001, 1.5899999999999999E+001, 16.0, 16.1, 16.2,
      1.6299999999999997E+001, 16.4, 16.5, 16.6, 16.7, 1.6799999999999997E+001,
      16.9, 17.0, 17.1, 17.2, 1.7299999999999997E+001, 17.4, 17.5, 17.6, 17.7,
      1.7799999999999997E+001, 17.9, 18.0, 18.1, 18.2, 1.8299999999999997E+001,
      18.4, 18.5, 18.6, 18.7, 1.8799999999999997E+001, 18.9, 19.0, 19.1, 19.2,
      1.9299999999999997E+001, 19.4, 19.5, 19.6, 19.7, 1.9799999999999997E+001,
      19.9, 20.0, 20.1, 20.2, 2.0299999999999997E+001, 20.4, 20.5, 20.6, 20.7,
      2.0799999999999997E+001, 20.9, 21.0, 21.1, 21.2, 2.1299999999999997E+001,
      21.4, 21.5, 21.6, 21.7, 2.1799999999999997E+001, 21.9, 22.0, 22.1, 22.2,
      22.3, 22.4, 22.5, 22.6, 22.7, 22.8, 22.9, 23.0, 23.1, 23.2, 23.3, 23.4,
      23.5, 23.6, 23.7, 23.8, 23.9, 24.0, 24.1, 24.2, 24.3, 24.4, 24.5, 24.6,
      24.7, 24.8, 24.9, 25.0, 25.1, 25.2, 25.3, 25.4, 25.5, 25.6, 25.7, 25.8,
      25.9, 26.0, 26.1, 26.2, 26.3, 26.4, 26.5, 26.6, 26.7, 26.8, 26.9, 27.0,
      27.1, 27.2, 27.3, 27.4, 27.5, 27.6, 27.7, 27.8, 27.9, 28.0, 28.1, 28.2,
      28.3, 28.4, 28.5, 28.6, 28.7, 28.8, 28.9, 29.0, 29.1, 29.2, 29.3, 29.4,
      29.5, 29.6, 29.7, 29.8, 29.9, 30.0 } ;

    static real_T pDataValues[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      2.3140015302908101E+000, -1.3884009181744863E+000, 5.2359877559830392E-001,
      5.2359877559830370E-001, 5.2359877559829127E-001, 5.2359877559829870E-001,
      5.2359877559830059E-001, 5.2359877559830026E-001, 5.2359877559830004E-001,
      5.2359877559830292E-001, 5.2359877559829238E-001, 5.2359877559830081E-001,
      5.2359877559829981E-001, 5.2359877559829715E-001, 5.2359877559830093E-001,
      5.2359877559829804E-001, 5.2359877559829526E-001, 5.2359877559829460E-001,
      5.2359877559830126E-001, 5.2359877559830448E-001, 5.2359877559829993E-001,
      5.2359877559829382E-001, 5.2359877559829882E-001, 5.2359877559830004E-001,
      5.2359877559829826E-001, 5.2359877559830215E-001, 5.2359877559830081E-001,
      5.2359877559830159E-001, 5.2359877559828605E-001, 5.2359877559830637E-001,
      5.2359877559830426E-001, 5.2359877559829271E-001, 5.2359877559830159E-001,
      5.2359877559830104E-001, 5.2359877559830059E-001, 5.2359877559829737E-001,
      5.2359877559829548E-001, 5.2359877559830115E-001, 5.2359877559829748E-001,
      5.2359877559830037E-001, 5.2359877559829726E-001, 5.2359877559829249E-001,
      5.2359877559829970E-001, 5.2359877559830492E-001, 5.2359877559829926E-001,
      5.2359877559829227E-001, 5.2359877559829548E-001, 5.2359877559830437E-001,
      5.2359877559829526E-001, 5.2359877559829704E-001, 5.2359877559829704E-001,
      5.2359877559830004E-001, 5.2359877559830137E-001, 5.2359877559830093E-001,
      5.2359877559830093E-001, 5.2359877559830303E-001, 5.2359877559829493E-001,
      5.2359877559829826E-001, 4.6351290144514401E-001, 4.8701100421075805E-001,
      4.8168782905231300E-001, 4.6554159552225116E-001, 4.4956042052366441E-001,
      4.3374289016473905E-001, 4.1808757078189968E-001, 4.0259300963339983E-001,
      3.8725773558843041E-001, 3.7208025980311216E-001, 3.5705907641049561E-001,
      3.4219266316557728E-001, 3.2747948212716693E-001, 3.1291798029961904E-001,
      2.9850659028008103E-001, 2.8424373091054395E-001, 2.7012780789064583E-001,
      2.5615721443080375E-001, 2.4233033185016212E-001, 2.2864553020564388E-001,
      2.1510116889381842E-001, 2.0169559725453559E-001, 1.8842715516502023E-001,
      1.7529417363544578E-001, 1.6229497538220009E-001, 1.4942787541332894E-001,
      1.3669118160363830E-001, 1.2408319524476076E-001, 1.1160221162701203E-001,
      9.9246520573912658E-002, 8.7014407009240846E-002, 7.4904151479360401E-002,
      6.2914030717541744E-002, 5.1042318150120884E-002, 3.9287284443344221E-002,
      2.7647198016290835E-002, 1.6120325558104737E-002, 4.7049325474519624E-003,
      -6.6007162539859802E-003, -1.7798356285677958E-002,
      -2.8889722696324401E-002, -3.9876549857215615E-002,
      -5.0760570863198747E-002, -6.1543517058522935E-002,
      -7.2227117553263318E-002, -8.2813098750433634E-002,
      -9.3303183879441434E-002, -1.0369909251877797E-001,
      -1.1400254015486255E-001, -1.2421523770273571E-001,
      -1.3433889106838898E-001, -1.4437520069171539E-001,
      -1.5432586110659857E-001, -1.6419256049228403E-001,
      -1.7397698024627503E-001, -1.8368079453470640E-001,
      -1.9330566988409029E-001, -2.0285326472729215E-001,
      -2.1232522900789941E-001, -2.2172320373744142E-001,
      -2.3104882059342208E-001, -2.4030370149978114E-001,
      -2.4948945821610366E-001, -2.5860769193496108E-001,
      -2.6765999287285303E-001, -2.7664793987317643E-001,
      -2.8557310000424413E-001, -2.9443702816865763E-001,
      -3.0324126670618390E-001, -3.1198734500631370E-001,
      -3.2067677912035686E-001, -3.2931107137688226E-001,
      -3.3789171000192247E-001, -3.4642016873088977E-001,
      -3.5489790644288738E-001, -3.6332636677835856E-001,
      -3.7170697776544021E-001, -3.8004115145286571E-001,
      -3.8833028354023152E-001, -3.9657575301068493E-001,
      -4.0477892176676106E-001, -4.1294113426646312E-001,
      -4.2106371716626029E-001, -4.2914797895641310E-001,
      -4.3719520960912833E-001, -4.4520668021888843E-001,
      -4.5318364265237737E-001, -4.6112732918972199E-001,
      -4.6903895218322539E-001, -4.7691970369893982E-001,
      -4.8477075517191931E-001, -4.9259325706312135E-001,
      -5.0038833850818820E-001, -5.0815710697596683E-001,
      -5.1590064792785084E-001, -5.1522842553131620E-001,
      -5.0237840977112214E-001, -5.2359877559829837E-001,
      -5.2359877559829660E-001, -5.2359877559830170E-001,
      -5.2359877559829693E-001, -5.2359877559829782E-001,
      -5.2359877559829537E-001, -5.2359877559829948E-001,
      -5.2359877559830059E-001, -5.2359877559829926E-001,
      -5.2359877559829771E-001, -5.2359877559829893E-001,
      -5.2359877559829582E-001, -5.2359877559829859E-001,
      -5.2359877559829548E-001, -5.2359877559830104E-001,
      -5.2359877559830126E-001, -5.2359877559830137E-001,
      -5.2359877559829682E-001, -5.2359877559829848E-001,
      -5.2359877559830170E-001, -5.2359877559830104E-001,
      -5.2359877559829460E-001, -5.2359877559830037E-001,
      -5.2359877559829859E-001, -5.2359877559829737E-001,
      -5.2359877559829771E-001, -5.2359877559830659E-001,
      -5.2359877559829837E-001, -5.2359877559829526E-001,
      -5.2359877559830037E-001, -5.2359877559829349E-001,
      -5.2359877559830159E-001, -5.2359877559829815E-001,
      -5.2359877559829415E-001, -5.2359877559830148E-001,
      -5.2359877559829859E-001, -5.2359877559829782E-001,
      -5.2359877559829660E-001, -5.2359877559829848E-001,
      -5.2359877559829848E-001, -5.2359877559829660E-001,
      -5.2359877559829782E-001, -5.2359877559829582E-001,
      1.7904027546925212E+000, -1.9119996937727928E+000,
      -1.9119996937727928E+000, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

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
