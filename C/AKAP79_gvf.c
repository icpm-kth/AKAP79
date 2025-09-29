#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

/* Enums will be used for indexing purposes.   */
enum stateVariable { _RiiP, _RiiP_cAMP, _RiiP_C, _RiiP_C_cAMP, _C, _Rii_cAMP, _Rii_C_cAMP, _RiiP_CaN, _RiiP_cAMP_CaN, _AKAR4_C, _AKAR4p, numStateVar };
enum param { _k5_1, _k1_2, _k3_2, _k2_3, _k3_4, _k4_3, _k4_1, _k1_4, _k8_7, _k7_8, _k5_6, _k6_5, _k8_5, _k7_6, _k6_7, _k6_2, _k5_8, _k3p_7off, _k4p_8off, _k3p_3off, _k4p_4off, _k3p_7on, _k4p_8on, _k3p_3on, _k4p_4on, _kf_C_AKAR4, _kb_C_AKAR4, _kcat_AKARp, _km3OFF, _km4OFF, _km3ON, _km4ON, _KD_T, _b_AKAP, _AKAR4_ConservedConst, _CaN_ConservedConst, _Rii_C_ConservedConst, _cAMP_ConservedConst, _Rii_ConservedConst, numParam };
enum func { _AKAR4pOUT, numFunc };

/* The error codes indicate how many values a function returns.                             */
/* Each function expects the output buffer to be allocated with at least that many values   */

/* ODE vector field: y' = f(t,y;p)   */
int AKAP79_vf(double t, const double y_[], double *f_, void *par){
	double *p_=par;
	if (!y_ || !f_) return 11;
/* 	constants   */
/* 	parameter values   */
	double k5_1 = p_[_k5_1];                                            /* [  0] */
	double k1_2 = p_[_k1_2];                                            /* [  1] */
	double k3_2 = p_[_k3_2];                                            /* [  2] */
	double k2_3 = p_[_k2_3];                                            /* [  3] */
	double k3_4 = p_[_k3_4];                                            /* [  4] */
	double k4_3 = p_[_k4_3];                                            /* [  5] */
	double k4_1 = p_[_k4_1];                                            /* [  6] */
	double k1_4 = p_[_k1_4];                                            /* [  7] */
	double k8_7 = p_[_k8_7];                                            /* [  8] */
	double k7_8 = p_[_k7_8];                                            /* [  9] */
	double k5_6 = p_[_k5_6];                                            /* [ 10] */
	double k6_5 = p_[_k6_5];                                            /* [ 11] */
	double k8_5 = p_[_k8_5];                                            /* [ 12] */
	double k7_6 = p_[_k7_6];                                            /* [ 13] */
	double k6_7 = p_[_k6_7];                                            /* [ 14] */
	double k6_2 = p_[_k6_2];                                            /* [ 15] */
	double k5_8 = p_[_k5_8];                                            /* [ 16] */
	double k3p_7off = p_[_k3p_7off];                                    /* [ 17] */
	double k4p_8off = p_[_k4p_8off];                                    /* [ 18] */
	double k3p_3off = p_[_k3p_3off];                                    /* [ 19] */
	double k4p_4off = p_[_k4p_4off];                                    /* [ 20] */
	double k3p_7on = p_[_k3p_7on];                                      /* [ 21] */
	double k4p_8on = p_[_k4p_8on];                                      /* [ 22] */
	double k3p_3on = p_[_k3p_3on];                                      /* [ 23] */
	double k4p_4on = p_[_k4p_4on];                                      /* [ 24] */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [ 25] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [ 26] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [ 27] */
	double km3OFF = p_[_km3OFF];                                        /* [ 28] */
	double km4OFF = p_[_km4OFF];                                        /* [ 29] */
	double km3ON = p_[_km3ON];                                          /* [ 30] */
	double km4ON = p_[_km4ON];                                          /* [ 31] */
	double KD_T = p_[_KD_T];                                            /* [ 32] */
	double b_AKAP = p_[_b_AKAP];                                        /* [ 33] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [ 34] */
	double CaN_ConservedConst = p_[_CaN_ConservedConst];                /* [ 35] */
	double Rii_C_ConservedConst = p_[_Rii_C_ConservedConst];            /* [ 36] */
	double cAMP_ConservedConst = p_[_cAMP_ConservedConst];              /* [ 37] */
	double Rii_ConservedConst = p_[_Rii_ConservedConst];                /* [ 38] */
/* 	state variables   */
	double RiiP = y_[_RiiP];                                            /* [  0] */
	double RiiP_cAMP = y_[_RiiP_cAMP];                                  /* [  1] */
	double RiiP_C = y_[_RiiP_C];                                        /* [  2] */
	double RiiP_C_cAMP = y_[_RiiP_C_cAMP];                              /* [  3] */
	double C = y_[_C];                                                  /* [  4] */
	double Rii_cAMP = y_[_Rii_cAMP];                                    /* [  5] */
	double Rii_C_cAMP = y_[_Rii_C_cAMP];                                /* [  6] */
	double RiiP_CaN = y_[_RiiP_CaN];                                    /* [  7] */
	double RiiP_cAMP_CaN = y_[_RiiP_cAMP_CaN];                          /* [  8] */
	double AKAR4_C = y_[_AKAR4_C];                                      /* [  9] */
	double AKAR4p = y_[_AKAR4p];                                        /* [ 10] */
/* 	expressions   */
	double k3p_7 = b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	double k4p_4 = b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	double k4p_8 = b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	double k3p_3 = b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
	double k4_4p = b_AKAP * ((k4p_4 + k4p_8)/km4ON ) + (1 - b_AKAP) * (k4p_4 + k4p_8) / km4OFF;
	double k3_3p = b_AKAP * ((k3p_3 + k3p_7)/km3ON) + (1 - b_AKAP) * (k3p_3 + k3p_7)/km3OFF;
	double k2_1 = k1_2 * KD_T;
	double AKAR4 = (AKAR4_ConservedConst - (AKAR4_C+AKAR4p));
	double CaN = (CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN));
	double Rii_C = (Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	double cAMP = (cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	double Rii = (Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	double reaction_51 = k5_1*Rii_C;
	double reaction_14 = k4_1*RiiP*C - k1_4*RiiP_C;
	double reaction_12 = k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	double reaction_43 = k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	double reaction_23 = k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	double reaction_78 = k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	double reaction_56 = k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	double reaction_76 = k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	double reaction_62 = k6_2*Rii_C_cAMP;
	double reaction_58 = k8_5*Rii*C - k5_8*Rii_C;
	double reaction_44 = k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	double reaction_33 = k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	double reaction_48 = k4p_8*RiiP_CaN;
	double reaction_37 = k3p_7*RiiP_cAMP_CaN;
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	memset(f_,0,sizeof(double)*11); /* initialize with 0.0 */
	f_[0] = -reaction_14-reaction_43-reaction_44;
	f_[1] = +reaction_43-reaction_23-reaction_33;
	f_[2] = +reaction_51+reaction_14-reaction_12;
	f_[3] = +reaction_12+reaction_23+reaction_62;
	f_[4] = -reaction_14-reaction_23-reaction_76-reaction_58-reaction_1+reaction_2;
	f_[5] = +reaction_78-reaction_76+reaction_37;
	f_[6] = +reaction_56+reaction_76-reaction_62;
	f_[7] = +reaction_44-reaction_48;
	f_[8] = +reaction_33-reaction_37;
	f_[9] = +reaction_1-reaction_2;
	f_[10] = +reaction_2;
	return GSL_SUCCESS;
}

/* ODE Jacobian: df(t,y;p)/dy   */
int AKAP79_jac(double t, const double y_[], double *jac_, double *dfdt_, void *par){
	double *p_=par;
	if (!y_ || !jac_) return 121;
/* 	constants   */
/* 	parameter values   */
	double k5_1 = p_[_k5_1];                                            /* [  0] */
	double k1_2 = p_[_k1_2];                                            /* [  1] */
	double k3_2 = p_[_k3_2];                                            /* [  2] */
	double k2_3 = p_[_k2_3];                                            /* [  3] */
	double k3_4 = p_[_k3_4];                                            /* [  4] */
	double k4_3 = p_[_k4_3];                                            /* [  5] */
	double k4_1 = p_[_k4_1];                                            /* [  6] */
	double k1_4 = p_[_k1_4];                                            /* [  7] */
	double k8_7 = p_[_k8_7];                                            /* [  8] */
	double k7_8 = p_[_k7_8];                                            /* [  9] */
	double k5_6 = p_[_k5_6];                                            /* [ 10] */
	double k6_5 = p_[_k6_5];                                            /* [ 11] */
	double k8_5 = p_[_k8_5];                                            /* [ 12] */
	double k7_6 = p_[_k7_6];                                            /* [ 13] */
	double k6_7 = p_[_k6_7];                                            /* [ 14] */
	double k6_2 = p_[_k6_2];                                            /* [ 15] */
	double k5_8 = p_[_k5_8];                                            /* [ 16] */
	double k3p_7off = p_[_k3p_7off];                                    /* [ 17] */
	double k4p_8off = p_[_k4p_8off];                                    /* [ 18] */
	double k3p_3off = p_[_k3p_3off];                                    /* [ 19] */
	double k4p_4off = p_[_k4p_4off];                                    /* [ 20] */
	double k3p_7on = p_[_k3p_7on];                                      /* [ 21] */
	double k4p_8on = p_[_k4p_8on];                                      /* [ 22] */
	double k3p_3on = p_[_k3p_3on];                                      /* [ 23] */
	double k4p_4on = p_[_k4p_4on];                                      /* [ 24] */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [ 25] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [ 26] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [ 27] */
	double km3OFF = p_[_km3OFF];                                        /* [ 28] */
	double km4OFF = p_[_km4OFF];                                        /* [ 29] */
	double km3ON = p_[_km3ON];                                          /* [ 30] */
	double km4ON = p_[_km4ON];                                          /* [ 31] */
	double KD_T = p_[_KD_T];                                            /* [ 32] */
	double b_AKAP = p_[_b_AKAP];                                        /* [ 33] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [ 34] */
	double CaN_ConservedConst = p_[_CaN_ConservedConst];                /* [ 35] */
	double Rii_C_ConservedConst = p_[_Rii_C_ConservedConst];            /* [ 36] */
	double cAMP_ConservedConst = p_[_cAMP_ConservedConst];              /* [ 37] */
	double Rii_ConservedConst = p_[_Rii_ConservedConst];                /* [ 38] */
/* 	state variables   */
	double RiiP = y_[_RiiP];                                            /* [  0] */
	double RiiP_cAMP = y_[_RiiP_cAMP];                                  /* [  1] */
	double RiiP_C = y_[_RiiP_C];                                        /* [  2] */
	double RiiP_C_cAMP = y_[_RiiP_C_cAMP];                              /* [  3] */
	double C = y_[_C];                                                  /* [  4] */
	double Rii_cAMP = y_[_Rii_cAMP];                                    /* [  5] */
	double Rii_C_cAMP = y_[_Rii_C_cAMP];                                /* [  6] */
	double RiiP_CaN = y_[_RiiP_CaN];                                    /* [  7] */
	double RiiP_cAMP_CaN = y_[_RiiP_cAMP_CaN];                          /* [  8] */
	double AKAR4_C = y_[_AKAR4_C];                                      /* [  9] */
	double AKAR4p = y_[_AKAR4p];                                        /* [ 10] */
/* 	expressions   */
	double k3p_7 = b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	double k4p_4 = b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	double k4p_8 = b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	double k3p_3 = b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
	double k4_4p = b_AKAP * ((k4p_4 + k4p_8)/km4ON ) + (1 - b_AKAP) * (k4p_4 + k4p_8) / km4OFF;
	double k3_3p = b_AKAP * ((k3p_3 + k3p_7)/km3ON) + (1 - b_AKAP) * (k3p_3 + k3p_7)/km3OFF;
	double k2_1 = k1_2 * KD_T;
	double AKAR4 = (AKAR4_ConservedConst - (AKAR4_C+AKAR4p));
	double CaN = (CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN));
	double Rii_C = (Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	double cAMP = (cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	double Rii = (Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	double reaction_51 = k5_1*Rii_C;
	double reaction_14 = k4_1*RiiP*C - k1_4*RiiP_C;
	double reaction_12 = k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	double reaction_43 = k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	double reaction_23 = k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	double reaction_78 = k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	double reaction_56 = k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	double reaction_76 = k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	double reaction_62 = k6_2*Rii_C_cAMP;
	double reaction_58 = k8_5*Rii*C - k5_8*Rii_C;
	double reaction_44 = k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	double reaction_33 = k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	double reaction_48 = k4p_8*RiiP_CaN;
	double reaction_37 = k3p_7*RiiP_cAMP_CaN;
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	memset(jac_,0,sizeof(double)*121); /* initialize with 0.0 */
	/*[ 0, 0]*/  jac_[0] = -(k4_1*C+k4_3*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))+((b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON+((1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN)));
	/*[ 0, 1]*/  jac_[1] = k4_3*RiiP+k3_4;
	/*[ 0, 2]*/  jac_[2] = k1_4;
	/*[ 0, 3]*/  jac_[3] = k4_3*RiiP;
	/*[ 0, 4]*/  jac_[4] = -k4_1*RiiP;
	/*[ 0, 5]*/  jac_[5] = k4_3*RiiP;
	/*[ 0, 6]*/  jac_[6] = k4_3*RiiP;
	/*[ 0, 7]*/  jac_[7] = ((b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON+((1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*RiiP+b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off;
	/*[ 0, 8]*/  jac_[8] = k4_3*RiiP+((b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON+((1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*RiiP;
	/*[ 1, 0]*/  jac_[11] = k4_3*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 1, 1]*/  jac_[12] = -(k4_3*RiiP+k3_4+k3_2*C+((b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON+((1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN)));
	/*[ 1, 3]*/  jac_[14] = k2_3-k4_3*RiiP;
	/*[ 1, 4]*/  jac_[15] = -k3_2*RiiP_cAMP;
	/*[ 1, 5]*/  jac_[16] = -k4_3*RiiP;
	/*[ 1, 6]*/  jac_[17] = -k4_3*RiiP;
	/*[ 1, 7]*/  jac_[18] = ((b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON+((1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*RiiP_cAMP;
	/*[ 1, 8]*/  jac_[19] = ((b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON+((1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*RiiP_cAMP+b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off-k4_3*RiiP;
	/*[ 2, 0]*/  jac_[22] = k4_1*C;
	/*[ 2, 1]*/  jac_[23] = k1_2*RiiP_C;
	/*[ 2, 2]*/  jac_[24] = -(k1_4+k5_1+k1_2*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN)));
	/*[ 2, 3]*/  jac_[25] = k1_2*RiiP_C+k1_2*KD_T-k5_1;
	/*[ 2, 4]*/  jac_[26] = k4_1*RiiP-k5_1;
	/*[ 2, 5]*/  jac_[27] = k1_2*RiiP_C;
	/*[ 2, 6]*/  jac_[28] = k1_2*RiiP_C-k5_1;
	/*[ 2, 8]*/  jac_[30] = k1_2*RiiP_C;
	/*[ 2, 9]*/  jac_[31] = -k5_1;
	/*[ 3, 1]*/  jac_[34] = k3_2*C-k1_2*RiiP_C;
	/*[ 3, 2]*/  jac_[35] = k1_2*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 3, 3]*/  jac_[36] = -(k2_3+k1_2*RiiP_C+k1_2*KD_T);
	/*[ 3, 4]*/  jac_[37] = k3_2*RiiP_cAMP;
	/*[ 3, 5]*/  jac_[38] = -k1_2*RiiP_C;
	/*[ 3, 6]*/  jac_[39] = k6_2-k1_2*RiiP_C;
	/*[ 3, 8]*/  jac_[41] = -k1_2*RiiP_C;
	/*[ 4, 0]*/  jac_[44] = k8_5*C-k4_1*C;
	/*[ 4, 1]*/  jac_[45] = k8_5*C-k3_2*C;
	/*[ 4, 2]*/  jac_[46] = k1_4-k5_8;
	/*[ 4, 3]*/  jac_[47] = k2_3-k5_8;
	/*[ 4, 4]*/  jac_[48] = -(k4_1*RiiP+k3_2*RiiP_cAMP+k7_6*Rii_cAMP+k8_5*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))+k8_5*C+k5_8+kf_C_AKAR4*(AKAR4_ConservedConst-(AKAR4_C+AKAR4p)));
	/*[ 4, 5]*/  jac_[49] = k8_5*C-k7_6*C;
	/*[ 4, 6]*/  jac_[50] = k6_7-k5_8;
	/*[ 4, 7]*/  jac_[51] = k8_5*C;
	/*[ 4, 8]*/  jac_[52] = k8_5*C;
	/*[ 4, 9]*/  jac_[53] = kf_C_AKAR4*C+kb_C_AKAR4-(k8_5*C+k5_8)+kcat_AKARp;
	/*[ 4,10]*/  jac_[54] = kf_C_AKAR4*C;
	/*[ 5, 0]*/  jac_[55] = -k8_7*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 5, 1]*/  jac_[56] = -(k8_7*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))+k8_7*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN)));
	/*[ 5, 3]*/  jac_[58] = -k8_7*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	/*[ 5, 4]*/  jac_[59] = k8_7*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))-k7_6*Rii_cAMP;
	/*[ 5, 5]*/  jac_[60] = -(k8_7*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))+k8_7*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))+k7_8+k7_6*C);
	/*[ 5, 6]*/  jac_[61] = k6_7-k8_7*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	/*[ 5, 7]*/  jac_[62] = -k8_7*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 5, 8]*/  jac_[63] = b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off-(k8_7*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))+k8_7*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN)));
	/*[ 5, 9]*/  jac_[64] = k8_7*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 6, 1]*/  jac_[67] = -k5_6*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	/*[ 6, 2]*/  jac_[68] = -k5_6*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 6, 3]*/  jac_[69] = -(k5_6*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))+k5_6*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C)));
	/*[ 6, 4]*/  jac_[70] = k7_6*Rii_cAMP-k5_6*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 6, 5]*/  jac_[71] = k7_6*C-k5_6*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	/*[ 6, 6]*/  jac_[72] = -(k6_7+k5_6*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))+k5_6*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))+k6_5+k6_2);
	/*[ 6, 8]*/  jac_[74] = -k5_6*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	/*[ 6, 9]*/  jac_[75] = -k5_6*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 7, 0]*/  jac_[77] = ((b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON+((1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN));
	/*[ 7, 7]*/  jac_[84] = -(((b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON+((1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*RiiP+b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off);
	/*[ 7, 8]*/  jac_[85] = -((b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON+((1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*RiiP;
	/*[ 8, 1]*/  jac_[89] = ((b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON+((1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN));
	/*[ 8, 7]*/  jac_[95] = -((b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON+((1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*RiiP_cAMP;
	/*[ 8, 8]*/  jac_[96] = -(((b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON+((1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*RiiP_cAMP+b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off);
	/*[ 9, 4]*/  jac_[103] = kf_C_AKAR4*(AKAR4_ConservedConst-(AKAR4_C+AKAR4p));
	/*[ 9, 9]*/  jac_[108] = -(kf_C_AKAR4*C+kb_C_AKAR4+kcat_AKARp);
	/*[ 9,10]*/  jac_[109] = -kf_C_AKAR4*C;
	/*[10, 9]*/  jac_[119] = kcat_AKARp;
	return GSL_SUCCESS;
}

/* ODE parameter Jacobian: df(t,y;p)/dp   */
int AKAP79_jacp(double t, const double y_[], double *jacp_, double *dfdt_, void *par){
	double *p_=par;
	if (!y_ || !jacp_) return 429;
/* 	constants   */
/* 	parameter values   */
	double k5_1 = p_[_k5_1];                                            /* [  0] */
	double k1_2 = p_[_k1_2];                                            /* [  1] */
	double k3_2 = p_[_k3_2];                                            /* [  2] */
	double k2_3 = p_[_k2_3];                                            /* [  3] */
	double k3_4 = p_[_k3_4];                                            /* [  4] */
	double k4_3 = p_[_k4_3];                                            /* [  5] */
	double k4_1 = p_[_k4_1];                                            /* [  6] */
	double k1_4 = p_[_k1_4];                                            /* [  7] */
	double k8_7 = p_[_k8_7];                                            /* [  8] */
	double k7_8 = p_[_k7_8];                                            /* [  9] */
	double k5_6 = p_[_k5_6];                                            /* [ 10] */
	double k6_5 = p_[_k6_5];                                            /* [ 11] */
	double k8_5 = p_[_k8_5];                                            /* [ 12] */
	double k7_6 = p_[_k7_6];                                            /* [ 13] */
	double k6_7 = p_[_k6_7];                                            /* [ 14] */
	double k6_2 = p_[_k6_2];                                            /* [ 15] */
	double k5_8 = p_[_k5_8];                                            /* [ 16] */
	double k3p_7off = p_[_k3p_7off];                                    /* [ 17] */
	double k4p_8off = p_[_k4p_8off];                                    /* [ 18] */
	double k3p_3off = p_[_k3p_3off];                                    /* [ 19] */
	double k4p_4off = p_[_k4p_4off];                                    /* [ 20] */
	double k3p_7on = p_[_k3p_7on];                                      /* [ 21] */
	double k4p_8on = p_[_k4p_8on];                                      /* [ 22] */
	double k3p_3on = p_[_k3p_3on];                                      /* [ 23] */
	double k4p_4on = p_[_k4p_4on];                                      /* [ 24] */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [ 25] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [ 26] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [ 27] */
	double km3OFF = p_[_km3OFF];                                        /* [ 28] */
	double km4OFF = p_[_km4OFF];                                        /* [ 29] */
	double km3ON = p_[_km3ON];                                          /* [ 30] */
	double km4ON = p_[_km4ON];                                          /* [ 31] */
	double KD_T = p_[_KD_T];                                            /* [ 32] */
	double b_AKAP = p_[_b_AKAP];                                        /* [ 33] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [ 34] */
	double CaN_ConservedConst = p_[_CaN_ConservedConst];                /* [ 35] */
	double Rii_C_ConservedConst = p_[_Rii_C_ConservedConst];            /* [ 36] */
	double cAMP_ConservedConst = p_[_cAMP_ConservedConst];              /* [ 37] */
	double Rii_ConservedConst = p_[_Rii_ConservedConst];                /* [ 38] */
/* 	state variables   */
	double RiiP = y_[_RiiP];                                            /* [  0] */
	double RiiP_cAMP = y_[_RiiP_cAMP];                                  /* [  1] */
	double RiiP_C = y_[_RiiP_C];                                        /* [  2] */
	double RiiP_C_cAMP = y_[_RiiP_C_cAMP];                              /* [  3] */
	double C = y_[_C];                                                  /* [  4] */
	double Rii_cAMP = y_[_Rii_cAMP];                                    /* [  5] */
	double Rii_C_cAMP = y_[_Rii_C_cAMP];                                /* [  6] */
	double RiiP_CaN = y_[_RiiP_CaN];                                    /* [  7] */
	double RiiP_cAMP_CaN = y_[_RiiP_cAMP_CaN];                          /* [  8] */
	double AKAR4_C = y_[_AKAR4_C];                                      /* [  9] */
	double AKAR4p = y_[_AKAR4p];                                        /* [ 10] */
/* 	expressions   */
	double k3p_7 = b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	double k4p_4 = b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	double k4p_8 = b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	double k3p_3 = b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
	double k4_4p = b_AKAP * ((k4p_4 + k4p_8)/km4ON ) + (1 - b_AKAP) * (k4p_4 + k4p_8) / km4OFF;
	double k3_3p = b_AKAP * ((k3p_3 + k3p_7)/km3ON) + (1 - b_AKAP) * (k3p_3 + k3p_7)/km3OFF;
	double k2_1 = k1_2 * KD_T;
	double AKAR4 = (AKAR4_ConservedConst - (AKAR4_C+AKAR4p));
	double CaN = (CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN));
	double Rii_C = (Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	double cAMP = (cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	double Rii = (Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	double reaction_51 = k5_1*Rii_C;
	double reaction_14 = k4_1*RiiP*C - k1_4*RiiP_C;
	double reaction_12 = k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	double reaction_43 = k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	double reaction_23 = k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	double reaction_78 = k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	double reaction_56 = k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	double reaction_76 = k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	double reaction_62 = k6_2*Rii_C_cAMP;
	double reaction_58 = k8_5*Rii*C - k5_8*Rii_C;
	double reaction_44 = k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	double reaction_33 = k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	double reaction_48 = k4p_8*RiiP_CaN;
	double reaction_37 = k3p_7*RiiP_cAMP_CaN;
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	memset(jacp_,0,sizeof(double)*429); /* initialize with 0.0 */
	/*[ 0, 4]*/  jacp_[4] = RiiP_cAMP;
	/*[ 0, 5]*/  jacp_[5] = -(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))*RiiP;
	/*[ 0, 6]*/  jacp_[6] = -RiiP*C;
	/*[ 0, 7]*/  jacp_[7] = RiiP_C;
	/*[ 0,18]*/  jacp_[18] = -((b_AKAP*(1-b_AKAP))/km4ON+gsl_pow_2(1-b_AKAP)/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN));
	/*[ 0,20]*/  jacp_[20] = (1-b_AKAP)*RiiP_CaN-((b_AKAP*(1-b_AKAP))/km4ON+gsl_pow_2(1-b_AKAP)/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN));
	/*[ 0,22]*/  jacp_[22] = -(gsl_pow_2(b_AKAP)/km4ON+((1-b_AKAP)*b_AKAP)/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN));
	/*[ 0,24]*/  jacp_[24] = b_AKAP*RiiP_CaN-(gsl_pow_2(b_AKAP)/km4ON+((1-b_AKAP)*b_AKAP)/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN));
	/*[ 0,29]*/  jacp_[29] = ((CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP*(1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/gsl_pow_2(km4OFF);
	/*[ 0,31]*/  jacp_[31] = ((CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP*b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/gsl_pow_2(km4ON);
	/*[ 0,33]*/  jacp_[33] = (k4p_4on-k4p_4off)*RiiP_CaN-((b_AKAP*(k4p_4on-k4p_4off+k4p_8on-k4p_8off)+b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off)/km4ON+((1-b_AKAP)*(k4p_4on-k4p_4off+k4p_8on-k4p_8off)-(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN));
	/*[ 0,35]*/  jacp_[35] = -((b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON+((1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*RiiP;
	/*[ 0,37]*/  jacp_[37] = -k4_3*RiiP;
	/*[ 1, 2]*/  jacp_[41] = -RiiP_cAMP*C;
	/*[ 1, 3]*/  jacp_[42] = RiiP_C_cAMP;
	/*[ 1, 4]*/  jacp_[43] = -RiiP_cAMP;
	/*[ 1, 5]*/  jacp_[44] = (cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))*RiiP;
	/*[ 1,17]*/  jacp_[56] = -((b_AKAP*(1-b_AKAP))/km3ON+gsl_pow_2(1-b_AKAP)/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP;
	/*[ 1,19]*/  jacp_[58] = (1-b_AKAP)*RiiP_cAMP_CaN-((b_AKAP*(1-b_AKAP))/km3ON+gsl_pow_2(1-b_AKAP)/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP;
	/*[ 1,21]*/  jacp_[60] = -(gsl_pow_2(b_AKAP)/km3ON+((1-b_AKAP)*b_AKAP)/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP;
	/*[ 1,23]*/  jacp_[62] = b_AKAP*RiiP_cAMP_CaN-(gsl_pow_2(b_AKAP)/km3ON+((1-b_AKAP)*b_AKAP)/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP;
	/*[ 1,28]*/  jacp_[67] = (RiiP_cAMP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*(1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/gsl_pow_2(km3OFF);
	/*[ 1,30]*/  jacp_[69] = (RiiP_cAMP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/gsl_pow_2(km3ON);
	/*[ 1,33]*/  jacp_[72] = (k3p_3on-k3p_3off)*RiiP_cAMP_CaN-((b_AKAP*(k3p_3on-k3p_3off+k3p_7on-k3p_7off)+b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off)/km3ON+((1-b_AKAP)*(k3p_3on-k3p_3off+k3p_7on-k3p_7off)-(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP;
	/*[ 1,35]*/  jacp_[74] = -((b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON+((1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*RiiP_cAMP;
	/*[ 1,37]*/  jacp_[76] = k4_3*RiiP;
	/*[ 2, 0]*/  jacp_[78] = Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C);
	/*[ 2, 1]*/  jacp_[79] = KD_T*RiiP_C_cAMP-RiiP_C*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 2, 6]*/  jacp_[84] = RiiP*C;
	/*[ 2, 7]*/  jacp_[85] = -RiiP_C;
	/*[ 2,32]*/  jacp_[110] = k1_2*RiiP_C_cAMP;
	/*[ 2,36]*/  jacp_[114] = k5_1;
	/*[ 2,37]*/  jacp_[115] = -k1_2*RiiP_C;
	/*[ 3, 1]*/  jacp_[118] = RiiP_C*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))-KD_T*RiiP_C_cAMP;
	/*[ 3, 2]*/  jacp_[119] = RiiP_cAMP*C;
	/*[ 3, 3]*/  jacp_[120] = -RiiP_C_cAMP;
	/*[ 3,15]*/  jacp_[132] = Rii_C_cAMP;
	/*[ 3,32]*/  jacp_[149] = -k1_2*RiiP_C_cAMP;
	/*[ 3,37]*/  jacp_[154] = k1_2*RiiP_C;
	/*[ 4, 2]*/  jacp_[158] = -RiiP_cAMP*C;
	/*[ 4, 3]*/  jacp_[159] = RiiP_C_cAMP;
	/*[ 4, 6]*/  jacp_[162] = -RiiP*C;
	/*[ 4, 7]*/  jacp_[163] = RiiP_C;
	/*[ 4,12]*/  jacp_[168] = -(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))*C;
	/*[ 4,13]*/  jacp_[169] = -Rii_cAMP*C;
	/*[ 4,14]*/  jacp_[170] = Rii_C_cAMP;
	/*[ 4,16]*/  jacp_[172] = Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C);
	/*[ 4,25]*/  jacp_[181] = -C*(AKAR4_ConservedConst-(AKAR4_C+AKAR4p));
	/*[ 4,26]*/  jacp_[182] = AKAR4_C;
	/*[ 4,27]*/  jacp_[183] = AKAR4_C;
	/*[ 4,34]*/  jacp_[190] = -kf_C_AKAR4*C;
	/*[ 4,36]*/  jacp_[192] = k5_8;
	/*[ 4,38]*/  jacp_[194] = -k8_5*C;
	/*[ 5, 8]*/  jacp_[203] = (cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	/*[ 5, 9]*/  jacp_[204] = -Rii_cAMP;
	/*[ 5,13]*/  jacp_[208] = -Rii_cAMP*C;
	/*[ 5,14]*/  jacp_[209] = Rii_C_cAMP;
	/*[ 5,17]*/  jacp_[212] = (1-b_AKAP)*RiiP_cAMP_CaN;
	/*[ 5,21]*/  jacp_[216] = b_AKAP*RiiP_cAMP_CaN;
	/*[ 5,33]*/  jacp_[228] = (k3p_7on-k3p_7off)*RiiP_cAMP_CaN;
	/*[ 5,37]*/  jacp_[232] = k8_7*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	/*[ 5,38]*/  jacp_[233] = k8_7*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 6,10]*/  jacp_[244] = (Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 6,11]*/  jacp_[245] = -Rii_C_cAMP;
	/*[ 6,13]*/  jacp_[247] = Rii_cAMP*C;
	/*[ 6,14]*/  jacp_[248] = -Rii_C_cAMP;
	/*[ 6,15]*/  jacp_[249] = -Rii_C_cAMP;
	/*[ 6,36]*/  jacp_[270] = k5_6*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 6,37]*/  jacp_[271] = k5_6*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	/*[ 7,18]*/  jacp_[291] = ((b_AKAP*(1-b_AKAP))/km4ON+gsl_pow_2(1-b_AKAP)/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))-(1-b_AKAP)*RiiP_CaN;
	/*[ 7,20]*/  jacp_[293] = ((b_AKAP*(1-b_AKAP))/km4ON+gsl_pow_2(1-b_AKAP)/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))-(1-b_AKAP)*RiiP_CaN;
	/*[ 7,22]*/  jacp_[295] = (gsl_pow_2(b_AKAP)/km4ON+((1-b_AKAP)*b_AKAP)/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))-b_AKAP*RiiP_CaN;
	/*[ 7,24]*/  jacp_[297] = (gsl_pow_2(b_AKAP)/km4ON+((1-b_AKAP)*b_AKAP)/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))-b_AKAP*RiiP_CaN;
	/*[ 7,29]*/  jacp_[302] = (-(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP*(1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/gsl_pow_2(km4OFF);
	/*[ 7,31]*/  jacp_[304] = (-(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP*b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/gsl_pow_2(km4ON);
	/*[ 7,33]*/  jacp_[306] = ((b_AKAP*(k4p_4on-k4p_4off+k4p_8on-k4p_8off)+b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off)/km4ON+((1-b_AKAP)*(k4p_4on-k4p_4off+k4p_8on-k4p_8off)-(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))-(k4p_4on-k4p_4off)*RiiP_CaN-(k4p_8on-k4p_8off)*RiiP_CaN;
	/*[ 7,35]*/  jacp_[308] = ((b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON+((1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*RiiP;
	/*[ 8,17]*/  jacp_[329] = ((b_AKAP*(1-b_AKAP))/km3ON+gsl_pow_2(1-b_AKAP)/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP-(1-b_AKAP)*RiiP_cAMP_CaN;
	/*[ 8,19]*/  jacp_[331] = ((b_AKAP*(1-b_AKAP))/km3ON+gsl_pow_2(1-b_AKAP)/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP-(1-b_AKAP)*RiiP_cAMP_CaN;
	/*[ 8,21]*/  jacp_[333] = (gsl_pow_2(b_AKAP)/km3ON+((1-b_AKAP)*b_AKAP)/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP-b_AKAP*RiiP_cAMP_CaN;
	/*[ 8,23]*/  jacp_[335] = (gsl_pow_2(b_AKAP)/km3ON+((1-b_AKAP)*b_AKAP)/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP-b_AKAP*RiiP_cAMP_CaN;
	/*[ 8,28]*/  jacp_[340] = (-RiiP_cAMP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*(1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/gsl_pow_2(km3OFF);
	/*[ 8,30]*/  jacp_[342] = (-RiiP_cAMP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/gsl_pow_2(km3ON);
	/*[ 8,33]*/  jacp_[345] = ((b_AKAP*(k3p_3on-k3p_3off+k3p_7on-k3p_7off)+b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off)/km3ON+((1-b_AKAP)*(k3p_3on-k3p_3off+k3p_7on-k3p_7off)-(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP-(k3p_3on-k3p_3off)*RiiP_cAMP_CaN-(k3p_7on-k3p_7off)*RiiP_cAMP_CaN;
	/*[ 8,35]*/  jacp_[347] = ((b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON+((1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*RiiP_cAMP;
	/*[ 9,25]*/  jacp_[376] = C*(AKAR4_ConservedConst-(AKAR4_C+AKAR4p));
	/*[ 9,26]*/  jacp_[377] = -AKAR4_C;
	/*[ 9,27]*/  jacp_[378] = -AKAR4_C;
	/*[ 9,34]*/  jacp_[385] = kf_C_AKAR4*C;
	/*[10,27]*/  jacp_[417] = AKAR4_C;
	return GSL_SUCCESS;
}

/* Output Function (Observables)   */
int AKAP79_func(double t, const double y_[], double *func_, void *par){
	double *p_=par;
	if (!y_ || !func_) return 1;
/* 	constants   */
/* 	parameter values   */
	double k5_1 = p_[_k5_1];                                            /* [  0] */
	double k1_2 = p_[_k1_2];                                            /* [  1] */
	double k3_2 = p_[_k3_2];                                            /* [  2] */
	double k2_3 = p_[_k2_3];                                            /* [  3] */
	double k3_4 = p_[_k3_4];                                            /* [  4] */
	double k4_3 = p_[_k4_3];                                            /* [  5] */
	double k4_1 = p_[_k4_1];                                            /* [  6] */
	double k1_4 = p_[_k1_4];                                            /* [  7] */
	double k8_7 = p_[_k8_7];                                            /* [  8] */
	double k7_8 = p_[_k7_8];                                            /* [  9] */
	double k5_6 = p_[_k5_6];                                            /* [ 10] */
	double k6_5 = p_[_k6_5];                                            /* [ 11] */
	double k8_5 = p_[_k8_5];                                            /* [ 12] */
	double k7_6 = p_[_k7_6];                                            /* [ 13] */
	double k6_7 = p_[_k6_7];                                            /* [ 14] */
	double k6_2 = p_[_k6_2];                                            /* [ 15] */
	double k5_8 = p_[_k5_8];                                            /* [ 16] */
	double k3p_7off = p_[_k3p_7off];                                    /* [ 17] */
	double k4p_8off = p_[_k4p_8off];                                    /* [ 18] */
	double k3p_3off = p_[_k3p_3off];                                    /* [ 19] */
	double k4p_4off = p_[_k4p_4off];                                    /* [ 20] */
	double k3p_7on = p_[_k3p_7on];                                      /* [ 21] */
	double k4p_8on = p_[_k4p_8on];                                      /* [ 22] */
	double k3p_3on = p_[_k3p_3on];                                      /* [ 23] */
	double k4p_4on = p_[_k4p_4on];                                      /* [ 24] */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [ 25] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [ 26] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [ 27] */
	double km3OFF = p_[_km3OFF];                                        /* [ 28] */
	double km4OFF = p_[_km4OFF];                                        /* [ 29] */
	double km3ON = p_[_km3ON];                                          /* [ 30] */
	double km4ON = p_[_km4ON];                                          /* [ 31] */
	double KD_T = p_[_KD_T];                                            /* [ 32] */
	double b_AKAP = p_[_b_AKAP];                                        /* [ 33] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [ 34] */
	double CaN_ConservedConst = p_[_CaN_ConservedConst];                /* [ 35] */
	double Rii_C_ConservedConst = p_[_Rii_C_ConservedConst];            /* [ 36] */
	double cAMP_ConservedConst = p_[_cAMP_ConservedConst];              /* [ 37] */
	double Rii_ConservedConst = p_[_Rii_ConservedConst];                /* [ 38] */
/* 	state variables   */
	double RiiP = y_[_RiiP];                                            /* [  0] */
	double RiiP_cAMP = y_[_RiiP_cAMP];                                  /* [  1] */
	double RiiP_C = y_[_RiiP_C];                                        /* [  2] */
	double RiiP_C_cAMP = y_[_RiiP_C_cAMP];                              /* [  3] */
	double C = y_[_C];                                                  /* [  4] */
	double Rii_cAMP = y_[_Rii_cAMP];                                    /* [  5] */
	double Rii_C_cAMP = y_[_Rii_C_cAMP];                                /* [  6] */
	double RiiP_CaN = y_[_RiiP_CaN];                                    /* [  7] */
	double RiiP_cAMP_CaN = y_[_RiiP_cAMP_CaN];                          /* [  8] */
	double AKAR4_C = y_[_AKAR4_C];                                      /* [  9] */
	double AKAR4p = y_[_AKAR4p];                                        /* [ 10] */
/* 	expressions   */
	double k3p_7 = b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	double k4p_4 = b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	double k4p_8 = b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	double k3p_3 = b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
	double k4_4p = b_AKAP * ((k4p_4 + k4p_8)/km4ON ) + (1 - b_AKAP) * (k4p_4 + k4p_8) / km4OFF;
	double k3_3p = b_AKAP * ((k3p_3 + k3p_7)/km3ON) + (1 - b_AKAP) * (k3p_3 + k3p_7)/km3OFF;
	double k2_1 = k1_2 * KD_T;
	double AKAR4 = (AKAR4_ConservedConst - (AKAR4_C+AKAR4p));
	double CaN = (CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN));
	double Rii_C = (Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	double cAMP = (cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	double Rii = (Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	double reaction_51 = k5_1*Rii_C;
	double reaction_14 = k4_1*RiiP*C - k1_4*RiiP_C;
	double reaction_12 = k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	double reaction_43 = k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	double reaction_23 = k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	double reaction_78 = k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	double reaction_56 = k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	double reaction_76 = k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	double reaction_62 = k6_2*Rii_C_cAMP;
	double reaction_58 = k8_5*Rii*C - k5_8*Rii_C;
	double reaction_44 = k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	double reaction_33 = k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	double reaction_48 = k4p_8*RiiP_CaN;
	double reaction_37 = k3p_7*RiiP_cAMP_CaN;
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	func_[_AKAR4pOUT] = (AKAR4p*5)*71.67+100;
	return GSL_SUCCESS;
}

/* Output function Jacobian: dF(t,y;p)/dx   */
int AKAP79_funcJac(double t, const double y_[], double *funcJac_, void *par){
	double *p_=par;
	if (!y_ || !funcJac_) return 11;
/* 	constants   */
/* 	parameter values   */
	double k5_1 = p_[_k5_1];                                            /* [  0] */
	double k1_2 = p_[_k1_2];                                            /* [  1] */
	double k3_2 = p_[_k3_2];                                            /* [  2] */
	double k2_3 = p_[_k2_3];                                            /* [  3] */
	double k3_4 = p_[_k3_4];                                            /* [  4] */
	double k4_3 = p_[_k4_3];                                            /* [  5] */
	double k4_1 = p_[_k4_1];                                            /* [  6] */
	double k1_4 = p_[_k1_4];                                            /* [  7] */
	double k8_7 = p_[_k8_7];                                            /* [  8] */
	double k7_8 = p_[_k7_8];                                            /* [  9] */
	double k5_6 = p_[_k5_6];                                            /* [ 10] */
	double k6_5 = p_[_k6_5];                                            /* [ 11] */
	double k8_5 = p_[_k8_5];                                            /* [ 12] */
	double k7_6 = p_[_k7_6];                                            /* [ 13] */
	double k6_7 = p_[_k6_7];                                            /* [ 14] */
	double k6_2 = p_[_k6_2];                                            /* [ 15] */
	double k5_8 = p_[_k5_8];                                            /* [ 16] */
	double k3p_7off = p_[_k3p_7off];                                    /* [ 17] */
	double k4p_8off = p_[_k4p_8off];                                    /* [ 18] */
	double k3p_3off = p_[_k3p_3off];                                    /* [ 19] */
	double k4p_4off = p_[_k4p_4off];                                    /* [ 20] */
	double k3p_7on = p_[_k3p_7on];                                      /* [ 21] */
	double k4p_8on = p_[_k4p_8on];                                      /* [ 22] */
	double k3p_3on = p_[_k3p_3on];                                      /* [ 23] */
	double k4p_4on = p_[_k4p_4on];                                      /* [ 24] */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [ 25] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [ 26] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [ 27] */
	double km3OFF = p_[_km3OFF];                                        /* [ 28] */
	double km4OFF = p_[_km4OFF];                                        /* [ 29] */
	double km3ON = p_[_km3ON];                                          /* [ 30] */
	double km4ON = p_[_km4ON];                                          /* [ 31] */
	double KD_T = p_[_KD_T];                                            /* [ 32] */
	double b_AKAP = p_[_b_AKAP];                                        /* [ 33] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [ 34] */
	double CaN_ConservedConst = p_[_CaN_ConservedConst];                /* [ 35] */
	double Rii_C_ConservedConst = p_[_Rii_C_ConservedConst];            /* [ 36] */
	double cAMP_ConservedConst = p_[_cAMP_ConservedConst];              /* [ 37] */
	double Rii_ConservedConst = p_[_Rii_ConservedConst];                /* [ 38] */
/* 	state variables   */
	double RiiP = y_[_RiiP];                                            /* [  0] */
	double RiiP_cAMP = y_[_RiiP_cAMP];                                  /* [  1] */
	double RiiP_C = y_[_RiiP_C];                                        /* [  2] */
	double RiiP_C_cAMP = y_[_RiiP_C_cAMP];                              /* [  3] */
	double C = y_[_C];                                                  /* [  4] */
	double Rii_cAMP = y_[_Rii_cAMP];                                    /* [  5] */
	double Rii_C_cAMP = y_[_Rii_C_cAMP];                                /* [  6] */
	double RiiP_CaN = y_[_RiiP_CaN];                                    /* [  7] */
	double RiiP_cAMP_CaN = y_[_RiiP_cAMP_CaN];                          /* [  8] */
	double AKAR4_C = y_[_AKAR4_C];                                      /* [  9] */
	double AKAR4p = y_[_AKAR4p];                                        /* [ 10] */
/* 	expressions   */
	double k3p_7 = b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	double k4p_4 = b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	double k4p_8 = b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	double k3p_3 = b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
	double k4_4p = b_AKAP * ((k4p_4 + k4p_8)/km4ON ) + (1 - b_AKAP) * (k4p_4 + k4p_8) / km4OFF;
	double k3_3p = b_AKAP * ((k3p_3 + k3p_7)/km3ON) + (1 - b_AKAP) * (k3p_3 + k3p_7)/km3OFF;
	double k2_1 = k1_2 * KD_T;
	double AKAR4 = (AKAR4_ConservedConst - (AKAR4_C+AKAR4p));
	double CaN = (CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN));
	double Rii_C = (Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	double cAMP = (cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	double Rii = (Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	double reaction_51 = k5_1*Rii_C;
	double reaction_14 = k4_1*RiiP*C - k1_4*RiiP_C;
	double reaction_12 = k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	double reaction_43 = k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	double reaction_23 = k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	double reaction_78 = k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	double reaction_56 = k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	double reaction_76 = k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	double reaction_62 = k6_2*Rii_C_cAMP;
	double reaction_58 = k8_5*Rii*C - k5_8*Rii_C;
	double reaction_44 = k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	double reaction_33 = k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	double reaction_48 = k4p_8*RiiP_CaN;
	double reaction_37 = k3p_7*RiiP_cAMP_CaN;
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	memset(funcJac_,0,sizeof(double)*11); /* initialize with 0.0 */
	/*[ 0,10]*/  funcJac_[10] = 358.35;
	return GSL_SUCCESS;
}

/* Output function parameter Jacobian: dF(t,y;p)/dp   */
int AKAP79_funcJacp(double t, const double y_[], double *funcJacp_, void *par){
	double *p_=par;
	if (!y_ || !funcJacp_) return 39;
/* 	constants   */
/* 	parameter values   */
	double k5_1 = p_[_k5_1];                                            /* [  0] */
	double k1_2 = p_[_k1_2];                                            /* [  1] */
	double k3_2 = p_[_k3_2];                                            /* [  2] */
	double k2_3 = p_[_k2_3];                                            /* [  3] */
	double k3_4 = p_[_k3_4];                                            /* [  4] */
	double k4_3 = p_[_k4_3];                                            /* [  5] */
	double k4_1 = p_[_k4_1];                                            /* [  6] */
	double k1_4 = p_[_k1_4];                                            /* [  7] */
	double k8_7 = p_[_k8_7];                                            /* [  8] */
	double k7_8 = p_[_k7_8];                                            /* [  9] */
	double k5_6 = p_[_k5_6];                                            /* [ 10] */
	double k6_5 = p_[_k6_5];                                            /* [ 11] */
	double k8_5 = p_[_k8_5];                                            /* [ 12] */
	double k7_6 = p_[_k7_6];                                            /* [ 13] */
	double k6_7 = p_[_k6_7];                                            /* [ 14] */
	double k6_2 = p_[_k6_2];                                            /* [ 15] */
	double k5_8 = p_[_k5_8];                                            /* [ 16] */
	double k3p_7off = p_[_k3p_7off];                                    /* [ 17] */
	double k4p_8off = p_[_k4p_8off];                                    /* [ 18] */
	double k3p_3off = p_[_k3p_3off];                                    /* [ 19] */
	double k4p_4off = p_[_k4p_4off];                                    /* [ 20] */
	double k3p_7on = p_[_k3p_7on];                                      /* [ 21] */
	double k4p_8on = p_[_k4p_8on];                                      /* [ 22] */
	double k3p_3on = p_[_k3p_3on];                                      /* [ 23] */
	double k4p_4on = p_[_k4p_4on];                                      /* [ 24] */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [ 25] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [ 26] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [ 27] */
	double km3OFF = p_[_km3OFF];                                        /* [ 28] */
	double km4OFF = p_[_km4OFF];                                        /* [ 29] */
	double km3ON = p_[_km3ON];                                          /* [ 30] */
	double km4ON = p_[_km4ON];                                          /* [ 31] */
	double KD_T = p_[_KD_T];                                            /* [ 32] */
	double b_AKAP = p_[_b_AKAP];                                        /* [ 33] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [ 34] */
	double CaN_ConservedConst = p_[_CaN_ConservedConst];                /* [ 35] */
	double Rii_C_ConservedConst = p_[_Rii_C_ConservedConst];            /* [ 36] */
	double cAMP_ConservedConst = p_[_cAMP_ConservedConst];              /* [ 37] */
	double Rii_ConservedConst = p_[_Rii_ConservedConst];                /* [ 38] */
/* 	state variables   */
	double RiiP = y_[_RiiP];                                            /* [  0] */
	double RiiP_cAMP = y_[_RiiP_cAMP];                                  /* [  1] */
	double RiiP_C = y_[_RiiP_C];                                        /* [  2] */
	double RiiP_C_cAMP = y_[_RiiP_C_cAMP];                              /* [  3] */
	double C = y_[_C];                                                  /* [  4] */
	double Rii_cAMP = y_[_Rii_cAMP];                                    /* [  5] */
	double Rii_C_cAMP = y_[_Rii_C_cAMP];                                /* [  6] */
	double RiiP_CaN = y_[_RiiP_CaN];                                    /* [  7] */
	double RiiP_cAMP_CaN = y_[_RiiP_cAMP_CaN];                          /* [  8] */
	double AKAR4_C = y_[_AKAR4_C];                                      /* [  9] */
	double AKAR4p = y_[_AKAR4p];                                        /* [ 10] */
/* 	expressions   */
	double k3p_7 = b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	double k4p_4 = b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	double k4p_8 = b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	double k3p_3 = b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
	double k4_4p = b_AKAP * ((k4p_4 + k4p_8)/km4ON ) + (1 - b_AKAP) * (k4p_4 + k4p_8) / km4OFF;
	double k3_3p = b_AKAP * ((k3p_3 + k3p_7)/km3ON) + (1 - b_AKAP) * (k3p_3 + k3p_7)/km3OFF;
	double k2_1 = k1_2 * KD_T;
	double AKAR4 = (AKAR4_ConservedConst - (AKAR4_C+AKAR4p));
	double CaN = (CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN));
	double Rii_C = (Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	double cAMP = (cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	double Rii = (Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	double reaction_51 = k5_1*Rii_C;
	double reaction_14 = k4_1*RiiP*C - k1_4*RiiP_C;
	double reaction_12 = k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	double reaction_43 = k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	double reaction_23 = k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	double reaction_78 = k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	double reaction_56 = k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	double reaction_76 = k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	double reaction_62 = k6_2*Rii_C_cAMP;
	double reaction_58 = k8_5*Rii*C - k5_8*Rii_C;
	double reaction_44 = k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	double reaction_33 = k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	double reaction_48 = k4p_8*RiiP_CaN;
	double reaction_37 = k3p_7*RiiP_cAMP_CaN;
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	memset(funcJacp_,0,sizeof(double)*39); /* initialize with 0.0 */
	return GSL_SUCCESS;
}

int AKAP79_default(double t, double *p_){
	if (!p_) return numParam;
/* 	constants   */
	memset(p_,0,sizeof(double)*39); /* initialize with 0.0 */
	p_[_k5_1] = 46.5411365826226;
	p_[_k1_2] = 1.24347737295428;
	p_[_k3_2] = 0.00328439415820813;
	p_[_k2_3] = 0.587853369688012;
	p_[_k3_4] = 0.0730805556930531;
	p_[_k4_3] = 0.00712225662720034;
	p_[_k4_1] = 0.00519774644439109;
	p_[_k1_4] = 0.00117992862984651;
	p_[_k8_7] = 0.00708844744872235;
	p_[_k7_8] = 0.0833086316093731;
	p_[_k5_6] = 0.665243673375593;
	p_[_k6_5] = 0.594408977425014;
	p_[_k8_5] = 0.100291450155288;
	p_[_k7_6] = 0.0979128290361247;
	p_[_k6_7] = 0.0222138094412742;
	p_[_k6_2] = 33.8372184382292;
	p_[_k5_8] = 0.000200941903978572;
	p_[_k3p_7off] = 0.294046396805502;
	p_[_k4p_8off] = 0.294046396805502;
	p_[_k3p_3off] = 14.9166596278814;
	p_[_k4p_4off] = 14.9166596278814;
	p_[_k3p_7on] = 0.43067428164868;
	p_[_k4p_8on] = 0.43067428164868;
	p_[_k3p_3on] = 1.72247710096425;
	p_[_k4p_4on] = 1.72247710096425;
	p_[_kf_C_AKAR4] = 0.0180933753586079;
	p_[_kb_C_AKAR4] = 0.104602112392241;
	p_[_kcat_AKARp] = 10.1811826795126;
	p_[_km3OFF] = 102.235674360709;
	p_[_km4OFF] = 102.235674360709;
	p_[_km3ON] = 0.986951983065571;
	p_[_km4ON] = 0.986951983065571;
	p_[_KD_T] = 0.667954721425184;
	p_[_AKAR4_ConservedConst] = 0.2;
	p_[_CaN_ConservedConst] = 1.5;
	p_[_Rii_C_ConservedConst] = 0.63;
	p_[_Rii_ConservedConst] = 6.3;
	return GSL_SUCCESS;
}

int AKAP79_init(double t, double *y_, void *par){
	double *p_=par;
	if (!y_ || !y_) return 11;
/* 	constants   */
/* 	parameter values   */
	double k5_1 = p_[_k5_1];                                            /* [  0] */
	double k1_2 = p_[_k1_2];                                            /* [  1] */
	double k3_2 = p_[_k3_2];                                            /* [  2] */
	double k2_3 = p_[_k2_3];                                            /* [  3] */
	double k3_4 = p_[_k3_4];                                            /* [  4] */
	double k4_3 = p_[_k4_3];                                            /* [  5] */
	double k4_1 = p_[_k4_1];                                            /* [  6] */
	double k1_4 = p_[_k1_4];                                            /* [  7] */
	double k8_7 = p_[_k8_7];                                            /* [  8] */
	double k7_8 = p_[_k7_8];                                            /* [  9] */
	double k5_6 = p_[_k5_6];                                            /* [ 10] */
	double k6_5 = p_[_k6_5];                                            /* [ 11] */
	double k8_5 = p_[_k8_5];                                            /* [ 12] */
	double k7_6 = p_[_k7_6];                                            /* [ 13] */
	double k6_7 = p_[_k6_7];                                            /* [ 14] */
	double k6_2 = p_[_k6_2];                                            /* [ 15] */
	double k5_8 = p_[_k5_8];                                            /* [ 16] */
	double k3p_7off = p_[_k3p_7off];                                    /* [ 17] */
	double k4p_8off = p_[_k4p_8off];                                    /* [ 18] */
	double k3p_3off = p_[_k3p_3off];                                    /* [ 19] */
	double k4p_4off = p_[_k4p_4off];                                    /* [ 20] */
	double k3p_7on = p_[_k3p_7on];                                      /* [ 21] */
	double k4p_8on = p_[_k4p_8on];                                      /* [ 22] */
	double k3p_3on = p_[_k3p_3on];                                      /* [ 23] */
	double k4p_4on = p_[_k4p_4on];                                      /* [ 24] */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [ 25] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [ 26] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [ 27] */
	double km3OFF = p_[_km3OFF];                                        /* [ 28] */
	double km4OFF = p_[_km4OFF];                                        /* [ 29] */
	double km3ON = p_[_km3ON];                                          /* [ 30] */
	double km4ON = p_[_km4ON];                                          /* [ 31] */
	double KD_T = p_[_KD_T];                                            /* [ 32] */
	double b_AKAP = p_[_b_AKAP];                                        /* [ 33] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [ 34] */
	double CaN_ConservedConst = p_[_CaN_ConservedConst];                /* [ 35] */
	double Rii_C_ConservedConst = p_[_Rii_C_ConservedConst];            /* [ 36] */
	double cAMP_ConservedConst = p_[_cAMP_ConservedConst];              /* [ 37] */
	double Rii_ConservedConst = p_[_Rii_ConservedConst];                /* [ 38] */
	memset(y_,0,sizeof(double)*11); /* initialize with 0.0 */
	return GSL_SUCCESS;
}

