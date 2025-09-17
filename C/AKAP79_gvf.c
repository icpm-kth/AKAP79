#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

/* string.h for memset() */
enum stateVariable { _Rii,_cAMP,_RiiP,_Rii_C,_RiiP_cAMP,_RiiP_C,_RiiP_C_cAMP,_C,_Rii_cAMP,_Rii_C_cAMP,_CaN,_RiiP_CaN,_RiiP_cAMP_CaN,_AKAR4,_AKAR4_C,_AKAR4p, numStateVar }; /* state variable indexes  */
enum param { _k5_1,_k1_2,_k3_2,_k2_3,_k3_4,_k4_3,_k4_1,_k1_4,_k8_7,_k7_8,_k5_6,_k6_5,_k8_5,_k7_6,_k6_7,_k6_2,_k5_8,_AKAPoff_1,_AKAPoff_3,_AKAPon_1,_AKAPon_3,_kf_C_AKAR4,_kb_C_AKAR4,_kcat_AKARp,_kmOFF,_kmON,_KD_T,_b_AKAP, numParam }; /* parameter indexes  */
enum func { _AKAR4pOUT, numFunc }; /* parameter indexes  */

/* The error code indicates how to pre-allocate memory
 * for output values such as `f_`. The _vf function returns
 * the number of state variables, if any of the args are `NULL`.
 * evaluation errors can be indicated by negative return values.
 * GSL_SUCCESS (0) is returned when no error occurred.
 */

/* ode vector field: y'=f(t,y;p) */
int AKAP79_vf(double t, const double y_[], double f_[], void *par)
{
	double *p_=par;
	if (!y_ || !f_) return 16;
	double k5_1=p_[0];
	double k1_2=p_[1];
	double k3_2=p_[2];
	double k2_3=p_[3];
	double k3_4=p_[4];
	double k4_3=p_[5];
	double k4_1=p_[6];
	double k1_4=p_[7];
	double k8_7=p_[8];
	double k7_8=p_[9];
	double k5_6=p_[10];
	double k6_5=p_[11];
	double k8_5=p_[12];
	double k7_6=p_[13];
	double k6_7=p_[14];
	double k6_2=p_[15];
	double k5_8=p_[16];
	double AKAPoff_1=p_[17];
	double AKAPoff_3=p_[18];
	double AKAPon_1=p_[19];
	double AKAPon_3=p_[20];
	double kf_C_AKAR4=p_[21];
	double kb_C_AKAR4=p_[22];
	double kcat_AKARp=p_[23];
	double kmOFF=p_[24];
	double kmON=p_[25];
	double KD_T=p_[26];
	double b_AKAP=p_[27];
	double Rii=y_[0];
	double cAMP=y_[1];
	double RiiP=y_[2];
	double Rii_C=y_[3];
	double RiiP_cAMP=y_[4];
	double RiiP_C=y_[5];
	double RiiP_C_cAMP=y_[6];
	double C=y_[7];
	double Rii_cAMP=y_[8];
	double Rii_C_cAMP=y_[9];
	double CaN=y_[10];
	double RiiP_CaN=y_[11];
	double RiiP_cAMP_CaN=y_[12];
	double AKAR4=y_[13];
	double AKAR4_C=y_[14];
	double AKAR4p=y_[15];
	double k3p_7=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k4p_4=b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4p_8=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k3p_3=b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4_4p=b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF;
	double k3_3p=b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF;
	double k2_1=k1_2 * KD_T;
	double reaction_51=k5_1*Rii_C;
	double reaction_14=k4_1*RiiP*C - k1_4*RiiP_C;
	double reaction_12=k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	double reaction_43=k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	double reaction_23=k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	double reaction_78=k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	double reaction_56=k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	double reaction_76=k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	double reaction_62=k6_2*Rii_C_cAMP;
	double reaction_58=k8_5*Rii*C - k5_8*Rii_C;
	double reaction_44=k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	double reaction_33=k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	double reaction_48=k4p_8*RiiP_CaN;
	double reaction_37=k3p_7*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	f_[_Rii] = -reaction_78-reaction_58+reaction_48; /* Rii */
	f_[_cAMP] = -reaction_12-reaction_43-reaction_78-reaction_56; /* cAMP */
	f_[_RiiP] = -reaction_14-reaction_43-reaction_44; /* RiiP */
	f_[_Rii_C] = -reaction_51-reaction_56+reaction_58; /* Rii_C */
	f_[_RiiP_cAMP] = +reaction_43-reaction_23-reaction_33; /* RiiP_cAMP */
	f_[_RiiP_C] = +reaction_51+reaction_14-reaction_12; /* RiiP_C */
	f_[_RiiP_C_cAMP] = +reaction_12+reaction_23+reaction_62; /* RiiP_C_cAMP */
	f_[_C] = -reaction_14-reaction_23-reaction_76-reaction_58-reaction_1+reaction_2; /* C */
	f_[_Rii_cAMP] = +reaction_78-reaction_76+reaction_37; /* Rii_cAMP */
	f_[_Rii_C_cAMP] = +reaction_56+reaction_76-reaction_62; /* Rii_C_cAMP */
	f_[_CaN] = -reaction_44-reaction_33+reaction_48+reaction_37; /* CaN */
	f_[_RiiP_CaN] = +reaction_44-reaction_48; /* RiiP_CaN */
	f_[_RiiP_cAMP_CaN] = +reaction_33-reaction_37; /* RiiP_cAMP_CaN */
	f_[_AKAR4] = -reaction_1; /* AKAR4 */
	f_[_AKAR4_C] = +reaction_1-reaction_2; /* AKAR4_C */
	f_[_AKAR4p] = +reaction_2; /* AKAR4p */
	return GSL_SUCCESS;
}
int AKAP79_netflux(double t, double y_[], double *flux, void *par){
	double *p_=par;
	if (!y_ || !flux) return 16;
	double k5_1 = p_[0];
	double k1_2 = p_[1];
	double k3_2 = p_[2];
	double k2_3 = p_[3];
	double k3_4 = p_[4];
	double k4_3 = p_[5];
	double k4_1 = p_[6];
	double k1_4 = p_[7];
	double k8_7 = p_[8];
	double k7_8 = p_[9];
	double k5_6 = p_[10];
	double k6_5 = p_[11];
	double k8_5 = p_[12];
	double k7_6 = p_[13];
	double k6_7 = p_[14];
	double k6_2 = p_[15];
	double k5_8 = p_[16];
	double AKAPoff_1 = p_[17];
	double AKAPoff_3 = p_[18];
	double AKAPon_1 = p_[19];
	double AKAPon_3 = p_[20];
	double kf_C_AKAR4 = p_[21];
	double kb_C_AKAR4 = p_[22];
	double kcat_AKARp = p_[23];
	double kmOFF = p_[24];
	double kmON = p_[25];
	double KD_T = p_[26];
	double b_AKAP = p_[27];
	double Rii = y_[0];
	double cAMP = y_[1];
	double RiiP = y_[2];
	double Rii_C = y_[3];
	double RiiP_cAMP = y_[4];
	double RiiP_C = y_[5];
	double RiiP_C_cAMP = y_[6];
	double C = y_[7];
	double Rii_cAMP = y_[8];
	double Rii_C_cAMP = y_[9];
	double CaN = y_[10];
	double RiiP_CaN = y_[11];
	double RiiP_cAMP_CaN = y_[12];
	double AKAR4 = y_[13];
	double AKAR4_C = y_[14];
	double AKAR4p = y_[15];
	double k3p_7 = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k4p_4 = b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4p_8 = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k3p_3 = b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4_4p = b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF;
	double k3_3p = b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF;
	double k2_1 = k1_2 * KD_T;
	flux[0] = k5_1*Rii_C;
	flux[1] = k4_1*RiiP*C - k1_4*RiiP_C;
	flux[2] = k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	flux[3] = k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	flux[4] = k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	flux[5] = k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	flux[6] = k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	flux[7] = k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	flux[8] = k6_2*Rii_C_cAMP;
	flux[9] = k8_5*Rii*C - k5_8*Rii_C;
	flux[10] = k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	flux[11] = k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	flux[12] = k4p_8*RiiP_CaN;
	flux[13] = k3p_7*RiiP_cAMP_CaN;
	flux[14] = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	flux[15] = kcat_AKARp*AKAR4_C;
	return GSL_SUCCESS;
}

int AKAP79_fwdflux(double t, double y_[], double *flux, void *par){
	double *p_=par;
	if (!y_ || !flux) return 16;
	double k5_1 = p_[0];
	double k1_2 = p_[1];
	double k3_2 = p_[2];
	double k2_3 = p_[3];
	double k3_4 = p_[4];
	double k4_3 = p_[5];
	double k4_1 = p_[6];
	double k1_4 = p_[7];
	double k8_7 = p_[8];
	double k7_8 = p_[9];
	double k5_6 = p_[10];
	double k6_5 = p_[11];
	double k8_5 = p_[12];
	double k7_6 = p_[13];
	double k6_7 = p_[14];
	double k6_2 = p_[15];
	double k5_8 = p_[16];
	double AKAPoff_1 = p_[17];
	double AKAPoff_3 = p_[18];
	double AKAPon_1 = p_[19];
	double AKAPon_3 = p_[20];
	double kf_C_AKAR4 = p_[21];
	double kb_C_AKAR4 = p_[22];
	double kcat_AKARp = p_[23];
	double kmOFF = p_[24];
	double kmON = p_[25];
	double KD_T = p_[26];
	double b_AKAP = p_[27];
	double Rii = y_[0];
	double cAMP = y_[1];
	double RiiP = y_[2];
	double Rii_C = y_[3];
	double RiiP_cAMP = y_[4];
	double RiiP_C = y_[5];
	double RiiP_C_cAMP = y_[6];
	double C = y_[7];
	double Rii_cAMP = y_[8];
	double Rii_C_cAMP = y_[9];
	double CaN = y_[10];
	double RiiP_CaN = y_[11];
	double RiiP_cAMP_CaN = y_[12];
	double AKAR4 = y_[13];
	double AKAR4_C = y_[14];
	double AKAR4p = y_[15];
	double k3p_7 = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k4p_4 = b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4p_8 = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k3p_3 = b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4_4p = b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF;
	double k3_3p = b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF;
	double k2_1 = k1_2 * KD_T;
	// k5_1*Rii_C
	flux[0] = k5_1*Rii_C;
	// k4_1*RiiP*C - k1_4*RiiP_C
	flux[1] = k4_1*RiiP*C ;
	// k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP
	flux[2] = k1_2*RiiP_C*cAMP ;
	// k4_3*cAMP*RiiP - k3_4*RiiP_cAMP
	flux[3] = k4_3*cAMP*RiiP ;
	// k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP
	flux[4] = k3_2*RiiP_cAMP*C ;
	// k8_7*cAMP*Rii - k7_8*Rii_cAMP
	flux[5] = k8_7*cAMP*Rii ;
	// k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP
	flux[6] = k5_6*Rii_C*cAMP ;
	// k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP
	flux[7] = k7_6*Rii_cAMP*C ;
	// k6_2*Rii_C_cAMP
	flux[8] = k6_2*Rii_C_cAMP;
	// k8_5*Rii*C - k5_8*Rii_C
	flux[9] = k8_5*Rii*C ;
	// k4_4p*RiiP*CaN - k4p_4*RiiP_CaN
	flux[10] = k4_4p*RiiP*CaN ;
	// k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN
	flux[11] = k3_3p*CaN*RiiP_cAMP ;
	// k4p_8*RiiP_CaN
	flux[12] = k4p_8*RiiP_CaN;
	// k3p_7*RiiP_cAMP_CaN
	flux[13] = k3p_7*RiiP_cAMP_CaN;
	// kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	flux[14] = kf_C_AKAR4*C*AKAR4 ;
	// kcat_AKARp*AKAR4_C
	flux[15] = kcat_AKARp*AKAR4_C;
	return GSL_SUCCESS;
}

int AKAP79_bwdflux(double t, double y_[], double *flux, void *par){
	double *p_=par;
	if (!y_ || !flux) return 16;
	double k5_1 = p_[0];
	double k1_2 = p_[1];
	double k3_2 = p_[2];
	double k2_3 = p_[3];
	double k3_4 = p_[4];
	double k4_3 = p_[5];
	double k4_1 = p_[6];
	double k1_4 = p_[7];
	double k8_7 = p_[8];
	double k7_8 = p_[9];
	double k5_6 = p_[10];
	double k6_5 = p_[11];
	double k8_5 = p_[12];
	double k7_6 = p_[13];
	double k6_7 = p_[14];
	double k6_2 = p_[15];
	double k5_8 = p_[16];
	double AKAPoff_1 = p_[17];
	double AKAPoff_3 = p_[18];
	double AKAPon_1 = p_[19];
	double AKAPon_3 = p_[20];
	double kf_C_AKAR4 = p_[21];
	double kb_C_AKAR4 = p_[22];
	double kcat_AKARp = p_[23];
	double kmOFF = p_[24];
	double kmON = p_[25];
	double KD_T = p_[26];
	double b_AKAP = p_[27];
	double Rii = y_[0];
	double cAMP = y_[1];
	double RiiP = y_[2];
	double Rii_C = y_[3];
	double RiiP_cAMP = y_[4];
	double RiiP_C = y_[5];
	double RiiP_C_cAMP = y_[6];
	double C = y_[7];
	double Rii_cAMP = y_[8];
	double Rii_C_cAMP = y_[9];
	double CaN = y_[10];
	double RiiP_CaN = y_[11];
	double RiiP_cAMP_CaN = y_[12];
	double AKAR4 = y_[13];
	double AKAR4_C = y_[14];
	double AKAR4p = y_[15];
	double k3p_7 = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k4p_4 = b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4p_8 = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k3p_3 = b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4_4p = b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF;
	double k3_3p = b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF;
	double k2_1 = k1_2 * KD_T;
	flux[0] = 0.0; // k5_1*Rii_C
	flux[1] =  k1_4*RiiP_C; // k4_1*RiiP*C - k1_4*RiiP_C
	flux[2] =  k2_1*RiiP_C_cAMP; // k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP
	flux[3] =  k3_4*RiiP_cAMP; // k4_3*cAMP*RiiP - k3_4*RiiP_cAMP
	flux[4] =  k2_3*RiiP_C_cAMP; // k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP
	flux[5] =  k7_8*Rii_cAMP; // k8_7*cAMP*Rii - k7_8*Rii_cAMP
	flux[6] =  k6_5*Rii_C_cAMP; // k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP
	flux[7] =  k6_7*Rii_C_cAMP; // k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP
	flux[8] = 0.0; // k6_2*Rii_C_cAMP
	flux[9] =  k5_8*Rii_C; // k8_5*Rii*C - k5_8*Rii_C
	flux[10] =  k4p_4*RiiP_CaN; // k4_4p*RiiP*CaN - k4p_4*RiiP_CaN
	flux[11] =  k3p_3*RiiP_cAMP_CaN; // k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN
	flux[12] = 0.0; // k4p_8*RiiP_CaN
	flux[13] = 0.0; // k3p_7*RiiP_cAMP_CaN
	flux[14] =  kb_C_AKAR4*AKAR4_C; // kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	flux[15] = 0.0; // kcat_AKARp*AKAR4_C
	return GSL_SUCCESS;
}

/* ode Jacobian df(t,y;p)/dy */
int AKAP79_jac(double t, const double y_[], double *jac_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jac_) return 16*16;
	double k5_1=p_[0];
	double k1_2=p_[1];
	double k3_2=p_[2];
	double k2_3=p_[3];
	double k3_4=p_[4];
	double k4_3=p_[5];
	double k4_1=p_[6];
	double k1_4=p_[7];
	double k8_7=p_[8];
	double k7_8=p_[9];
	double k5_6=p_[10];
	double k6_5=p_[11];
	double k8_5=p_[12];
	double k7_6=p_[13];
	double k6_7=p_[14];
	double k6_2=p_[15];
	double k5_8=p_[16];
	double AKAPoff_1=p_[17];
	double AKAPoff_3=p_[18];
	double AKAPon_1=p_[19];
	double AKAPon_3=p_[20];
	double kf_C_AKAR4=p_[21];
	double kb_C_AKAR4=p_[22];
	double kcat_AKARp=p_[23];
	double kmOFF=p_[24];
	double kmON=p_[25];
	double KD_T=p_[26];
	double b_AKAP=p_[27];
	double Rii=y_[0];
	double cAMP=y_[1];
	double RiiP=y_[2];
	double Rii_C=y_[3];
	double RiiP_cAMP=y_[4];
	double RiiP_C=y_[5];
	double RiiP_C_cAMP=y_[6];
	double C=y_[7];
	double Rii_cAMP=y_[8];
	double Rii_C_cAMP=y_[9];
	double CaN=y_[10];
	double RiiP_CaN=y_[11];
	double RiiP_cAMP_CaN=y_[12];
	double AKAR4=y_[13];
	double AKAR4_C=y_[14];
	double AKAR4p=y_[15];
	double k3p_7=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k4p_4=b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4p_8=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k3p_3=b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4_4p=b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF;
	double k3_3p=b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF;
	double k2_1=k1_2 * KD_T;
	double reaction_51=k5_1*Rii_C;
	double reaction_14=k4_1*RiiP*C - k1_4*RiiP_C;
	double reaction_12=k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	double reaction_43=k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	double reaction_23=k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	double reaction_78=k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	double reaction_56=k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	double reaction_76=k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	double reaction_62=k6_2*Rii_C_cAMP;
	double reaction_58=k8_5*Rii*C - k5_8*Rii_C;
	double reaction_44=k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	double reaction_33=k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	double reaction_48=k4p_8*RiiP_CaN;
	double reaction_37=k3p_7*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	memset(jac_,0,sizeof(double)*numStateVar*numStateVar); /* 256 */
/* column 1 (df/dy_0) */
	jac_[0] = (-cAMP*k8_7)-C*k8_5; /* [0, 0] */
	jac_[16] = -cAMP*k8_7; /* [1, 0] */
	jac_[48] = C*k8_5; /* [3, 0] */
	jac_[112] = -C*k8_5; /* [7, 0] */
	jac_[128] = cAMP*k8_7; /* [8, 0] */
/* column 2 (df/dy_1) */
	jac_[1] = -Rii*k8_7; /* [0, 1] */
	jac_[17] = (-Rii*k8_7)-Rii_C*k5_6-RiiP*k4_3-RiiP_C*k1_2; /* [1, 1] */
	jac_[33] = -RiiP*k4_3; /* [2, 1] */
	jac_[49] = -Rii_C*k5_6; /* [3, 1] */
	jac_[65] = RiiP*k4_3; /* [4, 1] */
	jac_[81] = -RiiP_C*k1_2; /* [5, 1] */
	jac_[97] = RiiP_C*k1_2; /* [6, 1] */
	jac_[129] = Rii*k8_7; /* [8, 1] */
	jac_[145] = Rii_C*k5_6; /* [9, 1] */
/* column 3 (df/dy_2) */
	jac_[18] = -cAMP*k4_3; /* [1, 2] */
	jac_[34] = (-CaN*k4_4p)-cAMP*k4_3-C*k4_1; /* [2, 2] */
	jac_[66] = cAMP*k4_3; /* [4, 2] */
	jac_[82] = C*k4_1; /* [5, 2] */
	jac_[114] = -C*k4_1; /* [7, 2] */
	jac_[162] = -CaN*k4_4p; /* [10, 2] */
	jac_[178] = CaN*k4_4p; /* [11, 2] */
/* column 4 (df/dy_3) */
	jac_[3] = k5_8; /* [0, 3] */
	jac_[19] = -cAMP*k5_6; /* [1, 3] */
	jac_[51] = (-k5_8)-cAMP*k5_6-k5_1; /* [3, 3] */
	jac_[83] = k5_1; /* [5, 3] */
	jac_[115] = k5_8; /* [7, 3] */
	jac_[147] = cAMP*k5_6; /* [9, 3] */
/* column 5 (df/dy_4) */
	jac_[20] = k3_4; /* [1, 4] */
	jac_[36] = k3_4; /* [2, 4] */
	jac_[68] = (-k3_4)-CaN*k3_3p-C*k3_2; /* [4, 4] */
	jac_[100] = C*k3_2; /* [6, 4] */
	jac_[116] = -C*k3_2; /* [7, 4] */
	jac_[164] = -CaN*k3_3p; /* [10, 4] */
	jac_[196] = CaN*k3_3p; /* [12, 4] */
/* column 6 (df/dy_5) */
	jac_[21] = -cAMP*k1_2; /* [1, 5] */
	jac_[37] = k1_4; /* [2, 5] */
	jac_[85] = (-k1_4)-cAMP*k1_2; /* [5, 5] */
	jac_[101] = cAMP*k1_2; /* [6, 5] */
	jac_[117] = k1_4; /* [7, 5] */
/* column 7 (df/dy_6) */
	jac_[22] = k2_1; /* [1, 6] */
	jac_[70] = k2_3; /* [4, 6] */
	jac_[86] = k2_1; /* [5, 6] */
	jac_[102] = (-k2_3)-k2_1; /* [6, 6] */
	jac_[118] = k2_3; /* [7, 6] */
/* column 8 (df/dy_7) */
	jac_[7] = -Rii*k8_5; /* [0, 7] */
	jac_[39] = -RiiP*k4_1; /* [2, 7] */
	jac_[55] = Rii*k8_5; /* [3, 7] */
	jac_[71] = -RiiP_cAMP*k3_2; /* [4, 7] */
	jac_[87] = RiiP*k4_1; /* [5, 7] */
	jac_[103] = RiiP_cAMP*k3_2; /* [6, 7] */
	jac_[119] = (-AKAR4*kf_C_AKAR4)-Rii*k8_5-Rii_cAMP*k7_6-RiiP*k4_1-RiiP_cAMP*k3_2; /* [7, 7] */
	jac_[135] = -Rii_cAMP*k7_6; /* [8, 7] */
	jac_[151] = Rii_cAMP*k7_6; /* [9, 7] */
	jac_[215] = -AKAR4*kf_C_AKAR4; /* [13, 7] */
	jac_[231] = AKAR4*kf_C_AKAR4; /* [14, 7] */
/* column 9 (df/dy_8) */
	jac_[8] = k7_8; /* [0, 8] */
	jac_[24] = k7_8; /* [1, 8] */
	jac_[120] = -C*k7_6; /* [7, 8] */
	jac_[136] = (-k7_8)-C*k7_6; /* [8, 8] */
	jac_[152] = C*k7_6; /* [9, 8] */
/* column 10 (df/dy_9) */
	jac_[25] = k6_5; /* [1, 9] */
	jac_[57] = k6_5; /* [3, 9] */
	jac_[105] = k6_2; /* [6, 9] */
	jac_[121] = k6_7; /* [7, 9] */
	jac_[137] = k6_7; /* [8, 9] */
	jac_[153] = (-k6_7)-k6_5-k6_2; /* [9, 9] */
/* column 11 (df/dy_10) */
	jac_[42] = -RiiP*k4_4p; /* [2, 10] */
	jac_[74] = -RiiP_cAMP*k3_3p; /* [4, 10] */
	jac_[170] = (-RiiP*k4_4p)-RiiP_cAMP*k3_3p; /* [10, 10] */
	jac_[186] = RiiP*k4_4p; /* [11, 10] */
	jac_[202] = RiiP_cAMP*k3_3p; /* [12, 10] */
/* column 12 (df/dy_11) */
	jac_[11] = k4p_8; /* [0, 11] */
	jac_[43] = k4p_4; /* [2, 11] */
	jac_[171] = k4p_8+k4p_4; /* [10, 11] */
	jac_[187] = (-k4p_8)-k4p_4; /* [11, 11] */
/* column 13 (df/dy_12) */
	jac_[76] = k3p_3; /* [4, 12] */
	jac_[140] = k3p_7; /* [8, 12] */
	jac_[172] = k3p_7+k3p_3; /* [10, 12] */
	jac_[204] = (-k3p_7)-k3p_3; /* [12, 12] */
/* column 14 (df/dy_13) */
	jac_[125] = -C*kf_C_AKAR4; /* [7, 13] */
	jac_[221] = -C*kf_C_AKAR4; /* [13, 13] */
	jac_[237] = C*kf_C_AKAR4; /* [14, 13] */
/* column 15 (df/dy_14) */
	jac_[126] = kcat_AKARp+kb_C_AKAR4; /* [7, 14] */
	jac_[222] = kb_C_AKAR4; /* [13, 14] */
	jac_[238] = (-kcat_AKARp)-kb_C_AKAR4; /* [14, 14] */
	jac_[254] = kcat_AKARp; /* [15, 14] */
/* column 16 (df/dy_15) */
	return GSL_SUCCESS;
}
/* ode parameter Jacobian df(t,y;p)/dp */
int AKAP79_jacp(double t, const double y_[], double *jacp_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jacp_) return 16*28;
	double k5_1=p_[0];
	double k1_2=p_[1];
	double k3_2=p_[2];
	double k2_3=p_[3];
	double k3_4=p_[4];
	double k4_3=p_[5];
	double k4_1=p_[6];
	double k1_4=p_[7];
	double k8_7=p_[8];
	double k7_8=p_[9];
	double k5_6=p_[10];
	double k6_5=p_[11];
	double k8_5=p_[12];
	double k7_6=p_[13];
	double k6_7=p_[14];
	double k6_2=p_[15];
	double k5_8=p_[16];
	double AKAPoff_1=p_[17];
	double AKAPoff_3=p_[18];
	double AKAPon_1=p_[19];
	double AKAPon_3=p_[20];
	double kf_C_AKAR4=p_[21];
	double kb_C_AKAR4=p_[22];
	double kcat_AKARp=p_[23];
	double kmOFF=p_[24];
	double kmON=p_[25];
	double KD_T=p_[26];
	double b_AKAP=p_[27];
	double Rii=y_[0];
	double cAMP=y_[1];
	double RiiP=y_[2];
	double Rii_C=y_[3];
	double RiiP_cAMP=y_[4];
	double RiiP_C=y_[5];
	double RiiP_C_cAMP=y_[6];
	double C=y_[7];
	double Rii_cAMP=y_[8];
	double Rii_C_cAMP=y_[9];
	double CaN=y_[10];
	double RiiP_CaN=y_[11];
	double RiiP_cAMP_CaN=y_[12];
	double AKAR4=y_[13];
	double AKAR4_C=y_[14];
	double AKAR4p=y_[15];
	double k3p_7=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k4p_4=b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4p_8=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k3p_3=b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4_4p=b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF;
	double k3_3p=b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF;
	double k2_1=k1_2 * KD_T;
	double reaction_51=k5_1*Rii_C;
	double reaction_14=k4_1*RiiP*C - k1_4*RiiP_C;
	double reaction_12=k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	double reaction_43=k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	double reaction_23=k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	double reaction_78=k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	double reaction_56=k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	double reaction_76=k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	double reaction_62=k6_2*Rii_C_cAMP;
	double reaction_58=k8_5*Rii*C - k5_8*Rii_C;
	double reaction_44=k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	double reaction_33=k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	double reaction_48=k4p_8*RiiP_CaN;
	double reaction_37=k3p_7*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	memset(jacp_,0,sizeof(double)*numStateVar*numParam); /* 448 */
/* column 1 (df/dp_0) */
	jacp_[84] = -Rii_C; /* [3, 0] */
	jacp_[140] = Rii_C; /* [5, 0] */
/* column 2 (df/dp_1) */
	jacp_[29] = -RiiP_C*cAMP; /* [1, 1] */
	jacp_[141] = -RiiP_C*cAMP; /* [5, 1] */
	jacp_[169] = RiiP_C*cAMP; /* [6, 1] */
/* column 3 (df/dp_2) */
	jacp_[114] = -C*RiiP_cAMP; /* [4, 2] */
	jacp_[170] = C*RiiP_cAMP; /* [6, 2] */
	jacp_[198] = -C*RiiP_cAMP; /* [7, 2] */
/* column 4 (df/dp_3) */
	jacp_[115] = RiiP_C_cAMP; /* [4, 3] */
	jacp_[171] = -RiiP_C_cAMP; /* [6, 3] */
	jacp_[199] = RiiP_C_cAMP; /* [7, 3] */
/* column 5 (df/dp_4) */
	jacp_[32] = RiiP_cAMP; /* [1, 4] */
	jacp_[60] = RiiP_cAMP; /* [2, 4] */
	jacp_[116] = -RiiP_cAMP; /* [4, 4] */
/* column 6 (df/dp_5) */
	jacp_[33] = -RiiP*cAMP; /* [1, 5] */
	jacp_[61] = -RiiP*cAMP; /* [2, 5] */
	jacp_[117] = RiiP*cAMP; /* [4, 5] */
/* column 7 (df/dp_6) */
	jacp_[62] = -C*RiiP; /* [2, 6] */
	jacp_[146] = C*RiiP; /* [5, 6] */
	jacp_[202] = -C*RiiP; /* [7, 6] */
/* column 8 (df/dp_7) */
	jacp_[63] = RiiP_C; /* [2, 7] */
	jacp_[147] = -RiiP_C; /* [5, 7] */
	jacp_[203] = RiiP_C; /* [7, 7] */
/* column 9 (df/dp_8) */
	jacp_[8] = -Rii*cAMP; /* [0, 8] */
	jacp_[36] = -Rii*cAMP; /* [1, 8] */
	jacp_[232] = Rii*cAMP; /* [8, 8] */
/* column 10 (df/dp_9) */
	jacp_[9] = Rii_cAMP; /* [0, 9] */
	jacp_[37] = Rii_cAMP; /* [1, 9] */
	jacp_[233] = -Rii_cAMP; /* [8, 9] */
/* column 11 (df/dp_10) */
	jacp_[38] = -Rii_C*cAMP; /* [1, 10] */
	jacp_[94] = -Rii_C*cAMP; /* [3, 10] */
	jacp_[262] = Rii_C*cAMP; /* [9, 10] */
/* column 12 (df/dp_11) */
	jacp_[39] = Rii_C_cAMP; /* [1, 11] */
	jacp_[95] = Rii_C_cAMP; /* [3, 11] */
	jacp_[263] = -Rii_C_cAMP; /* [9, 11] */
/* column 13 (df/dp_12) */
	jacp_[12] = -C*Rii; /* [0, 12] */
	jacp_[96] = C*Rii; /* [3, 12] */
	jacp_[208] = -C*Rii; /* [7, 12] */
/* column 14 (df/dp_13) */
	jacp_[209] = -C*Rii_cAMP; /* [7, 13] */
	jacp_[237] = -C*Rii_cAMP; /* [8, 13] */
	jacp_[265] = C*Rii_cAMP; /* [9, 13] */
/* column 15 (df/dp_14) */
	jacp_[210] = Rii_C_cAMP; /* [7, 14] */
	jacp_[238] = Rii_C_cAMP; /* [8, 14] */
	jacp_[266] = -Rii_C_cAMP; /* [9, 14] */
/* column 16 (df/dp_15) */
	jacp_[183] = Rii_C_cAMP; /* [6, 15] */
	jacp_[267] = -Rii_C_cAMP; /* [9, 15] */
/* column 17 (df/dp_16) */
	jacp_[16] = Rii_C; /* [0, 16] */
	jacp_[100] = -Rii_C; /* [3, 16] */
	jacp_[212] = Rii_C; /* [7, 16] */
/* column 18 (df/dp_17) */
/* column 19 (df/dp_18) */
/* column 20 (df/dp_19) */
/* column 21 (df/dp_20) */
/* column 22 (df/dp_21) */
	jacp_[217] = -AKAR4*C; /* [7, 21] */
	jacp_[385] = -AKAR4*C; /* [13, 21] */
	jacp_[413] = AKAR4*C; /* [14, 21] */
/* column 23 (df/dp_22) */
	jacp_[218] = AKAR4_C; /* [7, 22] */
	jacp_[386] = AKAR4_C; /* [13, 22] */
	jacp_[414] = -AKAR4_C; /* [14, 22] */
/* column 24 (df/dp_23) */
	jacp_[219] = AKAR4_C; /* [7, 23] */
	jacp_[415] = -AKAR4_C; /* [14, 23] */
	jacp_[443] = AKAR4_C; /* [15, 23] */
/* column 25 (df/dp_24) */
/* column 26 (df/dp_25) */
/* column 27 (df/dp_26) */
/* column 28 (df/dp_27) */
	return GSL_SUCCESS;
}
/* ode Functions F(t,y;p) */
int AKAP79_func(double t, const double y_[], double *func_, void *par)
{
	double *p_=par;
	if (!y_ || !func_) return 1;
	double k5_1=p_[0];
	double k1_2=p_[1];
	double k3_2=p_[2];
	double k2_3=p_[3];
	double k3_4=p_[4];
	double k4_3=p_[5];
	double k4_1=p_[6];
	double k1_4=p_[7];
	double k8_7=p_[8];
	double k7_8=p_[9];
	double k5_6=p_[10];
	double k6_5=p_[11];
	double k8_5=p_[12];
	double k7_6=p_[13];
	double k6_7=p_[14];
	double k6_2=p_[15];
	double k5_8=p_[16];
	double AKAPoff_1=p_[17];
	double AKAPoff_3=p_[18];
	double AKAPon_1=p_[19];
	double AKAPon_3=p_[20];
	double kf_C_AKAR4=p_[21];
	double kb_C_AKAR4=p_[22];
	double kcat_AKARp=p_[23];
	double kmOFF=p_[24];
	double kmON=p_[25];
	double KD_T=p_[26];
	double b_AKAP=p_[27];
	double Rii=y_[0];
	double cAMP=y_[1];
	double RiiP=y_[2];
	double Rii_C=y_[3];
	double RiiP_cAMP=y_[4];
	double RiiP_C=y_[5];
	double RiiP_C_cAMP=y_[6];
	double C=y_[7];
	double Rii_cAMP=y_[8];
	double Rii_C_cAMP=y_[9];
	double CaN=y_[10];
	double RiiP_CaN=y_[11];
	double RiiP_cAMP_CaN=y_[12];
	double AKAR4=y_[13];
	double AKAR4_C=y_[14];
	double AKAR4p=y_[15];
	double k3p_7=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k4p_4=b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4p_8=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k3p_3=b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4_4p=b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF;
	double k3_3p=b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF;
	double k2_1=k1_2 * KD_T;
	double reaction_51=k5_1*Rii_C;
	double reaction_14=k4_1*RiiP*C - k1_4*RiiP_C;
	double reaction_12=k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	double reaction_43=k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	double reaction_23=k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	double reaction_78=k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	double reaction_56=k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	double reaction_76=k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	double reaction_62=k6_2*Rii_C_cAMP;
	double reaction_58=k8_5*Rii*C - k5_8*Rii_C;
	double reaction_44=k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	double reaction_33=k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	double reaction_48=k4p_8*RiiP_CaN;
	double reaction_37=k3p_7*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	func_[_AKAR4pOUT] = (AKAR4p*5)*71.67+100; /* AKAR4pOUT */
	return GSL_SUCCESS;
}
/* Function Jacobian dF(t,y;p)/dy */
int AKAP79_funcJac(double t, const double y_[], double *funcJac_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJac_) return 16;
	double k5_1=p_[0];
	double k1_2=p_[1];
	double k3_2=p_[2];
	double k2_3=p_[3];
	double k3_4=p_[4];
	double k4_3=p_[5];
	double k4_1=p_[6];
	double k1_4=p_[7];
	double k8_7=p_[8];
	double k7_8=p_[9];
	double k5_6=p_[10];
	double k6_5=p_[11];
	double k8_5=p_[12];
	double k7_6=p_[13];
	double k6_7=p_[14];
	double k6_2=p_[15];
	double k5_8=p_[16];
	double AKAPoff_1=p_[17];
	double AKAPoff_3=p_[18];
	double AKAPon_1=p_[19];
	double AKAPon_3=p_[20];
	double kf_C_AKAR4=p_[21];
	double kb_C_AKAR4=p_[22];
	double kcat_AKARp=p_[23];
	double kmOFF=p_[24];
	double kmON=p_[25];
	double KD_T=p_[26];
	double b_AKAP=p_[27];
	double Rii=y_[0];
	double cAMP=y_[1];
	double RiiP=y_[2];
	double Rii_C=y_[3];
	double RiiP_cAMP=y_[4];
	double RiiP_C=y_[5];
	double RiiP_C_cAMP=y_[6];
	double C=y_[7];
	double Rii_cAMP=y_[8];
	double Rii_C_cAMP=y_[9];
	double CaN=y_[10];
	double RiiP_CaN=y_[11];
	double RiiP_cAMP_CaN=y_[12];
	double AKAR4=y_[13];
	double AKAR4_C=y_[14];
	double AKAR4p=y_[15];
	double k3p_7=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k4p_4=b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4p_8=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k3p_3=b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4_4p=b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF;
	double k3_3p=b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF;
	double k2_1=k1_2 * KD_T;
	double reaction_51=k5_1*Rii_C;
	double reaction_14=k4_1*RiiP*C - k1_4*RiiP_C;
	double reaction_12=k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	double reaction_43=k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	double reaction_23=k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	double reaction_78=k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	double reaction_56=k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	double reaction_76=k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	double reaction_62=k6_2*Rii_C_cAMP;
	double reaction_58=k8_5*Rii*C - k5_8*Rii_C;
	double reaction_44=k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	double reaction_33=k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	double reaction_48=k4p_8*RiiP_CaN;
	double reaction_37=k3p_7*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	memset(funcJac_,0,sizeof(double)*numFunc*numStateVar); /* 16 */
/* column 1 (dF/dy_0) */
/* column 2 (dF/dy_1) */
/* column 3 (dF/dy_2) */
/* column 4 (dF/dy_3) */
/* column 5 (dF/dy_4) */
/* column 6 (dF/dy_5) */
/* column 7 (dF/dy_6) */
/* column 8 (dF/dy_7) */
/* column 9 (dF/dy_8) */
/* column 10 (dF/dy_9) */
/* column 11 (dF/dy_10) */
/* column 12 (dF/dy_11) */
/* column 13 (dF/dy_12) */
/* column 14 (dF/dy_13) */
/* column 15 (dF/dy_14) */
/* column 16 (dF/dy_15) */
	funcJac_[15] = 358.35; /* [0, 15] */
	return GSL_SUCCESS;
}
/* Function parameter Jacobian dF(t,y;p)/dp */
int AKAP79_funcJacp(double t, const double y_[], double *funcJacp_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJacp_) return 28;
	double k5_1=p_[0];
	double k1_2=p_[1];
	double k3_2=p_[2];
	double k2_3=p_[3];
	double k3_4=p_[4];
	double k4_3=p_[5];
	double k4_1=p_[6];
	double k1_4=p_[7];
	double k8_7=p_[8];
	double k7_8=p_[9];
	double k5_6=p_[10];
	double k6_5=p_[11];
	double k8_5=p_[12];
	double k7_6=p_[13];
	double k6_7=p_[14];
	double k6_2=p_[15];
	double k5_8=p_[16];
	double AKAPoff_1=p_[17];
	double AKAPoff_3=p_[18];
	double AKAPon_1=p_[19];
	double AKAPon_3=p_[20];
	double kf_C_AKAR4=p_[21];
	double kb_C_AKAR4=p_[22];
	double kcat_AKARp=p_[23];
	double kmOFF=p_[24];
	double kmON=p_[25];
	double KD_T=p_[26];
	double b_AKAP=p_[27];
	double Rii=y_[0];
	double cAMP=y_[1];
	double RiiP=y_[2];
	double Rii_C=y_[3];
	double RiiP_cAMP=y_[4];
	double RiiP_C=y_[5];
	double RiiP_C_cAMP=y_[6];
	double C=y_[7];
	double Rii_cAMP=y_[8];
	double Rii_C_cAMP=y_[9];
	double CaN=y_[10];
	double RiiP_CaN=y_[11];
	double RiiP_cAMP_CaN=y_[12];
	double AKAR4=y_[13];
	double AKAR4_C=y_[14];
	double AKAR4p=y_[15];
	double k3p_7=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k4p_4=b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4p_8=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double k3p_3=b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double k4_4p=b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF;
	double k3_3p=b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF;
	double k2_1=k1_2 * KD_T;
	double reaction_51=k5_1*Rii_C;
	double reaction_14=k4_1*RiiP*C - k1_4*RiiP_C;
	double reaction_12=k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	double reaction_43=k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	double reaction_23=k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	double reaction_78=k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	double reaction_56=k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	double reaction_76=k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	double reaction_62=k6_2*Rii_C_cAMP;
	double reaction_58=k8_5*Rii*C - k5_8*Rii_C;
	double reaction_44=k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	double reaction_33=k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	double reaction_48=k4p_8*RiiP_CaN;
	double reaction_37=k3p_7*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	memset(funcJacp_,0,sizeof(double)*numFunc*numParam); /* 28 */
/* column 1 (dF/dp_0) */
/* column 2 (dF/dp_1) */
/* column 3 (dF/dp_2) */
/* column 4 (dF/dp_3) */
/* column 5 (dF/dp_4) */
/* column 6 (dF/dp_5) */
/* column 7 (dF/dp_6) */
/* column 8 (dF/dp_7) */
/* column 9 (dF/dp_8) */
/* column 10 (dF/dp_9) */
/* column 11 (dF/dp_10) */
/* column 12 (dF/dp_11) */
/* column 13 (dF/dp_12) */
/* column 14 (dF/dp_13) */
/* column 15 (dF/dp_14) */
/* column 16 (dF/dp_15) */
/* column 17 (dF/dp_16) */
/* column 18 (dF/dp_17) */
/* column 19 (dF/dp_18) */
/* column 20 (dF/dp_19) */
/* column 21 (dF/dp_20) */
/* column 22 (dF/dp_21) */
/* column 23 (dF/dp_22) */
/* column 24 (dF/dp_23) */
/* column 25 (dF/dp_24) */
/* column 26 (dF/dp_25) */
/* column 27 (dF/dp_26) */
/* column 28 (dF/dp_27) */
	return GSL_SUCCESS;
}
/* ode default parameters */
int AKAP79_default(double t, void *par)
{
	double *p_=par;
	if (!p_) return 28;
	memset(p_,0,sizeof(double)*numParam);
	p_[_k5_1] = 40.9050728070571;
	p_[_k1_2] = 1.15152189489404;
	p_[_k3_2] = 0.00411137537163148;
	p_[_k2_3] = 0.646565161402751;
	p_[_k3_4] = 0.0745642630670475;
	p_[_k4_3] = 0.0079505887529295;
	p_[_k4_1] = 0.0050721295784356;
	p_[_k1_4] = 0.0014740460115423;
	p_[_k8_7] = 0.00752117649276288;
	p_[_k7_8] = 0.0930377780667662;
	p_[_k5_6] = 0.56644582478049;
	p_[_k6_5] = 0.740455630685478;
	p_[_k8_5] = 0.108133945637089;
	p_[_k7_6] = 0.104354732702695;
	p_[_k6_7] = 0.0194383554779221;
	p_[_k6_2] = 35.3091749693311;
	p_[_k5_8] = 0.000150778070025973;
	p_[_AKAPoff_1] = 0.255654840526295;
	p_[_AKAPoff_3] = 14.8554290677037;
	p_[_AKAPon_1] = 0.451889374844688;
	p_[_AKAPon_3] = 1.67083963931664;
	p_[_kf_C_AKAR4] = 0.0180884373490592;
	p_[_kb_C_AKAR4] = 0.105853209380652;
	p_[_kcat_AKARp] = 10.2388587090355;
	p_[_kmOFF] = 100.472974563026;
	p_[_kmON] = 1.01497947500322;
	p_[_KD_T] = 0.695775064100595;
	return GSL_SUCCESS;
}
/* ode initial values */
int AKAP79_init(double t, double *y_, void *par)
{
	double *p_=par;
	if (!y_) return 16;
	double k5_1=p_[0];
	double k1_2=p_[1];
	double k3_2=p_[2];
	double k2_3=p_[3];
	double k3_4=p_[4];
	double k4_3=p_[5];
	double k4_1=p_[6];
	double k1_4=p_[7];
	double k8_7=p_[8];
	double k7_8=p_[9];
	double k5_6=p_[10];
	double k6_5=p_[11];
	double k8_5=p_[12];
	double k7_6=p_[13];
	double k6_7=p_[14];
	double k6_2=p_[15];
	double k5_8=p_[16];
	double AKAPoff_1=p_[17];
	double AKAPoff_3=p_[18];
	double AKAPon_1=p_[19];
	double AKAPon_3=p_[20];
	double kf_C_AKAR4=p_[21];
	double kb_C_AKAR4=p_[22];
	double kcat_AKARp=p_[23];
	double kmOFF=p_[24];
	double kmON=p_[25];
	double KD_T=p_[26];
	double b_AKAP=p_[27];
	/* the initial value of y may depend on the parameters. */
	memset(y_,0,sizeof(double)*numStateVar);
	y_[_Rii] = 6.3;
	y_[_Rii_C] = 0.63;
	y_[_CaN] = 1.5;
	y_[_AKAR4] = 0.2;
	return GSL_SUCCESS;
}
