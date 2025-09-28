#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

/* string.h for memset() */
enum stateVariable { _Rii,_cAMP,_RiiP,_Rii_C,_RiiP_cAMP,_RiiP_C,_RiiP_C_cAMP,_C,_Rii_cAMP,_Rii_C_cAMP,_CaN,_RiiP_CaN,_RiiP_cAMP_CaN,_AKAR4,_AKAR4_C,_AKAR4p, numStateVar }; /* state variable indexes  */
enum param { _k5_1,_k1_2,_k3_2,_k2_3,_k3_4,_k4_3,_k4_1,_k1_4,_k8_7,_k7_8,_k5_6,_k6_5,_k8_5,_k7_6,_k6_7,_k6_2,_k5_8,_k3p_7off,_k4p_8off,_k3p_3off,_k4p_4off,_k3p_7on,_k4p_8on,_k3p_3on,_k4p_4on,_kf_C_AKAR4,_kb_C_AKAR4,_kcat_AKARp,_kmOFF,_kmON,_KD_T,_b_AKAP, numParam }; /* parameter indexes  */
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
	double k3p_7off=p_[17];
	double k4p_8off=p_[18];
	double k3p_3off=p_[19];
	double k4p_4off=p_[20];
	double k3p_7on=p_[21];
	double k4p_8on=p_[22];
	double k3p_3on=p_[23];
	double k4p_4on=p_[24];
	double kf_C_AKAR4=p_[25];
	double kb_C_AKAR4=p_[26];
	double kcat_AKARp=p_[27];
	double kmOFF=p_[28];
	double kmON=p_[29];
	double KD_T=p_[30];
	double b_AKAP=p_[31];
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
	double k3p_7=b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	double k4p_4=b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	double k4p_8=b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	double k3p_3=b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
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
	double k3p_7off = p_[17];
	double k4p_8off = p_[18];
	double k3p_3off = p_[19];
	double k4p_4off = p_[20];
	double k3p_7on = p_[21];
	double k4p_8on = p_[22];
	double k3p_3on = p_[23];
	double k4p_4on = p_[24];
	double kf_C_AKAR4 = p_[25];
	double kb_C_AKAR4 = p_[26];
	double kcat_AKARp = p_[27];
	double kmOFF = p_[28];
	double kmON = p_[29];
	double KD_T = p_[30];
	double b_AKAP = p_[31];
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
	double k3p_7 = b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	double k4p_4 = b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	double k4p_8 = b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	double k3p_3 = b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
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
	double k3p_7off = p_[17];
	double k4p_8off = p_[18];
	double k3p_3off = p_[19];
	double k4p_4off = p_[20];
	double k3p_7on = p_[21];
	double k4p_8on = p_[22];
	double k3p_3on = p_[23];
	double k4p_4on = p_[24];
	double kf_C_AKAR4 = p_[25];
	double kb_C_AKAR4 = p_[26];
	double kcat_AKARp = p_[27];
	double kmOFF = p_[28];
	double kmON = p_[29];
	double KD_T = p_[30];
	double b_AKAP = p_[31];
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
	double k3p_7 = b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	double k4p_4 = b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	double k4p_8 = b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	double k3p_3 = b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
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
	double k3p_7off = p_[17];
	double k4p_8off = p_[18];
	double k3p_3off = p_[19];
	double k4p_4off = p_[20];
	double k3p_7on = p_[21];
	double k4p_8on = p_[22];
	double k3p_3on = p_[23];
	double k4p_4on = p_[24];
	double kf_C_AKAR4 = p_[25];
	double kb_C_AKAR4 = p_[26];
	double kcat_AKARp = p_[27];
	double kmOFF = p_[28];
	double kmON = p_[29];
	double KD_T = p_[30];
	double b_AKAP = p_[31];
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
	double k3p_7 = b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	double k4p_4 = b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	double k4p_8 = b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	double k3p_3 = b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
	double k4_4p = b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF;
	double k3_3p = b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF;
	double k2_1 = k1_2 * KD_T;
	flux[0] = 0.0; // k5_1*Rii_C
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
	double k3p_7off=p_[17];
	double k4p_8off=p_[18];
	double k3p_3off=p_[19];
	double k4p_4off=p_[20];
	double k3p_7on=p_[21];
	double k4p_8on=p_[22];
	double k3p_3on=p_[23];
	double k4p_4on=p_[24];
	double kf_C_AKAR4=p_[25];
	double kb_C_AKAR4=p_[26];
	double kcat_AKARp=p_[27];
	double kmOFF=p_[28];
	double kmON=p_[29];
	double KD_T=p_[30];
	double b_AKAP=p_[31];
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
	double k3p_7=b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	double k4p_4=b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	double k4p_8=b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	double k3p_3=b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
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
/* column 2 (df/dy_1) */
/* column 3 (df/dy_2) */
/* column 4 (df/dy_3) */
/* column 5 (df/dy_4) */
/* column 6 (df/dy_5) */
/* column 7 (df/dy_6) */
/* column 8 (df/dy_7) */
/* column 9 (df/dy_8) */
/* column 10 (df/dy_9) */
/* column 11 (df/dy_10) */
/* column 12 (df/dy_11) */
/* column 13 (df/dy_12) */
/* column 14 (df/dy_13) */
/* column 15 (df/dy_14) */
/* column 16 (df/dy_15) */
	return GSL_SUCCESS;
}
/* ode parameter Jacobian df(t,y;p)/dp */
int AKAP79_jacp(double t, const double y_[], double *jacp_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jacp_) return 16*32;
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
	double k3p_7off=p_[17];
	double k4p_8off=p_[18];
	double k3p_3off=p_[19];
	double k4p_4off=p_[20];
	double k3p_7on=p_[21];
	double k4p_8on=p_[22];
	double k3p_3on=p_[23];
	double k4p_4on=p_[24];
	double kf_C_AKAR4=p_[25];
	double kb_C_AKAR4=p_[26];
	double kcat_AKARp=p_[27];
	double kmOFF=p_[28];
	double kmON=p_[29];
	double KD_T=p_[30];
	double b_AKAP=p_[31];
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
	double k3p_7=b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	double k4p_4=b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	double k4p_8=b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	double k3p_3=b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
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
	memset(jacp_,0,sizeof(double)*numStateVar*numParam); /* 512 */
/* column 1 (df/dp_0) */
/* column 2 (df/dp_1) */
/* column 3 (df/dp_2) */
/* column 4 (df/dp_3) */
/* column 5 (df/dp_4) */
/* column 6 (df/dp_5) */
/* column 7 (df/dp_6) */
/* column 8 (df/dp_7) */
/* column 9 (df/dp_8) */
/* column 10 (df/dp_9) */
/* column 11 (df/dp_10) */
/* column 12 (df/dp_11) */
/* column 13 (df/dp_12) */
/* column 14 (df/dp_13) */
/* column 15 (df/dp_14) */
/* column 16 (df/dp_15) */
/* column 17 (df/dp_16) */
/* column 18 (df/dp_17) */
/* column 19 (df/dp_18) */
/* column 20 (df/dp_19) */
/* column 21 (df/dp_20) */
/* column 22 (df/dp_21) */
/* column 23 (df/dp_22) */
/* column 24 (df/dp_23) */
/* column 25 (df/dp_24) */
/* column 26 (df/dp_25) */
/* column 27 (df/dp_26) */
/* column 28 (df/dp_27) */
/* column 29 (df/dp_28) */
/* column 30 (df/dp_29) */
/* column 31 (df/dp_30) */
/* column 32 (df/dp_31) */
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
	double k3p_7off=p_[17];
	double k4p_8off=p_[18];
	double k3p_3off=p_[19];
	double k4p_4off=p_[20];
	double k3p_7on=p_[21];
	double k4p_8on=p_[22];
	double k3p_3on=p_[23];
	double k4p_4on=p_[24];
	double kf_C_AKAR4=p_[25];
	double kb_C_AKAR4=p_[26];
	double kcat_AKARp=p_[27];
	double kmOFF=p_[28];
	double kmON=p_[29];
	double KD_T=p_[30];
	double b_AKAP=p_[31];
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
	double k3p_7=b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	double k4p_4=b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	double k4p_8=b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	double k3p_3=b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
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
	double k3p_7off=p_[17];
	double k4p_8off=p_[18];
	double k3p_3off=p_[19];
	double k4p_4off=p_[20];
	double k3p_7on=p_[21];
	double k4p_8on=p_[22];
	double k3p_3on=p_[23];
	double k4p_4on=p_[24];
	double kf_C_AKAR4=p_[25];
	double kb_C_AKAR4=p_[26];
	double kcat_AKARp=p_[27];
	double kmOFF=p_[28];
	double kmON=p_[29];
	double KD_T=p_[30];
	double b_AKAP=p_[31];
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
	double k3p_7=b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	double k4p_4=b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	double k4p_8=b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	double k3p_3=b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
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
	return GSL_SUCCESS;
}
/* Function parameter Jacobian dF(t,y;p)/dp */
int AKAP79_funcJacp(double t, const double y_[], double *funcJacp_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJacp_) return 32;
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
	double k3p_7off=p_[17];
	double k4p_8off=p_[18];
	double k3p_3off=p_[19];
	double k4p_4off=p_[20];
	double k3p_7on=p_[21];
	double k4p_8on=p_[22];
	double k3p_3on=p_[23];
	double k4p_4on=p_[24];
	double kf_C_AKAR4=p_[25];
	double kb_C_AKAR4=p_[26];
	double kcat_AKARp=p_[27];
	double kmOFF=p_[28];
	double kmON=p_[29];
	double KD_T=p_[30];
	double b_AKAP=p_[31];
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
	double k3p_7=b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	double k4p_4=b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	double k4p_8=b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	double k3p_3=b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
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
	memset(funcJacp_,0,sizeof(double)*numFunc*numParam); /* 32 */
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
/* column 29 (dF/dp_28) */
/* column 30 (dF/dp_29) */
/* column 31 (dF/dp_30) */
/* column 32 (dF/dp_31) */
	return GSL_SUCCESS;
}
/* ode default parameters */
int AKAP79_default(double t, void *par)
{
	double *p_=par;
	if (!p_) return 32;
	memset(p_,0,sizeof(double)*numParam);
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
	p_[_kmOFF] = 102.235674360709;
	p_[_kmON] = 0.986951983065571;
	p_[_KD_T] = 0.667954721425184;
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
	double k3p_7off=p_[17];
	double k4p_8off=p_[18];
	double k3p_3off=p_[19];
	double k4p_4off=p_[20];
	double k3p_7on=p_[21];
	double k4p_8on=p_[22];
	double k3p_3on=p_[23];
	double k4p_4on=p_[24];
	double kf_C_AKAR4=p_[25];
	double kb_C_AKAR4=p_[26];
	double kcat_AKARp=p_[27];
	double kmOFF=p_[28];
	double kmON=p_[29];
	double KD_T=p_[30];
	double b_AKAP=p_[31];
	/* the initial value of y may depend on the parameters. */
	memset(y_,0,sizeof(double)*numStateVar);
	y_[_Rii] = 6.3;
	y_[_Rii_C] = 0.63;
	y_[_CaN] = 1.5;
	y_[_AKAR4] = 0.2;
	return GSL_SUCCESS;
}
