# ODE vector field: y' = f(t,y;p)
AKAP79_vf <- function(t, state, parameters) {
##	constants
##	parameter values
	k5_1 = parameters[1];
	k1_2 = parameters[2];
	k3_2 = parameters[3];
	k2_3 = parameters[4];
	k3_4 = parameters[5];
	k4_3 = parameters[6];
	k4_1 = parameters[7];
	k8_7 = parameters[8];
	k7_8 = parameters[9];
	k5_6 = parameters[10];
	k6_5 = parameters[11];
	k8_5 = parameters[12];
	k7_6 = parameters[13];
	k6_7 = parameters[14];
	k6_2 = parameters[15];
	k3p_7off = parameters[16];
	k4p_8off = parameters[17];
	k3p_3off = parameters[18];
	k4p_4off = parameters[19];
	k3p_7on = parameters[20];
	k4p_8on = parameters[21];
	k3p_3on = parameters[22];
	k4p_4on = parameters[23];
	kf_C_AKAR4 = parameters[24];
	kb_C_AKAR4 = parameters[25];
	kcat_AKARp = parameters[26];
	km3OFF = parameters[27];
	km4OFF = parameters[28];
	km3ON = parameters[29];
	km4ON = parameters[30];
	KD_T = parameters[31];
	b_AKAP = parameters[32];
	AKAR4_ConservedConst = parameters[33];
	CaN_ConservedConst = parameters[34];
	Rii_C_ConservedConst = parameters[35];
	cAMP_ConservedConst = parameters[36];
	Rii_ConservedConst = parameters[37];
##	state variables
	RiiP = state[1];
	RiiP_cAMP = state[2];
	RiiP_C = state[3];
	RiiP_C_cAMP = state[4];
	C = state[5];
	Rii_cAMP = state[6];
	Rii_C_cAMP = state[7];
	RiiP_CaN = state[8];
	RiiP_cAMP_CaN = state[9];
	AKAR4_C = state[10];
	AKAR4p = state[11];
##	expressions
	k3p_7 = b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	k4p_4 = b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	k4p_8 = b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	k3p_3 = b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
	k4_4p = b_AKAP * ((k4p_4 + k4p_8)/km4ON ) + (1 - b_AKAP) * (k4p_4 + k4p_8) / km4OFF;
	k3_3p = b_AKAP * ((k3p_3 + k3p_7)/km3ON) + (1 - b_AKAP) * (k3p_3 + k3p_7)/km3OFF;
	k2_1 = k1_2 * KD_T;
	k1_4 = k1_2 * k2_3 * k3_4 * k4_1 / (k4_3 * k3_2 * k2_1);
	k5_8 = k5_6 * k6_7 * k7_8 * k8_5 / (k8_7 * k7_6 * k6_5);
	AKAR4 = (AKAR4_ConservedConst - (AKAR4_C+AKAR4p));
	CaN = (CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN));
	Rii_C = (Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	cAMP = (cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	Rii = (Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	reaction_51 = k5_1*Rii_C;
	reaction_14 = k4_1*RiiP*C - k1_4*RiiP_C;
	reaction_12 = k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	reaction_43 = k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	reaction_23 = k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	reaction_78 = k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	reaction_56 = k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	reaction_76 = k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	reaction_62 = k6_2*Rii_C_cAMP;
	reaction_58 = k8_5*Rii*C - k5_8*Rii_C;
	reaction_44 = k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	reaction_33 = k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	reaction_48 = k4p_8*RiiP_CaN;
	reaction_37 = k3p_7*RiiP_cAMP_CaN;
	reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	reaction_2 = kcat_AKARp*AKAR4_C;
	f_ <- numeric(11)
	f_[1] <- -reaction_14-reaction_43-reaction_44
	f_[2] <- +reaction_43-reaction_23-reaction_33
	f_[3] <- +reaction_51+reaction_14-reaction_12
	f_[4] <- +reaction_12+reaction_23+reaction_62
	f_[5] <- -reaction_14-reaction_23-reaction_76-reaction_58-reaction_1+reaction_2
	f_[6] <- +reaction_78-reaction_76+reaction_37
	f_[7] <- +reaction_56+reaction_76-reaction_62
	f_[8] <- +reaction_44-reaction_48
	f_[9] <- +reaction_33-reaction_37
	f_[10] <- +reaction_1-reaction_2
	f_[11] <- +reaction_2
	names(f_) <- c("RiiP", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4_C", "AKAR4p")
	return(list(f_))
}

# ODE Jacobian: df(t,y;p)/dy
AKAP79_jac <- function(t, state, parameters) {
##	constants
##	parameter values
	k5_1 = parameters[1];
	k1_2 = parameters[2];
	k3_2 = parameters[3];
	k2_3 = parameters[4];
	k3_4 = parameters[5];
	k4_3 = parameters[6];
	k4_1 = parameters[7];
	k8_7 = parameters[8];
	k7_8 = parameters[9];
	k5_6 = parameters[10];
	k6_5 = parameters[11];
	k8_5 = parameters[12];
	k7_6 = parameters[13];
	k6_7 = parameters[14];
	k6_2 = parameters[15];
	k3p_7off = parameters[16];
	k4p_8off = parameters[17];
	k3p_3off = parameters[18];
	k4p_4off = parameters[19];
	k3p_7on = parameters[20];
	k4p_8on = parameters[21];
	k3p_3on = parameters[22];
	k4p_4on = parameters[23];
	kf_C_AKAR4 = parameters[24];
	kb_C_AKAR4 = parameters[25];
	kcat_AKARp = parameters[26];
	km3OFF = parameters[27];
	km4OFF = parameters[28];
	km3ON = parameters[29];
	km4ON = parameters[30];
	KD_T = parameters[31];
	b_AKAP = parameters[32];
	AKAR4_ConservedConst = parameters[33];
	CaN_ConservedConst = parameters[34];
	Rii_C_ConservedConst = parameters[35];
	cAMP_ConservedConst = parameters[36];
	Rii_ConservedConst = parameters[37];
##	state variables
	RiiP = state[1];
	RiiP_cAMP = state[2];
	RiiP_C = state[3];
	RiiP_C_cAMP = state[4];
	C = state[5];
	Rii_cAMP = state[6];
	Rii_C_cAMP = state[7];
	RiiP_CaN = state[8];
	RiiP_cAMP_CaN = state[9];
	AKAR4_C = state[10];
	AKAR4p = state[11];
##	expressions
	k3p_7 = b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	k4p_4 = b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	k4p_8 = b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	k3p_3 = b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
	k4_4p = b_AKAP * ((k4p_4 + k4p_8)/km4ON ) + (1 - b_AKAP) * (k4p_4 + k4p_8) / km4OFF;
	k3_3p = b_AKAP * ((k3p_3 + k3p_7)/km3ON) + (1 - b_AKAP) * (k3p_3 + k3p_7)/km3OFF;
	k2_1 = k1_2 * KD_T;
	k1_4 = k1_2 * k2_3 * k3_4 * k4_1 / (k4_3 * k3_2 * k2_1);
	k5_8 = k5_6 * k6_7 * k7_8 * k8_5 / (k8_7 * k7_6 * k6_5);
	AKAR4 = (AKAR4_ConservedConst - (AKAR4_C+AKAR4p));
	CaN = (CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN));
	Rii_C = (Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	cAMP = (cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	Rii = (Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	reaction_51 = k5_1*Rii_C;
	reaction_14 = k4_1*RiiP*C - k1_4*RiiP_C;
	reaction_12 = k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	reaction_43 = k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	reaction_23 = k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	reaction_78 = k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	reaction_56 = k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	reaction_76 = k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	reaction_62 = k6_2*Rii_C_cAMP;
	reaction_58 = k8_5*Rii*C - k5_8*Rii_C;
	reaction_44 = k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	reaction_33 = k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	reaction_48 = k4p_8*RiiP_CaN;
	reaction_37 = k3p_7*RiiP_cAMP_CaN;
	reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	reaction_2 = kcat_AKARp*AKAR4_C;
	jac_ <- matrix(0.0,11,11)
	jac_[1,1] <- -(k4_1*C+k4_3*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))+((b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON+((1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN)))
	jac_[2,1] <- k4_3*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))
	jac_[3,1] <- k4_1*C
	jac_[5,1] <- k8_5*C-k4_1*C
	jac_[6,1] <- -k8_7*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))
	jac_[8,1] <- ((b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON+((1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))
	jac_[1,2] <- k4_3*RiiP+k3_4
	jac_[2,2] <- -(k4_3*RiiP+k3_4+k3_2*C+((b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON+((1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN)))
	jac_[3,2] <- k1_2*RiiP_C
	jac_[4,2] <- k3_2*C-k1_2*RiiP_C
	jac_[5,2] <- k8_5*C-k3_2*C
	jac_[6,2] <- -(k8_7*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))+k8_7*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN)))
	jac_[7,2] <- -k5_6*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))
	jac_[9,2] <- ((b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON+((1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))
	jac_[1,3] <- (k1_2*k2_3*k3_4*k4_1)/(k4_3*k3_2*k1_2*KD_T)
	jac_[3,3] <- -((k1_2*k2_3*k3_4*k4_1)/(k4_3*k3_2*k1_2*KD_T)+k5_1+k1_2*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN)))
	jac_[4,3] <- k1_2*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))
	jac_[5,3] <- (k1_2*k2_3*k3_4*k4_1)/(k4_3*k3_2*k1_2*KD_T)-(k5_6*k6_7*k7_8*k8_5)/(k8_7*k7_6*k6_5)
	jac_[7,3] <- -k5_6*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))
	jac_[1,4] <- k4_3*RiiP
	jac_[2,4] <- k2_3-k4_3*RiiP
	jac_[3,4] <- k1_2*RiiP_C+k1_2*KD_T-k5_1
	jac_[4,4] <- -(k2_3+k1_2*RiiP_C+k1_2*KD_T)
	jac_[5,4] <- k2_3-(k5_6*k6_7*k7_8*k8_5)/(k8_7*k7_6*k6_5)
	jac_[6,4] <- -k8_7*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))
	jac_[7,4] <- -(k5_6*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))+k5_6*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C)))
	jac_[1,5] <- -k4_1*RiiP
	jac_[2,5] <- -k3_2*RiiP_cAMP
	jac_[3,5] <- k4_1*RiiP-k5_1
	jac_[4,5] <- k3_2*RiiP_cAMP
	jac_[5,5] <- -(k4_1*RiiP+k3_2*RiiP_cAMP+k7_6*Rii_cAMP+k8_5*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))+k8_5*C-(-k5_6*k6_7*k7_8*k8_5)/(k8_7*k7_6*k6_5)+kf_C_AKAR4*(AKAR4_ConservedConst-(AKAR4_C+AKAR4p)))
	jac_[6,5] <- k8_7*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))-k7_6*Rii_cAMP
	jac_[7,5] <- k7_6*Rii_cAMP-k5_6*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))
	jac_[10,5] <- kf_C_AKAR4*(AKAR4_ConservedConst-(AKAR4_C+AKAR4p))
	jac_[1,6] <- k4_3*RiiP
	jac_[2,6] <- -k4_3*RiiP
	jac_[3,6] <- k1_2*RiiP_C
	jac_[4,6] <- -k1_2*RiiP_C
	jac_[5,6] <- k8_5*C-k7_6*C
	jac_[6,6] <- -(k8_7*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))+k8_7*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))+k7_8+k7_6*C)
	jac_[7,6] <- k7_6*C-k5_6*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))
	jac_[1,7] <- k4_3*RiiP
	jac_[2,7] <- -k4_3*RiiP
	jac_[3,7] <- k1_2*RiiP_C-k5_1
	jac_[4,7] <- k6_2-k1_2*RiiP_C
	jac_[5,7] <- k6_7-(k5_6*k6_7*k7_8*k8_5)/(k8_7*k7_6*k6_5)
	jac_[6,7] <- k6_7-k8_7*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))
	jac_[7,7] <- -(k6_7+k5_6*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))+k5_6*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))+k6_5+k6_2)
	jac_[1,8] <- ((b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON+((1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*RiiP+b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off
	jac_[2,8] <- ((b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON+((1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*RiiP_cAMP
	jac_[5,8] <- k8_5*C
	jac_[6,8] <- -k8_7*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))
	jac_[8,8] <- -(((b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON+((1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*RiiP+b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off)
	jac_[9,8] <- -((b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON+((1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*RiiP_cAMP
	jac_[1,9] <- k4_3*RiiP+((b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON+((1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*RiiP
	jac_[2,9] <- ((b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON+((1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*RiiP_cAMP+b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off-k4_3*RiiP
	jac_[3,9] <- k1_2*RiiP_C
	jac_[4,9] <- -k1_2*RiiP_C
	jac_[5,9] <- k8_5*C
	jac_[6,9] <- b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off-(k8_7*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))+k8_7*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN)))
	jac_[7,9] <- -k5_6*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))
	jac_[8,9] <- -((b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON+((1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*RiiP
	jac_[9,9] <- -(((b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON+((1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*RiiP_cAMP+b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off)
	jac_[3,10] <- -k5_1
	jac_[5,10] <- (-k5_6*k6_7*k7_8*k8_5)/(k8_7*k7_6*k6_5)-k8_5*C+kf_C_AKAR4*C+kb_C_AKAR4+kcat_AKARp
	jac_[6,10] <- k8_7*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))
	jac_[7,10] <- -k5_6*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))
	jac_[10,10] <- -(kf_C_AKAR4*C+kb_C_AKAR4+kcat_AKARp)
	jac_[11,10] <- kcat_AKARp
	jac_[5,11] <- kf_C_AKAR4*C
	jac_[10,11] <- -kf_C_AKAR4*C
## no names for return value
	return(jac_)
}

# ODE parameter Jacobian: df(t,y;p)/dp
AKAP79_jacp <- function(t, state, parameters) {
##	constants
##	parameter values
	k5_1 = parameters[1];
	k1_2 = parameters[2];
	k3_2 = parameters[3];
	k2_3 = parameters[4];
	k3_4 = parameters[5];
	k4_3 = parameters[6];
	k4_1 = parameters[7];
	k8_7 = parameters[8];
	k7_8 = parameters[9];
	k5_6 = parameters[10];
	k6_5 = parameters[11];
	k8_5 = parameters[12];
	k7_6 = parameters[13];
	k6_7 = parameters[14];
	k6_2 = parameters[15];
	k3p_7off = parameters[16];
	k4p_8off = parameters[17];
	k3p_3off = parameters[18];
	k4p_4off = parameters[19];
	k3p_7on = parameters[20];
	k4p_8on = parameters[21];
	k3p_3on = parameters[22];
	k4p_4on = parameters[23];
	kf_C_AKAR4 = parameters[24];
	kb_C_AKAR4 = parameters[25];
	kcat_AKARp = parameters[26];
	km3OFF = parameters[27];
	km4OFF = parameters[28];
	km3ON = parameters[29];
	km4ON = parameters[30];
	KD_T = parameters[31];
	b_AKAP = parameters[32];
	AKAR4_ConservedConst = parameters[33];
	CaN_ConservedConst = parameters[34];
	Rii_C_ConservedConst = parameters[35];
	cAMP_ConservedConst = parameters[36];
	Rii_ConservedConst = parameters[37];
##	state variables
	RiiP = state[1];
	RiiP_cAMP = state[2];
	RiiP_C = state[3];
	RiiP_C_cAMP = state[4];
	C = state[5];
	Rii_cAMP = state[6];
	Rii_C_cAMP = state[7];
	RiiP_CaN = state[8];
	RiiP_cAMP_CaN = state[9];
	AKAR4_C = state[10];
	AKAR4p = state[11];
##	expressions
	k3p_7 = b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	k4p_4 = b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	k4p_8 = b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	k3p_3 = b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
	k4_4p = b_AKAP * ((k4p_4 + k4p_8)/km4ON ) + (1 - b_AKAP) * (k4p_4 + k4p_8) / km4OFF;
	k3_3p = b_AKAP * ((k3p_3 + k3p_7)/km3ON) + (1 - b_AKAP) * (k3p_3 + k3p_7)/km3OFF;
	k2_1 = k1_2 * KD_T;
	k1_4 = k1_2 * k2_3 * k3_4 * k4_1 / (k4_3 * k3_2 * k2_1);
	k5_8 = k5_6 * k6_7 * k7_8 * k8_5 / (k8_7 * k7_6 * k6_5);
	AKAR4 = (AKAR4_ConservedConst - (AKAR4_C+AKAR4p));
	CaN = (CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN));
	Rii_C = (Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	cAMP = (cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	Rii = (Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	reaction_51 = k5_1*Rii_C;
	reaction_14 = k4_1*RiiP*C - k1_4*RiiP_C;
	reaction_12 = k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	reaction_43 = k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	reaction_23 = k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	reaction_78 = k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	reaction_56 = k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	reaction_76 = k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	reaction_62 = k6_2*Rii_C_cAMP;
	reaction_58 = k8_5*Rii*C - k5_8*Rii_C;
	reaction_44 = k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	reaction_33 = k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	reaction_48 = k4p_8*RiiP_CaN;
	reaction_37 = k3p_7*RiiP_cAMP_CaN;
	reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	reaction_2 = kcat_AKARp*AKAR4_C;
	jacp_ <- matrix(0.0,11,37)
	jacp_[3,1] <- Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C)
	jacp_[1,2] <- (k4_3*k3_2*k1_2*KD_T*RiiP_C*k2_3*k3_4*k4_1-RiiP_C*k1_2*k2_3*k3_4*k4_1*k4_3*k3_2*KD_T)/(k4_3*k3_2*k1_2*KD_T)^2
	jacp_[3,2] <- -((k4_3*k3_2*k1_2*KD_T*RiiP_C*k2_3*k3_4*k4_1-RiiP_C*k1_2*k2_3*k3_4*k4_1*k4_3*k3_2*KD_T)/(k4_3*k3_2*k1_2*KD_T)^2+RiiP_C*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))-KD_T*RiiP_C_cAMP)
	jacp_[4,2] <- RiiP_C*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))-KD_T*RiiP_C_cAMP
	jacp_[5,2] <- (k4_3*k3_2*k1_2*KD_T*RiiP_C*k2_3*k3_4*k4_1-RiiP_C*k1_2*k2_3*k3_4*k4_1*k4_3*k3_2*KD_T)/(k4_3*k3_2*k1_2*KD_T)^2
	jacp_[1,3] <- (-RiiP_C*k1_2*k2_3*k3_4*k4_1*k4_3*k1_2*KD_T)/(k4_3*k3_2*k1_2*KD_T)^2
	jacp_[2,3] <- -RiiP_cAMP*C
	jacp_[3,3] <- (RiiP_C*k1_2*k2_3*k3_4*k4_1*k4_3*k1_2*KD_T)/(k4_3*k3_2*k1_2*KD_T)^2
	jacp_[4,3] <- RiiP_cAMP*C
	jacp_[5,3] <- (-RiiP_C*k1_2*k2_3*k3_4*k4_1*k4_3*k1_2*KD_T)/(k4_3*k3_2*k1_2*KD_T)^2-RiiP_cAMP*C
	jacp_[1,4] <- (RiiP_C*k1_2*k3_4*k4_1)/(k4_3*k3_2*k1_2*KD_T)
	jacp_[2,4] <- RiiP_C_cAMP
	jacp_[3,4] <- -(RiiP_C*k1_2*k3_4*k4_1)/(k4_3*k3_2*k1_2*KD_T)
	jacp_[4,4] <- -RiiP_C_cAMP
	jacp_[5,4] <- (RiiP_C*k1_2*k3_4*k4_1)/(k4_3*k3_2*k1_2*KD_T)+RiiP_C_cAMP
	jacp_[1,5] <- (RiiP_C*k1_2*k2_3*k4_1)/(k4_3*k3_2*k1_2*KD_T)+RiiP_cAMP
	jacp_[2,5] <- -RiiP_cAMP
	jacp_[3,5] <- -(RiiP_C*k1_2*k2_3*k4_1)/(k4_3*k3_2*k1_2*KD_T)
	jacp_[5,5] <- (RiiP_C*k1_2*k2_3*k4_1)/(k4_3*k3_2*k1_2*KD_T)
	jacp_[1,6] <- (-RiiP_C*k1_2*k2_3*k3_4*k4_1*k3_2*k1_2*KD_T)/(k4_3*k3_2*k1_2*KD_T)^2-(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))*RiiP
	jacp_[2,6] <- (cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))*RiiP
	jacp_[3,6] <- (RiiP_C*k1_2*k2_3*k3_4*k4_1*k3_2*k1_2*KD_T)/(k4_3*k3_2*k1_2*KD_T)^2
	jacp_[5,6] <- (-RiiP_C*k1_2*k2_3*k3_4*k4_1*k3_2*k1_2*KD_T)/(k4_3*k3_2*k1_2*KD_T)^2
	jacp_[1,7] <- (RiiP_C*k1_2*k2_3*k3_4)/(k4_3*k3_2*k1_2*KD_T)-RiiP*C
	jacp_[3,7] <- RiiP*C-(RiiP_C*k1_2*k2_3*k3_4)/(k4_3*k3_2*k1_2*KD_T)
	jacp_[5,7] <- (RiiP_C*k1_2*k2_3*k3_4)/(k4_3*k3_2*k1_2*KD_T)-RiiP*C
	jacp_[5,8] <- -((Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))*k5_6*k6_7*k7_8*k8_5*k7_6*k6_5)/(k8_7*k7_6*k6_5)^2
	jacp_[6,8] <- (cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))
	jacp_[5,9] <- ((Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))*k5_6*k6_7*k8_5)/(k8_7*k7_6*k6_5)
	jacp_[6,9] <- -Rii_cAMP
	jacp_[5,10] <- ((Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))*k6_7*k7_8*k8_5)/(k8_7*k7_6*k6_5)
	jacp_[7,10] <- (Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))
	jacp_[5,11] <- -((Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))*k5_6*k6_7*k7_8*k8_5*k8_7*k7_6)/(k8_7*k7_6*k6_5)^2
	jacp_[7,11] <- -Rii_C_cAMP
	jacp_[5,12] <- ((Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))*k5_6*k6_7*k7_8)/(k8_7*k7_6*k6_5)-(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))*C
	jacp_[5,13] <- -(Rii_cAMP*C+((Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))*k5_6*k6_7*k7_8*k8_5*k8_7*k6_5)/(k8_7*k7_6*k6_5)^2)
	jacp_[6,13] <- -Rii_cAMP*C
	jacp_[7,13] <- Rii_cAMP*C
	jacp_[5,14] <- Rii_C_cAMP+((Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))*k5_6*k7_8*k8_5)/(k8_7*k7_6*k6_5)
	jacp_[6,14] <- Rii_C_cAMP
	jacp_[7,14] <- -Rii_C_cAMP
	jacp_[4,15] <- Rii_C_cAMP
	jacp_[7,15] <- -Rii_C_cAMP
	jacp_[2,16] <- -((b_AKAP*(1-b_AKAP))/km3ON+(1-b_AKAP)^2/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP
	jacp_[6,16] <- (1-b_AKAP)*RiiP_cAMP_CaN
	jacp_[9,16] <- ((b_AKAP*(1-b_AKAP))/km3ON+(1-b_AKAP)^2/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP-(1-b_AKAP)*RiiP_cAMP_CaN
	jacp_[1,17] <- -((b_AKAP*(1-b_AKAP))/km4ON+(1-b_AKAP)^2/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))
	jacp_[8,17] <- ((b_AKAP*(1-b_AKAP))/km4ON+(1-b_AKAP)^2/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))-(1-b_AKAP)*RiiP_CaN
	jacp_[2,18] <- (1-b_AKAP)*RiiP_cAMP_CaN-((b_AKAP*(1-b_AKAP))/km3ON+(1-b_AKAP)^2/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP
	jacp_[9,18] <- ((b_AKAP*(1-b_AKAP))/km3ON+(1-b_AKAP)^2/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP-(1-b_AKAP)*RiiP_cAMP_CaN
	jacp_[1,19] <- (1-b_AKAP)*RiiP_CaN-((b_AKAP*(1-b_AKAP))/km4ON+(1-b_AKAP)^2/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))
	jacp_[8,19] <- ((b_AKAP*(1-b_AKAP))/km4ON+(1-b_AKAP)^2/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))-(1-b_AKAP)*RiiP_CaN
	jacp_[2,20] <- -(b_AKAP^2/km3ON+((1-b_AKAP)*b_AKAP)/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP
	jacp_[6,20] <- b_AKAP*RiiP_cAMP_CaN
	jacp_[9,20] <- (b_AKAP^2/km3ON+((1-b_AKAP)*b_AKAP)/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP-b_AKAP*RiiP_cAMP_CaN
	jacp_[1,21] <- -(b_AKAP^2/km4ON+((1-b_AKAP)*b_AKAP)/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))
	jacp_[8,21] <- (b_AKAP^2/km4ON+((1-b_AKAP)*b_AKAP)/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))-b_AKAP*RiiP_CaN
	jacp_[2,22] <- b_AKAP*RiiP_cAMP_CaN-(b_AKAP^2/km3ON+((1-b_AKAP)*b_AKAP)/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP
	jacp_[9,22] <- (b_AKAP^2/km3ON+((1-b_AKAP)*b_AKAP)/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP-b_AKAP*RiiP_cAMP_CaN
	jacp_[1,23] <- b_AKAP*RiiP_CaN-(b_AKAP^2/km4ON+((1-b_AKAP)*b_AKAP)/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))
	jacp_[8,23] <- (b_AKAP^2/km4ON+((1-b_AKAP)*b_AKAP)/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))-b_AKAP*RiiP_CaN
	jacp_[5,24] <- -C*(AKAR4_ConservedConst-(AKAR4_C+AKAR4p))
	jacp_[10,24] <- C*(AKAR4_ConservedConst-(AKAR4_C+AKAR4p))
	jacp_[5,25] <- AKAR4_C
	jacp_[10,25] <- -AKAR4_C
	jacp_[5,26] <- AKAR4_C
	jacp_[10,26] <- -AKAR4_C
	jacp_[11,26] <- AKAR4_C
	jacp_[2,27] <- (RiiP_cAMP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*(1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF^2
	jacp_[9,27] <- (-RiiP_cAMP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*(1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF^2
	jacp_[1,28] <- ((CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP*(1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF^2
	jacp_[8,28] <- (-(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP*(1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF^2
	jacp_[2,29] <- (RiiP_cAMP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON^2
	jacp_[9,29] <- (-RiiP_cAMP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON^2
	jacp_[1,30] <- ((CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP*b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON^2
	jacp_[8,30] <- (-(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP*b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON^2
	jacp_[1,31] <- (-RiiP_C*k1_2*k2_3*k3_4*k4_1*k4_3*k3_2*k1_2)/(k4_3*k3_2*k1_2*KD_T)^2
	jacp_[3,31] <- (RiiP_C*k1_2*k2_3*k3_4*k4_1*k4_3*k3_2*k1_2)/(k4_3*k3_2*k1_2*KD_T)^2+k1_2*RiiP_C_cAMP
	jacp_[4,31] <- -k1_2*RiiP_C_cAMP
	jacp_[5,31] <- (-RiiP_C*k1_2*k2_3*k3_4*k4_1*k4_3*k3_2*k1_2)/(k4_3*k3_2*k1_2*KD_T)^2
	jacp_[1,32] <- (k4p_4on-k4p_4off)*RiiP_CaN-((b_AKAP*(k4p_4on-k4p_4off+k4p_8on-k4p_8off)+b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off)/km4ON+((1-b_AKAP)*(k4p_4on-k4p_4off+k4p_8on-k4p_8off)-(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))
	jacp_[2,32] <- (k3p_3on-k3p_3off)*RiiP_cAMP_CaN-((b_AKAP*(k3p_3on-k3p_3off+k3p_7on-k3p_7off)+b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off)/km3ON+((1-b_AKAP)*(k3p_3on-k3p_3off+k3p_7on-k3p_7off)-(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP
	jacp_[6,32] <- (k3p_7on-k3p_7off)*RiiP_cAMP_CaN
	jacp_[8,32] <- ((b_AKAP*(k4p_4on-k4p_4off+k4p_8on-k4p_8off)+b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off)/km4ON+((1-b_AKAP)*(k4p_4on-k4p_4off+k4p_8on-k4p_8off)-(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))-(k4p_4on-k4p_4off)*RiiP_CaN-(k4p_8on-k4p_8off)*RiiP_CaN
	jacp_[9,32] <- ((b_AKAP*(k3p_3on-k3p_3off+k3p_7on-k3p_7off)+b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off)/km3ON+((1-b_AKAP)*(k3p_3on-k3p_3off+k3p_7on-k3p_7off)-(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP-(k3p_3on-k3p_3off)*RiiP_cAMP_CaN-(k3p_7on-k3p_7off)*RiiP_cAMP_CaN
	jacp_[5,33] <- -kf_C_AKAR4*C
	jacp_[10,33] <- kf_C_AKAR4*C
	jacp_[1,34] <- -((b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON+((1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*RiiP
	jacp_[2,34] <- -((b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON+((1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*RiiP_cAMP
	jacp_[8,34] <- ((b_AKAP*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4ON+((1-b_AKAP)*(b_AKAP*k4p_4on+(1-b_AKAP)*k4p_4off+b_AKAP*k4p_8on+(1-b_AKAP)*k4p_8off))/km4OFF)*RiiP
	jacp_[9,34] <- ((b_AKAP*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3ON+((1-b_AKAP)*(b_AKAP*k3p_3on+(1-b_AKAP)*k3p_3off+b_AKAP*k3p_7on+(1-b_AKAP)*k3p_7off))/km3OFF)*RiiP_cAMP
	jacp_[3,35] <- k5_1
	jacp_[5,35] <- (k5_6*k6_7*k7_8*k8_5)/(k8_7*k7_6*k6_5)
	jacp_[7,35] <- k5_6*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))
	jacp_[1,36] <- -k4_3*RiiP
	jacp_[2,36] <- k4_3*RiiP
	jacp_[3,36] <- -k1_2*RiiP_C
	jacp_[4,36] <- k1_2*RiiP_C
	jacp_[6,36] <- k8_7*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))
	jacp_[7,36] <- k5_6*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))
	jacp_[5,37] <- -k8_5*C
	jacp_[6,37] <- k8_7*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))
## no names for return value
	return(jacp_)
}

# Output Function (Observables)
AKAP79_func <- function(t, state, parameters) {
##	constants
##	parameter values
	k5_1 = parameters[1];
	k1_2 = parameters[2];
	k3_2 = parameters[3];
	k2_3 = parameters[4];
	k3_4 = parameters[5];
	k4_3 = parameters[6];
	k4_1 = parameters[7];
	k8_7 = parameters[8];
	k7_8 = parameters[9];
	k5_6 = parameters[10];
	k6_5 = parameters[11];
	k8_5 = parameters[12];
	k7_6 = parameters[13];
	k6_7 = parameters[14];
	k6_2 = parameters[15];
	k3p_7off = parameters[16];
	k4p_8off = parameters[17];
	k3p_3off = parameters[18];
	k4p_4off = parameters[19];
	k3p_7on = parameters[20];
	k4p_8on = parameters[21];
	k3p_3on = parameters[22];
	k4p_4on = parameters[23];
	kf_C_AKAR4 = parameters[24];
	kb_C_AKAR4 = parameters[25];
	kcat_AKARp = parameters[26];
	km3OFF = parameters[27];
	km4OFF = parameters[28];
	km3ON = parameters[29];
	km4ON = parameters[30];
	KD_T = parameters[31];
	b_AKAP = parameters[32];
	AKAR4_ConservedConst = parameters[33];
	CaN_ConservedConst = parameters[34];
	Rii_C_ConservedConst = parameters[35];
	cAMP_ConservedConst = parameters[36];
	Rii_ConservedConst = parameters[37];
##	state variables
	RiiP = state[1];
	RiiP_cAMP = state[2];
	RiiP_C = state[3];
	RiiP_C_cAMP = state[4];
	C = state[5];
	Rii_cAMP = state[6];
	Rii_C_cAMP = state[7];
	RiiP_CaN = state[8];
	RiiP_cAMP_CaN = state[9];
	AKAR4_C = state[10];
	AKAR4p = state[11];
##	expressions
	k3p_7 = b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off;
	k4p_4 = b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off;
	k4p_8 = b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off;
	k3p_3 = b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off;
	k4_4p = b_AKAP * ((k4p_4 + k4p_8)/km4ON ) + (1 - b_AKAP) * (k4p_4 + k4p_8) / km4OFF;
	k3_3p = b_AKAP * ((k3p_3 + k3p_7)/km3ON) + (1 - b_AKAP) * (k3p_3 + k3p_7)/km3OFF;
	k2_1 = k1_2 * KD_T;
	k1_4 = k1_2 * k2_3 * k3_4 * k4_1 / (k4_3 * k3_2 * k2_1);
	k5_8 = k5_6 * k6_7 * k7_8 * k8_5 / (k8_7 * k7_6 * k6_5);
	AKAR4 = (AKAR4_ConservedConst - (AKAR4_C+AKAR4p));
	CaN = (CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN));
	Rii_C = (Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	cAMP = (cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	Rii = (Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	reaction_51 = k5_1*Rii_C;
	reaction_14 = k4_1*RiiP*C - k1_4*RiiP_C;
	reaction_12 = k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP;
	reaction_43 = k4_3*cAMP*RiiP - k3_4*RiiP_cAMP;
	reaction_23 = k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP;
	reaction_78 = k8_7*cAMP*Rii - k7_8*Rii_cAMP;
	reaction_56 = k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP;
	reaction_76 = k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP;
	reaction_62 = k6_2*Rii_C_cAMP;
	reaction_58 = k8_5*Rii*C - k5_8*Rii_C;
	reaction_44 = k4_4p*RiiP*CaN - k4p_4*RiiP_CaN;
	reaction_33 = k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN;
	reaction_48 = k4p_8*RiiP_CaN;
	reaction_37 = k3p_7*RiiP_cAMP_CaN;
	reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	reaction_2 = kcat_AKARp*AKAR4_C;
	func_ <- numeric(1)
	func_[1] <- (AKAR4p*5)*71.67+100
	names(func_) <- c("AKAR4pOUT")
	return(func_)
}

AKAP79_default <- function(t) {
##	constants
	parameters <- numeric(37)
	parameters[1] <- 46.5411365826226
	parameters[2] <- 1.24347737295428
	parameters[3] <- 0.00328439415820813
	parameters[4] <- 0.587853369688012
	parameters[5] <- 0.0730805556930531
	parameters[6] <- 0.00712225662720034
	parameters[7] <- 0.00519774644439109
	parameters[8] <- 0.00708844744872235
	parameters[9] <- 0.0833086316093731
	parameters[10] <- 0.665243673375593
	parameters[11] <- 0.594408977425014
	parameters[12] <- 0.100291450155288
	parameters[13] <- 0.0979128290361247
	parameters[14] <- 0.0222138094412742
	parameters[15] <- 33.8372184382292
	parameters[16] <- 0.294046396805502
	parameters[17] <- 0.294046396805502
	parameters[18] <- 14.9166596278814
	parameters[19] <- 14.9166596278814
	parameters[20] <- 0.43067428164868
	parameters[21] <- 0.43067428164868
	parameters[22] <- 1.72247710096425
	parameters[23] <- 1.72247710096425
	parameters[24] <- 0.0180933753586079
	parameters[25] <- 0.104602112392241
	parameters[26] <- 10.1811826795126
	parameters[27] <- 102.235674360709
	parameters[28] <- 102.235674360709
	parameters[29] <- 0.986951983065571
	parameters[30] <- 0.986951983065571
	parameters[31] <- 0.667954721425184
	parameters[33] <- 0.2
	parameters[34] <- 1.5
	parameters[35] <- 0.63
	parameters[37] <- 6.3
	names(parameters) <- c("k5_1", "k1_2", "k3_2", "k2_3", "k3_4", "k4_3", "k4_1", "k8_7", "k7_8", "k5_6", "k6_5", "k8_5", "k7_6", "k6_7", "k6_2", "k3p_7off", "k4p_8off", "k3p_3off", "k4p_4off", "k3p_7on", "k4p_8on", "k3p_3on", "k4p_4on", "kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp", "km3OFF", "km4OFF", "km3ON", "km4ON", "KD_T", "b_AKAP", "AKAR4_ConservedConst", "CaN_ConservedConst", "Rii_C_ConservedConst", "cAMP_ConservedConst", "Rii_ConservedConst")
	return(parameters)
}

AKAP79_init <- function(t, parameters) {
##	constants
##	parameter values
	k5_1 = parameters[1];
	k1_2 = parameters[2];
	k3_2 = parameters[3];
	k2_3 = parameters[4];
	k3_4 = parameters[5];
	k4_3 = parameters[6];
	k4_1 = parameters[7];
	k8_7 = parameters[8];
	k7_8 = parameters[9];
	k5_6 = parameters[10];
	k6_5 = parameters[11];
	k8_5 = parameters[12];
	k7_6 = parameters[13];
	k6_7 = parameters[14];
	k6_2 = parameters[15];
	k3p_7off = parameters[16];
	k4p_8off = parameters[17];
	k3p_3off = parameters[18];
	k4p_4off = parameters[19];
	k3p_7on = parameters[20];
	k4p_8on = parameters[21];
	k3p_3on = parameters[22];
	k4p_4on = parameters[23];
	kf_C_AKAR4 = parameters[24];
	kb_C_AKAR4 = parameters[25];
	kcat_AKARp = parameters[26];
	km3OFF = parameters[27];
	km4OFF = parameters[28];
	km3ON = parameters[29];
	km4ON = parameters[30];
	KD_T = parameters[31];
	b_AKAP = parameters[32];
	AKAR4_ConservedConst = parameters[33];
	CaN_ConservedConst = parameters[34];
	Rii_C_ConservedConst = parameters[35];
	cAMP_ConservedConst = parameters[36];
	Rii_ConservedConst = parameters[37];
	state <- numeric(11)
	names(state) <- c("RiiP", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4_C", "AKAR4p")
	return(state)
}

# a variable that collects all functions into one list:
model <- list(vf=AKAP79_vf, jac=AKAP79_jac, jacp=AKAP79_jacp, default=AKAP79_default, init=AKAP79_init, func=AKAP79_func,name='AKAP79')
