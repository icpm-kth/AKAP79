# require("deSolve")

# ode vector field: y'=f(t,y;p)
AKAP79_vf <- function(t, state, parameters)
{
	k5_1 <- parameters[1]
	k1_2 <- parameters[2]
	k3_2 <- parameters[3]
	k2_3 <- parameters[4]
	k3_4 <- parameters[5]
	k4_3 <- parameters[6]
	k4_1 <- parameters[7]
	k1_4 <- parameters[8]
	k8_7 <- parameters[9]
	k7_8 <- parameters[10]
	k5_6 <- parameters[11]
	k6_5 <- parameters[12]
	k8_5 <- parameters[13]
	k7_6 <- parameters[14]
	k6_7 <- parameters[15]
	k6_2 <- parameters[16]
	k5_8 <- parameters[17]
	k3p_7off <- parameters[18]
	k4p_8off <- parameters[19]
	k3p_3off <- parameters[20]
	k4p_4off <- parameters[21]
	k3p_7on <- parameters[22]
	k4p_8on <- parameters[23]
	k3p_3on <- parameters[24]
	k4p_4on <- parameters[25]
	kf_C_AKAR4 <- parameters[26]
	kb_C_AKAR4 <- parameters[27]
	kcat_AKARp <- parameters[28]
	kmOFF <- parameters[29]
	kmON <- parameters[30]
	KD_T <- parameters[31]
	b_AKAP <- parameters[32]
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
	k3p_7 <- b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off
	k4p_4 <- b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off
	k4p_8 <- b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off
	k3p_3 <- b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off
	k4_4p <- b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF
	k3_3p <- b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF
	k2_1 <- k1_2 * KD_T
	reaction_51 <- k5_1*Rii_C
	reaction_14 <- k4_1*RiiP*C - k1_4*RiiP_C
	reaction_12 <- k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP
	reaction_43 <- k4_3*cAMP*RiiP - k3_4*RiiP_cAMP
	reaction_23 <- k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP
	reaction_78 <- k8_7*cAMP*Rii - k7_8*Rii_cAMP
	reaction_56 <- k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP
	reaction_76 <- k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP
	reaction_62 <- k6_2*Rii_C_cAMP
	reaction_58 <- k8_5*Rii*C - k5_8*Rii_C
	reaction_44 <- k4_4p*RiiP*CaN - k4p_4*RiiP_CaN
	reaction_33 <- k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN
	reaction_48 <- k4p_8*RiiP_CaN
	reaction_37 <- k3p_7*RiiP_cAMP_CaN
	reaction_1 <- kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2 <- kcat_AKARp*AKAR4_C
	f_<-vector(mode='numeric',len=16)
	f_[1] <- -reaction_78-reaction_58+reaction_48# Rii
	f_[2] <- -reaction_12-reaction_43-reaction_78-reaction_56# cAMP
	f_[3] <- -reaction_14-reaction_43-reaction_44# RiiP
	f_[4] <- -reaction_51-reaction_56+reaction_58# Rii_C
	f_[5] <- +reaction_43-reaction_23-reaction_33# RiiP_cAMP
	f_[6] <- +reaction_51+reaction_14-reaction_12# RiiP_C
	f_[7] <- +reaction_12+reaction_23+reaction_62# RiiP_C_cAMP
	f_[8] <- -reaction_14-reaction_23-reaction_76-reaction_58-reaction_1+reaction_2# C
	f_[9] <- +reaction_78-reaction_76+reaction_37# Rii_cAMP
	f_[10] <- +reaction_56+reaction_76-reaction_62# Rii_C_cAMP
	f_[11] <- -reaction_44-reaction_33+reaction_48+reaction_37# CaN
	f_[12] <- +reaction_44-reaction_48# RiiP_CaN
	f_[13] <- +reaction_33-reaction_37# RiiP_cAMP_CaN
	f_[14] <- -reaction_1# AKAR4
	f_[15] <- +reaction_1-reaction_2# AKAR4_C
	f_[16] <- +reaction_2# AKAR4p
	names(f_) <- c("Rii", "cAMP", "RiiP", "Rii_C", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "CaN", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4", "AKAR4_C", "AKAR4p")
## for some weird reason deSolve wants this to be a list:
	return(list(f_))
}
AKAP79_netflux <- function(t, state, parameters){
	k5_1 <- parameters[1]
	k1_2 <- parameters[2]
	k3_2 <- parameters[3]
	k2_3 <- parameters[4]
	k3_4 <- parameters[5]
	k4_3 <- parameters[6]
	k4_1 <- parameters[7]
	k1_4 <- parameters[8]
	k8_7 <- parameters[9]
	k7_8 <- parameters[10]
	k5_6 <- parameters[11]
	k6_5 <- parameters[12]
	k8_5 <- parameters[13]
	k7_6 <- parameters[14]
	k6_7 <- parameters[15]
	k6_2 <- parameters[16]
	k5_8 <- parameters[17]
	k3p_7off <- parameters[18]
	k4p_8off <- parameters[19]
	k3p_3off <- parameters[20]
	k4p_4off <- parameters[21]
	k3p_7on <- parameters[22]
	k4p_8on <- parameters[23]
	k3p_3on <- parameters[24]
	k4p_4on <- parameters[25]
	kf_C_AKAR4 <- parameters[26]
	kb_C_AKAR4 <- parameters[27]
	kcat_AKARp <- parameters[28]
	kmOFF <- parameters[29]
	kmON <- parameters[30]
	KD_T <- parameters[31]
	b_AKAP <- parameters[32]
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
	netFlux <- numeric(16)
	k3p_7 <- b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off
	k4p_4 <- b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off
	k4p_8 <- b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off
	k3p_3 <- b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off
	k4_4p <- b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF
	k3_3p <- b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF
	k2_1 <- k1_2 * KD_T
	netFlux[1] <- k5_1*Rii_C
	netFlux[2] <- k4_1*RiiP*C - k1_4*RiiP_C
	netFlux[3] <- k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP
	netFlux[4] <- k4_3*cAMP*RiiP - k3_4*RiiP_cAMP
	netFlux[5] <- k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP
	netFlux[6] <- k8_7*cAMP*Rii - k7_8*Rii_cAMP
	netFlux[7] <- k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP
	netFlux[8] <- k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP
	netFlux[9] <- k6_2*Rii_C_cAMP
	netFlux[10] <- k8_5*Rii*C - k5_8*Rii_C
	netFlux[11] <- k4_4p*RiiP*CaN - k4p_4*RiiP_CaN
	netFlux[12] <- k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN
	netFlux[13] <- k4p_8*RiiP_CaN
	netFlux[14] <- k3p_7*RiiP_cAMP_CaN
	netFlux[15] <- kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	netFlux[16] <- kcat_AKARp*AKAR4_C
	return(netFlux)
}

AKAP79_fwdflux <- function(t, state, parameters){
	k5_1 <- parameters[1]
	k1_2 <- parameters[2]
	k3_2 <- parameters[3]
	k2_3 <- parameters[4]
	k3_4 <- parameters[5]
	k4_3 <- parameters[6]
	k4_1 <- parameters[7]
	k1_4 <- parameters[8]
	k8_7 <- parameters[9]
	k7_8 <- parameters[10]
	k5_6 <- parameters[11]
	k6_5 <- parameters[12]
	k8_5 <- parameters[13]
	k7_6 <- parameters[14]
	k6_7 <- parameters[15]
	k6_2 <- parameters[16]
	k5_8 <- parameters[17]
	k3p_7off <- parameters[18]
	k4p_8off <- parameters[19]
	k3p_3off <- parameters[20]
	k4p_4off <- parameters[21]
	k3p_7on <- parameters[22]
	k4p_8on <- parameters[23]
	k3p_3on <- parameters[24]
	k4p_4on <- parameters[25]
	kf_C_AKAR4 <- parameters[26]
	kb_C_AKAR4 <- parameters[27]
	kcat_AKARp <- parameters[28]
	kmOFF <- parameters[29]
	kmON <- parameters[30]
	KD_T <- parameters[31]
	b_AKAP <- parameters[32]
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
	fFlux <- numeric(16)
	k3p_7 <- b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off
	k4p_4 <- b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off
	k4p_8 <- b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off
	k3p_3 <- b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off
	k4_4p <- b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF
	k3_3p <- b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF
	k2_1 <- k1_2 * KD_T
	# k5_1*Rii_C
	fFlux[1] <- k5_1*Rii_C
	# k4_1*RiiP*C - k1_4*RiiP_C
	fFlux[2] <- k4_1*RiiP*C  - 0
	# k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP
	fFlux[3] <- k1_2*RiiP_C*cAMP  - 0
	# k4_3*cAMP*RiiP - k3_4*RiiP_cAMP
	fFlux[4] <- k4_3*cAMP*RiiP  - 0
	# k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP
	fFlux[5] <- k3_2*RiiP_cAMP*C  - 0
	# k8_7*cAMP*Rii - k7_8*Rii_cAMP
	fFlux[6] <- k8_7*cAMP*Rii  - 0
	# k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP
	fFlux[7] <- k5_6*Rii_C*cAMP  - 0
	# k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP
	fFlux[8] <- k7_6*Rii_cAMP*C  - 0
	# k6_2*Rii_C_cAMP
	fFlux[9] <- k6_2*Rii_C_cAMP
	# k8_5*Rii*C - k5_8*Rii_C
	fFlux[10] <- k8_5*Rii*C  - 0
	# k4_4p*RiiP*CaN - k4p_4*RiiP_CaN
	fFlux[11] <- k4_4p*RiiP*CaN  - 0
	# k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN
	fFlux[12] <- k3_3p*CaN*RiiP_cAMP  - 0
	# k4p_8*RiiP_CaN
	fFlux[13] <- k4p_8*RiiP_CaN
	# k3p_7*RiiP_cAMP_CaN
	fFlux[14] <- k3p_7*RiiP_cAMP_CaN
	# kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	fFlux[15] <- kf_C_AKAR4*C*AKAR4  - 0
	# kcat_AKARp*AKAR4_C
	fFlux[16] <- kcat_AKARp*AKAR4_C
	return(fFlux)
}

AKAP79_bwdflux <- function(t, state, parameters){
	k5_1 <- parameters[1]
	k1_2 <- parameters[2]
	k3_2 <- parameters[3]
	k2_3 <- parameters[4]
	k3_4 <- parameters[5]
	k4_3 <- parameters[6]
	k4_1 <- parameters[7]
	k1_4 <- parameters[8]
	k8_7 <- parameters[9]
	k7_8 <- parameters[10]
	k5_6 <- parameters[11]
	k6_5 <- parameters[12]
	k8_5 <- parameters[13]
	k7_6 <- parameters[14]
	k6_7 <- parameters[15]
	k6_2 <- parameters[16]
	k5_8 <- parameters[17]
	k3p_7off <- parameters[18]
	k4p_8off <- parameters[19]
	k3p_3off <- parameters[20]
	k4p_4off <- parameters[21]
	k3p_7on <- parameters[22]
	k4p_8on <- parameters[23]
	k3p_3on <- parameters[24]
	k4p_4on <- parameters[25]
	kf_C_AKAR4 <- parameters[26]
	kb_C_AKAR4 <- parameters[27]
	kcat_AKARp <- parameters[28]
	kmOFF <- parameters[29]
	kmON <- parameters[30]
	KD_T <- parameters[31]
	b_AKAP <- parameters[32]
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
	bFlux <- numeric(16)
	k3p_7 <- b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off
	k4p_4 <- b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off
	k4p_8 <- b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off
	k3p_3 <- b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off
	k4_4p <- b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF
	k3_3p <- b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF
	k2_1 <- k1_2 * KD_T
	bFlux[1] <- 0# k5_1*Rii_C
	return(bFlux)
}

# ode Jacobian df(t,y;p)/dy
AKAP79_jac<-function(t, state, parameters)
{
	k5_1 <- parameters[1]
	k1_2 <- parameters[2]
	k3_2 <- parameters[3]
	k2_3 <- parameters[4]
	k3_4 <- parameters[5]
	k4_3 <- parameters[6]
	k4_1 <- parameters[7]
	k1_4 <- parameters[8]
	k8_7 <- parameters[9]
	k7_8 <- parameters[10]
	k5_6 <- parameters[11]
	k6_5 <- parameters[12]
	k8_5 <- parameters[13]
	k7_6 <- parameters[14]
	k6_7 <- parameters[15]
	k6_2 <- parameters[16]
	k5_8 <- parameters[17]
	k3p_7off <- parameters[18]
	k4p_8off <- parameters[19]
	k3p_3off <- parameters[20]
	k4p_4off <- parameters[21]
	k3p_7on <- parameters[22]
	k4p_8on <- parameters[23]
	k3p_3on <- parameters[24]
	k4p_4on <- parameters[25]
	kf_C_AKAR4 <- parameters[26]
	kb_C_AKAR4 <- parameters[27]
	kcat_AKARp <- parameters[28]
	kmOFF <- parameters[29]
	kmON <- parameters[30]
	KD_T <- parameters[31]
	b_AKAP <- parameters[32]
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
	k3p_7 <- b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off
	k4p_4 <- b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off
	k4p_8 <- b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off
	k3p_3 <- b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off
	k4_4p <- b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF
	k3_3p <- b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF
	k2_1 <- k1_2 * KD_T
	reaction_51 <- k5_1*Rii_C
	reaction_14 <- k4_1*RiiP*C - k1_4*RiiP_C
	reaction_12 <- k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP
	reaction_43 <- k4_3*cAMP*RiiP - k3_4*RiiP_cAMP
	reaction_23 <- k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP
	reaction_78 <- k8_7*cAMP*Rii - k7_8*Rii_cAMP
	reaction_56 <- k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP
	reaction_76 <- k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP
	reaction_62 <- k6_2*Rii_C_cAMP
	reaction_58 <- k8_5*Rii*C - k5_8*Rii_C
	reaction_44 <- k4_4p*RiiP*CaN - k4p_4*RiiP_CaN
	reaction_33 <- k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN
	reaction_48 <- k4p_8*RiiP_CaN
	reaction_37 <- k3p_7*RiiP_cAMP_CaN
	reaction_1 <- kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2 <- kcat_AKARp*AKAR4_C
	jac_ <- matrix(0.0,16,16)
# column 1
# column 2
# column 3
# column 4
# column 5
# column 6
# column 7
# column 8
# column 9
# column 10
# column 11
# column 12
# column 13
# column 14
# column 15
# column 16
	rownames(jac_) <- c("Rii", "cAMP", "RiiP", "Rii_C", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "CaN", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4", "AKAR4_C", "AKAR4p")
	colnames(jac_) <- c("Rii", "cAMP", "RiiP", "Rii_C", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "CaN", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4", "AKAR4_C", "AKAR4p")
	return(jac_)
}
# ode parameter Jacobian df(t,y;p)/dp
AKAP79_jacp<-function(t, state, parameters)
{
	k5_1 <- parameters[1]
	k1_2 <- parameters[2]
	k3_2 <- parameters[3]
	k2_3 <- parameters[4]
	k3_4 <- parameters[5]
	k4_3 <- parameters[6]
	k4_1 <- parameters[7]
	k1_4 <- parameters[8]
	k8_7 <- parameters[9]
	k7_8 <- parameters[10]
	k5_6 <- parameters[11]
	k6_5 <- parameters[12]
	k8_5 <- parameters[13]
	k7_6 <- parameters[14]
	k6_7 <- parameters[15]
	k6_2 <- parameters[16]
	k5_8 <- parameters[17]
	k3p_7off <- parameters[18]
	k4p_8off <- parameters[19]
	k3p_3off <- parameters[20]
	k4p_4off <- parameters[21]
	k3p_7on <- parameters[22]
	k4p_8on <- parameters[23]
	k3p_3on <- parameters[24]
	k4p_4on <- parameters[25]
	kf_C_AKAR4 <- parameters[26]
	kb_C_AKAR4 <- parameters[27]
	kcat_AKARp <- parameters[28]
	kmOFF <- parameters[29]
	kmON <- parameters[30]
	KD_T <- parameters[31]
	b_AKAP <- parameters[32]
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
	k3p_7<-b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off
	k4p_4<-b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off
	k4p_8<-b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off
	k3p_3<-b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off
	k4_4p<-b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF
	k3_3p<-b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF
	k2_1<-k1_2 * KD_T
	reaction_51<-k5_1*Rii_C
	reaction_14<-k4_1*RiiP*C - k1_4*RiiP_C
	reaction_12<-k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP
	reaction_43<-k4_3*cAMP*RiiP - k3_4*RiiP_cAMP
	reaction_23<-k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP
	reaction_78<-k8_7*cAMP*Rii - k7_8*Rii_cAMP
	reaction_56<-k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP
	reaction_76<-k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP
	reaction_62<-k6_2*Rii_C_cAMP
	reaction_58<-k8_5*Rii*C - k5_8*Rii_C
	reaction_44<-k4_4p*RiiP*CaN - k4p_4*RiiP_CaN
	reaction_33<-k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN
	reaction_48<-k4p_8*RiiP_CaN
	reaction_37<-k3p_7*RiiP_cAMP_CaN
	reaction_1<-kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2<-kcat_AKARp*AKAR4_C
	jacp_ <- matrix(0.0,16,32)
# column 1
# column 2
# column 3
# column 4
# column 5
# column 6
# column 7
# column 8
# column 9
# column 10
# column 11
# column 12
# column 13
# column 14
# column 15
# column 16
# column 17
# column 18
# column 19
# column 20
# column 21
# column 22
# column 23
# column 24
# column 25
# column 26
# column 27
# column 28
# column 29
# column 30
# column 31
# column 32
	rownames(jacp_) <- c("Rii", "cAMP", "RiiP", "Rii_C", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "CaN", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4", "AKAR4_C", "AKAR4p")
	colnames(jacp_) <- c("k5_1", "k1_2", "k3_2", "k2_3", "k3_4", "k4_3", "k4_1", "k1_4", "k8_7", "k7_8", "k5_6", "k6_5", "k8_5", "k7_6", "k6_7", "k6_2", "k5_8", "k3p_7off", "k4p_8off", "k3p_3off", "k4p_4off", "k3p_7on", "k4p_8on", "k3p_3on", "k4p_4on", "kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp", "kmOFF", "kmON", "KD_T", "b_AKAP")
	return(jacp_)
}
# ode Functions F(t,y;p)
AKAP79_func<-function(t, state, parameters)
{
	k5_1 <- parameters[1]
	k1_2 <- parameters[2]
	k3_2 <- parameters[3]
	k2_3 <- parameters[4]
	k3_4 <- parameters[5]
	k4_3 <- parameters[6]
	k4_1 <- parameters[7]
	k1_4 <- parameters[8]
	k8_7 <- parameters[9]
	k7_8 <- parameters[10]
	k5_6 <- parameters[11]
	k6_5 <- parameters[12]
	k8_5 <- parameters[13]
	k7_6 <- parameters[14]
	k6_7 <- parameters[15]
	k6_2 <- parameters[16]
	k5_8 <- parameters[17]
	k3p_7off <- parameters[18]
	k4p_8off <- parameters[19]
	k3p_3off <- parameters[20]
	k4p_4off <- parameters[21]
	k3p_7on <- parameters[22]
	k4p_8on <- parameters[23]
	k3p_3on <- parameters[24]
	k4p_4on <- parameters[25]
	kf_C_AKAR4 <- parameters[26]
	kb_C_AKAR4 <- parameters[27]
	kcat_AKARp <- parameters[28]
	kmOFF <- parameters[29]
	kmON <- parameters[30]
	KD_T <- parameters[31]
	b_AKAP <- parameters[32]
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
	k3p_7 <- b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off
	k4p_4 <- b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off
	k4p_8 <- b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off
	k3p_3 <- b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off
	k4_4p <- b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF
	k3_3p <- b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF
	k2_1 <- k1_2 * KD_T
	reaction_51 <- k5_1*Rii_C
	reaction_14 <- k4_1*RiiP*C - k1_4*RiiP_C
	reaction_12 <- k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP
	reaction_43 <- k4_3*cAMP*RiiP - k3_4*RiiP_cAMP
	reaction_23 <- k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP
	reaction_78 <- k8_7*cAMP*Rii - k7_8*Rii_cAMP
	reaction_56 <- k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP
	reaction_76 <- k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP
	reaction_62 <- k6_2*Rii_C_cAMP
	reaction_58 <- k8_5*Rii*C - k5_8*Rii_C
	reaction_44 <- k4_4p*RiiP*CaN - k4p_4*RiiP_CaN
	reaction_33 <- k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN
	reaction_48 <- k4p_8*RiiP_CaN
	reaction_37 <- k3p_7*RiiP_cAMP_CaN
	reaction_1 <- kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2 <- kcat_AKARp*AKAR4_C
	func_ <- vector(mode='numeric',len=1)
	func_[1] <- (AKAR4p*5)*71.67+100 # AKAR4pOUT
	names(func_) <- c("AKAR4pOUT")
	return(func_)
}
# output function Jacobian dF(t,y;p)/dp
AKAP79_funcJac<-function(t, state, parameters)
{
	k5_1 <- parameters[1]
	k1_2 <- parameters[2]
	k3_2 <- parameters[3]
	k2_3 <- parameters[4]
	k3_4 <- parameters[5]
	k4_3 <- parameters[6]
	k4_1 <- parameters[7]
	k1_4 <- parameters[8]
	k8_7 <- parameters[9]
	k7_8 <- parameters[10]
	k5_6 <- parameters[11]
	k6_5 <- parameters[12]
	k8_5 <- parameters[13]
	k7_6 <- parameters[14]
	k6_7 <- parameters[15]
	k6_2 <- parameters[16]
	k5_8 <- parameters[17]
	k3p_7off <- parameters[18]
	k4p_8off <- parameters[19]
	k3p_3off <- parameters[20]
	k4p_4off <- parameters[21]
	k3p_7on <- parameters[22]
	k4p_8on <- parameters[23]
	k3p_3on <- parameters[24]
	k4p_4on <- parameters[25]
	kf_C_AKAR4 <- parameters[26]
	kb_C_AKAR4 <- parameters[27]
	kcat_AKARp <- parameters[28]
	kmOFF <- parameters[29]
	kmON <- parameters[30]
	KD_T <- parameters[31]
	b_AKAP <- parameters[32]
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
	k3p_7<-b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off
	k4p_4<-b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off
	k4p_8<-b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off
	k3p_3<-b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off
	k4_4p<-b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF
	k3_3p<-b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF
	k2_1<-k1_2 * KD_T
	reaction_51<-k5_1*Rii_C
	reaction_14<-k4_1*RiiP*C - k1_4*RiiP_C
	reaction_12<-k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP
	reaction_43<-k4_3*cAMP*RiiP - k3_4*RiiP_cAMP
	reaction_23<-k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP
	reaction_78<-k8_7*cAMP*Rii - k7_8*Rii_cAMP
	reaction_56<-k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP
	reaction_76<-k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP
	reaction_62<-k6_2*Rii_C_cAMP
	reaction_58<-k8_5*Rii*C - k5_8*Rii_C
	reaction_44<-k4_4p*RiiP*CaN - k4p_4*RiiP_CaN
	reaction_33<-k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN
	reaction_48<-k4p_8*RiiP_CaN
	reaction_37<-k3p_7*RiiP_cAMP_CaN
	reaction_1<-kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2<-kcat_AKARp*AKAR4_C
	fjac_ <- matrix(0.0,1,16)
# column 1
# column 2
# column 3
# column 4
# column 5
# column 6
# column 7
# column 8
# column 9
# column 10
# column 11
# column 12
# column 13
# column 14
# column 15
# column 16
	rownames(fjac_) <- c("AKAR4pOUT")
	colnames(fjac_) <- c("Rii", "cAMP", "RiiP", "Rii_C", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "CaN", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4", "AKAR4_C", "AKAR4p")
	return(fjac_)
}
# output function parameter Jacobian dF(t,y;p)/dp
AKAP79_funcJacp<-function(t, state, parameters)
{
	k5_1 <- parameters[1]
	k1_2 <- parameters[2]
	k3_2 <- parameters[3]
	k2_3 <- parameters[4]
	k3_4 <- parameters[5]
	k4_3 <- parameters[6]
	k4_1 <- parameters[7]
	k1_4 <- parameters[8]
	k8_7 <- parameters[9]
	k7_8 <- parameters[10]
	k5_6 <- parameters[11]
	k6_5 <- parameters[12]
	k8_5 <- parameters[13]
	k7_6 <- parameters[14]
	k6_7 <- parameters[15]
	k6_2 <- parameters[16]
	k5_8 <- parameters[17]
	k3p_7off <- parameters[18]
	k4p_8off <- parameters[19]
	k3p_3off <- parameters[20]
	k4p_4off <- parameters[21]
	k3p_7on <- parameters[22]
	k4p_8on <- parameters[23]
	k3p_3on <- parameters[24]
	k4p_4on <- parameters[25]
	kf_C_AKAR4 <- parameters[26]
	kb_C_AKAR4 <- parameters[27]
	kcat_AKARp <- parameters[28]
	kmOFF <- parameters[29]
	kmON <- parameters[30]
	KD_T <- parameters[31]
	b_AKAP <- parameters[32]
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
	k3p_7<-b_AKAP * k3p_7on + (1 - b_AKAP) * k3p_7off
	k4p_4<-b_AKAP * k4p_4on  +  (1 - b_AKAP)* k4p_4off
	k4p_8<-b_AKAP * k4p_8on + (1 - b_AKAP) * k4p_8off
	k3p_3<-b_AKAP * k3p_3on  +  (1 - b_AKAP)* k3p_3off
	k4_4p<-b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF
	k3_3p<-b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF
	k2_1<-k1_2 * KD_T
	reaction_51<-k5_1*Rii_C
	reaction_14<-k4_1*RiiP*C - k1_4*RiiP_C
	reaction_12<-k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP
	reaction_43<-k4_3*cAMP*RiiP - k3_4*RiiP_cAMP
	reaction_23<-k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP
	reaction_78<-k8_7*cAMP*Rii - k7_8*Rii_cAMP
	reaction_56<-k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP
	reaction_76<-k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP
	reaction_62<-k6_2*Rii_C_cAMP
	reaction_58<-k8_5*Rii*C - k5_8*Rii_C
	reaction_44<-k4_4p*RiiP*CaN - k4p_4*RiiP_CaN
	reaction_33<-k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN
	reaction_48<-k4p_8*RiiP_CaN
	reaction_37<-k3p_7*RiiP_cAMP_CaN
	reaction_1<-kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2<-kcat_AKARp*AKAR4_C
	fjacp_ <- matrix(0.0,1,32)
# column 1
# column 2
# column 3
# column 4
# column 5
# column 6
# column 7
# column 8
# column 9
# column 10
# column 11
# column 12
# column 13
# column 14
# column 15
# column 16
# column 17
# column 18
# column 19
# column 20
# column 21
# column 22
# column 23
# column 24
# column 25
# column 26
# column 27
# column 28
# column 29
# column 30
# column 31
# column 32
	rownames(fjacp_) <- c("AKAR4pOUT")
	colnames(fjacp_) <- c("k5_1", "k1_2", "k3_2", "k2_3", "k3_4", "k4_3", "k4_1", "k1_4", "k8_7", "k7_8", "k5_6", "k6_5", "k8_5", "k7_6", "k6_7", "k6_2", "k5_8", "k3p_7off", "k4p_8off", "k3p_3off", "k4p_4off", "k3p_7on", "k4p_8on", "k3p_3on", "k4p_4on", "kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp", "kmOFF", "kmON", "KD_T", "b_AKAP")
	return(fjacp_)
}
# ode default parameters; can depend on constants, and time  of initialization
AKAP79_default<-function(t=0.0)
{
	parameters <- vector(mode='numeric',len=32)
	parameters[1] <- 46.5411365826226
	parameters[2] <- 1.24347737295428
	parameters[3] <- 0.00328439415820813
	parameters[4] <- 0.587853369688012
	parameters[5] <- 0.0730805556930531
	parameters[6] <- 0.00712225662720034
	parameters[7] <- 0.00519774644439109
	parameters[8] <- 0.00117992862984651
	parameters[9] <- 0.00708844744872235
	parameters[10] <- 0.0833086316093731
	parameters[11] <- 0.665243673375593
	parameters[12] <- 0.594408977425014
	parameters[13] <- 0.100291450155288
	parameters[14] <- 0.0979128290361247
	parameters[15] <- 0.0222138094412742
	parameters[16] <- 33.8372184382292
	parameters[17] <- 0.000200941903978572
	parameters[18] <- 0.294046396805502
	parameters[19] <- 0.294046396805502
	parameters[20] <- 14.9166596278814
	parameters[21] <- 14.9166596278814
	parameters[22] <- 0.43067428164868
	parameters[23] <- 0.43067428164868
	parameters[24] <- 1.72247710096425
	parameters[25] <- 1.72247710096425
	parameters[26] <- 0.0180933753586079
	parameters[27] <- 0.104602112392241
	parameters[28] <- 10.1811826795126
	parameters[29] <- 102.235674360709
	parameters[30] <- 0.986951983065571
	parameters[31] <- 0.667954721425184
	parameters[32] <- 0
	names(parameters) <- c("k5_1", "k1_2", "k3_2", "k2_3", "k3_4", "k4_3", "k4_1", "k1_4", "k8_7", "k7_8", "k5_6", "k6_5", "k8_5", "k7_6", "k6_7", "k6_2", "k5_8", "k3p_7off", "k4p_8off", "k3p_3off", "k4p_4off", "k3p_7on", "k4p_8on", "k3p_3on", "k4p_4on", "kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp", "kmOFF", "kmON", "KD_T", "b_AKAP")
	return(parameters);
}
# ode initial values
AKAP79_init<-function(t=0.0, parameters=NA)
{
	k5_1 <- parameters[1]
	k1_2 <- parameters[2]
	k3_2 <- parameters[3]
	k2_3 <- parameters[4]
	k3_4 <- parameters[5]
	k4_3 <- parameters[6]
	k4_1 <- parameters[7]
	k1_4 <- parameters[8]
	k8_7 <- parameters[9]
	k7_8 <- parameters[10]
	k5_6 <- parameters[11]
	k6_5 <- parameters[12]
	k8_5 <- parameters[13]
	k7_6 <- parameters[14]
	k6_7 <- parameters[15]
	k6_2 <- parameters[16]
	k5_8 <- parameters[17]
	k3p_7off <- parameters[18]
	k4p_8off <- parameters[19]
	k3p_3off <- parameters[20]
	k4p_4off <- parameters[21]
	k3p_7on <- parameters[22]
	k4p_8on <- parameters[23]
	k3p_3on <- parameters[24]
	k4p_4on <- parameters[25]
	kf_C_AKAR4 <- parameters[26]
	kb_C_AKAR4 <- parameters[27]
	kcat_AKARp <- parameters[28]
	kmOFF <- parameters[29]
	kmON <- parameters[30]
	KD_T <- parameters[31]
	b_AKAP <- parameters[32]
	# the initial value may depend on the parameters. 
	state<-vector(mode='numeric',len=16)
	state[1] <- 6.3
	state[2] <- 0
	state[3] <- 0
	state[4] <- 0.63
	state[5] <- 0
	state[6] <- 0
	state[7] <- 0
	state[8] <- 0
	state[9] <- 0
	state[10] <- 0
	state[11] <- 1.5
	state[12] <- 0
	state[13] <- 0
	state[14] <- 0.2
	state[15] <- 0
	state[16] <- 0
	names(state) <- c("Rii", "cAMP", "RiiP", "Rii_C", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "CaN", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4", "AKAR4_C", "AKAR4p")
	return(state)
}
model<-list(vf=AKAP79_vf, jac=AKAP79_jac, jacp=AKAP79_jacp, func=AKAP79_func, funcJac=AKAP79_funcJac, funcJacp=AKAP79_funcJacp, init=AKAP79_init, par=AKAP79_default, name="AKAP79")
