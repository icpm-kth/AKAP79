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
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
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
	k3p_7 <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k4p_4 <- b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	k4p_8 <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k3p_3 <- b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
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
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
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
	k3p_7 <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k4p_4 <- b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	k4p_8 <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k3p_3 <- b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
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
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
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
	k3p_7 <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k4p_4 <- b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	k4p_8 <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k3p_3 <- b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
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
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
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
	k3p_7 <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k4p_4 <- b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	k4p_8 <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k3p_3 <- b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	k4_4p <- b_AKAP * ((k4p_4 + k3p_7)/kmON ) + (1 - b_AKAP) * (k4p_4 + k3p_7) / kmOFF
	k3_3p <- b_AKAP * ((k3p_3 + k4p_8)/kmON) + (1 - b_AKAP) * (k3p_3 + k4p_8)/kmOFF
	k2_1 <- k1_2 * KD_T
	bFlux[1] <- 0# k5_1*Rii_C
	bFlux[2] <-  k1_4*RiiP_C# k4_1*RiiP*C - k1_4*RiiP_C
	bFlux[3] <-  k2_1*RiiP_C_cAMP# k1_2*RiiP_C*cAMP - k2_1*RiiP_C_cAMP
	bFlux[4] <-  k3_4*RiiP_cAMP# k4_3*cAMP*RiiP - k3_4*RiiP_cAMP
	bFlux[5] <-  k2_3*RiiP_C_cAMP# k3_2*RiiP_cAMP*C - k2_3*RiiP_C_cAMP
	bFlux[6] <-  k7_8*Rii_cAMP# k8_7*cAMP*Rii - k7_8*Rii_cAMP
	bFlux[7] <-  k6_5*Rii_C_cAMP# k5_6*Rii_C*cAMP - k6_5*Rii_C_cAMP
	bFlux[8] <-  k6_7*Rii_C_cAMP# k7_6*Rii_cAMP*C - k6_7*Rii_C_cAMP
	bFlux[9] <- 0# k6_2*Rii_C_cAMP
	bFlux[10] <-  k5_8*Rii_C# k8_5*Rii*C - k5_8*Rii_C
	bFlux[11] <-  k4p_4*RiiP_CaN# k4_4p*RiiP*CaN - k4p_4*RiiP_CaN
	bFlux[12] <-  k3p_3*RiiP_cAMP_CaN# k3_3p*CaN*RiiP_cAMP - k3p_3*RiiP_cAMP_CaN
	bFlux[13] <- 0# k4p_8*RiiP_CaN
	bFlux[14] <- 0# k3p_7*RiiP_cAMP_CaN
	bFlux[15] <-  kb_C_AKAR4*AKAR4_C# kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	bFlux[16] <- 0# kcat_AKARp*AKAR4_C
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
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
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
	k3p_7 <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k4p_4 <- b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	k4p_8 <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k3p_3 <- b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
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
	jac_[1,1] <- (-cAMP*k8_7)-C*k8_5
	jac_[2,1] <- -cAMP*k8_7
	jac_[4,1] <- C*k8_5
	jac_[8,1] <- -C*k8_5
	jac_[9,1] <- cAMP*k8_7
# column 2
	jac_[1,2] <- -Rii*k8_7
	jac_[2,2] <- (-Rii*k8_7)-Rii_C*k5_6-RiiP*k4_3-RiiP_C*k1_2
	jac_[3,2] <- -RiiP*k4_3
	jac_[4,2] <- -Rii_C*k5_6
	jac_[5,2] <- RiiP*k4_3
	jac_[6,2] <- -RiiP_C*k1_2
	jac_[7,2] <- RiiP_C*k1_2
	jac_[9,2] <- Rii*k8_7
	jac_[10,2] <- Rii_C*k5_6
# column 3
	jac_[2,3] <- -cAMP*k4_3
	jac_[3,3] <- (-CaN*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF))-cAMP*k4_3-C*k4_1
	jac_[5,3] <- cAMP*k4_3
	jac_[6,3] <- C*k4_1
	jac_[8,3] <- -C*k4_1
	jac_[11,3] <- -CaN*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
	jac_[12,3] <- CaN*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
# column 4
	jac_[1,4] <- k5_8
	jac_[2,4] <- -cAMP*k5_6
	jac_[4,4] <- (-k5_8)-cAMP*k5_6-k5_1
	jac_[6,4] <- k5_1
	jac_[8,4] <- k5_8
	jac_[10,4] <- cAMP*k5_6
# column 5
	jac_[2,5] <- k3_4
	jac_[3,5] <- k3_4
	jac_[5,5] <- (-CaN*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF))-k3_4-C*k3_2
	jac_[7,5] <- C*k3_2
	jac_[8,5] <- -C*k3_2
	jac_[11,5] <- -CaN*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
	jac_[13,5] <- CaN*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
# column 6
	jac_[2,6] <- -cAMP*k1_2
	jac_[3,6] <- k1_4
	jac_[6,6] <- (-k1_4)-cAMP*k1_2
	jac_[7,6] <- cAMP*k1_2
	jac_[8,6] <- k1_4
# column 7
	jac_[2,7] <- KD_T*k1_2
	jac_[5,7] <- k2_3
	jac_[6,7] <- KD_T*k1_2
	jac_[7,7] <- (-k2_3)-KD_T*k1_2
	jac_[8,7] <- k2_3
# column 8
	jac_[1,8] <- -Rii*k8_5
	jac_[3,8] <- -RiiP*k4_1
	jac_[4,8] <- Rii*k8_5
	jac_[5,8] <- -RiiP_cAMP*k3_2
	jac_[6,8] <- RiiP*k4_1
	jac_[7,8] <- RiiP_cAMP*k3_2
	jac_[8,8] <- (-AKAR4*kf_C_AKAR4)-Rii*k8_5-Rii_cAMP*k7_6-RiiP*k4_1-RiiP_cAMP*k3_2
	jac_[9,8] <- -Rii_cAMP*k7_6
	jac_[10,8] <- Rii_cAMP*k7_6
	jac_[14,8] <- -AKAR4*kf_C_AKAR4
	jac_[15,8] <- AKAR4*kf_C_AKAR4
# column 9
	jac_[1,9] <- k7_8
	jac_[2,9] <- k7_8
	jac_[8,9] <- -C*k7_6
	jac_[9,9] <- (-k7_8)-C*k7_6
	jac_[10,9] <- C*k7_6
# column 10
	jac_[2,10] <- k6_5
	jac_[4,10] <- k6_5
	jac_[7,10] <- k6_2
	jac_[8,10] <- k6_7
	jac_[9,10] <- k6_7
	jac_[10,10] <- (-k6_7)-k6_5-k6_2
# column 11
	jac_[3,11] <- -RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
	jac_[5,11] <- -RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
	jac_[11,11] <- (-RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF))-RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
	jac_[12,11] <- RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
	jac_[13,11] <- RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
# column 12
	jac_[1,12] <- AKAPon_1*b_AKAP+AKAPoff_1*(1-b_AKAP)
	jac_[3,12] <- AKAPon_3*b_AKAP+AKAPoff_3*(1-b_AKAP)
	jac_[11,12] <- AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)
	jac_[12,12] <- (-AKAPon_3*b_AKAP)-AKAPon_1*b_AKAP-AKAPoff_3*(1-b_AKAP)-AKAPoff_1*(1-b_AKAP)
# column 13
	jac_[5,13] <- AKAPon_3*b_AKAP+AKAPoff_3*(1-b_AKAP)
	jac_[9,13] <- AKAPon_1*b_AKAP+AKAPoff_1*(1-b_AKAP)
	jac_[11,13] <- AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)
	jac_[13,13] <- (-AKAPon_3*b_AKAP)-AKAPon_1*b_AKAP-AKAPoff_3*(1-b_AKAP)-AKAPoff_1*(1-b_AKAP)
# column 14
	jac_[8,14] <- -C*kf_C_AKAR4
	jac_[14,14] <- -C*kf_C_AKAR4
	jac_[15,14] <- C*kf_C_AKAR4
# column 15
	jac_[8,15] <- kcat_AKARp+kb_C_AKAR4
	jac_[14,15] <- kb_C_AKAR4
	jac_[15,15] <- (-kcat_AKARp)-kb_C_AKAR4
	jac_[16,15] <- kcat_AKARp
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
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
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
	k3p_7<-b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k4p_4<-b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	k4p_8<-b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k3p_3<-b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
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
	jacp_ <- matrix(0.0,16,28)
# column 1
	jacp_[4,1] <- -Rii_C
	jacp_[6,1] <- Rii_C
# column 2
	jacp_[2,2] <- KD_T*RiiP_C_cAMP-RiiP_C*cAMP
	jacp_[6,2] <- KD_T*RiiP_C_cAMP-RiiP_C*cAMP
	jacp_[7,2] <- RiiP_C*cAMP-KD_T*RiiP_C_cAMP
# column 3
	jacp_[5,3] <- -C*RiiP_cAMP
	jacp_[7,3] <- C*RiiP_cAMP
	jacp_[8,3] <- -C*RiiP_cAMP
# column 4
	jacp_[5,4] <- RiiP_C_cAMP
	jacp_[7,4] <- -RiiP_C_cAMP
	jacp_[8,4] <- RiiP_C_cAMP
# column 5
	jacp_[2,5] <- RiiP_cAMP
	jacp_[3,5] <- RiiP_cAMP
	jacp_[5,5] <- -RiiP_cAMP
# column 6
	jacp_[2,6] <- -RiiP*cAMP
	jacp_[3,6] <- -RiiP*cAMP
	jacp_[5,6] <- RiiP*cAMP
# column 7
	jacp_[3,7] <- -C*RiiP
	jacp_[6,7] <- C*RiiP
	jacp_[8,7] <- -C*RiiP
# column 8
	jacp_[3,8] <- RiiP_C
	jacp_[6,8] <- -RiiP_C
	jacp_[8,8] <- RiiP_C
# column 9
	jacp_[1,9] <- -Rii*cAMP
	jacp_[2,9] <- -Rii*cAMP
	jacp_[9,9] <- Rii*cAMP
# column 10
	jacp_[1,10] <- Rii_cAMP
	jacp_[2,10] <- Rii_cAMP
	jacp_[9,10] <- -Rii_cAMP
# column 11
	jacp_[2,11] <- -Rii_C*cAMP
	jacp_[4,11] <- -Rii_C*cAMP
	jacp_[10,11] <- Rii_C*cAMP
# column 12
	jacp_[2,12] <- Rii_C_cAMP
	jacp_[4,12] <- Rii_C_cAMP
	jacp_[10,12] <- -Rii_C_cAMP
# column 13
	jacp_[1,13] <- -C*Rii
	jacp_[4,13] <- C*Rii
	jacp_[8,13] <- -C*Rii
# column 14
	jacp_[8,14] <- -C*Rii_cAMP
	jacp_[9,14] <- -C*Rii_cAMP
	jacp_[10,14] <- C*Rii_cAMP
# column 15
	jacp_[8,15] <- Rii_C_cAMP
	jacp_[9,15] <- Rii_C_cAMP
	jacp_[10,15] <- -Rii_C_cAMP
# column 16
	jacp_[7,16] <- Rii_C_cAMP
	jacp_[10,16] <- -Rii_C_cAMP
# column 17
	jacp_[1,17] <- Rii_C
	jacp_[4,17] <- -Rii_C
	jacp_[8,17] <- Rii_C
# column 18
	jacp_[1,18] <- RiiP_CaN*(1-b_AKAP)
	jacp_[3,18] <- -CaN*RiiP*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)
	jacp_[5,18] <- -CaN*RiiP_cAMP*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)
	jacp_[9,18] <- RiiP_cAMP_CaN*(1-b_AKAP)
	jacp_[11,18] <- (-CaN*RiiP_cAMP*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF))-CaN*RiiP*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)+RiiP_cAMP_CaN*(1-b_AKAP)+RiiP_CaN*(1-b_AKAP)
	jacp_[12,18] <- CaN*RiiP*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)-RiiP_CaN*(1-b_AKAP)
	jacp_[13,18] <- CaN*RiiP_cAMP*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)-RiiP_cAMP_CaN*(1-b_AKAP)
# column 19
	jacp_[3,19] <- RiiP_CaN*(1-b_AKAP)-CaN*RiiP*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)
	jacp_[5,19] <- RiiP_cAMP_CaN*(1-b_AKAP)-CaN*RiiP_cAMP*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)
	jacp_[11,19] <- (-CaN*RiiP_cAMP*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF))-CaN*RiiP*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)+RiiP_cAMP_CaN*(1-b_AKAP)+RiiP_CaN*(1-b_AKAP)
	jacp_[12,19] <- CaN*RiiP*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)-RiiP_CaN*(1-b_AKAP)
	jacp_[13,19] <- CaN*RiiP_cAMP*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)-RiiP_cAMP_CaN*(1-b_AKAP)
# column 20
	jacp_[1,20] <- RiiP_CaN*b_AKAP
	jacp_[3,20] <- -CaN*RiiP*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)
	jacp_[5,20] <- -CaN*RiiP_cAMP*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)
	jacp_[9,20] <- RiiP_cAMP_CaN*b_AKAP
	jacp_[11,20] <- (-CaN*RiiP_cAMP*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF))-CaN*RiiP*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)+RiiP_cAMP_CaN*b_AKAP+RiiP_CaN*b_AKAP
	jacp_[12,20] <- CaN*RiiP*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)-RiiP_CaN*b_AKAP
	jacp_[13,20] <- CaN*RiiP_cAMP*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)-RiiP_cAMP_CaN*b_AKAP
# column 21
	jacp_[3,21] <- RiiP_CaN*b_AKAP-CaN*RiiP*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)
	jacp_[5,21] <- RiiP_cAMP_CaN*b_AKAP-CaN*RiiP_cAMP*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)
	jacp_[11,21] <- (-CaN*RiiP_cAMP*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF))-CaN*RiiP*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)+RiiP_cAMP_CaN*b_AKAP+RiiP_CaN*b_AKAP
	jacp_[12,21] <- CaN*RiiP*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)-RiiP_CaN*b_AKAP
	jacp_[13,21] <- CaN*RiiP_cAMP*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)-RiiP_cAMP_CaN*b_AKAP
# column 22
	jacp_[8,22] <- -AKAR4*C
	jacp_[14,22] <- -AKAR4*C
	jacp_[15,22] <- AKAR4*C
# column 23
	jacp_[8,23] <- AKAR4_C
	jacp_[14,23] <- AKAR4_C
	jacp_[15,23] <- -AKAR4_C
# column 24
	jacp_[8,24] <- AKAR4_C
	jacp_[15,24] <- -AKAR4_C
	jacp_[16,24] <- AKAR4_C
# column 25
	jacp_[3,25] <- (CaN*RiiP*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF^2
	jacp_[5,25] <- (CaN*RiiP_cAMP*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF^2
	jacp_[11,25] <- (CaN*RiiP_cAMP*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF^2+(CaN*RiiP*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF^2
	jacp_[12,25] <- -(CaN*RiiP*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF^2
	jacp_[13,25] <- -(CaN*RiiP_cAMP*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF^2
# column 26
	jacp_[3,26] <- (CaN*RiiP*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON^2
	jacp_[5,26] <- (CaN*RiiP_cAMP*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON^2
	jacp_[11,26] <- (CaN*RiiP_cAMP*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON^2+(CaN*RiiP*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON^2
	jacp_[12,26] <- -(CaN*RiiP*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON^2
	jacp_[13,26] <- -(CaN*RiiP_cAMP*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON^2
# column 27
	jacp_[2,27] <- RiiP_C_cAMP*k1_2
	jacp_[6,27] <- RiiP_C_cAMP*k1_2
	jacp_[7,27] <- -RiiP_C_cAMP*k1_2
# column 28
	jacp_[1,28] <- (AKAPon_1-AKAPoff_1)*RiiP_CaN
	jacp_[3,28] <- (AKAPon_3-AKAPoff_3)*RiiP_CaN-CaN*RiiP*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF)
	jacp_[5,28] <- (AKAPon_3-AKAPoff_3)*RiiP_cAMP_CaN-CaN*RiiP_cAMP*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF)
	jacp_[9,28] <- (AKAPon_1-AKAPoff_1)*RiiP_cAMP_CaN
	jacp_[11,28] <- (-CaN*RiiP_cAMP*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF))-CaN*RiiP*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF)+(AKAPon_3-AKAPoff_3)*RiiP_cAMP_CaN+(AKAPon_1-AKAPoff_1)*RiiP_cAMP_CaN+(AKAPon_3-AKAPoff_3)*RiiP_CaN+(AKAPon_1-AKAPoff_1)*RiiP_CaN
	jacp_[12,28] <- CaN*RiiP*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF)-(AKAPon_3-AKAPoff_3)*RiiP_CaN-(AKAPon_1-AKAPoff_1)*RiiP_CaN
	jacp_[13,28] <- CaN*RiiP_cAMP*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF)-(AKAPon_3-AKAPoff_3)*RiiP_cAMP_CaN-(AKAPon_1-AKAPoff_1)*RiiP_cAMP_CaN
	rownames(jacp_) <- c("Rii", "cAMP", "RiiP", "Rii_C", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "CaN", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4", "AKAR4_C", "AKAR4p")
	colnames(jacp_) <- c("k5_1", "k1_2", "k3_2", "k2_3", "k3_4", "k4_3", "k4_1", "k1_4", "k8_7", "k7_8", "k5_6", "k6_5", "k8_5", "k7_6", "k6_7", "k6_2", "k5_8", "AKAPoff_1", "AKAPoff_3", "AKAPon_1", "AKAPon_3", "kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp", "kmOFF", "kmON", "KD_T", "b_AKAP")
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
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
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
	k3p_7 <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k4p_4 <- b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	k4p_8 <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k3p_3 <- b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
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
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
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
	k3p_7<-b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k4p_4<-b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	k4p_8<-b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k3p_3<-b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
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
	fjac_[1,16] <- 358.35
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
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
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
	k3p_7<-b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k4p_4<-b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	k4p_8<-b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	k3p_3<-b_AKAP * AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
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
	fjacp_ <- matrix(0.0,1,28)
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
	rownames(fjacp_) <- c("AKAR4pOUT")
	colnames(fjacp_) <- c("k5_1", "k1_2", "k3_2", "k2_3", "k3_4", "k4_3", "k4_1", "k1_4", "k8_7", "k7_8", "k5_6", "k6_5", "k8_5", "k7_6", "k6_7", "k6_2", "k5_8", "AKAPoff_1", "AKAPoff_3", "AKAPon_1", "AKAPon_3", "kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp", "kmOFF", "kmON", "KD_T", "b_AKAP")
	return(fjacp_)
}
# ode default parameters; can depend on constants, and time  of initialization
AKAP79_default<-function(t=0.0)
{
	parameters <- vector(mode='numeric',len=28)
	parameters[1] <- 33.0369541036815
	parameters[2] <- 0.494310686986835
	parameters[3] <- 0.00544502652842421
	parameters[4] <- 0.0155955250282695
	parameters[5] <- 0.00159955802861467
	parameters[6] <- 0.0149968483550237
	parameters[7] <- 0.0380189396320561
	parameters[8] <- 0.00260015956316527
	parameters[9] <- 0.0149968483550237
	parameters[10] <- 0.00159955802861467
	parameters[11] <- 0.49545019080479
	parameters[12] <- 1.41253754462275
	parameters[13] <- 2.09893988362352
	parameters[14] <- 0.298538261891796
	parameters[15] <- 0.0179887091512879
	parameters[16] <- 33.0369541036815
	parameters[17] <- 0.000299916251898765
	parameters[18] <- 2.60015956316527
	parameters[19] <- 19.9986186963274
	parameters[20] <- 0.449779854893288
	parameters[21] <- 1.99986186963274
	parameters[22] <- 0.0179887091512879
	parameters[23] <- 0.105925372517729
	parameters[24] <- 10.2093948370768
	parameters[25] <- 100
	parameters[26] <- 1
	parameters[27] <- 0.699841996002273
	parameters[28] <- 0
	names(parameters) <- c("k5_1", "k1_2", "k3_2", "k2_3", "k3_4", "k4_3", "k4_1", "k1_4", "k8_7", "k7_8", "k5_6", "k6_5", "k8_5", "k7_6", "k6_7", "k6_2", "k5_8", "AKAPoff_1", "AKAPoff_3", "AKAPon_1", "AKAPon_3", "kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp", "kmOFF", "kmON", "KD_T", "b_AKAP")
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
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
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
