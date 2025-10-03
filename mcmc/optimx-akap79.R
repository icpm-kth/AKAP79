library(uqsa)
library(parallel)
library(rgsl)
library(SBtabVFGEN)
library(optimx)

modelFiles <- dir("..",pattern="[.]tsv$",full.names=TRUE)
sb <- SBtabVFGEN::sbtab_from_tsv(modelFiles,verbose=FALSE)
suppressMessages(
	modelName <- checkModel(
		comment(sb),
		sprintf("./%s.so",comment(sb))
	)
)
stopifnot(file.exists("experiments.RDS"))
experiments <- readRDS("experiments.RDS")
options(mc.cores = length(experiments))
stopifnot(all(startsWith(trimws(sb$Parameter[["!Scale"]]),"log10")))
parMCMC <- sb$Parameter[["!DefaultValue"]]  # this already is in log10-space
names(parMCMC) <- rownames(sb$Parameter)
Median <- sb$Parameter[["!Median"]]
stdv <- sb$Parameter[["!Std"]]
dprior <- dNormalPrior(mean=Median,sd=stdv)
sim <- simulator.c(experiments,modelName,parMap=log10ParMap)
llf <- logLikelihoodFunc(experiments)

MINIMIZE_THIS <- function(parMCMC){
	attr(parMCMC,"simulations") <- sim(parMCMC); 
	return(-1.0 * (llf(parMCMC)+log(dprior(parMCMC))))
}

O <- optimx(parMCMC,MINIMIZE_THIS, method="BFGS",itnmax=1e5)
p <- as.numeric(O[1,seq_along(parMCMC)])
y <- sim(p)

dev.new()
par(mfrow=c(3,6))
for (i in seq_along(y)){
	tm <- experiments[[i]]$outputTime
	o <- experiments[[i]]$outputValues[,1]
	e <- experiments[[i]]$errorValues[,1]
	plot(tm,o,ylim=c(0,200),xlab="time", ylab="AKAR4p",main=rownames(O)[1])
	arrows(tm,o,tm,o+e,length=0.01,angle=90)
	arrows(tm,o,tm,o-e,length=0.01,angle=90)
	lines(tm,y[[i]]$func[1,,1],col="blue4")
}

