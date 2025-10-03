library(uqsa)
library(parallel)
library(rgsl)
library(SBtabVFGEN)
library(optimx)

source("post.R")

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
sim <- simcf(experiments,modelName,parMap=log10ParMap)
llf <- logLikelihoodFunc(experiments)

MINIMIZE_THIS <- function(parMCMC){
	attr(parMCMC,"simulations") <- sim(parMCMC); 
	return(-1.0 * (llf(parMCMC) + log(dprior(parMCMC))))
}

O <- optimx(parMCMC,MINIMIZE_THIS, method="Nelder-Mead",itnmax=1e4)

optimxResult <- function(res, initialValue){
	message(sprintf("           value: %g",res$value[1]))
	message(sprintf("      func evals: %g",res$fevals[1]))
	message(sprintf("      grad evals: %g",res$gevals[1]))
	message(sprintf("  num iterations: %g",res$niter[1]))
	message(sprintf("convergence code: %g",res$convcode[1]))
	v <- res[1,seq_along(initialValue)]
	d <- norm(v-initialValue,type="2")
	message(sprintf("      difference: %g (%g %%)",d,100*norm(d,type="2")/norm(initialValue,type="2")))
	return(v)
}

p <- optimxResult(O,parMCMC)
y <- sim(p)

pdf(file="optimx-llf.pdf",width=cm(40) %as% "inches", height=cm(30) %as% "inches")
par(mfrow=c(3,3))
for (i in seq_along(y)){
	tm <- experiments[[i]]$outputTime
	o <- experiments[[i]]$outputValues[,1]
	e <- experiments[[i]]$errorValues[,1]
	plot(tm,o,ylim=c(0,200),xlab="time", ylab="AKAR4p",main=rownames(O)[1])
	errorBars(tm,o,e)
	lines(tm,y[[i]]$func[1,,1],col="blue4",lwd=3)
}
dev.off()
