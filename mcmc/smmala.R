#!/usr/bin/env Rscript
## ----setup--------------------------------------------------------------------
library(uqsa)
library(parallel)
library(rgsl)
library(SBtabVFGEN)
library(pbdMPI)

Branch <- system2("git",c("branch","--show-current"),stdout=TRUE)
start_time <- Sys.time()                         # measure sampling time
comm  <- 0
r <- pbdMPI::comm.rank(comm=comm)
cs <- pbdMPI::comm.size(comm=comm)
pbdMPI::init()
attr(comm,"rank") <- r
attr(comm,"size") <- cs

scale <- \(x) as.integer(round(log10(x))+1)

# random number seed:

set.seed(137*r + 1337*cs)
R <- round(cs/2)

N <- 10000
h <- 1e-2
cycl <- 7
PREFIX <- "main"
## possibly adjust sampling-settings
a <- commandArgs(trailingOnly=TRUE)
if (!is.null(a) && length(a)>0 && length(a) %% 2 == 0){
	for (i in seq(1,length(a),by=2)){
		key = a[i]
		val = a[i+1]
		switch(key,
			"-N"={N <- as.integer(val)},
			"-h"={h <- as.double(val)},
			"-c"={cycl <- as.integer(val)},
			"-R"={R <- as.integer(val)},
			"--prefix"={PREFIX=as.character(val)},
			{warning(sprintf("unknown option %s %s",key,val))}
		)
	}
}
beta <- (1.0 - ((r %% R)/R))^2
scale <- \(x) return(round(log10(x))+1)

Hash <- digest::sha1(x=c(N,R,cycl,r,cs,beta))
shortHash <- substr(Hash,1,8)

## ----label="SBtab content"----------------------------------------------------
modelFiles <- dir("..",pattern="[.]tsv$",full.names=TRUE)
sb <- SBtabVFGEN::sbtab_from_tsv(modelFiles,verbose=FALSE)

suppressMessages(
	modelName <- checkModel(
		comment(sb),
		sprintf("./%s.so",comment(sb))
	)
)

experiments <- sbtab.data(sb)

for (i in seq_along(experiments)){
	t_ <- experiments[[i]]$outputTimes
	nt <- length(t_)
	D <- t(experiments[[i]]$outputValues)
	D[is.na(D)] <- 0.0
	SD <- t(experiments[[i]]$errorValues) ## or: D*0.05+apply(D,1,FUN=max,na.rm=TRUE)*0.05
	SD[is.na(SD)] <- Inf
	experiments[[i]]$time <- t_
	experiments[[i]]$data <- D
	experiments[[i]]$stdv <- SD
}

options(mc.cores = length(experiments))
## ----default------------------------------------------------------------------
n <- length(experiments[[1]]$input)

stopifnot(all(startsWith(trimws(sb$Parameter[["!Scale"]]),"log10")))
Median <- sb$Parameter[["!Median"]]
parMCMC <- sb$Parameter[["!DefaultValue"]]  # this already is in log10-space
stdv <- sb$Parameter[["!Std"]]              # this as well
if (is.null(stdv) || any(!is.numeric(stdv)) || any(is.na(stdv))) {
	warning("no standard error («!Std» field) in 'SBtab$Parameter'")
	stdv <- parMCMC*0.5 + 0.5 + 0.5*max(parMCMC)
}
dprior <- dNormalPrior(mean=Median,sd=stdv)
rprior <- rNormalPrior(mean=Median,sd=stdv)
gprior <- gradLog_NormalPrior(mean=Median,sd=stdv)

## log-likelihood
llf <- function(parMCMC){
	return(sum(sapply(attr(parMCMC,"simulations"),\(y) y$logLikelihood)))
}

## gradient of the log-likelihood
gllf <- function(parMCMC){
	i <- seq_along(parMCMC)
	J <- log10ParMapJac(parMCMC)
	return(as.numeric(rowSums(sapply(attr(parMCMC,"simulations"),\(y) y$gradLogLikelihood))[i] %*% J))
}

## Fisher Information
fi <- function(parMCMC){
	i <- seq_along(parMCMC)
	J <- log10ParMapJac(parMCMC)
	fi <- rowSums(sapply(attr(parMCMC,"simulations"),\(y) y$FisherInformation[i,i,1],simplify="array"),dims=2) # sums over 3rd dimension
	return(t(J) %*% fi %*% J)
}

fiPrior <- solve(diag(stdv, length(parMCMC)))

## ----simulate-----------------------------------------------------------------
sim <- simfi(experiments,modelName,log10ParMap)

X <- NULL # this is to suppress one intial warning
sampleSize <- round(exp(seq(log(1000),log(N),length.out=cycl)))

## ----sample-------------------------------------------------------------------
smmala <- mcmcUpdate(
	simulate=sim,
	experiments=experiments,
	logLikelihood=llf,
	dprior=dprior,
	gradLogLikelihood=gllf,
	gprior=gprior,
	fisherInformation=fi,
	fisherInformationPrior=fiPrior
)
serialSMMALA <- mcmc(smmala)
ptSMMALA <- mcmc_mpi(smmala,comm=comm,swapFunc=pbdMPI_bcast_reduce_temperatures)

for (i in seq(cycl)){
	rm(X)
	x <- mcmcInit(
		beta,
		parMCMC,
		simulate=sim,
		logLikelihood=llf,
		dprior,
		gllf,
		gprior,
		fi
	)
	if (i < round(cycl/2)){
		for (k in seq(5)){
			s <- serialSMMALA(x,100,h)
			ar <- attr(s,"acceptanceRate")
			ml <- max(attr(s,"logLikelihood"))
			message(
				sprintf("[smmala-r%0*ic%0*ik%i] h=%6g, a=%.2f, β=%7g, ml=%g",
					as.integer(round(log10(cs)))+1,as.integer(r),
					as.integer(round(log10(cycl)))+1,as.integer(i),
					k,
					h,
					ar,
					beta,
					ml
				)
			)
			h <- h * (0.5 + ar^4/(0.5^4 + ar^4)) # 0.5 is the target acceptance value
			x <- attr(s,"lastPoint")
		}
	}
	pbdMPI::barrier()
	## This is where the main amount of work is done:
	s <- ptSMMALA(x,sampleSize[i],0.25*h)
	pbdMPI::barrier()
	commonName <- sprintf("%s-%s-branch-Sample",
		PREFIX,
		Branch
	)
	sid <- round(i*cs+r) # a sequential id number that increases with iteration i and rank
	maxsid <- round(cycl*cs+r)
	sampleFile=sprintf("%s-%0*x%s-cycle-%0*i-rank-%0*i-of-%0*i.RDS",
		commonName,
		scale(maxsid),sid,
		shortHash,
		scale(cycl),i,
		scale(cs),r,
		scale(cs),cs)
	saveRDS(s,file=sampleFile)
	saveRDS(s,file=sampleFile)
	## free up memory for next sample
	ar <- as.integer(round(100*attr(s,"acceptanceRate")))
	sr <- as.integer(round(100*attr(s,"swapRate")))
	rm(s)
	gc()
	pbdMPI::barrier()
	f <- dir(pattern=sprintf("^%s-.*-rank-[0-9]+-of-%0*i[.]RDS$",commonName,scale(cs),cs))
	X <- uqsa::gatherSample(f,beta)
	effar <- as.integer(round(100*mean(diff(sort(attr(X,"logLikelihood")))>1e-14)))
	attr(X,"beta") <- beta
	parMCMC <- as.numeric(tail(X,1))
	cat(
		sprintf("rank %*i/%i finished %*i iterations on cycle %*i with acc. rate of %3i %% and swap rate of %3i %% (log10(step-size) is %5.4f) and acc. rate %3i %% for β=%g.\n",
			scale(cs),r,
			cs,
			scale(tail(sampleSize,1)),sampleSize[i],
			scale(cycl),i,
			ar,
			sr,
			log10(h),
			effar,
			beta
		)
	)
}
saveRDS(X,file=sprintf("%s-final-sample-%i-%i-lb100-%i.RDS",PREFIX,r,cs,round(log(beta)*100)))
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)
finalize()
print(warnings())
