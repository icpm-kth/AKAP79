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

# random number seed:

set.seed(137*r + 1337*cs)
R <- round(cs/2)

N <- 10000
h <- 2e-2
cycl <- 7
PREFIX <- system2("date",args=c("--iso-8601"),stdout=TRUE)
## possibly adjust sampling-settings
a <- commandArgs(trailingOnly=TRUE)
if (!is.null(a) && length(a)>0){
	if (length(a)%%2 != 0) {
		warning("an option (key-value-pairs) is missing the value.")
	}
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
options(mc.cores = length(experiments))
## ----default------------------------------------------------------------------
n <- length(experiments[[1]]$input)

stopifnot(all(startsWith(trimws(sb$Parameter[["!Scale"]]),"log10")))
Median <- sb$Parameter[["!Median"]]
parMCMC <- sb$Parameter[["!DefaultValue"]]  # this already is in log10-space
names(parMCMC) <- rownames(sb$Parameter)
stdv <- sb$Parameter[["!Std"]]              # this as well
if (is.null(stdv) || any(!is.numeric(stdv)) || any(is.na(stdv))) {
	warning("no standard error («!Std» field) in 'SBtab$Parameter'")
	stdv <- parMCMC*0.5 + 0.5 + 0.5*max(parMCMC)
}
dprior <- dNormalPrior(mean=Median,sd=stdv)
rprior <- rNormalPrior(mean=Median,sd=stdv)
## ----simulate-----------------------------------------------------------------
sim <- simcf(experiments,modelName,log10ParMap)

llf <- logLikelihoodFunc(experiments)

X <- NULL # this is to suppress one intial warning
sampleSize <- round(exp(seq(log(100),log(N),length.out=cycl)))

## save the chosen settings for this run
save(
	R,dprior,rprior,Median,parMCMC,stdv,experiments,sim,llf,sb,sampleSize,Branch,a,
	file=sprintf("%s-%s-branch-%s-settings.Rdata",
		PREFIX,
		Branch,
		shortHash)
)

## ----sample-------------------------------------------------------------------
metropolis <- mcmcUpdate(
	simulate=sim,
	experiments=experiments,
	logLikelihood=llf,
	dprior=dprior
)
mhmcmc <- mcmc(metropolis) # sequential version for step-size tuning
ptMetropolis <- mcmc_mpi(
	metropolis,
	comm=comm,
	swapDelay=0,
	swapFunc=pbdMPI_bcast_reduce_temperatures
)
for (i in seq(cycl)){
	rm(X)
	x <- mcmcInit(
		beta,
		parMCMC,
		simulate=sim,
		logLikelihood=llf,
		dprior
	)
	names(x) <- rownames(sb$Parameter)
	if (i < round(cycl/2)){
		for (k in seq(5)){
			s <- mhmcmc(x,100,h)
			ar <- attr(s,"acceptanceRate")
			ml <- max(attr(s,"logLikelihood"))
			message(
				sprintf("[metropolis-r%0*ic%0*ik%i] h=%6g, a=%.2f, β=%7g, ml=%g",
					as.integer(round(log10(cs)))+1,as.integer(r),
					as.integer(round(log10(cycl)))+1,as.integer(i),
					k,
					h,
					ar,
					beta,
					ml
				)
			)
			h <- h * (0.5 + ar^4/(0.25^4 + ar^4)) # 0.25 is the target acceptance value
		}
	}
	## This is where the main amount of work is done:
	s <- ptMetropolis(x,sampleSize[i],h/2)
	##
	commonName <- sprintf("%s-%s-branch-Sample",
		PREFIX,
		Branch
	)
	sid <- round(i*cs+r)
	maxsid <- round(cycl*cs+r)
	sampleFile=sprintf("%s-%0*x%s-cycle-%0*i-rank-%0*i-of-%0*i.RDS",
		commonName,
		scale(maxsid),sid,
		shortHash,
		scale(cycl),i,
		scale(cs),r,
		scale(cs),cs
		)
	saveRDS(s,file=sampleFile)
	## free up memory for next sample
	ar <- as.integer(round(100*attr(s,"acceptanceRate")))
	sr <- as.integer(round(100*attr(s,"swapRate")))
	rm(s)
	gc()
	pbdMPI::barrier()
	f <- dir(pattern=sprintf(
		"^%s-.*-cycle-%0*i-rank-[0-9]+-of-%0*i[.]RDS$",
		commonName,
		scale(cycl),i,
		scale(cs),cs)
	)
	X <- uqsa::gatherSample(f,beta)
	effar <- as.integer(round(100*mean(diff(sort(attr(X,"logLikelihood")))>1e-14)))
	attr(X,"beta") <- beta
	parMCMC <- as.numeric(tail(X,1))
	names(parMCMC) <- rownames(sb$Parameter)
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
time_ <- difftime(Sys.time(),start_time,units="min")
saveRDS(X,
	file=sprintf("%s-final-sample-%i-%i-lb100-%i-duration-%i-minutes.RDS",
		PREFIX,
		r,cs,
		round(log(beta)*100),
		round(time_)
	)
)
print(time_)
finalize()
print(warnings())
