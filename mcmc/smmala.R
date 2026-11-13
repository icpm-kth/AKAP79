#!/usr/bin/env Rscript
## ----setup--------------------------------------------------------------------
library(uqsa)
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
R <- round(cs)

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
m <- model_from_tsv("..")
o <- as_ode(m) |> write_and_compile()
ex<- experiments(m,o)

options(mc.cores = length(ex))
## ----default------------------------------------------------------------------
n <- length(ex[[1]]$input)

median <- m$Parameter$median
p0 <- values(m$Parameter)        # this already is in log10-space
stdv <- uncertainty(m$Parameter) # this as well
if (is.null(stdv) || any(!is.numeric(stdv)) || any(is.na(stdv))) {
	warning("no standard error («stdv» field) in model files")
	stdv <- p0*0.5 + 0.5 + 0.5*max(p0)
}
dprior <- dNormalPrior(mean=median,sd=stdv)
rprior <- rNormalPrior(mean=median,sd=stdv)
gprior <- gNormalPrior(mean=median,sd=stdv)

fiPrior <- solve(diag(stdv, length(p0)))

## ----simulate-----------------------------------------------------------------
sim <- simfi(ex,o,log10ParMap)

X <- NULL # this is to suppress one intial warning
sampleSize <- round(exp(seq(log(1000),log(N),length.out=cycl)))

## ----sample-------------------------------------------------------------------
smmala <- smmala_update(
	simulate=sim,
	dprior=dprior,
	gprior=gprior,
	fisherInformationPrior=fiPrior
)

serial_smmala <- mcmc(smmala)
ptSMMALA <- mcmc_mpi(smmala,comm)
parMCMC <- p0
for (i in seq(cycl)){
	rm(X)
	x <- mcmc_init(
		beta,
		parMCMC,
		simulate=sim,
		dprior=dprior,
		gradLogLikelihood = uqsa::gllf(log10ParMapJac),
		gprior=gprior,
		fisherInformation= uqsa::fi(log10ParMapJac)
	)
	if (i < round(cycl/2)){
		h <- tune_step_size(serial_smmala,x,h=1e-7)
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
	sampleFile=sprintf("%s-%0*x%s-cycle-%0*i-rank-%0*i-of-%0*i.rds",
		commonName,
		scale(maxsid),sid,
		shortHash,
		scale(cycl),i,
		scale(cs),r,
		scale(cs),cs)
	saveRDS(s,file=sampleFile)
	## free up memory for next sample
	ar <- as.integer(round(100*attr(s,"acceptanceRate")))
	sr <- as.integer(round(100*attr(s,"swapRate")))
	rm(s)
	gc()
	pbdMPI::barrier()
	f <- dir(pattern=sprintf("^%s-.*-rank-[0-9]+-of-%0*i[.]rds$",commonName,scale(cs),cs))
	X <- gatherSample(f,beta)
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
saveRDS(X,file=sprintf("%s-final-sample-%i-%i-lb100-%i.rds",PREFIX,r,cs,round(log(beta)*100)))
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)
finalize()
print(warnings())
