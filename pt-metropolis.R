#!/usr/bin/env Rscript
library(uqsa)
library(pbdMPI)

start_time <- Sys.time()
comm  <- 0
pbdMPI::init()
r <- pbdMPI::comm.rank(comm=comm)
cs <- pbdMPI::comm.size(comm=comm)
attr(comm,"rank") <- r
attr(comm,"size") <- cs
TMPDIR <- ifelse(dir.exists('/dev/shm'),'/dev/shm',Sys.getenv("TMPDIR",unset='/tmp'))
beta <- (1.0 - (r/cs))^2           # Parallel Tempering: inverse temperature

a <- commandArgs(trailingOnly=TRUE)
if (length(a)>0){
	N <- as.integer(a[1])
} else {
	N <- 50000                     # default sample size
}

m <- model_from_tsv()
o <- as_ode(m,cla=FALSE) |> write_and_compile()
ex <- experiments(m,o)
## -----------------------------------------------------------------
RWM <- high_level_metropolis(m,o,ex)                               # (serial) Random Walk Metropolis
pt <- high_level_metropolis(m,o,ex,beta=beta,comm=comm)            # parallel Random Walk Metropolis
##------------------------------------------------------------------
h <- 1e-4 # initially
x <- attr(RWM,"init")
for (j in seq(2)){
	h <- tune_step_size(RWM,x,h=h) # adjust using test chain
	pbdMPI::barrier()
	rank_file.rds <- tempfile(
		pattern=sprintf("AKAP79-pt-smmala-sample-%i-rank-%i",j,r),
		tmpdir=TMPDIR,
		fileext='.rds'
	)
	X <- pt(x,N,h)                 # the main amount of work is done here
	saveRDS(X,file=rank_file.rds)
	pbdMPI::barrier()
	f <- dir(TMPDIR,pattern=sprintf("^AKAP79-pt-smmala-sample-%i-rank-.*rds$",j),full.names=TRUE)
	Z <- uqsa::gatherSample(f,beta) # <- we only want this value of beta, in this rank's workspace
	ordered_sample.rds <- file.path(TMPDIR,sprintf("AKAP79-temperature-ordered-pt-smmala-sample-%i-for-rank-%i.rds",j,r))
	saveRDS(Z,file=ordered_sample.rds)
	x <- mcmc_init(beta,as.numeric(tail(Z,1)),attr(RWM,"sim"),dprior=dprior)
	## clean up this loop's intermediate work.
	file.remove(rank_file.rds)
}
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)
cli::cli_alert_info("final sample can be found here: {ordered_sample.rds}")

cat(
	sprintf(
		"rank %02i of %02i finished with an acceptance rate of %f and swap rate of %f.\n",
		round(r),
		round(cs),
		attr(Z,"acceptanceRate"),
		attr(Z,"swapRate")
	)
)
finalize()
