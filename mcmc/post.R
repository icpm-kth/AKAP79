loadNamespace("SBtabVFGEN")
loadNamespace("rgsl")
loadNamespace("errors")

## this gets one sample matrix:
sampleQuality <- function(x){
	stopifnot(NROW(x)>NCOL(x))
	## load sbtab and experiments
	sb <- SBtabVFGEN::sbtab_from_tsv(dir("..",pattern="tsv$",full.names=TRUE))
	ex <- SBtabVFGEN::sbtab.data(sb)
	if (is.null(colnames(x))) colnames(x) <- rownames(sb$Parameter)
	modelName <- uqsa::checkModel(comment(sb),sprintf("./%s.so",comment(sb)))
	## load the sample of the requested temperature
	x <- uqsa::gatherSample(files,beta=1.0)
	colnames(x) <- rownames(sb$Parameter)
	## get the log-likelihood values
	l <- attr(x,"logLikelihood")
	plot(l,type="l")
	N <- length(l)
	## use a sub-sample of reduced size, omitting the burn-in
	## and calculate the auto-correlation length
	if (requireNamespace("hadron")){
		sub <- seq(ceiling(N*0.75),N)
		res <- hadron::uwerr(data=l[sub],pl=TRUE)
		cat("tau: ",res$tauint,"±",res$dtauint,"\n")
		tau <- errors::set_errors(res$tauint,res$dtauint)
	} else {
		a <- 1-mean(diff(sort(l))<1e-15+1e-15*max(l))
		tau <- 1.0/a
		tau <- errors::set_errors(tau,0.1*tau + 5) ## a wild guess
	}
	a <- mean(as.numeric(diff(sort(l))>1e-13))

	s <- uqsa::simulator.c(ex,modelName,parMap=uqsa::log10ParMap)
	i <- seq(ceiling(N/4),N,by=ceiling(as.numeric(tau)+1.0/a))
	X <- x[i,]
	L <- l[i]

	o <- order(L,decreasing=TRUE)
	Xo <- X[o,] # ordered
	attr(Xo,"logLikelihood") <- L[o]

	y <- s(t(Xo))
	for (i in seq_along(y)){
		rownames(y[[i]]$state) <- rownames(sb$Compound)
	}
	return(list(sample=Xo,tau=tau,simulations=y,sbtab=sb,experiments=ex))
}


## this gets many files
manyFilesSample <- function(files, beta=1.0){
	if (length(files)>20){
		cat(head(files),sep="\n")
		cat("[...]\n")
		cat(tail(files),sep="\n")
	} else {
		cat(files,sep="\n")
	}
	x <- uqsa::gatherSample(files,beta=1.0)
	q <- sampleQuality(x)
	comment(q) <- uqsa::determinePrefix(files)
	return(q)
}
