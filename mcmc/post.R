loadNamespace("SBtabVFGEN")
loadNamespace("rgsl")
loadNamespace("errors")

## this gets one sample matrix:
sampleQuality <- function(x, sb, ex){
	stopifnot(NROW(x)>NCOL(x))
	## load sbtab and experiments
	if (is.null(colnames(x))) colnames(x) <- rownames(sb$Parameter)
	modelName <- uqsa::checkModel(comment(sb),sprintf("./%s.so",comment(sb)))
	l <- attr(x,"logLikelihood")
	plot(l,type="l")
	N <- length(l)
	## use a sub-sample of reduced size, omitting the burn-in
	## and calculate the auto-correlation length
	sub <- seq(ceiling(N*0.75),N)
	if (requireNamespace("hadron")){
		res <- hadron::uwerr(data=l[sub],pl=TRUE)
		cat("tau: ",res$tauint,"±",res$dtauint,"\n")
		tau <- errors::set_errors(res$tauint,res$dtauint)
	} else  {
		A <- acf(l[sub],lag.max=round(length(sub)/2))
		a <- as.numeric(A$acf)
		tau <- sum(a[a>0.1])
		tau <- errors::set_errors(tau,1)
		## TODO: find out what the error is in this case
	}

	a <- 1-mean(diff(sort(l))<1e-15+1e-15*max(l))
	message(sprintf("acceptance^-1: %g",solve(a)))

	s <- uqsa::simulator.c(ex,modelName,parMap=uqsa::log10ParMap)
	i <- seq(ceiling(N/4),N,by=ceiling(as.numeric(tau)+1.0/a))
	X <- x[i,]
	L <- l[i]

	o <- order(L,decreasing=TRUE)
	Xo <- X[o,] # ordered
	attr(Xo,"logLikelihood") <- L[o]

	y <- s(t(Xo))
	for (i in seq_along(y)){
		rownames(y[[i]]$state) <- names(ex[[i]]$initialState)
	}
	return(list(sample=Xo,tau=tau,simulations=y,sbtab=sb,experiments=ex))
}


## this gets many files
manyFilesSample <- function(files, sb=NULL, ex=NULL, beta=1){
	if (is.null(sb)){
		sb <- SBtabVFGEN::sbtab_from_tsv(dir("..",pattern="tsv$",full.names=TRUE))
	}
	if (is.null(ex) && file.exists("experiments.RDS")){
		ex <- readRDS("experiments.RDS")
	} else if (is.null(ex) && file.exists("AKAP79.RDS")) {
		m <- SBtabVFGEN::sbtab_to_vfgen(sb)
		cl <- m$conservationLaws
		ex <- SBtabVFGEN::sbtab.data(sb,cl)
	} else {
		ex <- SBtabVFGEN::sbtab.data(sb)
	}
	if (length(files)>20){
		cat(head(files),sep="\n")
		cat("[...]\n")
		cat(tail(files),sep="\n")
	} else {
		cat(files,sep="\n")
	}
	x <- uqsa::gatherSample(files,beta=1.0)
	q <- sampleQuality(x, sb, ex)
	comment(q) <- uqsa::determinePrefix(files)
	return(q)
}

errorBars <- function(x,y,e,...){
	arrows(x,y,x,y-e,angle=90,...)
	arrows(x,y,x,y+e,angle=90,...)
	points(x,y,...)
}

#' this function tags the number as being in centimeters
#' @param x a numer
#' @return a number with attributes describing the unit
cm <- function(x){
	attr(x,"kind") <- "m"
	attr(x,"scale") <- -2
	attr(x,"exponent") <- 1
	return(x)
}

#' unit conversion helper function
#'
#' currently only converts between meters and inches, more to add
#' later.
#'
#' @param value some value with unit attributes
#' @param unit desired unit
#' @return a value in changed units
`%as%` <- function(value,unit){
	y <- value
	if (unit == "inches") {
		stopifnot(attr(value,"kind")=="m")
		y <- value*10^attr(value,"scale") * 39.37
		attr(y,"kind") <- "inches"
		attr(y,"scale") <- 1
		attr(y,"exponent") <- 1
	}
	return(y)
}

makePlots <- function(Q,yl=c(90,200),...){
	ex <- Q$experiments
	sb <- Q$sbtab
	s_ <- Q$simulations
	par(...)
	for (i in seq_along(ex)){
		t_ <- ex[[i]]$outputTimes
		y_ <- t(ex[[i]]$outputValues)
		e_ <- t(ex[[i]]$errorValues)
		h_ <- s_[[i]]$func['AKAR4pOUT',,]
		matplot(t_,h_,
			type="l",
			main=names(ex)[i],
			xlab="t",
			ylab="AKAR4p",
			lty=1,lwd=2,
			col=rgb(0.2,0,0.8,0.5),
			ylim=yl
		)
		errorBars(t_,y_,e_,length=0.01)
		lines(t_,h_[,1],col="red3")
	}
}
