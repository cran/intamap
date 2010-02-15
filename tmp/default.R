estimateParameters.default = function(object, ...) {
	stop("there is no default parameter estimation method")
}

spatialPredict.default = function(object, ...) {
	stop("there is no default prediction method")
}

postProcess.default = function(object, ...) {
	# smooth over boundaries?

	# spatial aggregation?

	# find out what to output
  	object$outputTable = getOutputTable(object)

	# write to data base

	return(object)
}

#		blockFat=TRUE,??

getOutputTable = function(object) {
	what = object$outputWhat
	pred = object$predictions
	nCols = length(what)
	out = matrix(NA, length(pred[[1]]), nCols)
	ns = c("x", "y")
	for (i in 1:nCols) {
		if (names(what)[i] == "mean") {
			out[,i] = pred[[1]]
			ns = c(ns, "mean")
		} else if (names(what)[i] == "variance") {
			out[,i] = pred[[2]]
			ns = c(ns, "variance")
		} else if (names(what)[i] == "quantile") {
			out[,i] = qnorm(what[[i]], pred[[1]], sqrt(pred[[2]]))
			ns = c(ns, paste("quantile", what[[i]], sep=""))
		} else if (names(what)[i] == "cumdistr") {
			out[,i] = pnorm(what[[i]], pred[[1]], sqrt(pred[[2]]))
			ns = c(ns, paste("cumdistr", what[[i]], sep=""))
		} else if (names(what)[i] == "excprob") {
			out[,i] = 1 - pnorm(what[[i]], pred[[1]], sqrt(pred[[2]]))
			ns = c(ns, paste("excprob", what[[i]], sep=""))
		} else
			stop(paste("unknown request: ", names(what)[i]))
	}
  	ret = cbind(coordinates(object$predictions), out)
	names(ret)[[2]] = ns
	ret
}
