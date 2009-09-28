preProcess.idw = function(object, ...) {
	# perhaps first do some method-specific stuff, then
	# call the default method for this here:
	NextMethod()
}

estimateParameters.idw = function(object, ..., idpRange = seq(0.1, 2.9, 0.1), nfold = 5) {
	# add parameter estimate
	mse = rep(NA, length(idpRange))
  if ("formulaString" %in% names(object)) formulaString = object$formulaString else formulaString = as.formula("value ~ 1")
  dots = list(...)
  if ("nmax" %in% names(dots)) {
    nmax = dots$nmax
  } else nmax = object$params$nmax
 	for (i in seq(along = idpRange)) {
  	mse[i] = mean(krige.cv(formulaString, object$observations, nfold = nfold, 
  	nmax = nmax, set = list(idp = idpRange[i]))$residual ** 2)	
	}
  best = which(mse == min(mse))[1]
	object$inverseDistancePower = idpRange[best]
	print(paste("best idp value found is", object$inverseDistancePower, "rmse", sqrt(mse[best])))
	return(object)
}

spatialPredict.idw = function(object, ...) {
  dots = list(...)
  if (!all(names(object$outputWhat) == "mean"))
    stop(paste("It is not possible to request other prediction types than mean for method idw",
          "requested",names(object$outputWhat))) 
  if ("formulaString" %in% names(object)) formulaString = object$formulaString else formulaString = as.formula("value ~ 1")
  if ("nmax" %in% names(dots)) {
    nmax = dots$nmax
  } else nmax = object$params$nmax
 	object$predictions = idw(formulaString, object$observations, object$predictionLocations, 
		nmax = nmax, idp = object$inverseDistancePower, debug.level = object$params$debug.level)
	return(object)
}

postProcess.idw = function(object, ...) {
	# smooth over boundaries?

	# spatial aggregation?

	# find out what to output

	# write to data base
  object$predictions = object$predictions[,-2]
  object = NextMethod(object,...)
	return(object)
}
