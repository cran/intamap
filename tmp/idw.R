preProcess.idw = function(object, ...) {
	# perhaps first do some method-specific stuff, then
	# call the default method for this here:
	NextMethod()
}

estimateParameters.idw = function(object,formulaString, ..., idpRange = seq(0.1, 2.9, 0.1), nfold = 5) {
	# add parameter estimate
	mse = rep(NA, length(idpRange))
	for (i in seq(along = idpRange)) {
	  if (missing(formulaString))
      if ("formulaString" %in% names(object)) formulaString = object$formulaString else formulaString = as.formula("value ~ 1")
  	mse[i] = mean(krige.cv(formulaString, object$pointData, nfold = nfold, 
			set = list(idp = idpRange[i]))$residual ** 2)	
	}
  best = which(mse == min(mse))[1]
	object$inverseDistancePower = idpRange[best]
	print(paste("best idp value found is", object$inverseDistancePower, "rmse", sqrt(mse[best])))
	return(object)
}

spatialPredict.idw = function(object, ...) {
	if ("formulaString" %in% names(object)) formulaString = object$formulaString else formulaString = as.formula("value ~ 1")
	object$predictions = idw(formulaString, object$pointData, object$predictionLocations, 
		idp = object$inverseDistancePower,...)
	return(object)
}

postProcess.idw = function(object, ...) {
	# smooth over boundaries?

	# spatial aggregation?

	# find out what to output

	# write to data base
	return(object)
}
