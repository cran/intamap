estimateParameters.linearVariogram = function(object, ...) {
	# no parameters to be estimated...
	return(object)
}

spatialPredict.linearVariogram = function(object, formulaString, ...) {
  if (missing(formulaString))
      if ("formulaString" %in% names(object)) formulaString = object$formulaString else formulaString = as.formula("value ~ 1")
	object$predictions = krige(formulaString, object$pointData, object$predictionLocations, vgm(1, "Lin", 0))
	return(object)
}
