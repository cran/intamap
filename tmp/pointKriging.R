estimateParameters.automap = function(object,...) {

  observations = object$pointData
  vario = autofitVariogram(object$formulaString,observations)
  if (inherits(vario,"autofitVariogram")) # i.e.,coming from package automap
  	vario = vario$var_model
  object$variogramModel = vario
  return(object)
}

spatialPredict.point = function(object,...) {
  
}


spatialPredict.automap = function(object,...) {
# This function does not take the clusters properly into account at the moment. 
# Variograms should be estimated separately, prediction locations needs to 
# be associated with a cluster and we need to figure out what to do in the case 
# of 
# 
  nclus = ifelse(object$params$doSegmentation & !is.null(object$clusters$clusterNumber),
                 object$clusters$clusterNumber,1) 
  for (iclus in 1:nclus) {
    if ("variogramModel" %in% names(object)) {
      pred = krige(object$formulaString,object$pointData, object$predictionLocations,object$variogramModel,nmax = 50)
    } else {
      predObj = autoKrige(object$formulaString,object$pointData,object$predictionLocations )
      pred = predObj$krige_output
      object$variogramModel = predObj$var_model
    }
    predictType = object$params$predictType

    object$predictions = pred
    if (!all(is.na(object$params$predictType$threshCor)))
      object$predictions = unbiasedKrige(object,...)$predictions
  }
  object
}
