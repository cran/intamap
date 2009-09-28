

estimateParameters.transGaussian = function(object, lambda, significant = TRUE, ...) {
  observations = object$observations
  formulaString = object$formulaString
  dataObs = observations[[as.character(formulaString[[2]]) ]]
  if (missing(lambda)) {
    test = doNonGauss(dataObs)
#    if (min(dataObs <=0)) {
#      pcor = sqrt(var(dataObs))/10000
#      if (min(dataObs) + pcor > 0 & length(dataObs <= 0) < length(dataObs)/4) {
#        dataObs = dataObs+pcor
#        object$TGcorrection = pcor
#      } 
#    }
    if (test || !significant)  lambda = bcFit(dataObs) else lambda = 1
#    if (lambda == 1 && !is.null(object$TGcorrection)) {
#       object$TGcorrection = 0
#       dataObs = observations[[as.character(formulaString[[2]]) ]]
#    } 
  }
  dataObsBC = bcTrans(dataObs,lambda)
  object$observations[[as.character(formulaString[[2]]) ]] = dataObsBC
  object$lambda = lambda
  object = estimateParameters.automap(object,...)
  object$observations = observations
  object
}


spatialPredict.transGaussian = function(object, nsim = 0, ...) {
  dots = list(...)
  if ("nmax" %in% names(dots)) {
    nmax = dots$nmax
  } else nmax = object$params$nmax
  if ("debug.level" %in% names(dots)) debug.level = dots$debug.level else
    debug.level = object$params$debug.level
    if (! "variogramModel" %in% names(object)) object = estimateParameters(object,...)
    if ("lambda" %in% names(object)) lambda = object$lambda else lambda = 1
    observations = object$observations
    formulaString = object$formulaString
#    if (!is.null(object$TGcorrection)) 
#      observations[[as.character(formulaString[[2]])]] = 
#             observations[[as.character(formulaString[[2]])]]+ object$TGcorrection
    pred = krigeTg(formulaString,observations,
           object$predictionLocations,object$variogramModel,nmax = nmax,
           debug.level = debug.level, lambda = lambda)
    pred = pred[c("var1TG.pred","var1TG.var")]
    names(pred) = c("var1.pred","var1.var")
    if (nsim >0) {
      pred2 = krigeTg(formulaString,observations,
           object$predictionLocations,object$variogramModel,nmax = nmax,
           debug.level = debug.level, nsim = nsim, lambda = lambda)
      pred@data = cbind(pred@data,pred2@data)
    }
#    if (!is.null(object$TGcorrection)) pred$var1.pred = pred$var1.pred - object$TGcorrection
    object$predictions = pred
    if ("MOK" %in% names(object$outputWhat) | "IWQSEL" %in% names(object$outputWhat))
      object$predictions = unbiasedKrige(object,debug.level = debug.level,...)$predictions
  object
}




bcFit = function(z, lambda = seq(-3,3,1/100), eps = 1/50) {
	bc = boxcox(z~1, lambda = lambda, plotit = FALSE)
	m = length(bc$x)
	lambda.index = (1:m)[bc$y == max(bc$y)][1]
	if (lambda.index == 1 || lambda.index == m)
		warning("optimal lambda found at the edge of search range")
	bc$x[lambda.index]
}


bcTrans = function(z, lambda) {
	if (lambda == 1.0)
		zt = z
	else if (abs(lambda) > 0.0)
		zt <- (z^lambda - 1)/lambda
	else zt <- log(z)
	zt
}
