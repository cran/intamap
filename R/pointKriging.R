
estimateParameters.automap = function(object,...) {

  dots = list(...)
  if ("debug.level" %in% names(dots)) 
  	debug.level = dots$debug.level 
  else 
    debug.level = object$params$debug.level
  observations = object$observations
  depVar = as.character(object$formulaString[[2]])

    
#estimate Anisotropy
  if (object$params$doAnisotropy) {
    object = estimateAnisotropy(object) 
    if (object$anisPar$doRotation && all(as.character(object$formulaString[[3]])=="1")){
			#rotate Data
				objTemp = object
        objTemp$observations=rotateAnisotropicData(objTemp$observations,objTemp$anisPar)
				#Estimate Variogram Model
        if (dim(coordinates(object$predictionLocations))[1] > 100000) 
            model = c("Sph", "Exp", "Gau") else model = c("Sph", "Exp", "Gau", "Ste")
				afv = autofitVariogram(objTemp$formulaString, objTemp$observations,
                   verbose = (debug.level >=2), model = model, ...)
			  vario = afv$var_model				
				ovar = var(observations[,depVar]@data)
      	if ((vario$model[2]  == "Gau" | (vario$model[2] == "Ste" && vario$kappa[2] > 2)) 
            && vario$psill[1] <= ovar/1e5 ) vario$psill[1] = ovar/1e5   
        #Combine the isotropic Model with anisotropy parameters
				vario$anis1[2]=1/objTemp$anisPar$ratio
				vario$ang1[2]=90-objTemp$anisPar$direction
        if (vario$ang1[2] < 0) vario$ang1[2] = vario$ang1[2] + 180				
				object$variogramModel=vario
#				vario$range=vario$range/objTemp$anisPar$ratio
  	} else {
		afv = autofitVariogram(object$formulaString,observations,verbose=(debug.level >=2),...)
    	object$variogramModel = afv$var_model
    }
  } else { 
  	afv = autofitVariogram(object$formulaString,observations,verbose=(debug.level >=2),...)
  	object$variogramModel = afv$var_model
  }
  if (debug.level >=2) print(object$variogramModel)
  object$sampleVariogram = afv$exp_var
  return(object)
}

spatialPredict.automap = function(object, nsim = 0, ...) {
# This function does not take the clusters properly into account at the moment. 
# Variograms should be estimated separately, prediction locations needs to 
# be associated with a cluster and we need to figure out what to do in the case 
# of 
# 
  dots = list(...)
  if ("nmax" %in% names(dots)) {
    nmax = dots$nmax
  } else nmax = object$params$nmax
  if ("debug.level" %in% names(dots)) debug.level = dots$debug.level else 
    debug.level = object$params$debug.level
  
    if (! "variogramModel" %in% names(object)) object = estimateParameters(object,...)
    
    pred = krige(object$formulaString,object$observations, 
           object$predictionLocations,object$variogramModel,nsim=nsim,nmax = nmax,debug.level = debug.level)
    if (nsim >0) {
      pred2 = krige(object$formulaString,object$observations, 
           object$predictionLocations,object$variogramModel,nmax = nmax,debug.level = debug.level)
      pred@data = cbind(pred2@data,pred@data)
    }

    object$predictions = pred
    if ("MOK" %in% names(object$outputWhat) | "IWQSEL" %in% names(object$outputWhat))
      object$predictions = unbiasedKrige(object,debug.level = debug.level,...)$predictions
  object
}


estimateParameters.yamamoto = function(object,...) {
  estimateParameters.automap(object,...)
}




spatialPredict.yamamoto = function(object, nsim = 0, ...) {
# 
  dots = list(...)
  if ("nmax" %in% names(dots)) {
    nmax = dots$nmax
  } else nmax = object$params$nmax
  formulaString = object$formulaString
  if ("debug.level" %in% names(dots)) debug.level = dots$debug.level else 
    debug.level = object$params$debug.level

  if (!"variogramModel" %in% names(object)) {
    afv = autofitVariogram(object$formulaString,object$observations,object$predictionLocations)
  	object$variogramModel = afv$var_model
	object$sampleVariogram = afv$exp_var
  }
                     
  predictions = yamamotoKrige(formulaString,object$observations, 
            object$predictionLocations,object$variogramModel,nsim=nsim,nmax = nmax,...)
  object$predictions = predictions
    if ("MOK" %in% names(object$outputWhat) | "IWQSEL" %in% names(object$outputWhat))
      object$predictions = unbiasedKrige(object,debug.level = debug.level,nsim = nsim, nmax = nmax,...)$predictions
  object
}

