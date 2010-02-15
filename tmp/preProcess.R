#####################################
#
# getIntamapParams - function for setting intamap parameters
# Input - parameters to be set
# Output - parameter list containing:
#
# formulaString = formula string for parameter estimation and 
#                 prediction functions
# doAnisotropy = defining whether anisitropy should be calculated
# removeBias = definition which biases to remove
# addBias = Biases that can be added again (if possible, e.g. regional biases 
#           might be added, localBias cannot be added
# biasRemovalMethod = Which method to use for bias removal, 
#                    "UK" (universal kriging using whole dataset) 
#                     or "LM" (local methods)
# doCluster = if the prediction region should be divided into a group of clusters
# maxCluster = The maximum number of clusters if clusters will be used
# numberOfClusters = If a fixed number of clusters is to be used (still 
#                    to be decided exactly how to use this
# nmax = maximum number of neighbours to use for kriging
# predictType = List of different prediction types (all TRUE/FALSE)
#     threshCor = Prediction corrected for shrinkage in kriging predictions
#                 MOK Modified ordinary kriging predictor
#                 IWQSEL IWQSEL-predictor
#     block = Predictions for block 
#     blockFat = Estimated fraction above threshold (if block kriging)
# thresh = Threshold for Fraction above threshold (blockFat) and exceedance probability (exc)
# isEmergency = parameter to disable certain functions, e.g. bias correction
# confProj = if Projections should be conformed, setting intCRS as 
#            interpolation projection
# processType = gaussian, nonGaussian, logNormal
########################################



getIntamapParams = function(oldPar, newPar,...){
  dots = list(...)
  twoLists = FALSE
  if (!missing(oldPar) && !inherits(oldPar,"IntamapParams")) {
    if (!missing(newPar)) {
      newPar2 = newPar 
      newPar = oldPar 
      twoLists = TRUE
    }  else newPar = oldPar
    oldPar = getIntamapDefaultParams()
  } else if (missing(oldPar)) {
    oldPar = getIntamapDefaultParams()
  }
#  if (!missing(newPar) & "thresh" %in% names(newPar) & !is.list(newPar$thresh)) newPar$thresh = as.list(newPar$thresh)
#  if (!missing(oldPar) & "thresh" %in% names(oldPar) & !is.list(oldPar$thresh)) oldPar$thresh = as.list(oldPar$thresh)
#  if (!missing(newPar) & "quantiles" %in% names(newPar$predictType) & !is.list(newPar$predictType$quantiles)) 
#           newPar$predictType$quantiles = as.list(newPar$predictType$quantiles)
#  if (!missing(oldPar) & "quantiles" %in% names(oldPar$predictType) & !is.list(oldPar$predictType$quantiles)) 
#           oldPar$predictType$quantiles = as.list(oldPar$predictType$quantiles)
  if (!missing(newPar)) oldPar = modifyList(oldPar,newPar)
  if (twoLists) oldPar = modifyList(oldPar,newPar2)
  if (length(dots) >0) oldPar = modifyList(oldPar,dots)
  class(oldPar) = "IntamapParams"
  return(oldPar)
}


getIntamapDefaultParams = function(doAnisotropy = FALSE, 
  removeBias= NA,  addBias = NA, biasRemovalMethod = "LM", doSegmentation = FALSE, maxCluster = 0, 
  numberOfClusters = 1,nmax = Inf, ngrid = 100, nsim = 100, block=numeric(0), processType="gaussian",
  predictType = list(threshCor=NA, yamamoto = FALSE), 
  isEmergency = FALSE,confProj = FALSE,... ) {
return(list(doAnisotropy = doAnisotropy, removeBias = removeBias, addBias = addBias,
  biasRemovalMethod = biasRemovalMethod, doSegmentation = doSegmentation, 
  maxCluster = maxCluster,
  numberOfClusters = numberOfClusters, nmax = nmax, ngrid = ngrid, nsim = nsim, processType = processType,
  predictType = predictType, 
  isEmergency = isEmergency, confProj = confProj,... ))
}

createIntamapObject = function(pointData, blockData, formulaString, predictionLocations=100,
  targetCRS,boundaries,boundaryLines,intCRS, params=list(),boundFile,lineFile,classes="idw",
  outputWhat = list(mean=TRUE, variance=TRUE), blockWhat = "none",...) {
  object = list()
  dots = list(...)
  if (!missing(pointData) && !extends(class(pointData),"Spatial")) 
  	stop("pointData not object of class Spatial*")
  if (!missing(blockData) && !extends(class(blockData),"Spatial")) 
  	stop("blockData not object of class Spatial*")
  if (!missing(pointData) & !missing(blockData)) 
  	stop("Warning: only few methods can handle both point observations and observations with support")
  if (missing(pointData) & missing(blockData)) 
  	stop("Neither pointData or blockData submitted, cannot perform interpolation without data")
  if (!missing(pointData)) 
  	object$pointData = pointData 
  if (!missing(blockData))
  	object$blockData = blockData 
  if (missing(formulaString)) {
    if (!missing(blockData)) {
      if ("value" %in% names(blockData)) {
        formulaString = "value~1"
      } else formulaString = paste(names(blockData)[1],"~1")
    }        
    if (!missing(pointData)) {
      if ("value" %in% names(pointData)) {
        formulaString = "value~1"
      } else formulaString = paste(names(pointData)[1],"~1")
    }        
    print(paste("createIntamapProject: formulaString is missing, using: ",formulaString))
  }
  if (!inherits(formulaString,"formula")) 
  	formulaString = as.formula(formulaString)
  object$formulaString = formulaString
  if (!is.numeric(predictionLocations)) {
    if (extends(class(predictionLocations),"Spatial")) {
      object$predictionLocations = predictionLocations
    } else stop("predictionLocations not spatial object or number of samples")    
  } else {
    if (!missing(boundaries)) {
      warning("createIntamapObject: No prediction locations submitted - sampling from boundaries")
      object$predictionLocations = spsample(boundaries,predictionLocations,"regular") 
    } else if (!missing(pointData)) {
      warning("createIntamapObject: No prediction locations submitted - sampling from bbox of pointData")
      object$predictionLocations = spsample(bbox(pointData),predictionLocations,"regular") 
    } else if (!missing(blockData)) {
      warning("createIntamapObject: No prediction locations submitted - finding point prediction locations")
      object$predictionLocations = spsample(bbox(blockData),predictionLocations,"regular")
    }
  }
  if (!missing(boundaries)) {
    object$boundaries = boundaries
    if (!is.na(proj4string(boundaries))) object$boundCRS = proj4string(boundaries)
  }
  if (!missing(targetCRS)) object$targetCRS = targetCRS
  if (!missing(intCRS)) object$targetCRS = intCRS
  if (!missing(pointData) && "regCode" %in% names(pointData)) 
              object$regCode = unique(pointData$regCode)
  if (!missing(blockData) && "regCode" %in% names(blockData)) 
              object$regCode = unique(blockData$regCode)

  if (missing(params)) object$params = getIntamapParams() else 
    object$params = getIntamapParams(params) 
  if (!missing(boundaries)) {
    objectboundaries = boundaries
  } else if (!missing(boundFile)) {
  	# EJP:
    #if (require(maptools)) object$boundaries = readShapePoly(boundFile) else
    #  warning("maptools not installed, not able to read boundaries")
	object$boundaries = readOGR(".", boundFile)
  }
  if (!missing(boundaryLines)) {
    object$boundaryLines = boundaryLines
  } else if (!missing(lineFile)) {
    object$lineFile = lineFile
    load(lineFile)
    object[[boundaryLines]] = lineFile
  }
  if (length(names(dots))>0)
  	object = modifyList(object,dots)
  object$outputWhat = outputWhat
  object$blockWhat = blockWhat
  class(object) = classes
  return(object)
}




#########################################################
# ConformProjection
#
# Input: Intamap object
#
# Output: Intamap object with projections of pointData and
#         predictionLocations equal to
#         (a) targetCRS - if not longlat
#         (b) observationCRS - if not longlat
#         (c) predictionCRS - if not longlat
#         (d) intCRS = "+init=epsg:3035" if both above are longlat
#
############################################################

conformProjections = function(object) {
  observations = object$pointData
  predictionLocations = object$predictionLocations
  obsCRS = proj4string(observations)
  predCRS = proj4string(predictionLocations)
  targetCRS = object$targetCRS
  if ("intCRS"%in% names(object)) {
    intCRS = object$intCRS
  } else {
      if ("longlat" %in% CRSargs(CRS(targetCRS))) {
          intCRS = targetCRS
      } else {
          if ("longlat" %in% obsCRS) {
              intCRS = obsCRS
              object$intCRS = intCRS
          } else {
              if ("longlat" %in% predCRS) {
                  intCRS = predCRS
                  object$intCRS = intCRS
              } else {
                  intCRS = "+init=epsg:3035"
                  object$intCRS = intCRS
              }
          }
      }
  }
  if (CRSargs(CRS(obsCRS)) != CRSargs(CRS(intCRS))) 
    object$pointData = spTransform(observations,CRS(intCRS))
  if (CRSargs(CRS(predCRS)) != CRSargs(CRS(intCRS))) 
    object$predictionLocations = spTransform(predictionLocations,CRS(intCRS))
  if (!is.null(object$boundaries)) {
  	boundaries = object$boundaries
    boundCRS = proj4string(object$boundaries)
  	if (CRSargs(CRS(boundCRS)) != CRSargs(CRS(intCRS))) 
    	object$boundaries = spTransform(boundaries,CRS(intCRS))
  }

	return(object)
}





######################################
# cleanData - function for removing troublesome parts of the data
#
# Input: Intamap object
# Output: Intamap object with data cleaned
#         This version only removes duplicates
#
############################################

cleanData = function(object) {
# What to do with multiple observations?
object2 = object[!duplicated(coordinates(object)),]
return(object2)
}

findElevation = function(object) {
  if (!is.null(object$elevation)) {
    return(object$elevation)
  } else {
    return(rep(NA,dim(object)[1]))
  }
}


setLocalGroup.default = function(object,...) {
#  print(paste("setLocalGroup.default",class(object)))
  return(object)
}

setLocalGroup.SpatialPointsDataFrame = function(object,...) {
#  print(paste("setLocalGroup.SpatialPointsDataFrame",class(object)))
  object = NextMethod()
  return(object)
}

setLocalGroup.automap = function(object,...) {
#  print(paste("setLocalGroup.automap",class(object)))
  object = NextMethod()
  return(object)
}

cleanRegions.default = function(regions,...) {
  return(regions)
}

cleanRegions.SpatialPointsDataFrame = function(regions,...) {
#  print(paste("cleanRegions.SpatialPointsDataFrame",class(object)))
  regions = NextMethod()
  return(regions)
}

cleanRegions.automap = function(regions,...) {
#  print(paste("cleanRegions.automap",class(object)))
  setLocalGroup.SpatialPointsDataFrame(regions) #EJP: is this supposed to happen?
  regions = NextMethod()
  return(regions)
}



#############################################
# Default preProcessing function
#
# Input: intamap object
#        lgFUN - function for grouping data locally
#        cid - country ID - or other regional grouping factor
#
# Output: intamap object with the following added
#         localBias - data frame with local biases
#         regionalBias - data frame with biases between countries
#         Modifications of pointData
#         Elevations added (still not properly implemented
#         Duplicated data observations deleted
#         Projections are conformed and an interpolation projection
#                     set if they do not conform, or dont have a projection
#
###########################################

preProcess.default = function(object,...) {
  dots = list(...)
  params = object$params
 	if ("formulaString" %in% names(object)) formulaString = object$formulaString else formulaString = as.formula("value ~ 1") 
  depVar = as.character(formulaString[[2]])
  if (params$confProj) 
  	object = conformProjections(object)
  observations = object$pointData
  intCRS = ifelse ("intCRS" %in% names(object), object$intCRS,NA)  
  if (!is.null(object$boundaries) && is.projected(object$boundaries) && is.null(object$boundCRS)) {
    boundCRS = proj4string(object$boundaries)
    object$boundCRS = boundCRS
  } else
    boundCRS = ifelse (is.null(object$boundCRS), NA, object$boundCRS)
  observations = cleanData(observations) # practically unimplemented
  observations$elev = findElevation(observations)  # practically unimplemented
  removeBias = params$removeBias
#  
  if ("regCode" %in% names(object)) {
    regCode = object$regCode
  } else if ("cid" %in% names(dots)) {
    cid = dots[["cid"]]
    icid = which(names(observations) == cid)
    regCode = observations[,icid]
    regCode = unique(regCode)
    object$regCode = regCode
  } else {
    regCode = NULL
  }
  if (!is.na(removeBias[[1]]) & !params$isEmergency) {
    if ("localBias" %in% removeBias ) {
      if ("localBias" %in% names(object)) {
        localBias = object$localBias
      } else if("gid" %in% names(dots)) {
        gid = dots[["gid"]]
      } else if ("lgFUN" %in% names(dots)) {
        FUN = match.fun(dots$lgFUN)
        observations = FUN(observations)
        gid = "group"               
      } else if ("group" %in% names(observations)) {
        gid = "group"      
      } else {
        observations = setLocalGroup(observations)
        gid = "group"               
      }
    }
    if (params$biasRemovalMethod == "UK") {
      if ("localBias" %in% removeBias ) {
        localBias = findLocalBias(observations,regCode,formulaString=formulaString,...)
        object$localBias = localBias
        observations = removeLocalBias(observations,localBias,depVar)
      }
      class(observations) = class(observations)[[1]]
      if ("regionalBias" %in% removeBias  & !params$isEmergency) {
        regionalBiasUK = findBiasUK(object)
        object$regionalBias = regionalBiasUK
      }
    } else if (params$biasRemovalMethod == "LM") {
      if ("localBias" %in% removeBias ) {
        localBias = findLocalBias(observations,regCode,formulaString=formulaString,...)
        object$localBias = localBias
        observations = removeLocalBias(observations,localBias,depVar)
      }
      class(observations) = class(observations)[[1]]
      if ("regionalBias" %in% removeBias  & !params$isEmergency) {
        if ("regionalBias" %in% names(object)) {
          regionalBias = object$regionalBias
        } else {
          if ("boundaryLines" %in% names(object)){ 
            boundaryLines = object$boundaryLines
          } else if ("boundaries" %in% names(object)) {
            boundaryLines = findBoundaryLines(regions=object$boundaries,
                projOrig = ifelse(params$confProj, intCRS,boundCRS), 
                projNew = intCRS)
                object$boundaryLines = boundaryLines
          } else warning("No boundaryLines or boundaries in object")
#
#          projOrig = "+proj=laea +lat_0=48 +lon_0=9 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m", 
#          projNew = "+init=epsg:3035" )
#  countryBoundaries = findCountryBoundaries("d:/svn/intamap/rdata/countryLim", fileType="Shape",
#      projOrig = "+proj=laea +lat_0=48 +lon_0=9 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m", 
#      projNew = "+init=epsg:3035" )
#      save(countryBoundaries,file = "cuntryBoundaries.rda")
          if (exists("boundaryLines")) regionalBias = findRegionalBias(observations,boundaryLines,formulaString=formulaString)
        }
        if (exists("boundaryLines")) {
          observations = removeRegionalBias(observations,regionalBias$regionalBias,depVar)
          object$regionalBias = regionalBias
        }
      }
    }
  }
  object$pointData = observations
  return(object)
}

# kobj = preProcess(krigingObject,lgFUN = setLocalGroup.eurdep)
