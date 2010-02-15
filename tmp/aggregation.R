
#getIntamapParams = function(formulaString= as.formula(value~1),doAnisotropy = FALSE, 
#  removeBias= NA,  addBias = NA, biasRemovalMethod = "LM", doSegmentation = FALSE, maxCluster = 0, 
#  numberOfClusters = 1,
#  predictType = data.frame(pred = TRUE, var=TRUE,exc=TRUE,threshCor=NA,yamamoto = FALSE,
#                  blockMean=TRUE,blockFat=TRUE),
#  thresh=200, isEmergency = FALSE,confProj = FALSE )
  
blockPredict = function(object,...) {
  object = spatialPredict.block(object,...)
  return(object)
} 
 

spatialPredict.block = function(object,...) {
# Define which methods can actually handle block kriging directly
  blockMethods = c("automap","idw")
# Extract some of the parameters and variables from the object
# - for easier access
  params = object$params
  dots = list(...)
  block = params$block
  if ("block" %in% names(dots)) block = dots$block
  nmax = params$nmax
  predictType = params$predictType
  threshCor = params$predictType$threshCor
  thresh = params$thresh
  yamamoto = predictType$yamamoto
  observations = object$pointData
  formul = object$formulaString
  predictionLocations = object$predictionLocations
  if (!inherits(predictionLocations,"SpatialPolygonsDataFrame") &
      !inherits(predictionLocations,"SpatialLinesDataFrame") &
      !length(block) >0 ) {
      warning("prediction locations are not blocks, polygons or lines")
      warning("calling point prediction")
      return(spatialPredict(object,...))
  }
  
  if (!("variogram" %in% names(object)) &
      !inherits(object,"idw")) object = estimateParameters(object)
  if (params$numberOfClusters == 1 & 
      params$processType %in% c("gaussian","logNormal") &
      sum(class(object) %in% blockMethods) >= 1) {
# Check if we might be able to use an analytical solution - 
# Possible if we want a prediction for Gaussian or logNormal distribution, only 1 cluster
# Only possible for certain methods - now hard coded into blockMethods
    if (inherits(object,"idw")) return(spatialPredict(object,...))
    form = object$formulastring
    depVar = as.character(form[1])
# We have to discuss whether it makes sense to have the test of logNormal here, 
# as it can easily create some inconsistencies.
# We have to be sure that the variogram is also found from the logarithmised data
    if (params$processType == "logNormal") observations[[depVar]] = log(observations[[depVar]])
      predictions = krige(object$formulaString,observations, predictionLocations,object$variogramModel,nmax = nmax)
    if (params$processType == "logNormal") {
      predictions$var1.pred = exp(predictions$var1.pred + predictions$var1.var/2) 
      warning("Note that the back-calculated lognormal predictor is only approximate, as the Lagrange parameter is not included")
    }
    if (predictType$exc & length(thresh)>0) {
      for (it in 1:length(thresh)){
        threshi = thresh[it]
        iname = paste("exc",thresh,sep="")
        predictions[[iname]] = 
        apply(predictions,MARGIN=1,thresh = threshi,FUN=function(arr,threshi) sum(I(arr>threshi))/length(arr))    
      }
    }
    print("performed ordinary block kriging")
  }
  
# Point predictions/simulations are necessary
  if (!is.na(threshCor) | (params$numberOfClusters > 1) | 
     !(all(!c(predictType$yamamoto,predictType$blockFat))) |
     sum(class(object) %in% blockMethods) == 0){

# Create a new object for point simulation
    pointObject = object
# I am not sure why it at some stage made sense to recalculate variogram at this place
#    if ((ivar = which("variogramModel" == names(object)))>0) {
#      pointObject = pointObject[-ivar]
#      class(pointObject) = class(object)
#    }
#
#   Create a grid
    ngrid = params$ngrid
    if ("boundaries" %in% names(object)) {
      boundaries = object$boundaries
      pointObject$predictionLocations = findGrid(newdata = boundaries,n = ngrid,
          sampleSubregions=TRUE,ID="regCode",...)
    } else {
      bb = bbox(observations)
      if ("cellsize" %in% names(dots)) {
        pointObject$predictionLocations = findGrid(newdata = observations,
            dots$cellsize,sampleSubregions=FALSE,...)
      } else {
        pointObject$predictionLocations = findGrid(newdata = observations,ngrid,
            sampleSubregions=FALSE,...)
      }
    }

# Do we need point predictions or simulations?
# Prediction only if we dont need variance or exceedance probability
    if (!predictType$var & !predictType$exc) {
      pointObject = spatialPredict(pointObject)
    } else {
      pointObject = spatialPredict(pointObject,nsim=params$nSim)
    }
    if (yamamoto) {
      variogram = object$variogram
      predPoint = yamamotoKrige(formul,observations, pointObject$predictionLocations,variogram,...)
    } else {
      predPoint = krige(object$formulaString,observations, pointObject$predictionLocations,variogram,nmax = nmax)
    }
    print(class(predPoint))
    object$predPoint = predPoint    
    if (!is.na(threshCor)) {
      print("from aggregation")
      print(summary(object$predPoint))
      object = unbiasedKrige(object,...)
    }
    print(class(predPoint))
    object$predictionLocations = predictionLocations
    object = mergePredictions(object)
  }
  
  return(object)

}

#      predictType = data.frame(pred = TRUE, var=TRUE,exc=TRUE,threshCor=NA,yamamoto = FALSE,
#                  blockFat=TRUE,quantiles = NA),
  
  
mergePredictions = function(object) {
  predPoint = object$predPoint
  coor = coordinates(predPoint)
  predictions = predPoint@data
  predLocations = object@predictionLocations
  params = object$params
  thresh = params$thresh
  predType = params$predType
# MOK or IWQSEL are not relevant for simulations

  if ("sim" %in% names(predictions)) {
    if (predType$pred) pred = matrix(0,nrow=dim(predLocations@data)[1],ncol=dim(predPoint)[2])
    if (predType$blockFat) blockFat = matrix(0,nrow=dim(predLocations@data)[1],ncol=dim(predPoint)[2])
    if (predType$yamamoto) yamamoto = matrix(0,nrow=dim(predLocations@data)[1],ncol=dim(predPoint)[2])
     
    for (isim in 1:dim(predictions)[2]) {
      input = predictions[,isim]
      if (predType$pred) pred[,isim] = aggregatePredictions(input,coor,predLocations,"mean")
      if (predType$blockFat) blockFat[,isim] = aggregatePredictions(input,coor,predLocations,"fat")
      if (predType$yamamoto) yamamoto[,isim] = aggregatePredictions(input,coor,predLocations,"fat")
      
    }
    # If variance or exceedance probabilities are needed from the simulations:
    if (predType$var) predLocations$var1.var = apply(pred,MARGIN=1,FUN=var)
    if (predType$exc & !is.na(thresh)) predLocations$exc = apply(pred,MARGIN=1,thresh = thresh,FUN=function(arr,thresh) sum(I(arr>thresh))/length(arr))    
    if (predType$pred) predLocations$sim = pred
    if (predType$blockFat) predLocations$blockFat = blockFat
  }  
  object$predictions = predLocations
return(object)
}  

#  predictType = data.frame(pred = TRUE, var=TRUE,exc=TRUE,threshCor=NA,yamamoto = FALSE,
#                  blockFat=TRUE),
  
  
#  aggregatePredictions = function( input, coordinates, areaAggreg, action, threshold, percentile,...)

# Aggregates results from a vector of interpolation or simulation results: argument "input".
# "coordinates" is the corresponding locations of class Spatial Points of "input". 
# The area of aggregation "areaAggreg" is a spatial polygon data frame.
# 
# Depending on the action argument, the function computes max, min, mean, acdf, 
# probability above a threshold, and percentile (between 0 and 1).

# If action = "acdf", any argument for function hist can be included.
# If action = "threshold", argument threshold must be informed.
# If action = "percentile", argument percentile must be informed.

  

exProb = function(arr,pred,unc,thresh) {
if (!missing(arr)) {
  ep = 1-mean(I(arr<thresh))
} else {
  ep = 1-pnorm(thresh,pred,sqrt(unc))
}
return(ep)
}





findGrid = function(newdata, cellsize, n=100, sampleSubregions=FALSE, 
                    cellmin=0, sMin = 4, block, regCode = "regCode",...) {
#  3 options:
#  1) If sampleSubregions is false then sample the whole region of newdata
#  2) If sampleSubregions is false subregions denominated the string "regCode" 
#     are sampled
#  3) If block is given, 
#  If both n and cellsize is given, cellsize is used, unless the sampling
#  will be done in subregions, where cellsize does not make sense
#  If cellmin gives empty polygons, sMin gives a minimum number of samples from
#  each polygon   
  if (!missing(block)) {

#supx = c(-5, -2,  2,  5,  5,  2, -2, -5,-5)
#supy = c(-2, -5, -5, -2,  2,  5,  5,  2,-2)
#block = data.frame(x=supx,y=supy)
     
     
    if (!missing(block)) {
      if (is.null(dim(block))) {
        block = expand.grid(x=c(-block[1]/2,block[1]/2),y=c(-block[2]/2,block[2]/2)) 
        block = block[c(1,2,4,3),]
        block[5,] = block[1,]
      }
      aR = coordinates(spsample(Polygon(block),n,"regular",offset=c(0.5,0.5)))
      naR = dim(aR)[1]
      coord = coordinates(newdata)
      for (i in 1:dim(newdata)[1]) {
        coor = t(matrix(rep(coord[i,],naR),ncol=naR))
        aR = coor+aR
        idNew=as.character(rep(i,naR))
        predGridNew = SpatialPointsDataFrame(aR,data=data.frame(id = idNew))
        if (i == 1) {
          id = idNew
          predGrid = predGridNew
        } else {
          id = c(id,idNew)
          predGridNew = SpatialPointsDataFrame(aR,data=data.frame(id = idNew))                     
          predGrid = rbind(predGrid,predGridNew)
        }
      }
    } 
  } else if (!sampleSubregions) {
    if (!missing(cellsize)) return(spsample(newdata,type="regular",cellsize=cellsize))
    return(spsample(newdata,n,type="regular"))
  } else {  
    if (!missing(cellsize)) n = bbArea(bbox(newdata))/(cellsize*cellsize)

    for (i in 1:length(newdata@polygons)) {
      nl = length(newdata@polygons[i][[1]]@Polygons)
      cArea = 0
    # Summing up areas of all subpolygons in a polygon
      for (ip in 1:nl) cArea = cArea + newdata@polygons[i][[1]]@Polygons[[ip]]@area
      if (sqrt(cArea/n)>cellmin) {
        aR = spsample(newdata[1,],type="regular", n=n)
      } else {
        aR = spsample(newdata[1,],type="regular",n=max(sMin,sqrt(cArea)/cellmin))
      }
      idNew=as.character(rep(newdata@data[i,][[regCode]],length(coordinates(aR)[,1])))
      predGridNew = SpatialPointsDataFrame(aR,data=data.frame(id = idNew))
# 
      if (i == 1) {
        id = idNew
        predGrid = predGridNew
      } else {
        id = c(id,idNew)
        predGridNew = SpatialPointsDataFrame(aR,data=data.frame(id = idNew))                     
        predGrid = rbind(predGrid,predGridNew)
      }
    }
    return(predGrid)   
    }
  }

# Example of use

# countries = readShapePoly("...")

# f1=findGrid(region=countries, gridSize=NA, samplesSubregions=2000, cellmin=2000, ID="COUNTRY")
# plot(f1[,1:2],pch=".")

# f2=findGrid(region=countries, gridSize=5000, samplesSubregions=2000, cellmin=2000, ID="COUNTRY")
# plot(f2[,1:2],pch=".")
 
# source("c:/svn/intamap/R/zzz.r")
#   krigingObject = spatialPredict(krigingObject)  # Does not yet take anisotropy into account






## Author: O. Baume
## Date: 21-08-2008

aggregatePredictions = function(predictions, coordinates, areaAggreg, action, threshold, 
	percentile,...) {

# Aggregates results from a vector of interpolation or simulation results: argument "input".
# "coordinates" is the corresponding locations of class Spatial Points of "input". 
# The area of aggregation "areaAggreg" is a spatial polygon data frame.
# 
# Depending on the action argument, the function computes max, min, mean, acdf, 
# probability above a threshold, and percentile (between 0 and 1).

# If action = "acdf", any argument for function hist can be included.
# If action = "threshold", argument threshold must be informed.
# If action = "percentile", argument percentile must be informed.

# Action "fat" has not yet been implemented
  areaGrid = overlay( coordinates, areaAggreg)
  val=predictions[!is.na(predictions[areaGrid])]
  if (action=="max")
    {
    return( max(val))
    }
  if (action=="min")
    {
    return( min(val))
    }
  if (action=="mean")
    {
    return( mean(val))
    }
  if (action=="acdf")
    {
    h=hist(val,...)
    acdf=cumsum(h$counts)/sum(h$counts)
    return(data.frame(breaks=h$breaks[2:length(h$breaks)],acdf))
    } 
  if (action=="threshold")
    {
    if ((threshold>max(val)) | (threshold<min(val)))
      {
      return(NULL)
      } else {
      h=hist(val,breaks=c(min(val),threshold,max(val)))
      return(1-cumsum(h$counts)[1]/sum(h$counts))
      }
    }
  if (action=="percentile")
    {
    if ((percentile>1) | (percentile<0))
      {
      return(NULL)
      } else {
      return(sort(val)[round(percentile*length(val))])
      }
    }
}
