  
blockPredict = function(object,...) {
  object = spatialPredict.block(object,...)
  return(object)
} 
 

spatialPredict.block = function(object,...) {
#spablock = function(object,...) {
# Define which methods can actually handle block kriging directly
  blockMethods = c("automap","idw")
# Extract some of the parameters and variables from the object
# - for easier access
  params = object$params
  blockWhat = object$blockWhat
  dots = list(...)
  block = params$block
  if ("block" %in% names(dots)) block = dots$block
  nmax = params$nmax

  observations = object$observations
  formul = object$formulaString
  predictionLocations = object$predictionLocations
  outputWhat = object$outputWhat
  if (!is(predictionLocations,"SpatialPolygons") &
      !is(predictionLocations,"SpatialLinesDataFrame") &
      !length(block) >0 ) {
      warning("prediction locations are not blocks, polygons or lines")
      warning("calling point prediction")
      return(spatialPredict(object,...))
  }
  
  if (!("variogramModel" %in% names(object)) &
      !inherits(object,"idw") & !inherits(object,"copula")) object = estimateParameters(object)
  if (params$processType %in% c("gaussian","logNormal") &
      sum(class(object) %in% blockMethods) >= 1) {
# Check if we might be able to use an analytical solution - 
# Possible if we want a prediction for Gaussian or logNormal distribution, only 1 cluster
# Only possible for certain methods - now hard coded into blockMethods
    if (inherits(object,"idw")) return(spatialPredict(object,...))
    form = object$formulastring
    depVar = as.character(form[2])
# We have to discuss whether it makes sense to have the test of logNormal here, 
# as it can easily create some inconsistencies.
# We have to be sure that the variogram is also found from the logarithmised data
    if (params$processType == "logNormal") observations[[depVar]] = log(observations[[depVar]])
      predictions = krige(object$formulaString,observations, predictionLocations,object$variogramModel,block = block,nmax = nmax)
    if (params$processType == "logNormal") {
      predictions$var1.pred = exp(predictions$var1.pred + predictions$var1.var/2) 
      warning("Note that the back-calculated lognormal predictor is only approximate, as the Lagrange parameter is not included")
    }
    print("performed ordinary block kriging")
    object$predictions = predictions
  }
  
# Point predictions/simulations are necessary
    if ("MOK" %in% names(object$outputWhat) | "IWQSEL" %in% names(object$outputWhat)
      | inherits(object,"yamamoto") | length(names(blockWhat)) > 0
      | !inherits(object,blockMethods)) {
# Create a new object for point simulation
    pointObject = object
#   Create a grid
    ngrid = params$ngrid
    if (is(predictionLocations,"SpatialPolygons")) {
      pointObject$predictionLocations = findGrid(predictionLocations = predictionLocations,
          n = ngrid, sampleSubregions=TRUE,...)
    } else if (length(block) >0 ) {
      pointObject$predictionLocations = findGrid(predictionLocations = predictionLocations,
          n = ngrid, sampleSubregions=TRUE, block = block,...)
    } else {
      bb = bbox(observations)
      if ("cellsize" %in% names(dots)) {
        pointObject$predictionLocations = findGrid(predictionLocations = observations,
            cellsize = dots$cellsize, sampleSubregions=FALSE,...)
      } else {
        pointObject$predictionLocations = findGrid(predictionLocations = observations,
          n = ngrid, sampleSubregions=FALSE,...)
      }
    }

# Do we need point predictions or simulations?
# Prediction only if only need mean
    
    nsim = ifelse(all(names(outputWhat)=="mean") & blockWhat == "none",0, params$nsim) 
    vmod = object$variogramModel
    nmax = ifelse(params$nmax == Inf & dim(coordinates(pointObject$predictionLocations))[1] > 200,20,params$nmax)
    pointObject = spatialPredict(pointObject,nsim=params$nsim,nmax = nmax)
    pointObject$predictions = pointObject$predictions[,grep("sim",names(pointObject$predictions))]
    object$pointPredictions = pointObject$predictions
    object$pointLocations = pointObject$predictionLocations
    object = spatialAggregate(object)
    if (object$params$debug.level < 2) object$pointPredictions = 
                  "pointPredictions deleted from object, debug.level < 2"
  }
  return(object)

}


  
spatialAggregate = function(object) {
  predictionLocations = object$predictionLocations
  pointLocations = object$pointLocations  
  ids = pointLocations$id
  id = unique(ids)
  if ("predictions"%in% names(object)) {
    predictions = object$predictions
  } else {
    predictions = predictionLocations
    if (!length(grep("DataFrame",class(predictions)))>0) 
      predictions = SpatialDataFrame(predictions,data=data.frame(id = id))
  }
  pointPredictions = object$pointPredictions
  coor = SpatialPoints(pointPredictions)
  pointPred = pointPredictions@data
  params = object$params
  thresh = params$thresh
  outputWhat = object$outputWhat
  blockWhat = object$blockWhat
  sims = pointPred[,grep("sim",names(pointPred))>0]

  predAggr = aggregate(sims,by=list(id=pointLocations$id),mean)
  predictions@data = data.frame(predictions@data,predAggr)
  
  if (length(blockWhat) > 0 && blockWhat != "none") {
    for (ib in 1:length(blockWhat)) {
      what = blockWhat[ib]
      if (names(what) == "fat") {
        thresh = what[[1]]
        fatf = function(arr,thresh) sum(I(arr>thresh))/length(arr)  
        fatx = aggregate(sims,by=list(id = ids),FUN = fatf,thresh=thresh)
        vmean = rowMeans(fatx[,-1])
        vvar = apply(fatx[,-1],MARGIN=1,FUN=function(arr) var(arr))
        vname = paste("fat",what[[1]],sep="")
        vnamevar = paste("fatVar",what[[1]],sep="")
      } else if (names(what) == "blockMax" && what[[1]]) {
        bmax = aggregate(sims,by=list(id=ids),FUN = max)
        vmean = rowMeans(bmax[,-1])
        vvar = apply(bmax[,-1],MARGIN=1,FUN=function(arr) var(arr))
        vname = "blockMax"
        vnamevar = "blockMaxVar"
      } else if (names(what) == "blockMin" && what[[1]]) {
        bmin = aggregate(sims,by=list(id=ids),FUN = min)
        vmean = rowMeans(bmin[,-1])
        vvar = apply(bmin[,-1],MARGIN=1,FUN=function(arr) var(arr))
        vname = "blockMin"
        vnamevar = "blockMinVar"
      }
      predictions@data[vname] = vmean        
      predictions@data[vnamevar] = vvar             
    }
  }
  object$predictions = predictions
  object
}  

  
  








findGrid = function(predictionLocations, cellsize, n=100, sampleSubregions=FALSE, 
                    block, sMin = 4, ...) {
#  newdata needs to be SpatialPointsDataFrame.
#  3 options:
#  1) If sampleSubregions is false then sample the whole region of newdata
#  2) If sampleSubregions is false subregions denominated the string "regCode" 
#     are sampled
#  3) If block is given, 
#  If both n and cellsize is given, cellsize is used, unless the sampling
#  will be done in subregions, where cellsize does not make sense
#  If cellmin gives empty polygons, sMin gives a minimum number of samples from
#  each polygon   
if (is(predictionLocations,"SpatialPolygons")) {
  if (sampleSubregions) {
    if (!missing(cellsize)) n = bbArea(bbox(predictionLocations))/(cellsize*cellsize)
    if (n < sMin) n = sMin
    ids = sapply(slot(predictionLocations, "polygons"), function(i) slot(i, "ID"))
    for (i in 1:length(predictionLocations@polygons)) {
      ldata = predictionLocations@polygons[i][[1]]
      nl = length(ldata@Polygons)
      cArea = 0
    # Summing up areas of all subpolygons in a polygon
      for (ip in 1:nl) cArea = cArea + ldata@Polygons[[ip]]@area
      aR = spsample(ldata,type="regular", n=n)
      naR = length(coordinates(aR)[,1])      
      if (is(predictionLocations,"SpatialPolygonsDataFrame")) {
        idl = ids[i]
        idNew = as.character(rep(idl,naR))
      } else idNew = as.character(rep(i,naR))
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
  } else {  
    if (!missing(cellsize)) predLoc = spsample(predictionLocations,type="regular",cellsize=cellsize) else 
      predGrid = spsample(predictionLocations,n,type="regular")
    predGrid = SpatialPointsDataFrame(predGrid,data = data.frame(id = overlay(predGrid,predictionLocations)))
  }
} else if (!missing(block)) {     
  if (is.null(dim(block))) {
		if (length(block) == 1)	block = rep(block,2)
    block = expand.grid(x=c(-block[1]/2,block[1]/2),y=c(-block[2]/2,block[2]/2)) 
    block = block[c(1,2,4,3),]
    block[5,] = block[1,]
  }
  aR = coordinates(spsample(Polygon(block),n,"regular",offset=c(0.5,0.5)))
  naR = dim(aR)[1]
  coord = coordinates(predictionLocations)
  for (i in 1:dim(coordinates(predictionLocations))[1]) {
    coor = t(matrix(rep(coord[i,],naR),ncol=naR))
    iaR = coor+aR
    idNew=as.character(rep(i,naR))
    predGridNew = SpatialPointsDataFrame(iaR,data=data.frame(id = idNew))
    if (i == 1) {
      id = idNew
      predGrid = predGridNew
    } else {
      id = c(id,idNew)
      predGridNew = SpatialPointsDataFrame(iaR,data=data.frame(id = idNew))                     
      predGrid = rbind(predGrid,predGridNew)
    } 
  }
} else stop ("Polygons or block size not given")
if (!is.na(proj4string(predictionLocations))) proj4string(predGrid) = proj4string(predictionLocations) 
predGrid
}

