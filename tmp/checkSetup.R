checkSetup = function(object, quiet = FALSE) {
	if (!quiet)
		cat("Checking object ... ")
# check PointData
	if (is.null(object$pointData))
		stop("no pointData provided")
	if (!is(object$pointData, "SpatialPointsDataFrame"))
		stop("pointData not of class SpatialPointsDataFrame")
	if (!"formulaString" %in% names(object)) 
	  stop("object does not contain a formulaString")
  depVar=as.character(object$formulaString[[2]])
	if (is.null(object$pointData[[depVar]]))
		stop("pointData has no attribute with name value")

# check PredictionLocations
	if (is.null(object$predictionLocations))
		stop("no predictionsLocations provided")
	if (!is(object$predictionLocations, "Spatial"))
		stop("predictionLocations should be of some Spatial class")

# Check if targetCRS was set
  if (!is.null(object$targetCRS)) {
# check if PointData CRS was set
	  if (is.na(proj4string(object$pointData)))
		  stop("can't reproject pointData when its CRS is not set")
# check if PredictionLocations CRS was set:
	  if (is.na(proj4string(object$predictionLocations)))
		  stop("can't reproject predictionLocations when CRS is not given")
  }
# check that the biases to be added is among the ones removed:
  if (!is.na(object$params$addBias) & 
                sum(object$params$addBias %in% object$params$removeBias) 
                                        < length(object$params$addBias))
    stop("Cannot add biases that have not been removed")

	if (!quiet)
		cat("OK\n")
	invisible(TRUE);
}



