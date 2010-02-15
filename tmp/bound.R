
###################################
#
# Bound.R included different functions for finding boundaries between regions (countries)
#
###################################

findBoundaryLines = function(regions, projOrig, projNew) {
  regions = cleanRegions(regions)
#    writePolyShape(regions,"regionLim")
  if (!missing(projOrig)) {
    proj4string(regions) = CRS(projOrig)
    if (!missing(projNew) & projOrig != projNew) {
      regions = spTransform(regions, CRS(projNew))
    }
  }
  boundaryLines = findBoundaries(regions)
#    save(regionalBoundaries,file="regionalBoundaries.r")

  if (!extends(class(boundaryLines),"Spatial")) coordinates(boundaryLines) = ~x+y
  if (!missing(projNew) & (is.na(is.projected(boundaryLines)) | 
                          !is.projected(boundaryLines)))
     proj4string(boundaryLines) = CRS(projNew)

  return(boundaryLines)
}





findBoundaries = function(regions) {
# Function which takes a shape file of regional boundaries as input,
# and gives back a list of the points that define each single border
regCode = unique(regions$regCode)
nRegCode = length(regCode)
if (exists("cabound")) rm(cabound)
cList = list()
polyList = list()
for (ic in 1:nRegCode) {
  rci = as.character(regCode[ic])
  c1 = regions["regCode" %in% names(regions) & regions$regCode == rci,]
  c1p = spdf2list2(c1)
  c1kor = c1p$poly
  c1df = as.data.frame(c1kor)
  c1n = c1df[!is.na(c1df$V1),]
  c1k = unique(c1n)
  cList[[ic]] = c1k
  polyList[[ic]] = c1

}
for (i in 1:(nRegCode-1)) {
  rci = as.character(regCode[i])
  c1k = cList[[i]]
  c1 = polyList[[i]]
  for (j in (i+1):nRegCode) {
    rcj = as.character(regCode[j])
    c2k = cList[[j]]
    c2 = polyList[[j]]
#    cat(paste("Checking",i,j,rci,rcj,"\n"))
    if ( commonArea(c1,c2)[[1]] > 0.001) {
      c21 = rbind(c1k,c2k)
      lbound = c21[duplicated(c21),]
      ldim = dim(lbound)[1]
      if (ldim > 0) {
        caz = data.frame(c1 = rci,c2 = rcj, lbound)
        if (exists("cabound")) {
          cabound = rbind(cabound,caz)
        } else {
          cabound = caz
        }
      }
    }
  }
}
names(cabound) = c("c1","c2","x","y")
return (cabound)
}


