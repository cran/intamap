commonArea = function(objecti,objectj) {
  bi = bbox(objecti)
  bj = bbox(objectj)
  iarea = bbArea(bi)
  jarea = bbArea(bj)
  bl = list()
  for (i in 1:2) bl[[i]] =  max(bi[[i]],bj[[i]])
  for (i in 3:4) bl[[i]] =  min(bi[[i]],bj[[i]])
  if (bl[[3]] >= bl[[1]] & bl[[4]] >= bl[[2]]) {
    larea = bbArea(bl)
    if (larea < 0.0001*max(iarea,jarea)) larea = 0.0001*max(iarea,jarea)
  } else {
    larea = 0
  }
#
  if (jarea< (0.0001*iarea)) jarea = 0.0001*iarea
  if (iarea< (0.0001*jarea)) iarea = 0.0001*jarea
  ilarea = larea/iarea
  jlarea = larea/jarea

#  print(bbox(objecti))
#  print(bbox(objectj))
#  print(bl)
#  cat(paste(iarea, jarea, larea, ilarea,jlarea,"\n"))

  return(list(ilarea,jlarea))
}


bbArea = function(bb) {
  xd = bb[[3]]-bb[[1]]
  yd = bb[[4]]-bb[[2]]
  xdd = xd*xd
  ydd = yd*yd
  hdd = xdd+ydd
  bbArea = sqrt(hdd)
}








 spdf2list2 = function (data) {
# Modified function to convert a SpatialPolygonsDataFrame into a list of polygons
# This is a new version, also checking if there are several polygons in each Polygons
# Original function taken from GeoXp
    if (class(data)[1] == "SpatialPolygonsDataFrame") {
        poly <- data@polygons
        n <- length(poly)
        for (i in 1:n) {
          np = length(poly[[i]]@Polygons)
          for (j in 1:np) {
            if (i ==1 & j ==1) {
              X <- poly[1][[1]]@Polygons[[1]]@labpt[1]
              Y <- poly[1][[1]]@Polygons[[1]]@labpt[2]
              contours <- rbind(NA, NA, NA, poly[1][[1]]@Polygons[[1]]@coords)
            } else {
              X <- rbind(X, poly[i][[1]]@Polygons[[j]]@labpt[1])
              Y <- rbind(Y, poly[i][[1]]@Polygons[[j]]@labpt[2])
              contours = rbind(contours, NA, NA, NA, poly[i][[1]]@Polygons[[j]]@coords)
            }
          }
        }
        contours = rbind(contours, NA, NA, NA)
        return(list(X = X, Y = Y, poly = contours))
    }
    else {
        tkmessageBox(message = "Must be a SpatialPolygonsDataFrame object",
            icon = "warning", type = "ok")
    }
}




dSolve = function(diffs) {
  d = as.data.frame(diffs$d)
  xnam <- paste("x", 1:dim(d)[2], sep="")
  names(d) = xnam
  fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"),"-1"))
  y = diffs$y
  v = diffs$v
  dat = cbind(d,y=y)
  ols = as.data.frame(summary(lm(fmla,dat))$coefficients)
  wls = as.data.frame(summary(lm(fmla,dat,weights = 1/v))$coefficients)
  names(ols) = c("ols","ols.std","ols t value", "ols.p")
  names(wls) = c("wls","wls.std","wls t value", "wls.p")
  return(cbind(ols,wls))
}

zip.file.extract.dir <- function (file, tmpdir, zipname = "R.zip", unzip = getOption("unzip"),zipdir)
{
    if (!is.character(unzip) || length(unzip) != 1)
        stop("'unzip' must be a single character string")
    if (!nzchar(unzip))
        unzip <- "internal"
    path <- dirname(file)
    topic <- basename(file)
    if (file.exists(file.path(path, zipname))) {
#        if (missing(zipdir)) {
#          tmpd <- tempdir()
#        } else {
#          tmpd <- zipdir
#        }
        tmpd = ifelse(missing(zipdir),tempdir(),zipdir)
        if (unzip != "internal") {
            cmd <- paste(unzip, "-oq", shQuote(file.path(path,
                zipname)), topic, " -d ", tmpd)
            res <- if (.Platform$OS.type == "windows")
                system(cmd, invisible = TRUE)
            else system(paste(cmd, "> /dev/null"))
            if (!res)
                file <- file.path(tmpd, topic)
        }
        else {
            rc <- .Internal(int.unzip(file.path(path, zipname),
                topic, tmpd))
            if (rc == 0)
                file <- file.path(tmpd, topic)
        }
    }
    file
}
