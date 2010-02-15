

getNuts = function(zname = "NUTS2006_20M") {

# Just to test the encoding of the escape characters
#  cat(paste("DE: \251 EuroGeographics bez\374glich der Verwaltungsgrenzen \n \n"))
  sname = basename(zname)
  sname = sub("_","-",sname)
  path = dirname(zname)
  sList = list.files(path = path, pattern=sname)
  sname = paste(ifelse(path != ".",paste(path,"/",sep=""),""),sname,sep="")
  print(sList)
  print(length(sList))

  if (length(sList) > 0) {
    if ((length(grep("dbf",sList)) > 0 & length(grep("shx",sList)) > 0 &
         length(grep("shp",sList)) > 0 & length(grep("prj",sList)) > 0 )) {
#      cat(paste("Shape-files",zname,"already appears to exist. Please delete \n"))
#      cat(paste("one of the files if you want to reload the shape. \n"))
#      cat(paste("These data can also be accessed by:\n"))
#      cat(paste("readShapePoly(\"",sname,"\") \n",sep=""))
#      calling readShapePoly is a problem because of UK/GB
      #EJP: if (require(maptools)) nuts = readShapePoly(sname) else
        #warning("maptools not installed, not able to read boundaries")
	  #EJP changed into:
	  nuts = readOGR(".", sname)
      return(nuts)
    }
  }

  cat(paste("You are now downloading a file from Eurostat \n"))
  cat(paste("http://epp.eurostat.ec.europa.eu/portal/page?_pageid=2254,64099847,2254_64185160&_dad=portal&_schema=PORTAL \n \n"))

  cat(paste("In addition to the general copyright and licence policy \n"))
  cat(paste("applicable to the whole Eurostat web site (available from \n"))
  cat(paste("http://epp.eurostat.ec.europa.eu/pls/portal/url/page/PGP_DS_LICENCE) \n"))
  cat(paste("the following specific provisions apply to the datasets you \n"))
  cat(paste("are downloading. The download and usage of these data is \n"))
  cat(paste("subject to the acceptance of the following clauses: \n \n"))

  cat(paste("1. The Commission agrees to grant the non-exclusive and not \n"))
  cat(paste("   transferable right to use and process the Eurostat/GISCO \n"))
  cat(paste("   geographical data downloaded from this page (the \"data\"). \n \n"))

  cat(paste("2. The permission to use the data is granted on condition that: \n"))
  cat(paste("    a) the data will not be used for commercial purposes; \n"))
  cat(paste("    b) the source will be acknowledged. A copyright notice, \n"))
  cat(paste("       as specified below, will have to be visible on any \n"))
  cat(paste("       printed or electronic publication using the data \n"))
  cat(paste("       downloaded from this page. \n \n"))

  cat(paste("             Copyright notice \n \n"))
  cat(paste("When data downloaded from this page is used in any printed or \n"))
  cat(paste("electronic publication, in addition to any other provisions \n"))
  cat(paste("applicable to the whole Eurostat web site, data source will \n"))
  cat(paste("have to be acknowledged in the legend of the map and in the \n"))
  cat(paste("introductory page of the publication with the following \n"))
  cat(paste("copyright notice: \n \n"))

  cat(paste("EN: \251 EuroGeographics for the administrative boundaries \n"))
  cat(paste("FR: \251 EuroGeographics pour les limites administratives \n"))
  cat(paste("DE: \251 EuroGeographics bez\374glich der Verwaltungsgrenzen \n \n"))

  cat(paste("For publications in languages other than English, French or \n"))
  cat(paste("German, the translation of the copyright notice in the language \n"))
  cat(paste("of the publication shall be used. \n \n"))

  cat(paste("Do you accept these conditions? (Y/N) \n"))
  ans = readLines(con=stdin(),n = 1)
  if (ans == "Y" | ans == "y") {
    zfile = paste(zname,".zip",sep="")
    download.file(paste("http://epp.eurostat.ec.europa.eu/pls/portal/docs/PAGE/PGP_DS_GISCO/PGE_DS_GISCO/TAB_GEO/",
             zname,"_shapefile.ZIP",sep=""),zfile,mode="wb")

    sfile = paste(sname,".dbf",sep="")
    zip.file.extract.dir(sfile,zipdir = getwd(),zipname=zfile)
    sfile = paste(sname,".prj",sep="")
    zip.file.extract.dir(sfile,zipdir = getwd(),zipname=zfile)
    sfile = paste(sname,".shp",sep="")
    zip.file.extract.dir(sfile,zipdir = getwd(),zipname=zfile)
    sfile = paste(sname,".shx",sep="")
    zip.file.extract.dir(sfile,zipdir = getwd(),zipname=zfile)
	#EJP:
    #if (require(maptools)) nuts = readShapePoly(sname) else
    #    warning("maptools not installed, not able to read boundaries")
	nuts = readOGR(".", sname)
    cat(paste("\nThe data has been downloaded and stored in the directory \n"))
    cat(paste(getwd(),"\n"))
#    cat(paste("Next time you can access these data calling:\n"))
#    cat(paste("readShapePoly(\"",sname,"\") \n",sep=""))
    return(nuts)
  } else {
    cat(paste("You did not agree to the conditions. Returning without data. \n"))
  }
}


nuts2Boundaries = function(nuts, zname) {
#if (!require(maptools)) error("maptools not installed, cannot read shape-file")
	if (missing(nuts)) nuts = getNuts(zname)
	nuts = nuts[nchar(as.character(nuts$NUTS_ID)) == 2,]

	iCountry = which(names(nuts)=="CNTR_CODE")
	names(nuts)[iCountry] = "regCode"
	if ("UK" %in% nuts$regCode) {
  		levels(nuts$regCode) = c(levels(nuts$regCode),"GB")
  		nuts$regCode[nuts$regCode == "UK"] = "GB"
	}
	return(nuts)
}
