set.seed(15331)
library(intamap)
#loadMeuse()

crs = CRS("+init=epsg:28992 +proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.417,50.3319,465.552,-0.398957,0.343988,-1.8774,4.0725 +units=m +no_defs")
data("meuse")
coordinates(meuse) <- ~x+y
proj4string(meuse) <- crs
data("meuse.grid")
coordinates(meuse.grid) <- ~x+y
gridded(meuse.grid) <- TRUE
proj4string(meuse.grid) <- crs

meuse$value = log(meuse$zinc)
meuse.grid = meuse.grid[sample(1:dim(meuse.grid)[1], 100),]
output = interpolate(meuse, meuse.grid, list(mean=T, variance=T, nsim = 5),methodName = "automap")
summary(t(output$outputTable))

output = interpolate(meuse, meuse.grid,
    optList = list(idpRange = seq(0.1, 2.9, 0.5), nfold = 3), 
    outputWhat = list(mean=TRUE),methodName = "idw")
summary(t(output$outputTable))


output = interpolate(meuse, meuse.grid, list(mean=T, variance=T),methodName = "transGaussian")
summary(t(output$outputTable))


data(meuse)
meuse = meuse[sample(dim(meuse)[1],30),]
meuse$value=meuse$zinc
coordinates(meuse) = ~x+y
mgrid = spsample(meuse,25,"regular")
gridded(mgrid) = TRUE
output = interpolate(meuse, mgrid, list(mean=T, variance=T, excprob = 1000,quantile = 0.5),
                     methodName = "copula")

output2 = interpolate(meuse, mgrid, list(mean=T, variance=T, excprob = 1000,quantile = 0.5),
                     methodName = "copula",optList = list(methodParameters = output$methodParameters))

output3 = interpolate(meuse, mgrid, list(mean=T, variance=T, excprob = 1000,quantile = 0.5),
                     methodName = "copula",optList = list(methodParameters = output$methodParameters, nclus = 5))

output4 = interpolate(meuse, mgrid, list(mean=T, variance=T, excprob = 1000,quantile = 0.5),
                     methodName = "automap",optList = list(methodParameters = output$methodParameters, nclus = 5))

output5 = interpolate(meuse, mgrid, list(mean=T, variance=T, excprob = 1000,quantile = 0.5),
                     methodName = "automap",optList = list(model = "Sph", nclus = 5))

summary(t(output2$outputTable))

output2$outputTable - output$outputTable
output3$outputTable - output$outputTable
output4$outputTable - output$outputTable
output5$outputTable - output$outputTable


