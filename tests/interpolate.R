set.seed(15331)
library(intamap)
data(meuse)
meuse$value=meuse$zinc
coordinates(meuse) = ~x+y
data(meuse.grid)
gridded(meuse.grid)=~x+y
output = interpolate(meuse, meuse.grid, list(mean=T, variance=T, nsim = 5),methodName = "automap")
summary(t(output$outputTable))


library(intamap)
data(meuse)
meuse$value=meuse$zinc
coordinates(meuse) = ~x+y
data(meuse.grid)
gridded(meuse.grid)=~x+y
output = interpolate(meuse, meuse.grid, list(mean=TRUE),methodName = "idw")
summary(t(output$outputTable))


library(intamap)
data(meuse)
data(meuse.grid)
coordinates(meuse) = ~x+y
coordinates(meuse.grid) = ~x+y
meuse$value=meuse$zinc
output = interpolate(meuse, meuse.grid, list(mean=T, variance=T),methodName = "transGaussian")
summary(t(output$outputTable))



library(intamap)
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

summary(t(output2$outputTable))

output2$outputTable - output$outputTable


