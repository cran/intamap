set.seed(15331)
library(intamap)
loadMeuse()
meuse$value = log(meuse$zinc)
meuse.grid = meuse.grid[sample(1:dim(meuse.grid)[1], 100),]
gridded(meuse.grid)=TRUE
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

summary(t(output2$outputTable))

output2$outputTable - output$outputTable


