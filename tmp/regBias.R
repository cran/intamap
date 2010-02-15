removeRegionalBias = function(object,regionalBias,depVar="value"){
  regCode = unique(object$regCode)
  for (ic in 1:length(regCode)){
    rci = regCode[ic]
    cBias = regionalBias$wls[regCode == rci]
    cVar =  regionalBias$wlsvar[regCode == rci]
    cNew = object[[depVar]][object$regCode == rci] - cBias
    object[[depVar]][object$regCode == rci] = cNew
  }
  return(object)
}



    

findRegionalBias = function(object,boundaryLines,formulaString=as.formula("value~1"),minKrige = 5) {
  regCode = unique(object$regCode)
  gdat = regionalVariograms(object,regCode,formulaString)
  cDiff = regionalDiff(object,regCode,gdat,boundaryLines,formulaString,minKrige=5)
  cDiff = unBias(cDiff,regCode)
  cBias = dSolve(cDiff)
  cBias = list(regionalBias = data.frame(regCode,cBias),cDiff=cDiff)
  return(cBias)
}


unBias = function(cDiff,regCode) {
  d = cDiff$d
  y = cDiff$y
  v = cDiff$v
#  
  d = rbind(d,0)
  cd = dim(d)[1]
  ino = which(regCode == "NO")
  iis = which(regCode == "IS")
  iuk = which(regCode == "UK" | regCode == "GB")
  if (length(iis) > 0) {
    d[cd,ino] = .5
    d[cd,iuk] = .5
    d[cd,iis]= -1
    y[cd] = 0
    v[cd] = max(v)
    cd = cd+1
    d = rbind(d,0)
  }
  ibe = which(regCode == "BE")
  ide = which(regCode == "DE")
  inl = which(regCode == "NL")
  ifr = which(regCode == "FR")
  iie = which(regCode == "IE")
  d[cd,ibe] = .25
  d[cd,ide] = .25
  d[cd,inl] = .25
  d[cd,ifr] = .25
  d[cd,iuk] = -0.5
  d[cd,iie] = -0.5
  y[cd] = 0
  v[cd] = max(v)

  icy = which(regCode == "CY")
  if (length(icy > 0)) {
    igr = which(regCode == "GR")
    cd = cd + 1
    d = rbind(d,0)
    d[cd,icy] = 1
    d[cd,igr] = -1
    y[cd] = 0
    v[cd] = max(v)
  }
  iru = which(regCode == "RU")
  if (length(iru) == 0) {
    cd = cd + 1
    d = rbind(d,0)
    ise = which(regCode == "SE")
    ifi = which(regCode == "FI")
    idk = which(regCode == "DK")
    iee = which(regCode == "EE")
    ilt = which(regCode == "LT")
    ilv = which(regCode == "LV")
    d[cd,ino] = 0.333
    d[cd,ise] = 0.333
    d[cd,ifi] = 0.334
    d[cd,iuk] = 0.2
    d[cd,ide] = 0.2
    d[cd,idk] = 0.2
    d[cd,iee] = 0.1
    d[cd,ilt] = 0.1
    d[cd,ilv] = 0.1
    d[cd,iis] = 0.1
    y[cd] = 0
    v[cd] = mean(v)
  }
  imt = which(regCode == "MT")
  if (length(imt) > 0) {
    ies = which(regCode == "ES")
    iit = which(regCode == "IT")
    cd = cd + 1
    d = rbind(d,0)
    d[cd,imt] = 1
    d[cd,ifr] = -0.333
    d[cd,ies] = -0.333
    d[cd,iit] = -0.334
    y[cd] = 0
    v[cd] = max(v)
  }
  ibg = which(regCode == "BG")
  if (length(ibg) == 0) {
    igr = which(regCode == "GR")
    cd = cd+1
    d = rbind(d,0)
    d[cd,igr] = 1
    d[cd,iit] = -1
    y[cd] = 0
    v[cd] = max(v)
  }
    

  cd = cd + 1
  d = rbind(d,0)
  d[cd,] = 1
  y[cd] = 0
  v[cd] = min(v)
  cDiff$d = d
  cDiff$y = y
  cDiff$v = v
  return(cDiff)
}


regionalDiff = function(object,regCode,gdat,boundaryLines,formulaString,minKrige=5) {
# Function to estimate the differences between all neighbouring countries.
# Function gives a kriging estimate on the border of each of the countries,
# and uses the difference between the estimates to define the bias.
ic = 0
nRegCode = length(regCode)
res  = data.frame(ic = c(1),i = c(1),j=c(1), c1 = c("A"),c2 = c("A"),a1=c(1),v1=c(1),a2=c(1),v2=c(1),v3=c(0),adiff = c(0),ldim = c(0))
d = matrix(0,ncol = nRegCode,nrow = 1)
y = array(c(0))
v = array(c(0))
for (i in 1:(nRegCode-1)) {
  rci = as.character(regCode[i])
  localDatai = object[object$regCode == rci,]
  ndat = dim(localDatai)[1]
  if (ndat > minKrige) imod = gdat[rci][[2]][[1]]
  for (j in (i+1):nRegCode) {
    rcj = as.character(regCode[j])
    boundaries = boundaryLines
    boundVec = (boundaries$c1 == rci & boundaries$c2 == rcj) |
                        (boundaries$c1 == rcj & boundaries$c2 == rci)
 #    cat(paste(i,j,rci,rcj,dim(lbound)[1],"\n"))
    if (sum(boundVec) > 0){
      lbound = boundaries[boundVec,]
#      cat(paste(rci,rcj,dim(lbound)[1],(dim(lbound)[1] > 0),"\n"))
      localDataj = object[object$regCode == rcj,]
      mdat = dim(localDataj)[1]
      if (mdat > minKrige) {
        jmod = gdat[rcj][[2]][[1]]
        bord = paste(rci,rcj,sep="")
#        cat(paste(rci,rcj,dim(lbound)[1],dim(lbound),"\n"))
        xlinNew = SpatialLines(list(Lines(list(Line(coordinates(lbound))),ID = bord)))
        proj4string(xlinNew) = CRS(proj4string(boundaries))
        cat(paste("interpolating border between ",rci,rcj,"\n"))
        stri = krige(formulaString, localDatai,xlinNew,imod)
        strj = krige(formulaString, localDataj,xlinNew,jmod)
      }
      a1 = stri$var1.pred
      v1 = sqrt(stri$var1.var)
      a2 = strj$var1.pred
      v2 = sqrt(strj$var1.var)
      v3 = sqrt(stri$var1.var+strj$var1.var)
      adiff = abs(stri$var1.pred-strj$var1.pred)
      ic = ic+1
      if (ic > 1) d = rbind(d,0)
      d[ic,i] = 1
      d[ic,j] = -1
      y[ic] = a1-a2
      v[ic] = v3
      res$ic = ic
      res$i = i
      res$j = j
      res$c1 = factor(res$c1,levels = c(levels(res$c1), rci))
      res$c2 = factor(res$c2,levels = c(levels(res$c2), rcj))
      res$c1 = rci
      res$c2 = rcj
      res$a1 = a1
      res$v1 = v1
      res$a2 = a2
      res$v2 = v2
      res$v3 = v3
      res$adiff = adiff
      res$ldim = dim(lbound)[1]
      if (ic == 1) {
        rest = res
      } else {
        rest = rbind(rest,res)
      }
    }
  }
}
#print(acdf)
return(list(d=d,y=y,v=v,rest=rest))
}

regionalVariograms = function(object,regCode,formulaString,minKrige = 5){
  nRegCode = length(regCode)
  ig = 0
  for (ic in 1:nRegCode) {
    rci = as.character(regCode[ic])
    localData = object[object$regCode == rci,]
    ndat = dim(localData)[1]
    if (ndat > minKrige) {
      regionalVariogram = autofitVariogram(formulaString,localData)$var_model
      if (ig == 0) {
        gdat = gstat(id = rci,formula = formulaString,model = regionalVariogram,data = localData)
        ig = 1
      } else {
        gdat = gstat(gdat,id = rci,formula = formulaString,model = regionalVariogram,data = localData)
      }
    }
  }
  return(gdat)
}
