sphervar = function(sill,nugg,rangev,xd) {
  if (xd < 0.000001) {
    0
  } else if (xd > rangev) {
    sill+nugg
  } else {                                
    fte = 1.5*sill*xd/rangev
    ste = 0.5*sill*(xd/rangev)^3
    nugg+fte-ste
  }
}

sYamKrige = function(newCor,cinv,Obs,val="value",sill,nugg,rangev,ikri){
# Function for predicting at one location
  if (!missing(ikri) && ikri %% 100 == 0) print(ikri)
  c0dist = spDistsN1(coordinates(Obs),newCor)
  c0arr = mapply(FUN = sphervar,xd = c0dist,MoreArgs = list(sill = sill,nugg = nugg, rangev = rangev))
  clen = length(c0arr)
  c0arr[(clen+1)] = 1

  ww = cinv %*% c0arr
  if (min(ww[1:clen]) < 0) {
    ww[1:clen] = ww[1:clen]-min(ww[1:clen])
    ww = ww/sum(ww[1:clen])
  }
  var1.pred = ww[1:clen] %*% as.matrix(Obs[val]@data)
  var1.ok = t(ww) %*% c0arr
  var1.var = t(ww[1:clen]) %*% ((as.matrix(Obs[val]@data)-as.numeric(var1.pred))^2)
  return(c(var1.pred=var1.pred,var1.ok=var1.ok,var1.var = var1.var))
}



yamamotoKrige = function(formula,Obs, newPoints,vmod,...) {
  depVar=as.character(formula[[2]])
  sill = vmod$psill[2]
  nugg = vmod$psill[1]
  rangev = vmod$range[2]
  cObs = coordinates(Obs)
  if (inherits(newPoints,"Spatial")) {
    cNew = coordinates(newPoints)
  } else {
    cNew = newPoints
  }
  if (is.null(dim(cNew))) cNew = as.data.frame(cNew)
  if (dim(cNew)[2] == 1) cNew = t(cNew)
  dvar = as.matrix(dist(cObs))
  cmat = apply(dvar,MARGIN=c(1,2),FUN = sphervar,sill = sill,nugg = nugg, rangev = rangev)
  dl = dim(cmat)[1]
  cmat = rbind(cmat,rep(1,dl))
  cmat = cbind(cmat,rep(1,dl+1) )
  dl = dl+1
  cmat[dl,dl] = 0
  cinv = solve(cmat)
  ikri = c(1:dl) 
  preds = t(apply(cNew,MARGIN = 1,FUN = sYamKrige,Obs = Obs, depVar=depVar,cinv = cinv,sill = sill, nugg = nugg, rangev = rangev,ikri=ikri))
  preds = as.data.frame(cbind(preds,cNew))
  coordinates(preds) = as.formula(paste("~",dimnames(cNew)[2][[1]][1],"+",dimnames(cNew)[2][[1]][2]))
#  print(preds)
  return(preds)
}



condSimYama = function(Obs,newPoints,isim=1,vmod,depVar="value",nmax = 25) {
  if (length(names(Obs)) == 1) depVar = names(Obs)
  sill = vmod$psill[2]
  nugg = vmod$psill[1]
  rangev = vmod$range[2]
  ccObs = coordinates(Obs)
  nObs = dim(ccObs)[1]
#  dataObs = matrix(nrow = nObs,ncol = ssim)
#  dataObs = t(mapply(Obs@data[[val]],FUN = function(X,ssim) c(rep(X,ssim)),MoreArgs = list(ssim)))
  sim = Obs@data[[depVar]]
  coords = coordinates(newPoints)
  nPred = dim(coords)[1]
  iord = sample(c(1:nPred))
#  nPred = 100

  cObs = matrix(ncol = 2,nrow=(nObs+nPred))
  cObs[1:nObs,] = ccObs
  iObs = nObs
  for (inew in 1:nPred) {
    if (iObs %% 100 == 0) cat(paste("Simulation number",isim,"Simulating point",iObs,"of",nPred,"\n"))
    cNew = coords[iord[inew],]
    c0dist = spDistsN1(cObs[1:iObs,],cNew)
    clen = length(c0dist)
    if (nmax < iObs) {
      jord = order(c0dist)
      c0dist = c0dist[jord[1:nmax]]
      rObs = cObs[jord[1:nmax],]
      dObs = sim[jord[1:nmax]]
    } else {
      rObs = cObs[1:nObs]
      dObs = sim
    }
    c0arr = mapply(FUN = sphervar,xd = c0dist,MoreArgs = list(sill = sill,nugg = nugg, rangev = rangev))
    rvar = as.matrix(dist(rObs))
    cmat = apply(rvar,MARGIN=c(1,2),FUN = sphervar,sill = sill,nugg = nugg, rangev = rangev)
    dl = dim(cmat)[1]
    clen = dl
    cmat = rbind(cmat,rep(1,dl))
    cmat = cbind(cmat,rep(1,dl+1) )
    dl = dl+1
    cmat[dl,dl] = 0
    cinv = solve(cmat)
    c0arr[dl] = 1

    ww = cinv %*% c0arr
    if (min(ww[1:clen]) < 0) {
      ww[1:clen] = ww[1:clen]-min(ww[1:clen])
      ww = ww/sum(ww[1:clen])
    }
    var1.pred = ww[1:clen] %*% as.matrix(dObs)
    pred.var = t(ww) %*% c0arr
    var1.var = t(ww[1:clen]) %*% ((as.matrix(dObs)-as.numeric(var1.pred))^2)
    iObs = iObs + 1
    cObs[iObs,] = cNew
#    cObs = rbind(cObs,cNew)
    sim[iObs] = rnorm(1,var1.pred,sqrt(var1.var))

  }

  sim = sim[(nObs+1):(nObs+nPred)]
  sim = sim[order(iord)]

  preds = as.data.frame(cbind(coords,sim = sim))
  coordinates(preds) = as.formula(paste("~",dimnames(as.data.frame(coords))[2][[1]][1],"+",dimnames(as.data.frame(coords))[2][[1]][2]))
#  print(preds)
  return(preds)
}

yamamotoSim = function(formulaString,Obs,newPoints,nsim,vmod,...) {
depVar=as.character(formulaString[[2]])
mSim = newPoints
for (i in 1:nsim) {
  cat(paste("Conditional simulation ",i,"\n"))
#  cSim = condSim(Obs,newPoints,isim=i,...)
  cSim = condSimYama(Obs,newPoints,isim=i,vmod = vmod,depVar=depVar,...)
  cat(paste("Conditional simulation ",i," done \n"))
  colname = paste("sim",i,sep="")
  names(cSim) = colname
  if (i == 1) {
    mSim = cSim
  } else {
    mSim[[colname]] = cSim[[colname]]
  }
}
return(mSim)
}
