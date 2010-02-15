#######################################################################
#Function	:EstimateAnisotropy function 								
#Description	:Receives a KrigingObject calls estimateAnisotropySc 	
#				function and returns a new KrigingObject which 			
#				includes anisotropy parameters estimation.													
#																		
#Input		:object		& The original Kriging Object				
#Output		:object		&including 									
#																		
#			&object$anisPar$(list(ratio=R,direction=theta.deg))		
#																		
#																		
#######################################################################
estimateAnisotropy<-function(object){
	
	# params = object$params
 	if ("formulaString" %in% names(object)) formulaString = object$formulaString else formulaString = as.formula("value ~ 1")	
	depVar=as.character(formulaString[[2]])

	xy<-as.matrix(coordinates(object$pointData))
	
	object$anisPar<-estimateAnisotropySc(xy[,1],xy[,2],object$pointData[[depVar]],method="linear",pl=FALSE)
	return(object)
}

#######################################################################
#Function	:rotateAnisotropicData function 								
#Description	:Receives an Intamap Object containing at least a node pointData
#		 containing data coordinates and a anisPar node containing 
#		 the anisotropy parameters that will be used for the data rotation. 		 	
#Input		:object		& The original Kriging Object.				
#Output		:object		&including rotated data.											
#				
#																		
#																		
#######################################################################
rotateAnisotropicData<-function(object){
	if ("formulaString" %in% names(object)) formulaString = object$formulaString else formulaString = as.formula("value ~ 1")	
	depVar=as.character(formulaString[[2]])
	
	xy<-as.matrix(coordinates(object$pointData))
	x<-as.matrix(xy[,1])
	y<-as.matrix(xy[,2])
	
	theta=(object$anisPar$direction)*pi/180
	R=object$anisPar$ratio

	#rotate data
	x_new=(x)*cos(theta)+y*sin(theta)
	y_new=R*(-x*sin(theta)+y*cos(theta))

	object$pointData = SpatialPointsDataFrame(cbind(x_new,y_new),data=object$pointData@data)
	
	return(object)
}




######################### Internal functions ##########################
#
#######################################################################
#Function	estimateAnisotropySc										
#													
#Description : 	This is the main fucntion used in order to calculate 
#		anisotropy parameters from scattered data
#		
#		Given a Cartesian coordinates system  of axes x and y and an ellipsoid	
#		of correlation isolevels has principal directions  a and b. The anisotropy 
#		rotation angle(theta.deg) express the  angle between principal direction a and axes x.
#			
#												
#Inputs		&x			&The x-axis coordinate of field			
# 			&y			&The y-axis coordinate of field			
#			&r			&The field value in (xi,yi) point		
#			&len		&The length of the interpolated field		
#			&method		&The method that we will use for interpolation 
#			&			&"cubic", "linear", "v4"		
#			&min.x		&minimum x value 	
#			&max.x		&maximum x value		
#			&min.y		&minimum y value		
#			&max.y		&maximum x value		
#			&deb		&toggles debugging mode	
#			&pl			&toggles plot			
#			&br			&borders coordinates		
#Outputs	&R			&Anisotropy ratio		
#			&theta.deg	&Anisotropy rotation angle in degrees 
#																	
#Packages : akima package 									 
#													
#######################################################################
estimateAnisotropySc<-function(x, y, r, len=length(x), method="linear", min.x=min(x), max.x=max(x), min.y=min(y), max.y=max(y),deb=FALSE,pl=FALSE,br){
            
	plot.borders=TRUE
#	if(missing(deb)==TRUE){ deb=FALSE}
#        if(missing(min.x)==TRUE){ min.x=min(x)}
#        if(missing(min.y)==TRUE){ min.y=min(y)}
#        if(missing(max.x)==TRUE){ max.x=max(x)}
#        if(missing(max.y)==TRUE){ max.y=max(y)}
#        if(missing(method)==TRUE){method="linear"}
#        if(missing(len)==TRUE)   {len=length(x)}
#	if(missing(pl))	{pl=FALSE}
	if(missing(br))	{plot.borders=FALSE}
            

	#mesh creation
	xn<-seq(min.x,max.x,length=round(sqrt(len)))
	yn<-seq(min.y,max.y,length=round(sqrt(len)))
        mesh<-meshgrid(xn,yn)

        #selection of interpolation method
        if (method=="cubic")   {  ri<-interp(x,y,r,xn,yn,linear=FALSE,extrap=TRUE,duplicate="mean")
					}	#akima package
        if (method=="linear")  {  ri<-interp(x,y,r,xn,yn,linear=TRUE,duplicate="mean")	
				 ri$z<-t(ri$z)	}				#akima package
        if  (method=="v4")     {  ri<-biharmonicSplineAnisotropy(x,y,r,mesh$x,mesh$y)}     
       		
	#calculate anisotropy parameters over regular grid 
	res=estimateAnisotropyGrid(mesh$x,mesh$y,ri$z)
	
		
	if((pl==TRUE & plot.borders==TRUE)){
	        image(mesh$x[1,],mesh$y[,1],t(ri$z),col=rainbow(20), #color.palette=rainbow,
	        plot.title=title(main=method),xlim=range(br[,1]),ylim=range(br[,2]))
		points(br,pch=".",xlim=(br[,1]),ylim=range(br[,2]))
	}else if ((pl==T & plot.borders==FALSE)){
		 image(mesh$x[1,],mesh$y[,1],t(ri$z),col=rainbow(20), #color.palette=rainbow,
	        plot.title=title(main=method))

	}


	dump<-jpdf(seq(1,3,by=0.01),seq(-90,90,by=1),res$R,res$theta.deg,len)

	return(list(ratio=res$R, direction=res$theta.deg,Q=res$Q,doRotation=dump$doRotation))
}
#######################################################################                              
# function [R , theta_deg ]= estimateAnisotropyGrid(xi ,yi , ri)  		
#     theta_deg , R :variables to b calculated              			
#     xi ,yi    : Coordinates on Grid                        			
#         ri    : Random field at (xi,yi)                    			
#######################################################################

estimateAnisotropyGrid<-function(xi,yi,ri){

	#Gradient function
	mgradient<-function(k,stx,sty){

    		gg<-function(m,st){

     		m<-as.matrix(m)
     		n=dim(m)[1]                     #rows
   	  	p=dim(m)[2]                     #columns
	
  	    	g<-matrix(0,n,p)
 	     	g[1,]<-(m[2,]- m[1,])/st        
 	     	g[n,]<-(m[n,]-m[n-1,])/st       
	      	g[2:(n-1),]=(m[3:n,]- m[1:(n-2),])/(2*st)                        
	      	return (g)
	    	}
	

	    y=gg(k,sty)                       
	    x=t(gg(t(k),stx))                 
	    return(list(x=x,y=y))             
	}




  stepx<-xi[1,2]-xi[1,1]
  stepy<-yi[2,1]-yi[1,1]
  divZ<-mgradient(ri,stepx,stepy)

  #get rid of Na
  ##############################
 # divZ$x[which(is.na(divZ$x))]=0
 # divZ$y[which(is.na(divZ$y))]=0
  ##############################

  divX2<-as.vector(divZ$x*divZ$x)
  divY2<-as.vector(divZ$y*divZ$y)
  divXY<-as.vector(divZ$x*divZ$y)

  Q11<-mean(divX2[is.na(divX2)==FALSE])
  Q22<-mean(divY2[is.na(divY2)==FALSE])
  Q12<-mean(divXY[is.na(divXY)==FALSE])
	
	if(Q11<10^-31){
	 return (list(R=1,theta.deg=0,Q=cbind(Q11,Q22,Q12)))
	}else{

  Zdiag<-Q22/Q11
  Zoff<-Q12/Q11
  theta<-0.5*atan(2*Zoff/(1-Zdiag))

   R<-sqrt(1+((1-Zdiag)/(Zdiag-(1+Zdiag)*(cos(theta))^2)))
   #Make sure that all results come in the same interval
	if (R<1){
		R=1/R
	
		if ((theta+pi/2>-pi/2) & (theta+pi/2)<(pi/2)){
			theta=theta+pi/2
		}else{
			theta=theta-pi/2
		}
		

	}
   		theta.deg<-theta*180/pi
	}
  return (list(R=R,theta.deg=theta.deg,Q=cbind(Q11,Q22,Q12)))

}
########################################################
#Function   : jpdf function 
#
#Description: calculates the joint pdf function for R
#		and theta over a given region 
#Arguments  : 
#	R 	: a mesh type matrix of the R region 	
#	theta 	: a mesh type matrix of the theta region
#	Rest	: the given ratio value
#	thetaEst: the given theta value
#######################################################
jpdf<-function(Rf,thetaf,Rest,thetaEst,N){
	
	
	erfc <- function(x){ 2 * pnorm(x * sqrt(2), lower = FALSE)}
		
	RRest=Rest
	tthetaEst=thetaEst
	
	if(N>1200){N=1200}

	Rest=1
	thetaEst=0
	
	
	#create mesh grids
	mesh<-meshgrid(Rf,thetaf)
	R<-mesh$x
	theta<-mesh$y

	#threshold value to decide if the given set is in 95% interval of Jpdf
	threshold=2.79548^2



	#first convert degrees to Radians
	theta=theta*pi/180
	thetaEst=thetaEst*pi/180

	
	#Covariance matrix elements -- constants 
	s1=2*(cos(thetaEst)^2 +Rest^2 * sin(thetaEst)^2 )
	s2=2*(sin(thetaEst)^2 +Rest^2 * cos(thetaEst)^2 )
	s12=2*sin(thetaEst)* cos(thetaEst)*(1-Rest^2)


	#Determinant
	d=s1*s2-s12^2
	
	qd = (R^2 +tan(theta)^2) /(1 + R^2 *tan(theta)^2)
	qo = tan(theta) * (1-R^2)/ (1+R^2 *tan(theta)^2)

	#Jacobian (qd,qo) ->(R,theta)
	J1= ((2*R*abs(R^2 -1)*(1/cos(theta))^6))
	J2=(1+R^2 *tan(theta)^2)^3
	J=J1/J2

	a = N/(2*d^2) * (qd^2 * s1^2  + 2 * qd * s12^2 + s2^2 - 4 * qo * s12 * (qd * s1 + s2) - 2 * qo^2 * (d - 2 * s1*s2))	
	b = -(N/d) * (qd * s1 - 2 * qo * s12 + s2);
	k = (N/d)^(3/2) / (4*sqrt(2)* pi^(3/2));

	#contour plot
	cc= a*s1^2 + b*s1+N

	p = (pi/180)*k*J* (exp(-N/2) * (-4 * sqrt(a)* b + (4*a + b^2)* exp(b^2/(8*a))*sqrt(2*pi)*erfc(b/(2*sqrt(2)*sqrt(a)))))/(8*a^(5/2))
	


	booleanC<-(cc<threshold)

	maxV<-max(mesh$x[which(booleanC)])
	print(paste("jpdf 95% interval max ratio value for statistical isotropy: ",maxV))
	doRotation=(max(mesh$x[which(booleanC)])<RRest)
		
	
	
	
#	contour(as.matrix(Rf),as.matrix(thetaf),t(p))
#	contour(Rf,thetaf,t(p))

	return(list(c=cc,p=p,bc=booleanC,doRotation=doRotation))




}




#######################################################################
#biharmonics Spline interpolation			
#   scattered data coordinates -> x, y			
#   data value @(x,y)          -> z			
#   grid                       ->xi, yi			
#######################################################################
biharmonicSplineAnisotropy<-function(x,y,z,xi,yi){
	#sortrows function
	SortMat <- function(Mat, Sort)
	{
 	       m <- do.call("order", as.data.frame(Mat[, Sort]))
	        Mat[m, ]
	}

        #sorting input
        Mat<-cbind(x,y,z)
        Mat<-SortMat(Mat,c(2,1))
        x<-Mat[,1]
        y<-Mat[,2]
        z<-Mat[,3]

        #Move to Complex field

        xy<-as.matrix(x+y*1i)
        d<-matrix(xy,length(xy),length(xy))

        #Calculate distances
        d<-abs(d-t(d))

        #remove zeros so we can use log
	ind=which(d==0)
	if(length(ind)>0){
        d[ind]=1
	}
        d=d^2*(log(d)-1)

        #return replace zeros to diag elements
        d[seq(1,length(d),by=(dim(d)[1]+1))]=0
	
	ind=which(is.na(d)==TRUE)
	if(length(ind)>0){
        	d[ind]=1
	}
        
	#calculate weights
        wweights=solve(d)%*%z
        rm(d)
        m=dim(xi)[1]
        n=dim(xi)[2]
        zi<-matrix(0,m,n)
        xy=t(xy)


        for(i in 1:m){
         for(j in 1:n){
             d = abs(xi[i,j]+(1i*yi[i,j])-xy)
             mask=which(d==0)

             if (length(mask)>0) {d[mask]=1}
             g<-(d)^2*(log(d)-1)
             if(length(mask)>0) { g[mask]=0}
              zi[i,j]=g%*% wweights
         }
        }

  return(list(x=xi,y=yi,z=zi))

}


	#Meshgrid function
	meshgrid <- function(a,b) {
	  list(
	       x=outer(b*0,a,FUN="+"),
	       y=outer(b,a*0,FUN="+")
	       )
	}


########################End of internal functions######################



