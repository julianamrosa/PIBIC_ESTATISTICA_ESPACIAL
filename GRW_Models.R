########################## GRW MODELS ##########################

# Macro for searching the optimum Bandwidth
# 
# REQUIRED PARAMETERS
# DATA = the name of the SAS data set to be used
# YVAR = the name of the dependent or response variable
# XVAR = the name of the independent or explicative variables. 
# A blank space should separate the names. Note: an intercept variable must not be created in advance
# DCOORD = the name of the SAS data set with the geographic coordinates
# OUTPUT = the name of the SAS data set to be used as output results
# MINV = the minimum distance between two locations i and k to be consider
# MIDDLEV = the middle distance between two locations i and k to be consider
# MAXV = the maximum distance between two locations i and k to be consider
# METHOD = there are three choices:
# FIXED_G asks the program to compute the bandwidth as fixed gaussian;
# FIXED_BSQ to compute the bandwidth as fixed bi-square; 
# ADAPTIVEN to compute the bandwidth as adaptive bi-square ($n$ values) and
# ADAPTIVE_BSQ to compute the bandwidth as adaptive bi-square ($one$ value)
# DISTANCEKM = if the distances between two locations will be computed in Km using the Basic formulae for calculating spherical distance. 
# The default value is NO, so the distance is computed using euclidian distance.

golden2 <- function(y, x, xvarglobal, weight, lat, long, output, method, model=GAUSSIAN, bandwidth,offset=NULL,distancekm=NO){
  E <- 10
  n <- length(y)
  wt <- matrix(1,nrow=n,ncol=1)
  if (is.null(offset)){
    offset <- matrix(0, nrow=n, ncol=1)
  }
  x <- cbind(matrix(1, nrow=n, ncol=1), x)
  nvar <- ncol(x)
  yhat <- matrix(0, nrow=n, ncol=1)
  alphai <- matrix(0, nrow=n, ncol=1)
  S <- matrix(0, nrow=n, ncol=1)
  
  # global estimates #
  
  if (model=="poisson" | model=="negbin"){
    uj <- (y+mean(y))/2
    nj <- log(uj)
    parg <- sum((y-uj)^2/uj)/(n-nvar)
    ddpar <- 1
    cont <- 1
    while (abs(ddpar)>0.000001 & cont<100){
      dpar <- 1
      parold <- parg
      cont1 <- 1
    }
    if (model=="poisson"){
      alphag <- E^-6
      parg <- 1/alphag
    }
    if (model=="negbin"){
      if (cont>1){
        parg <- 1/(sum((y-uj)^2/uj)/(n-nvar))
      }
      while (abs(dpar)>0.000001 & cont1<200){
        parg <- ifelse(parg<E^-10,E^-10,parg)
        g <- sum(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj))
        hess <- sum(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/((parg+uj)^2))
        hess <- ifelse(hess==0, E^-23,hess)
        par0 <- parg
        parg <- par0-solve(hess)*g
        dpar <- parg-par0
        cont <- cont1+1
        if(parg>E^6){
          parg <- E^6
          dpar <- 0
        }
        alphag <- 1/parg
      }
    }
    devg <- 0
    ddev <- 1
    cont2 <- 0
    while (abs(ddev)>0.000001 & cont2<100){
      uj <-ifelse(uj>E^100,E^100,uj)
      Ai <- (uj/(1+alphag%*%uj))+(y-uj)*(alphag%*%uj/1+2%*%alphag%*%uj+alphag^2%*%uj*uj)
      Ai <- ifelse(Ai<=0,E^-5,Ai)
      zj <- nj+(y-uj)/(Ai*(1+alphag%*%uj))-offset
      if (det(t(x)%*%(Ai*x))==0){
        bg <- matrix(0,nrow=nvar,ncol=1)
      }
      else{
        bg <- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
      }
      nj <- ifelse(nj>E^2,E^2,nj)
      uj <- exp(nj)
      olddev <- devg
      uj <- ifelse(uj<E^-150,E^-150,uj)
      tt <- y/yj
      tt <- ifelse(tt=0,E^-10,tt)
      if(model=="poisson"){
        devg <- 2*sum(y*log(tt)-(y-uj))
      }
      if(model=="negbin"){
        devg <- 2*sum(y%*%log(tt)-(y+1/alphag)*log((1+alphag*y)/(1+alphag*uj)))
      }
      if (cont2>100){
        ddev <- 0.0000001
      }
      else{
        ddev <- devg-olddev
        cont2 <- cont2+1
      }
      ujg <- uj
      cont <- cont+1
      ddpar <- parg-parold
    }
  }
  distance <- dist(COORD,"euclidean")
  seq <- seq(1,n)
  
  cv <- function(n, wt, x, xa, y, ujg, yhat, nvar,hv, coord, distance, seq, offset, alphag, alphai, S, parg){
    if(method=="adaptiven"){
      hv <- matrix(0,1,1)
      yhat <- matrix(0,1,1)
      output <<- hv[,"h"]
      # aqui faltam os laços sobre os quais eu perguntei ao professor
    }
    if(method=="fixed_g"|method=="fixed_bsq"|method="adaptive_bsq"){
    }
    for(j in i:n){
      seqi <- matrix(i,n,1)
      distan <- cbind(seqi, sequ, as.matrix(distance)[,i])
      if(distancekm=="yes"){
        dist[,3] <- dist[,3]*111
      }
    }
  u <- nrow(distance)
  w <- matrix(0,u,1)
  for(jj in 1:u){
    w[jj] <- exp(-0.5*(dist[jj,3]/h)^2)
    if(method=="fixed_bsq") {
      w[jj] <- (1-(dist[jj,3]/h)^2)^2
    } #nao precisa adptiven
    if(bandwidth=="cv"){
      w[i] <-0
    }
    if(method=="fixed_bsq"){
      posit <- which(!distance[,3]>h==0)
      posit <- distance[which(!distance[,3]>3==0)]
      w[posit] <- 0
    }
    if(method=="adaptive_bsq"){
      distance <- 
    }
  }
    
    #parei na linha 148
    
  } #fecha CV
  
  
} #fecha golden

# trocas ----
# position ---- posit

teste <- which(!exemplo==0)           #retorna apenas indices
teste <- exemplo[which(!exemplo==0)]  #retorna os valores




