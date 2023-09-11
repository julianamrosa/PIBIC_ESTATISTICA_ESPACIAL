library(readr)
library(dplyr)
library(tidyverse)

golden <- function(DATA,YVAR, XVAR, XVARGLOBAL=NULL, WEIGHT=NULL, LAT, LONG, 
                   GLOBALMIN="YES", METHOD, MODEL="NEGBIN", BANDWIDTH="CV",
                   OFFSET=NULL,DISTANCEKM="NO"){
  # distancekm, model e offset = default
  E <- 10
  y <<- DATA[,YVAR] #the name of the dependent or response variable
  x <<- DATA[XVAR] # the name of the independent or explicative variables.
  N <<- length(y)
  wt <<- rep(1, N)
  if (!is.null(WEIGHT)){
    wt <<- DATA[,WEIGHT]
    wt <<- as.matrix(wt)
  }
  Offset <<- rep(0, N)
  if (!is.null(OFFSET)){
    Offset <<- DATA[,OFFSET]
    Offset <<- as.matrix(Offset)
  }
  x <<- as.matrix(cbind(rep(1,N),x))
  nvar <<- ncol(x)
  if (!is.null(XVARGLOBAL)){
    xa <<- DATA[,XVARGLOBAL]
    xa <<- as.matrix(xa)
  }
  yhat <<- rep(0,N) 
  alphai <<- rep(0,N)   
  S <<- rep(0,N)
  
  # global estimates #
  y <<- as.numeric(y) 
  uj <- (y+mean(y))/2
  print(uj)
  nj <- log(uj)
  parg <<- sum((y-uj)^2/uj)/(N-nvar)
  ddpar <- 1
  cont <- 1
  while (abs(ddpar)>0.000001 & cont<100){
    dpar <- 1
    parold <- parg
    cont1 <- 1
    if (toupper(MODEL)=="POISSON"){
      alphag <<- E^(-6)
      parg <<- 1/alphag
    }
    else if (toupper(MODEL)=="NEGBIN"){
      if (cont>1){
        parg <<- 1/(sum((y-uj)^2/uj)/(N-nvar))
      }
      while (abs(dpar)>0.000001 & cont1<200){
        parg <<- ifelse(parg < E^(-10),E^(-10),parg)
        g <- sum(digamma(parg+y)- digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj))
        hess <- sum(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/((parg+uj)^2))
        hess <- ifelse(hess==0, E^(-23),hess)
        par0 <- parg
        parg <<- par0-solve(hess)%*%g
        dpar <- parg-par0
        cont1 <- cont1+1
        if(parg > E^6){
          parg <<- E^6
          dpar <- 0
        }
        parg <<- as.numeric(parg)
      }
      alphag <<- as.numeric(1/parg)
    }
    devg <- 0
    ddev <- 1
    cont2 <- 0
    while (abs(ddev)>0.000001 & cont2<100){
      uj <- ifelse(uj>E^100,E^100,uj)
      Ai <<- as.numeric((uj/(1+alphag*uj))+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj*uj)))
      Ai <<- ifelse(Ai<=0,E^(-5),Ai)
      # Ai <- as.numeric(Ai)
      zj <- nj+(y-uj)/(Ai*(1+alphag*uj))-Offset
      if (det(t(x)%*%(Ai*x))==0){ 
        bg <- rep(0,nvar)
      }
      else{
        bg <- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
      }
      nj <- x%*%bg+Offset
      nj <- ifelse(nj>E^2,E^2,nj)
      uj <- as.numeric(exp(nj))
      olddev <- devg
      uj <- ifelse(uj<E^(-150),E^(-150),uj)
      tt <- y/uj
      tt <- ifelse(tt==0,E^(-10),tt)
      if(toupper(MODEL)=="POISSON"){
        devg <- 2*sum(y*log(tt)-(y-uj))
      }
      if(toupper(MODEL)=="NEGBIN"){ #else if? 
        devg <- 2*sum(y*log(tt)-(y+1/alphag)*log((1+alphag*y)/(1+alphag*uj)))
      }
      if (cont2>100){
        ddev <- 0.0000001
      }
      else{
        ddev <- devg-olddev
      }
      cont2 <- cont2+1
    }
    ujg <<- uj
    cont <- cont+1
    ddpar <- parg-parold
  } #fecha while linha 39
  
  ###################################
  
  LONG <- DATA[, LONG]
  LAT  <- DATA[, LAT]
  COORD <<- matrix(c(LONG, LAT), ncol=2)
  distance <<- dist(COORD,"euclidean")
  sequ <<- 1:N
  cv <- function(h){ #linha 115 do SAS
    for (i in 1:N){ 
      for(j in 1:N){
        seqi <- rep(i,N)
        distan <- cbind(seqi, t(sequ), as.matrix(distance)[,i])
        if(toupper(DISTANCEKM)=="YES"){
          distan[,3] <- distan[,3]*111
        }
      }
      u <- nrow(distan)
      w <- rep(0,u)
      for(jj in 1:u){
        w[jj] <- exp(-0.5*(distan[jj,3]/h)^2)
        if(toupper(METHOD)=="FIXED_BSQ"){
          w[jj] <- (1-(distan[jj,3]/h)^2)^2  
        }
        if(toupper(BANDWIDTH)=="CV"){
          w[i] <-0
        }
      }
      if(toupper(METHOD)=="FIXED_BSQ"){
        position <- which(distan[,3]>h) 
        w[position] <- 0
      }
      else if(toupper(METHOD)=="ADAPTIVE_BSQ"){
        distan <- distan[order(distan[,3]),] 
        distan <- cbind(distan,1:nrow(distan)) #se necessario, usar t(1:nrow(distan))
        w <- matrix(0,N,2)
        hn <- distan[h,3]
        for(jj in 1:N){
          if(distan[jj,4]<=h){
            w[jj,1] <- (1-(distan[jj,3]/hn)^2)^2
          }
          else{
            w[jj,1] <- 0
          }
          w[jj,2] <- distan[jj,2]
        }
        if(toupper(BANDWIDTH)=="CV"){
          w[which(w[,2]==i)] <- 0
        }
        w <- w[order(w[,2]),]
        w <- w[,1]
      } #fecha else if
      uj <- ujg
      nj <- log(uj)
      par <- parg
      ddpar <- 1
      cont <- 1
      while (abs(ddpar) > 0.000001 & cont < 100) {
        dpar <- 1
        parold <- par
        cont1 <- 1
        if (toupper(MODEL) == "POISSON") {
          alpha <- E^(-6)
          par <- 1/alpha
        }
        if (toupper(MODEL) == "NEGBIN") {
          if (par <= E^(-5)) {
            if (i>1) {
              par <- 1/alphai[i - 1]
            }
          }
          while (abs(dpar) > 0.000001 & cont1 < 200) {
            par <- ifelse(par < E^(-10), E^(-10), par)
            g <- sum(w * wt * (digamma(par + y) - digamma(par) + log(par) + 1 - log(par + uj) - (par + y) / (par + uj)))
            hess <- sum(w * wt * (trigamma(par + y) - trigamma(par) + 1 / par - 2 / (par + uj) + (y + par) / (par + uj)^2))
            hess <- ifelse(hess == 0, E^(-23), hess)
            par0 <- par
            par <- par0 - solve(hess) %*% g
            dpar <- par - par0
            cont1 <- cont1 + 1
            if (par > E^6) {
              par <- E^6
              dpar <- 0
            }
          }
          alpha <- 1 / par
          alpha <- as.numeric(alpha)
        }
        dev <- 0
        ddev <- 1
        cont2 <- 1
        while (abs(ddev) > 0.000001 & cont2 < 100) {
          uj <- ifelse(uj > E^100, E^100, uj)
          Ai <<- as.numeric((uj/(1 + alpha * uj)) + (y-uj) * (alpha * uj / (1 + 2 * alpha * uj + alpha^2 *uj*uj)))
          Ai <<- ifelse(Ai <= 0, E^(-5), Ai)
          zj <- nj + (y - uj) / (Ai * (1 + alpha * uj)) - Offset
          if (det(t(x) %*% (w * Ai * x * wt)) == 0) {
            b <- matrix(0, nvar, 1)
          } 
          else {
            b <- solve(t(x)%*%(w*Ai*x*wt)) %*% t(x) %*% (w * Ai * wt * zj)
          }
          nj <- x %*% b + Offset
          nj <- ifelse(nj > E^2, E^2, nj)
          uj <- exp(nj)
          olddev <- dev
          uj <- ifelse(uj < E^(-150), E^(-150), uj)
          tt <- y / uj
          tt <- ifelse(tt == 0, E^(-10), tt)
          if (toupper(MODEL) == "POISSON") {
            dev <- 2 * sum(y * log(tt) - (y - uj))
          }
          if (toupper(MODEL) == "NEGBIN") {
            dev <- 2*sum(y*log(tt)-(y+1/alpha)*log((1+alpha*y)/(1+alpha*uj)))
          }
          if (cont2 > 100) {
            ddev <- 0.0000001
          } 
          else {
            ddev <- dev - olddev
          }
          cont2 <- cont2 + 1
        }
        cont <- cont + 1
        ddpar <- par - parold
      }
      if(toupper(METHOD)=="FIXED_G"|toupper(METHOD)=="FIXED_BSQ"|toupper(METHOD)=="ADAPTIVE_BSQ"){
        yhat[i] <<- uj[i]
        alphai[i] <<- alpha
        if(det(t(x)%*%as.matrix((w*Ai*x*wt)))==0){
          S[i] <<- 0
        }
        else{
          S[i] <<- (x[i,]%*%solve(t(x)%*%(w*Ai*x*wt))%*%t(x*w*Ai*wt))[i]
        }
        if(!is.null(XVARGLOBAL)){
          if(i==1){
            W_f <- cbind(matrix(i,N,1),w,t(seq(1,nrow(w)))) 
          }
          else{
            W_f <- rbind(W_f,c(cbind(matrix(i,N,1),w,t(seq(1,nrow(w)))))) 
          }
        }
      }
    }
    if(toupper(METHOD)=="FIXED_G"|toupper(METHOD)=="FIXED_BSQ"|toupper(METHOD)=="ADAPTIVE_BSQ"){
      if(!is.null(XVARGLOBAL)){
        uj <- (y+mean(y))/2
        nj <- log(uj)
        if(toupper(MODEL)=="POISSON"){
          alphag <<- E^(-6)
        }
        dev <- 0
        ddev <- 1
        cont2 <- 0
        while(abs(ddev)>0.000001 & cont2 < 100){
          uj <- ifelse(uj > E^100, E^100,uj)
          Aa <- (uj/(1+alphag*uj)) + (y-uj) * (alphag*uj/(1+2*alphag*uj+alphag^2*uj*uj))
          Aa <- ifelse(Aa <= 0, E^(-5), Aa)
          za <- nj+(y-uj)/(Aa*(1+alphag*uj)) - Offset
          if(det(t(xa)%*%(Aa*xa))==0){
            ba <- matrix(0,nvar,1)
          } else{
            ba <- solve(t(xa)%*%(Aa*xa))%*%t(xa)%*%(Aa*za)
          }
          nj <-xa %*% ba + Offset
          nj <- ifelse(nj > E^2, E^2, nj)
          uj <- exp(nj)
          olddev <- devga
          uj <- ifelse(uj < E^(-150), E^(-150), uj)
          tt <- y/uj
          tt <- ifelse(tt == 0, E^(-10),tt)
          if(toupper(MODEL)=="POISSON"){
            devga <- 2*sum(y*log(tt)-(y-uj)) 
          }
          if(toupper(MODEL)=="NEGBIN"){
            devga <- 2*sum(y*log(tt)-(y+1/alphag)*log((1+alphag*y)/(1+alphag*uj)))
          }
          ddev <- devga - olddev
          cont2 <- cont2 + 1
        }
        diffba <- 1
        contb <- 1
        while(abs(diffba) > 0.00001 & contb < 50){
          oldba <- ba
          hat_matrix <- matrix(0,N,N) 
          for(i in 1:N){
            m1 <- (i-1)*n+1
            m2 <- m1 + (n-1)
            w <- w_f[m1:m2,2]
            if(toupper(MODEL)=="POISSON"){
              alpha <- E^(-6)
            }
            if(toupper(MODEL)=="NEGBIN"){
              par <- 1/alphai[i]
              if(par <= E^(-5)){
                if(i>1){
                  par <- 1/alphai[i-1]
                }
              }
              while(abs(dpar)> 0.000001 & cont1 < 200){
                par <- ifelse(par < E^(-10),E^(-10),par)
                g <- sum(w*wt*(digamma(par+y)-digamma(par)+log(par)+1-log(par+uj)-(par+y)/(par+uj)))
                hess <- sum(w * wt * (trigamma(par + y) - trigamma(par) + 1 / par - 2 / (par + uj) + (y + par) / (par + uj)^2))
                hess <- ifelse(hess == 0, E^(-23), hess)
                par0 <- par
                par <- par0 - solve(hess) %*% g
                dpar <- par - par0
                cont1 <- cont1 + 1
                if(par > E^6){
                  par <- E^6
                  dpar <- 0
                }
              }
              alpha <- 1/par
            }
            print(alpha)
            dev <- 0
            ddev <- 1
            cont2 <- 0
            while (abs(ddev) > 0.000001 & cont2 < 100) {
              uj <- ifelse(uj > E^100, E^100, uj)
              Ai <- (uj / (1 + alpha * uj)) + (y - uj) * (alpha * uj / (1 + 2 * alpha * uj + alpha^2 * uj*uj))
              Ai <- ifelse(Ai <= 0, E^(-5), Ai)
              zj <- nj + (y - uj) / (Ai * (1 + alpha * uj)) - Offset - xa %*% ba
              if (det(t(x) %*% (w * Ai * x * wt)) == 0) {
                b <- matrix(0, nvar, 1)
              } else {
                b <- solve(t(x) %*% (w * Ai * x * wt)) %*% t(x) %*% (w * Ai * wt) %*% zj
              }
              nj <- x %*% b + Offset + xa %*% ba
              nj <- ifelse(nj > E^2, E^2, nj)
              uj <- exp(nj)
              olddev <- dev
              uj <- ifelse(uj < E^(-150), E^(-150), uj)
              tt <- y / uj
              tt <- ifelse(tt == 0, E^(-10), tt)
              
              if (toupper(MODEL) == "POISSON") {
                dev <- 2 * sum(y * log(tt) - (y - uj))
              }
              if (toupper(MODEL) == "NEGBIN") {
                dev <- 2 * sum(y * log(tt) - (y + 1 / alpha) * log((1 + alpha * y) / (1 + alpha * uj)))
              }
              
              ddev <- dev - olddev
              cont2 <- cont2 + 1
            }
            C <- solve(t(x) %*% (w * Ai * x * wt)) %*% t(x) %*% (w * Ai * wt)
            R_matrix <- rbind(R_matrix, (x[i,] %*% C))
            Z_matrix <- cbind(Z_matrix, (zj + Offset + xa %*% ba))
            yhat[i] <<- uj
          }
          hat_matrix <- R_matrix %*% Z_matrix / diag(Z_matrix)
          uj <- yhat
          nj <- log(uj)
          Aa <- (uj / (1 + alphag * uj)) + (y - uj) * (alphag * uj / (1 + 2 * alphag * uj + alphag^2 * uj*uj))
          Aa <- ifelse(Aa <= 0, E^(-5), Aa)
          za <- nj + (y - uj) / (Aa * (1 + alphag * uj)) - Offset
          ba <- solve((xa %*% Aa) %*% t(xa) %*% (I(n) - hat_matrix) %*% (xa %*% Aa %*% wt) %*% (I(n) - hat_matrix) + hat_matrix) %*% ((xa %*% Aa) %*% t(xa) %*% (I(n) - hat_matrix) %*% (za %*% wt))
          nj <- xa %*% ba + Offset
          nj <- ifelse(nj > E^2, E^2, nj)
          uj <- exp(nj)
          diffba <- sum(oldba - ba)
          contb <- contb + 1
        }
        S <<- diag((I(n) - hat_matrix) %*% xa %*% solve((xa %*% Aa) %*% t(xa) %*% (I(n) - hat_matrix) %*% xa %*% wt) %*% (xa %*% Aa %*% wt) %*% (I(n) - hat_matrix) + hat_matrix)
      }
      # CV <- sum((y - yhat) %*% wt %*% (y - yhat)) pode ser que tenha que transpor. Ver linha 328 do SAS
      CV <- t((y-yhat) %*% wt %*% (y - yhat))
      if (toupper(MODEL) == "POISSON") {
        ll <- sum(-yhat + y * log(yhat) - lgamma(y + 1))
        npar <- sum(S)
      }
      if (toupper(MODEL) == "NEGBIN") {
        ll <- sum(y * log(alphai * yhat) - (y + 1 / alphai) * log(1 + alphai * yhat) + lgamma(y + 1 / alphai) - lgamma(1 / alphai) - lgamma(y + 1))
        npar <- sum(S) + sum(S) / nvar
      }
      AIC <- 2 * npar - 2 * ll
      AICC <- AIC + (2 * npar * (npar + 1)) / (N - npar - 1)
      if (toupper(BANDWIDTH) == "AIC") { 
        CV <- AICC 
      }
    }
    res <- cbind(CV, npar)
    return(res)
  } #fecha CV!!!!
  # DEFINING GOLDEN SECTION SEARCH PARAMETERS #
  
  if(toupper(METHOD)=="FIXED_G"|toupper(METHOD)=="FIXED_BSQ"){
    ax <- 0
    bx <- floor(max(distance)+1)
    if(toupper(DISTANCEKM)=="YES"){
      bx <- bx*111
    } 
  }
  if(toupper(METHOD)=="ADAPTIVE_BSQ"){
    ax <- 5
    bx <- N
  }
  r <- 0.61803399
  tol <- 0.1
  if(toupper(GLOBALMIN)=="NO"){
    lower <- ax
    upper <- bx
    xmin <- matrix(0,1,2)
    GMY <- 1
    ax1 <- lower[GMY]
    bx1 <- upper[GMY]
  } else{
    lower <- cbind(ax,(1-r)*bx,r*bx)
    upper <- cbind((1-r)*bx,r*bx,bx)
    xmin <- matrix(0,3,2)
    for(GMY in 1:3){
      ax1 <- lower[GMY]
      bx1 <- upper[GMY]
    }
  }
  h0 <- ax1
  h3 <- bx1
  h1 <- bx1-r*(bx1-ax1)
  h2 <- ax1+r*(bx1-ax1)
  print(rbind(c("h0", "h1", "h2", "h3"), c(h0, h1, h2, h3)))
  
  # /***************************************/ #
  
  res1 <- cv(h1)
  CV1 <- res1[1]
  res2 <- cv(h2)
  CV2 <- res2[1]
  
  if(toupper(METHOD)=="FIXED_G"|toupper(METHOD)=="FIXED_BSQ"|toupper(METHOD)=="ADAPTIVE_BSQ"){
    var_ <- data.frame() #linha 403 do SAS
    if(GMY==1){
      var_ <- c(GMY,h1,CV1,h2,CV2) #revisar aqui
      names(var_) <- c('GMY','h1','CV1','h2','CV2')
      View(var_)
    }
  }  
  int <- 1
  print(cv(h1))
  while(abs(h3-h0) > tol*(abs(h1)+abs(h2)) & int<200){
    if(CV2<CV1){
      h0 <- h1
      h1 <- h3-r*(h3-h0)
      h2 <- h0+r*(h3-h0)
      CV1 <- CV2
      res2 <- cv(h2)
      CV2 <- res2[1]
    }
    else{
      h3 <- h2
      h1 <- h3-r*(h3-h0)
      h2 <- h0+r*(h3-h0)
      CV2 <- CV1
      res1 <- cv(h1)
      CV1 <- res1[1]
    }
    int <- int+1  
  } #fecha while
  if(CV1<CV2){
    golden <- CV1
    xmin[GMY,1] <- golden
    xmin[GMY,2] <- h1
    npar <- res1[1]
    if(toupper(METHOD)=="ADAPTIVE_BSQ"){
      xmin[GMY,2] <- floor(h1)
    }
  }
  else{
    golden <- CV2
    xmin[GMY,1] <- golden
    xmin[GMY,2] <- h2
    npar <- res2[1]
    if(toupper(METHOD)=='ADAPTIVE_BSQ'){
      xmin[GMY,2] <- floor(h2)
    }
  }
  if (toupper(METHOD)=="FIXED_G" | toupper(METHOD)=="FIXED_BSQ" | toupper(METHOD)=="ADAPTIVE_BSQ"){
    print(golden)
    print(c("xmin", xmin[GMY, 2]))
    if (toupper(BANDWIDTH=="AIC")){
      print(npar)
    }
  }
  min_bandwidth <- as.data.frame(xmin)
  names(min_bandwidth) <- c("golden", "bandwidth")
  if (GLOBALMIN=="YES"){
    print(min_bandwidth)
    xming <- xmin[, 2][xmin[, 1]==min(xmin[, 1])]
    print("Global Minimum")
    print(c(xming, "Da Silva and Mendes, 2018"))
  }
  if (toupper(METHOD)=="FIXEG_G"|toupper(METHOD)=="FIXED_BSQ"|toupper(METHOD)=="ADAPTIVE_BSQ"){
    h_<- min(min_bandwidth[,2]) #proc sql
  }
} #fecha golden

# TESTES #
example2 <- read.csv("example2.csv")
example2 <- example2 %>% 
  mutate(Le = log(Loe))

system.time(golden(example2,YVAR="Mort2564", XVAR=c('Professl','Elderly','OwnHome','Unemply'), 
                   LAT="Y", LONG="X",MODEL="NEGBIN",OFFSET="Le",METHOD="FIXED_G",BANDWIDTH="AIC",GLOBALMIN = 'NO')) #25 minutos

traceback()
# #golden(DATA=nakaya,YVAR=Mort2564,XVAR=c('Professl' 'Elderly' 'OwnHome' 'Unemply'),
#   LONG=x,LAT=Y,OUTPUT=band,MODEL=POISSON,OFFSET=Le,METHOD=FIXED_G,
#   BANDWIDTH=AIC)

# /*******************************************************************************/
#   /* Macro for estimating GWR Model */
#   /* REQUIRED PARAMETERS
# /*    DATA = the name of the SAS data set to be used
# /*    YVAR = the name of the dependent or response variable
# /*    XVAR = the name of the independent or explicative variables. A blank space
# /*           should separate the names. Note: an intercept variable must not be
# /*           created in advance
# /*   WEIGHT = the name of the sample weight variable 
# /*  DCOORD = the name of the SAS data set with the geographic coordinates
# /*    GRID = the name of the SAS data set with the grid of geographic coordinates
# /*           the standard errors of complex data
# /*     DHV = the name of the SAS data set with the bandwidth adaptive ($n$ values),
# /*           which must have an unique variable
# /*       H = A pre-defined bandwidth value for METHOD equal to FIXED or ADAPTIVE1
# /*    MAXV = the maximum distance between two locations i and k to be consider
# /*  METHOD = there are three choices:
#   /*           FIXED_G asks the program to compute the bandwidth as fixed gaussian;
# /*           FIXED_BSQ to compute the bandwidth as fixed bi-square; 
# /*           ADAPTIVEN to compute the bandwidth as adaptive bi-square ($n$ values) and
# /*			 ADAPTIVE_BSQ to compute the bandwidth as adaptive bi-square ($one$ value)
# /*  DISTANCEKM = if the distances between two locations will be computed in Km
# /*               using the Basic formulae for calculating spherical distance. The 
# /*               default value is NO, so the distance is computed using euclidian
# /*               distance.
# /********************************************************************************/

GWNBR <- function(DATA, YVAR, XVAR, XVARGLOBAL=NULL, WEIGHT=NULL, LAT, LONG,
                  GRID=NULL, DHV, METHOD, MODEL="NEGBIN",OFFSET=NULL,
                  DISTANCEKM="NO", H=NULL){
  E <- 10
  y <- DATA[,YVAR]
  x <- DATA[,which(names(DATA) %in% XVAR)]
  N <- length(y)
  xx <- as.matrix(cbind(rep(1,N),x))
  Yhat <- rep(0,N)
  Nvar <- ncol(x)
  if(!is.null(XVARGLOBAL)){
    xa <- as.matrix(XVARGLOBAL)
  }
  Wt <-rep(1, N)
  if(!is.null(WEIGHT)){
    Wt <- as.matrix(WEIGHT)
  }
  Offset <- rep(0, N)
  if(!is.null(OFFSET)){
    Offset <- as.numeric(OFFSET)
  }
  # global estimates # 
  x1 <- x
  if(!is.null(XVARGLOBAL)){
    x1 <<- cbind(x,xa)
  }
  nvar1 <- ncol(x1)
  uj <- (y+mean(y))/2
  nj <- log(uj)
  parg <- sum((y-uj)^2/uj)/(N-nvar1)
  ddpar <- 1
  cont <- 1
  while(abs(ddpar)>0.00001|cont<100){
    dpar <- 1
    parold <- parg
    cont1 <- 1
    if(toupper(MODEL)=="POISSON"){
      alphag <- E^(-6)
      parg <- 1/(sum((y-uj)^2/uj)/(N-nvar1))
    }
    if(toupper(MODEL)=="NEGBIN"){
      if(cont>1){
        parg <- 1/(sum((y-uj)^2/uj)/(N-nvar1))
      }
      while(abs(dpar)>0.000001|cont1<200){
        parg <- ifelse(parg<E^(-10),E^(-10),parg)
        g <- sum(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj))
        hess <- sum(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/((parg+uj)^2))
        hess <- ifelse(hess==0, E^-23,hess)
        par0 <- parg
        parg <- par0-solve(hess)%*%g
        dpar <- parg-par0
        cont <- cont1+1
        if(parg>E^6){
          parg <- E^6
          dpar <- 0
        }
      }  
      alphag <- 1/parg
    }
    devg <- 0
    ddev <- 1
    cont2 <- 0  
    while(abs(ddev)>0.000001|cont2<100){
      uj <- ifelse(uj>E^100,E^100,uj)
      Ai <<- (uj/(1+alphag*uj))+(y-uj)*(alphag*uj/1+2*alphag*uj+alphag^2*uj*uj)
      Ai <<- ifelse(Ai<=0,E^-5,Ai)
      zj <- nj+(y-uj)/(Ai*(1+alphag*uj))-Offset
      if (det(t(x)%*%(Ai*x))==0){
        bg <- rep(0,Nvar)
      }
      else{
        bg <- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
      }
      nj <- x%*%bg+Offset
      nj <- ifelse(nj>E^2,E^2,nj)
      uj <- exp(nj)
      olddev <- devg
      uj <- ifelse(uj<E^-150,E^-150,uj)
      tt <- y/uj
      tt <- ifelse(tt==0,E^-10,tt)
      if(toupper(MODEL)=="POISSON"){
        devg <- 2*sum(y*log(tt)-(y-uj))
      }
      if(toupper(MODEL)=="NEGBIN"){
        devg <- 2*sum(y*log(tt)-(y+1/alphag)*log((1+alphag*y)/(1+alphag*uj)))
      }
      if (cont2>100){
        ddev <- 0.0000001
      }
      else{
        ddev <- devg-olddev
        cont2 <- cont2+1
      }
    }
    Ujg <- uj
    cont <- cont+1
    ddpar <- parg-parold
  }
  varg <- diag(solve(t(x1*wt*Ai)%*%x1))
  if(toupper(MODEL)=="NEGBIN"){
    sealphag <- sqrt(1/abs(hess))/(parg^2)
  }
  if (!is.null(H)){
    hh <- H
    print(c("Bandwidth: ", hh))
  }
  else if (toupper(METHOD)=='FIXED_G'|toupper(METHOD)=='FIXED_BSQ'|toupper(METHOD)=='ADAPTIVE_BSQ'){
    hh <- h_
    print(c("Bandwidth: ",h))
  } # linha 596 SAS
  LONG <- DATA[, LONG]
  LAT  <- DATA[, LAT]
  COORD <<- matrix(c(LONG, LAT), ncol=2)
  POINTS <- matrix(c(x, y), ncol=2, byrow=F) # no SAS, usa-se uma condicao
  m <- nrow(POINTS)
  bi <- matrix(0,nvar*m,4)
  alphai<- matrix(0,m,3)
  BB <- matrix(0,nvar*n,n)
  rsqri <- matrix(0,m,1)
  sumwi <- matrix(0,m,1)
  varbi <- matrix(0,nvar*m,1)
  S <- matrix(0,m,1)
  S2 <- matrix(0,m,1)
  biT <-matrix(0,m,nvar+1)
  ym <- y - mean(y)
  
  # /******** calculating distance **********/ 
  distance <<- dist(COORD,"euclidean")
  sequ <<- 1:N
  distance <- as.data.frame(distance)
  for(i in 1:m){
    for(j in 1:N){
      seqi <<- matrix(i,n,1)
      distan <<- cbind(seqi, sequ, as.matrix(distance)[,i])
      if (toupper(DISTANCEKM)=="YES"){
        distan[,3] <<- distan[,3]*111
      }
    }
    u <<- nrow(distan)
    w <<- matrix(0,u,1)
    for(jj in 1:u){
      if(toupper(METHOD=="FIXED_G")){
        w[jj] <<- exp(-0.5*(distan[jj,3]/h)^2)
      }
      else if(toupper(METHOD=="FIXED_BSQ")){
        w[jj] <<- (1-(distan[jj,3]/h)^2)^2
      }
      if(is.null(GRID)){
        if(toupper(METHOD)=="ADAPTIVE_BSQ"){
          distan <<- distan[prder(distan[,3]),]
          distan <<- cbind(distan,1:nrow(distan))
          w <<- matrix(0,1)
          hn <<- distan[h,3]
          for(jj in 1:n){
            if(distan[jj,4] <= h){
              w[jj,1] <<- (1-(distan[jj,3]/hn)^2)^2 
            }
            else{
              w[jj,1] <<- 0
            }
            w <<- w[order(w[,2]),1]   
          }
          w <<-w[order(w[,2]),]
          w <<- w[,1]
        }
      }
    }
    # /******** model selection **********/ 
    uj <- ujg
    nj <- log(uj)
    par <- parg
    ddpar <- 1
    cont <- 1
    while(abs(ddpar) > 0.000001 & cont < 100) {
      dpar <- 1
      parold <- par
      cont1 <- 1
      if(toupper(MODEL)=="POISSON"){
        alpha <- E^(-6)
        par <- 1/alpha
      }
      if(touuper(MODEL)=="NEGBIN"){
        if(par <= E^(-5)){
          if(i > 1){
            par <- 1/alphai[i-1,2]
          }
        }
        while(abs(dpar) > 0.000001 & cont1 < 200){
          par <- choose(par < 1E-10, 1E-10, par)
          g <- sum(w * wt * (digamma(par + y) - digamma(par) + log(par) + 1 - log(par + uj) - (par + y) / (par + uj)))
          hess <- sum(w * wt * (trigamma(par + y) - trigamma(par) + 1 / par - 2 / (par + uj) + (y + par) / (par + uj))^2)
          hess <- choose(hess == 0, 1E-23, hess)
          par0 <- par
          par <- par0 - inv(hess) * g
          dpar <- par - par0
          cont1 <- cont1 + 1
          if(par > E^6){
            par <- E^6
            dpar <- 0
          }
        }
        alpha <- 1/par
      }
      dev <- 0
      ddev <- 1
      cont2 <- 0
      while (abs(ddev) > 0.000001 & cont2 < 100) {
        uj <- choose(uj > E^100, E^100, uj)
        Ai <<- (uj / (1 + alpha * uj)) + (y - uj) * (alpha * uj / (1 + 2 * alpha * uj + alpha^2 * uj * uj))
        Ai <<- choose(Ai <= 0, 1E-5, Ai)
        zj <- nj + (y - uj) / (Ai * (1 + alpha * uj)) - offset
        
        if (det(t(x) %*% (w * Ai * x * wt)) == 0) {
          b <- matrix(0, ncol = 1)
        } else {
          b <- solve(t(x) %*% (w * Ai * x * wt)) %*% (t(x) %*% (w * Ai * wt * zj))
        }
        
        nj <- x %*% b + offset
        nj <- choose(nj > 1E2, 1E2, nj)
        uj <- exp(nj)
        olddev <- dev
        uj <- choose(uj < 1E-150, 1E-150, uj)
        tt <- y / uj
        tt <- choose(tt == 0, 1E-10, tt)
        
        if (toupper(MODEL) == "POISSON") {
          dev <- 2 * sum(y * log(tt) - (y - uj))
        }
        
        if (toupper(MODEL) == "NEGBIN") {
          dev <- 2 * sum(y * log(tt) - (y + 1 / alpha) * log((1 + alpha * y) / (1 + alpha * uj)))
        }
        
        cont2 <- cont2 + 1
      }
      cont <- cont+1
      ddpar <- par - parold
    } #fecha linha 919
    if (det(t(x) %*% (w*Ai*x*wt)) == 0) {
      C <- matrix(0, nrow = nvar, ncol = 1)
    } else {
      C <- solve(t(x) %*% (w*Ai*x*wt)) %*% t(x) %*% (w*Ai*wt)
      varbs <- solve(t(x) %*% (sqrt(w)*Ai*x*wt))
      varb <- C %*%diag(1/Ai)%*% t(C)
    }
    if (toupper(MODEL)=="NEGBIN") {
      sealpha <- sqrt(1/abs(hess))/(par^2)
      alphai[i, 1] <- i
      alphai[i, 2] <- alpha
      alphai[i, 3] <- sealpha
    }
    # /*** standard GWR variance ****/ 
    if (!is.null(WEIGHT)) {
      varbgg <- varbs
    } else {
      varbgg <- varb
    }
    # /*******************************/
    m1 <- (i-1)*nvar+1
    m2 <- m1+(nvar-1)
    bi[m1:m2, 1] <- i
    bi[m1:m2, 2] <- b
    bi[m1:m2, 3] <- POINTS[i, 1]
    bi[m1:m2, 4] <- POINTS[i, 2]
    varbi[m1:m2, 1] <- diag(varb)
    if(is.null(GRID)){
      r <- x[i,]%*%C
      yhat[i] <- uj[i]
      m1 <-(i-1)*nvar+1
      m2 <- m1+(nvar-1)
      bi[m1:m2,1] <- i
      bi[m1:m2,2] <- b
      # /** creating non-stationarity matrix **/
      if (toupper(METHOD) != "ADAPTIVE_BSQ") {
        CCC <- cbind(x, w, wt)
        BB[m1:m2, ] <- solve(t(CCC[, 1:nvar]) %*% (CCC[, ncol(CCC)-1]*CCC[, 1:nvar] %*% CCC[, ncol(CCC)])) %*% t(CCC[, 1:nvar]) %*% t(CCC[, ncol(CCC)-1] %*% CCC[, ncol(CCC)]) #revisado!
      }
      varbi[m1:m2, 1] <- diag(varb)
      S[i] <- r[i]
      S2[i] <- r %*% t(r)
      biT[i, 1] <- i
      biT[i, 2:(nvar+1)] <- t(b)
      yhati <- uj
      TSS <- sum(y * w * t(wt) * y) - (sum(y * w * t(wt)))^2 / sum(w * wt)
      RSS <- sum((y - yhati) * wt) %*% (w * (y - yhati))
      rsqri[i] <- 1 - RSS / TSS
      w_ <- w
      w_ <- w_[order(w_[, 1]), ]
      sumwi[i] <- sum(w_[1:int(nrow(w_) * 1), 1])
      w_ <- w_[order(w_[,1]),]
      sumwi[i] <- sum(w_[1:int(nrow(w_)*1),1])
    }
    if (i == 1) { 
      W_f <- cbind(matrix(i, n, 1), W, t(1:nrow(w)))
    } else {
      W_f <- rbind(W_f, cbind(matrix(i, n, 1), W, t(1:nrow(w))))
    }
  } #fecha laco da linha 876
  w_f <- W_f
  W_f <- rbind(W_f, w_f)
  if (!is.null(XVARGLOBAL)) {
    uj <- (y + mean(y)) / 2
    nj <- log(uj)
    if (toupper(MODEL) == "POISSON") {
      alphag <- E^(-6)
    }
    devga <- 0
    ddev <- 1
    cont2 <- 0
    while (abs(ddev) > 0.000001 & cont2 < 100) {
      uj <- ifelse(uj > E^100, E^100, uj)
      Aa <- (uj / (1 + alphag * uj)) + (y - uj) * (alphag * uj / (1 + 2 * alphag * uj + alphag^2 * uj^2)) #revisar
      Aa <- ifelse(Aa <= 0, E^(-5), Aa)
      za <- nj+(y - uj)/ (Aa*(1 + alphag*uj)) - Offset
      if (det(t(xa) %*% (Aa * xa)) == 0) {
        ba <- matrix(0, nrow = nvar, ncol = 1)
      } else {
        ba <- solve(t(xa) %*% (Aa * xa)) %*% t(xa) %*% (Aa * za)
      }
      nj <- xa %*% ba + Offset
      nj <- ifelse(nj > E^2, E^2, nj)
      uj <- exp(nj)
      olddev <- devga
      uj <- ifelse(uj < E^(-150), E^(-150), uj)
      tt <- y / uj
      tt <- ifelse(tt == 0, E^(-10), tt)
      if (toupper(MODEL) == "POISSON") {
        devga <- 2 * sum(y * log(tt) - (y - uj))
      }
      if (toupper(MODEL) == "NEGBIN") {
        devga <- 2 * sum(y * log(tt) - (y + 1 / alphag) * log((1 + alphag * y) / (1 + alphag * uj)))
      }
      ddev <- devga - olddev
      cont2 <- cont2 + 1
    }
    diffba <-1
    contb <- 1
    while (abs(diffba) > 0.00001 & contb < 50) {
      oldba <- ba
      rm(hat_matrix, R_matrix, Z_matrix)
      for (i in 1:m) {
        m1 <- (i - 1) * m + 1
        m2 <- m1 + (m - 1)
        w <- w_f[m1:m2, 2]
        if (toupper(MODEL) == "POISSON") {
          alpha <- E^(-6)
        }
        if (toupper(MODEL) == "NEGBIN") {
          par <- 1 / alphai[i, 2]
          dev <- 0
          ddev <- 1
          cont2 <- 0
          while (abs(ddev) > 0.000001 & cont2 < 100) {
            uj <- ifelse(uj > E^100, E^100, uj)
            Ai <<- (uj / (1 + alpha * uj)) + (y - uj) * (alpha * uj / (1 + 2 * alpha * uj + alpha^2 * uj^2))
            Ai <<- ifelse(Ai <= 0, E^(-5), Ai)
            zj <- nj + (y - uj) / (Ai * (1 + alpha * uj)) - offset - xa %*% ba
            if (det(t(x) %*% (w * Ai %*% x %*% wt)) == 0) {
              b <- matrix(0, nrow = nvar, ncol = 1)
            } else {
              b <- solve(t(x) %*% (w * Ai %*% x %*% wt)) %*% t(x) %*% (w * Ai %*% wt %*% zj)
            }
            nj <- x %*% b + offset + xa %*% ba
            nj <- ifelse(nj > E^2, E^2, nj)
            uj <- exp(nj)
            olddev <- dev
            uj <- ifelse(uj < E^(-150), E^(-150), uj)
            tt <- y / uj
            tt <- ifelse(tt == 0, E^(-10), tt)
            if (toupper(MODEL) == "POISSON") {
              dev <- 2 * sum(y * log(tt) - (y - uj))
            }
            if (toupper(MODEL) == "NEGBIN") {
              dev <- 2 * sum(y * log(tt) - (y + 1 / alpha) * log((1 + alpha * y) / (1 + alpha * uj)))
            }
            ddev <- dev - olddev
            cont2 <- cont2 + 1
          }
          m1 <- (i - 1) * nvar + 1
          m2 <- m1 + (nvar - 1)
          bi[m1:m2, 1] <- i
          bi[m1:m2, 2] <- b
          C <- solve(t(x) %*% (w * Ai %*% x %*% wt)) %*% t(x) %*% (w * Ai %*% wt)
          varb <- C %*% diag(1 / Ai) %*% t(C)
          varbi[m1:m2, 1] <- diag(varb)
          if (toupper(MODEL) == "NEGBIN") {
            sealpha <- sqrt(1 / abs(hess)) / (par^2)
            alphai[i, 1] <- i
            alphai[i, 2] <- alpha
            alphai[i, 3] <- sealpha
          }
          R_matrix <- rbind(R_matrix, x[i,] %*% C)
          Z_matrix <- cbind(Z_matrix, zj + Offset + xa %*% ba)
          yhat[i] <- uj[i]
        }
      } #fecha laco linha 1087
      hat_matrix <- R_matrix %*% Z_matrix / diag(Z_matrix)
      uj <- yhat
      nj <- log(uj)
      Aa <- (uj / (1 + alphag * uj)) + (y - uj) * (alphag * uj / (1 + 2 * alphag * uj + alphag^2 * uj^2))
      Aa <- ifelse(Aa <= 0, 1E-5, Aa)
      za <- nj + (y - uj) / (Aa * (1 + alphag * uj)) - offset
      
      ba <- solve(t(xa) %*% Aa %*% (diag(n) - hat_matrix) %*% xa %*% wt) %*% t(xa) %*% Aa %*% (diag(n) - hat_matrix) %*% za %*% wt
      nj <- xa %*% ba + offset
      nj <- ifelse(nj > 1E2, 1E2, nj)
      uj <- exp(nj)
      diffba <- oldba - ba
      contb <- contb+1
    }
    Ca <- solve(t(xa) %*% Aa %*% (diag(n) - hat_matrix) %*% xa %*% wt) %*% t(xa) %*% Aa %*% wt %*% solve(t(xa) %*% Aa %*% (diag(n) - hat_matrix))
    varba <- diag(Ca %*% diag(1 / Aa) %*% t(Ca))
    S <- diag((diag(n) - hat_matrix) %*% xa %*% solve(t(xa) %*% Aa %*% (diag(n) - hat_matrix) %*% xa %*% wt) %*% t(xa) %*% Aa %*% wt %*% solve(t(xa) %*% Aa %*% (diag(n) - hat_matrix)) + hat_matrix)
    rm(hat_matrix, R_matrix, Z_matrix)
  }
  v11 <- sum(S)
  v2 <- sum(S2)
  v1 <- 2 * v11 - v2
  if (!is.null(XVARGLOBAL)) {
    v1 <- sum(S)
  }
  nparmodel <- n - v11
  if (v11 < v2) {
    v1 <- v11
  }
  v1 <- sum(S)
  if (!is.null(GRID)) {
    res <- y - yhat
    rsqr1 <- sum(res * wt) %*% res
    ym <- sum(y * wt) %*% y
    rsqr2 <- ym - sum(y * wt) %*% y^2 / sum(wt)
    rsqr <- 1 - rsqr1 / rsqr2
    rsqradj <- 1 - ((n - 1) / (n - v1)) * (1 - rsqr)
    sigma2 <- n * rsqr1 / ((n - v1) * sum(wt))
    root_mse <- sqrt(sigma2)
    cat("sigma2e:", sigma2, "Root MSE:", root_mse, "#GWR parameters:", v1, "#GWR parameters (model):", nparmodel, "#GWR parameters (variance):", v2, "\n")
    influence <- S
    resstd <- res / (sqrt(sigma2) * sqrt(abs(1 - influence)))
    CooksD <- resstd^2 * influence / (v1 * (1 - influence))
    df <- n - nvar
    stdbi <- sqrt(varbi)
    tstat <- bi[, 2] / stdbi
    probt <- 2 * (1 - pnorm(abs(tstat), df))
    malpha <- 0.05 * (nvar / v1)
    t_critical <- abs(qt(malpha / 2, df))
    ll <- -n * log(rsqr1 / n) / 2 - n * log(2 * acos(-1)) / 2 - sum((y - yhat) * (y - yhat)) / (2 * (rsqr1 / n))
    AIC <- 2 * v1 - 2 * ll
    AICc <- AIC + 2 * (v1 * (v1 + 1) / (n - v1 - 1))
    cat("AIC:", AIC, "AICc:", AICc, "\n")
    if (toupper(MODEL) == "POISSON") {
      tt <- y / yhat
      tt[tt == 0] <- E^(-10)
      dev <- 2 * sum(y * log(tt) - (y - yhat))
      ll <- sum(-yhat + y * log(yhat) - lgamma(y + 1))
      AIC <- 2 * v1 - 2 * ll
      AICc <- AIC + 2 * (v1 * (v1 + 1) / (n - v1 - 1))
      tt2 <- y / mean(y)
      tt2[tt2 == 0] <- 1E-10
      devnull <- 2 * sum(y * log(tt2) - (y - mean(y)))
      pctdev <- 1 - dev / devnull
      adjpctdev <- 1 - ((n - 1) / (n - v1)) * (1 - pctdev)
      cat("Deviance:", dev, "Full Log Likelihood:", ll, "PctDev:", pctdev, "Adj PctDev:", adjpctdev, "AIC:", AIC, "AICc:", AICc, "\n")
    }
    if (toupper(MODEL) == "NEGBIN") {
      tt <- y / yhat
      tt[tt == 0] <- E^(-10)
      dev <- 2 * sum(y * log(tt) - (y + 1 / alphai[, 2]) * log((1 + alphai[, 2] * y) / (1 + alphai[, 2] * yhat)))
      ll <- sum(y * log(alphai[, 2] * yhat) - (y + 1 / alphai[, 2]) * log(1 + alphai[, 2] * yhat) + lgamma(y + 1 / alphai[, 2]) - lgamma(1 / alphai[, 2]) - lgamma(y + 1))
      AIC <- 2 * (v1 + v1 / nvar) - 2 * ll
      AICc <- AIC + 2 * (v1 + v1 / nvar) * (v1 + v1 / nvar + 1) / (n - (v1 + v1 / nvar) - 1)
      tt2 <- y / mean(y)
      tt2[tt2 == 0] <- E^(-10)
      devnull <- 2 * sum(y * log(tt2) - (y + 1 / alphai[, 2]) * log((1 + alphai[, 2] * y) / (1 + alphai[, 2] * mean(y))))
      pctdev <- 1 - dev / devnull
      adjpctdev <- 1 - ((n - 1) / (n - (v1 + v1 / nvar))) * (1 - pctdev)
      cat("Deviance:", dev, "Full Log Likelihood:", ll, "PctDev:", pctdev, "Adj PctDev:", adjpctdev, "AIC:", AIC, "AICc:", AICc, "\n")
    }
    beta <- matrix(c(bi[, 1], bi[, 2]), nrow = n)
    beta2 <- beta
    if (toupper(MODEL) == "NEGBIN") {
      alpha <- matrix(c(alphai[, 1], alphai[, 2]), nrow = n)
      beta2 <- cbind(beta, alpha)
    }
    i <- seq(2, ncol(beta), 2)
    beta <- beta[, i]
    i <- seq(2, ncol(beta2), 2)
    beta2 <- beta2[, i]
    qntl <- quantile(beta2)
    qntl <- rbind(qntl, qntl[3, ] - qntl[1, ])
    descriptb <- cbind(beta2[, ], beta2[, ], beta2[, ])
  }
  if (!is.null(XVARGLOBAL)) {
    stdba <- sqrt(varba)
    tstatba <- ba / stdba
    probtba <- 2 * (1 - pt(abs(tstatba), df))
    ba_stda <- cbind(ba, stdba)
    
    rownames(ba_stda) <- XVARGLOBAL
    
    colnames(ba_stda) <- c("Par. Est.", "Std Error")
    
    print(ba_stda)
    print(data.frame(tstatba, probtba), row.names = XVARGLOBAL)
  }
  print(qntl, 
        label = "Quantiles of GWR Parameter Estimates",
        rownames = c("P25", "P50", "P75", "IQR"),
        colnames = c("Intercept", "xvar", if (toupper(MODEL) == "NEGBIN") "alpha" else NULL)
  )
  
  print(descriptb,
        label = "Descriptive Statistics",
        rownames = c("Mean", "Min", "Max"),
        colnames = c("Intercept", "xvar", if (toupper(MODEL) == "NEGBIN") "alpha" else NULL)
  )
  stdbeta <- matrix(stdbi, ncol = n)
  stdbeta2 <- stdbeta
  
  if (toupper(MODEL) == "NEGBIN") {
    stdalpha <- alphai[,3]
    stdbeta2 <- cbind(stdbeta, stdalpha)
  }
  
  qntls <- apply(stdbeta2, 2, function(col) quantile(col, probs = c(0.25, 0.5, 0.75)))
  qntls <- rbind(qntls, qntls[3,] - qntls[1,])
  descripts <- rbind(stdbeta2, stdbeta2, stdbeta2)
  
  malpha <- 0.05 * (nvar / v1)
  t_critical <- abs(qt(0.975, df))
  
  print(data.frame(malpha = malpha, t_critical = t_critical, df = df))
  colnames(qntls) <- c("Intercept", "xvar", if (toupper(MODEL) == "NEGBIN") "alpha" else NULL)
  print(qntls, 
        label = "Quantiles of GWR Standard Errors",
        rownames = c("P25", "P50", "P75", "IQR")
  )
  
  colnames(descripts) <- c("Intercept", "xvar", if (toupper(MODEL) == "NEGBIN") "alpha" else NULL)
  print(descripts,
        label = "Descriptive Statistics of Standard Errors",
        rownames = c("Mean", "Min", "Max")
  )
  # /****** Non-Stationarity Test *****************/
  if (is.null(GRID)) {
    if (toupper(METHOD) != "ADAPTIVE_BSQ") {
      BBk <- matrix(0, n, n)
      Vk <- matrix(0, nvar, 1)
      df1k <- matrix(0, nvar, 1)
      df2k <- matrix(0, nvar, 1)
      for (k in 1:nvar) {
        ek <- matrix(0, nvar, 1)
        ek[k] <- 1
        for (i in 1:n) {
          m1 <- (i - 1) * nvar + 1
          m2 <- m1 + (nvar - 1)
          BBk[i,] <- t(ek) %*% BB[m1:m2,]
        }
        Vk[k] <- t(y) %*% (1/n) * t(BBk) %*% (diag(n) - (1/n) * matrix(1, n, n)) %*% BBk %*% y
        df1k[k] <- sum((1/n) * t(BBk) %*% (diag(n) - (1/n) * matrix(1, n, n)) %*% BBk)
        df2k[k] <- sum(((1/n) * t(BBk) %*% (diag(n) - (1/n) * matrix(1, n, n)) %*% BBk) ** 2)
      }
      Vk <- ifelse(abs(Vk) <= E^(-8), 0, Vk)
      Fk <- (Vk / df1k) / sigma2
      ndf <- df1k ** 2 / df2k
      ddf <- n - v1
      ddf <- rep(ddf, nvar)
      probf <- 1 - pf(Fk, ndf, ddf)
      cat("\n\nNon-Stationarity Test (Leung et al., 2000)\n")
      print(data.frame(Vk = Vk, Fk = Fk, ndf = ndf, ddf = ddf, probf = probf))
    }
  }
  # /***** global estimates ***************/
  if (is.na(WEIGHT)) { #verificar se é NULL, NA ou alguma outra opcao
    vargd <- varg
    dfg <- n - nrow(bg)
    stdg <- sqrt(vargd)
    if (toupper(MODEL) == "NEGBIN") {
      alphag <- cbind(alphag)
      stdg <- cbind(stdg, sealphag)
      dfg <- dfg - 1
    }
    tg <- bg / stdg
    probtg <- 2 * (1 - pt(abs(tg), df = dfg))
    bg_stdg <- cbind(bg, stdg)
    cat("\nGlobal Parameter Estimates\n")
    print(data.frame(bg_stdg))
    cat("\n")
    print(data.frame(tg = tg, probtg = probtg))
    cat("\nNOTE: The denominator degrees of freedom for the t tests is", dfg, ".\n")
  }
  if (toupper(MODEL) == "POISSON") {
    yhatg <- exp(x1 %*% bg + offset)
    ll <- sum(-yhatg + y * log(yhatg) - lgamma(y + 1))
    AIC <- -2 * ll + 2 * nvar1
    AICc <- -2 * ll + 2 * nvar1 * (n / (n - nvar1 - 1))
    tt2 <- y / mean(y)
    tt2 <- ifelse(tt2 == 0, E^(-10), tt2)
    devnullg <- 2 * sum(y * log(tt2) - (y - mean(y)))
    pctdevg <- 1 - devg / devnullg
    adjpctdevg <- 1 - ((n - 1) / (n - nvar1)) * (1 - pctdevg)
    cat("\n")
    print(data.frame(
      devg = devg,
      ll = ll,
      pctdevg = pctdevg,
      adjpctdevg = adjpctdevg,
      AIC = AIC,
      AICc = AICc
    ))
  }
  if (toupper(MODEL) == "NEGBIN") {
    yhatg <- exp(x1 %*% bg[1:(nrow(bg) - 1)] + offset)
    ll <- sum(y * log(alphag * yhatg) - (y + 1 / alphag) * log(1 + alphag * yhatg) + lgamma(y + 1) - lgamma(1 / alphag) - lgamma(y + 1))
    AIC <- -2 * ll + 2 * (nvar1 + 1)
    AICc <- -2 * ll + 2 * (nvar1 + 1) * (n / (n - (nvar1 + 1) - 1))
    tt2 <- y / mean(y)
    tt2 <- ifelse(tt2 == 0, E^(-10), tt2)
    devnullg <- 2 * sum(y * log(tt2) - (y + 1 / alphag) * log((1 + alphag * y) / (1 + alphag * mean(y))))
    pctdevg <- 1 - devg / devnullg
    adjpctdevg <- 1 - ((n - 1) / (n - nvar1)) * (1 - pctdevg)
    cat("\n")
    print(data.frame(
      devg = devg,
      ll = ll,
      pctdevg = pctdevg,
      adjpctdevg = adjpctdevg,
      AIC = AIC,
      AICc = AICc
    ))
  }
  # /****************************************/
  if (is.null(GRID)) {
    res_ <- data.frame(var = c(wt, y, yhat, res, resstd, rsqri, influence, cooksD, sumwi))
    beta_ <- data.frame(id = bi[,1], B = bi[,2], x = bi[,3], y = bi[,4])
    parameters_ <- data.frame(id = bi[,1], B = bi[,2], x = bi[,3], y = bi[,4],
                              stdbi = stdbi[,1], tstat = tstat[,1], probt = probt[,1])
    
    if (toupper(MODEL) == "NEGBIN") {
      atstat <- alphai[,2] / alphai[,3]
      aprobtstat <- 2 * (1 - pnorm(abs(atstat)))
      siga_ <- ifelse(aprobtstat < 0.01 * (nvar / v1), "significant at 99%",
                      ifelse(aprobtstat < 0.05 * (nvar / v1), "significant at 95%",
                             ifelse(aprobtstat < 0.1 * (nvar / v1), "significant at 90%", "not significant at 90%")))
      
      alphai <- cbind(alphai, atstat, aprobtstat)
      alpha_ <- data.frame(x = alphai[,1], y = alphai[,2], id = alphai[,3],
                           alpha = alphai[,4], std = alphai[,5], tstat = alphai[,6], probt = alphai[,7])
      
      write.csv(alpha_, file = "_alpha_.csv", row.names = FALSE)
    }
  } else {
    beta_ <- data.frame(id = bi[,1], B = bi[,2], x = bi[,3], y = bi[,4])
    stdbi <- sqrt(varbi)
    tstat <- bi[,2] / stdbi
    tstat[is.na(tstat)] <- 0
    parameters_ <- data.frame(id = bi[,1], B = bi[,2], x = bi[,3], y = bi[,4],
                              stdbi = stdbi[,1], tstat = tstat)
    if (toupper(MODEL) == "NEGBIN") {
      atstat <- alphai[,2] / alphai[,3]
      atstat[alphai[,3] == 0] <- 0
      aprobtstat <- 2 * (1 - pnorm(abs(atstat)))
      alphai <- cbind(POINTS, alphai, atstat, aprobtstat)
      
      alpha_ <- data.frame(x = alphai[,1], y = alphai[,2], id = alphai[,3],
                           alpha = alphai[,4], std = alphai[,5], tstat = alphai[,6], probt = alphai[,7])
      
      write.csv(alpha_, file = "alpha_.csv", row.names = FALSE)
    }
  }
} #fecha GWNBR


# SUBSTITUICOES SAS -> R
## _h_ = h_
## h   = hh 
## _dist_ = distance
## seq = sequ
## dist = distan
## _w_ = w_
## _res_ = res_ 
## _beta_ = beta_
## _parameters_ = parameters_
## _siga_ = siga_
## _alpha_ = alpha_
## _beta_ = beta_
## _parameters_ = parameters_ 

# /*****************************************/ 


setwd('C:/Users/jehhv/OneDrive/Documentos/UnB/PIBIC/GWNBR')
library(readr)
GeorgiaData <- read_delim("GeorgiaData.csv", 
                          delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(GeorgiaData)
example2 <- read_csv("example2.csv")
example2 <- example2 %>% 
  mutate(Le = log(Loe))
View(example2)

system.time(golden(example2,YVAR="Mort2564", XVAR=c('Professl','Elderly','OwnHome','Unemply'), 
                   LAT="Y", LONG="X",MODEL="NEGBIN",OFFSET="Le",METHOD="FIXED_G",BANDWIDTH="AIC",GLOBALMIN='NO'))

GWNBR(example2,YVAR="Mort2564",XVAR=c('Professl','Elderly','OwnHome','Unemply'), 
      LAT="Y", LONG="X",MODEL="NEGBIN",OFFSET="Le",METHOD="FIXED_G", H=16780.349)

# --------------------------------------------------------------- # 

# Rodando com Georgia Data
georgia_data <- read_delim("GeorgiaData.csv", 
                          delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(georgia_data)

#Discretizando y e padronizando x
georgia_nb_std <- georgia_data
georgia_nb_std$PctBach <- as.integer(georgia_nb_std$PctBach)

#TotPop90
georgia_nb_std$TotPop90 <- (georgia_nb_std$TotPop90-mean(georgia_nb_std$TotPop90))/sd(georgia_nb_std$TotPop90)

#PctRural
georgia_nb_std$PctRural <- (georgia_nb_std$PctRural-mean(georgia_nb_std$PctRural))/sd(georgia_nb_std$PctRural)

#PctEld
georgia_nb_std$PctEld <- (georgia_nb_std$PctEld-mean(georgia_nb_std$PctEld))/sd(georgia_nb_std$PctEld)

#PctFB
georgia_nb_std$PctFB <- (georgia_nb_std$PctFB-mean(georgia_nb_std$PctFB))/sd(georgia_nb_std$PctFB)

#PctPov
georgia_nb_std$PctPov <- (georgia_nb_std$PctPov-mean(georgia_nb_std$PctPov))/sd(georgia_nb_std$PctPov)

#PctBlack
georgia_nb_std$PctBlack <- (georgia_nb_std$PctBlack-mean(georgia_nb_std$PctBlack))/sd(georgia_nb_std$PctBlack)


system.time(golden(georgia_nb_std,YVAR="PctBach", 
                   XVAR=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov", "PctBlack"), 
                   LAT="Y", LONG="X", MODEL="NEGBIN", METHOD="FIXED_G",
                   BANDWIDTH="AIC", GLOBALMIN='NO'))
