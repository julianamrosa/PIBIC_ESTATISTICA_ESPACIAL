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
# DISTANCEKM = if the distances between two locations will be computed in Km using the Basic formula for calculating spherical distance. 
# The default value is NO, so the distance is computed using euclidean distance.

golden2 <- function(y, x, xvarglobal, weight, lat, long, output, method, model=GAUSSIAN, bandwidth,offset=NULL,distancekm=NO){
  # distancekm, model e offset = default
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
      # resposta: remover adptiven
    }
    if(method=="fixed_g"|method=="fixed_bsq"|method=="adaptive_bsq"){
    }
    for(j in 1:n){
      seqi <- matrix(j,n,1)
      distan <- cbind(seqi, sequ, as.matrix(distance)[,j])
      if(distancekm=="yes"){
        distan[,3] <- distan[,3]*111
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
        posit <- which((!distance[,3]>h)==0)
        posit <- distance[which((!distance[,3]>3)==0)]
        w[posit] <- 0
      }
      if(method=="adaptive_bsq"){
        distance <- distance[order(distance[,3]),]
        distance <- cbind(distance,1:n)
        w <- matrix(0,n,2)
        hn <- distance[h,3]
        for(jj in 1:n){
          if(distance[jj,4]<=h){
            w[jj,1] <- (1-(distante[jj,3]/hn)^2)^2
          }
          else{
            w[jj,1] <- 0
            w[jj,2] <- distance[jj,2]
          }
          if(bandwidth=="cv"){
            w <- which(!w==0)
            w <- w[which(!w==0)]
            w <- 0
          }
          w <- w[order(w[,2]),]
          w <- w[,1]
        }
      } # fecha if linha 146
      if(model=="gaussian"){
        if (det(t(x)%*%(w*x*wt)*x)==0){
          b <- matrix(0,nvar,1)
        }
        else{
          b <- solve(t(x)%*%(w*x*wt)%*%t(x)%*%(w*y*wt))
        }
        if(method=="fixed_g"|method=="fixed_bsq"|method=="adaptive_bsq"){
          yhat[i] <- x[i,]%*%b
          if(det(t(x)%*%(w*x*wt)==0)){
            S[i] <- 0
          }
          else{
            S[i] <- (x[i,]%*%solve(t(x)%*%(w*x*wt))*t((x*w*wt)))[i]
          }
        if(is.null(xvarglobal)){
          hat_matrix <- rbind(hat_matrix,(x[i,]%*%solve(t(x)%*%(w*x*wt))*t((x*w*wt)))) #criando hat_matrix aqui, entao supostamente nao haveria problema com dimensão
          if(i==1){
            W_f <- cbind(matrix(i,n,1),w,seq(1,nrow(w)))
          }
          else{
            W_f <- cbind(W_f,c(rbind(matrix(i,n,1),w,t(seq(1,nrow(w)))))) #revisar linha
          }
          if(is.null(xvarglobal)){
            ba <- solve(t(t(xa)%*%diag(1,n,n)-hat_matrix))%*%(diag(1,n,n)-hat_matrix)%*%(xa*wt)%*%t(xa)%*%(t(diag(1,n,n)-hat_matrix))%*%(diag(1,n,n)-hat_matrix)*(y*wt) 
            ya <- y-xa%*%ba
            for(i in 1:n){
              m1 <- (i-1)%*%n+1
              m2 <- m1+(n-1)
              w <- W_f[m1:m2,2]
              if(det(t(x)%*%(w*x*wt))==0){
                b <- matrix(0,nvar,1)
              }
              else{
                b <- solve(t(x)%*%(w*x*wt))%*%t(x)*(w*ya*wt)
              }
              uj <- cbind(x,xa)*rbind(b,ba)
              yhat[i] <- uj[i]
            }
            S <- diag(xa%*%solve(t(xa)%*%t(diag(1,n,n)-hat_matrix)%*%(diag(1,n,n)-hat_matrix)%*%(xa*wt))%*%t((xa*wt))%*%t(diag(1,n,n)-hat_matrix)%*%(diag(1,n,n)-hat_matrix)+hat_matrix)
          }
          CV <- t((y-yhat)*wt)%*%(y-yhat)
          if(model=='poisson'){
            ll <- sum((-1)*yhat+y*log(yhat)-lgamma(y+1)) #confirmar se esse log nao é ln
            npar <- sum(S)
          }
          if(model=='negbin'){
            ll <- sum(y*log(alphai*yhat)-(y+1/alphai)*log(1+alphai*yhat)+lgamma(y+1/alphai)-lgamma(1/alphai)-lgamma(y+1))
            npar <- sum(S)+sum(S)/nvar
          }
          AIC <- 2*npar-2*11
          AICC <- AIC+(2*npar*(npar+1))/(n-npar-1)
          if(bandwidth=='AIC'){
            CV<- AICC
          }
        } #frcha linha 183 
      } #fecha linha 175 
    } #fecha if == gaussian linha 168
  } #fecha laço for na linha 133
  if(model=='logistic'){
    uj <- (y+y[,])/2 #revisar linha
    nj <- log(uj/(1-uj))
    dev <- 0; ddev <- 1; cont <- 0
    while(abs(ddev)>0.000001 & cont<100){
      cont <- cont+1
      uj <- ifelse(uj>E^100,E^100,uj)
      Ai <- uj*(1-uj)
      Ai <- ifelse(Ai<=0,E*(-5),Ai)
      zj <- nj+(y-uj)/Ai
      if (det(t(x)%*%(w*Ai*x*wt))==0){
        b <- matrix(0,nvar,1)
      }
      else{
        b <- solve(t(x)%*%(w*Ai*x*wt))*t(x)%*%(w*Ai*wt*zj)
      }
      nj <- x%*%b
      nj <- ifelse(nj>E^2,E^2,nj)
      uj <- exp(nj)/(1+exp(nj))
      olddev <- dev
      uj <- ifelse(uj<E^(-150),E^(-150),uj)
      tt <- y/uj
      tt < ifelse(tt==0,E^(-10),tt)
      uj <- ifelse(uj==1,0.99999,uj)
      tt2 <- ifelse(tt2==0,E^-10,tt2)
      dev <- 2*sum(y*log(tt))+(1-y)*log(tt2)
      if(cont>100){
        ddev <- 0.0000001
      }
      else{
        ddev <- dev-olddev
      }
    } #fecha while
    if(method=="fixed_g"|method=="fixed_bsq"|method=="adaptive_bsq"){
      yhat[i] <- uj[i]
      if(det(t(x)%*%(w*Ai*x*wt)==0)){
        S[i] <- 0
      }
      else{
        S[i] <- (x[i,]%*%solve(t(x)%*%(w*Ai*x*wt))*t((x*w*wt*Ai)))[i]
      }
      if(is.null(xvarglobal)){
        if(i==1){
          W_f <- cbind(matrix(i,n,1),w,t(seq(1,nrow(w))))
        }
        else{
          W_f <- rbind(W_f,c(cbind(matrix(i,n,1),w,t(seq(1,nrow(w))))))
        }
      }
      if(is.null(xvarglobal)){
        uj <- (y+y[,])/2
        nj <- log(uj/(1-uj))
        dev <- 0; ddev <- 1; cont <-1
      }
      while(abs(ddev)>0.000001 & cont<100){
        cont <- cont+1
        uj <- ifelse(uj>E^100,E^100,uj)
        Aa <- uj*(1-uj)
        Aa <- ifelse(Aa <= 0,E^-5,Aa)
        za <- nj+(y-uj)/Aa
        if(det(t(xa)%*%(Aa*xa*wt))==0){
          ba <- matrix(0,nvar,1)
        }
        else{
          ba <- solve(t(xa)%*%(Aa*xa*wt))%*%t(xa)*(Aa*za*wt)
        }
        nj <- xa%*%ba
        nj <- ifelse(nj>E^2,E^2,nj)
        uj <- exp(nj)/(1+exp(nj))
        olddev <- dev
        uj <- ifelse(uj<E^-150,E^-150,uj)
        tt <- y/uj
        tt <- ifelse(tt=0,E^-10,tt)
        uj <- ifelse(uj==1,0.99999,uj)
        tt2 <- (1-y)/(1-uj)
        tt2 <- ifelse(tt2==0,E^-10,tt2)
        dev <- 2*sum((y*log(tt))+(1-y)*log(tt2))
        ddev <- dev-olddev
      } #fecha while
      diffba <- 1; contb <- 1
    } #fecha method
    while(abs(diffba)>0.00001 & contb<50){
      oldba <- ba
    for(i in 1:n) {
      m1 <- (i-1)*n+1
      m2 <- m1+(n-1)
      w <- w_f[m1:m2,2] 
      dev <- 0; ddev <- 1; cont2<- 0
      while(abs(ddev)>0.000001 & cont2<100){
        uj <- ifelse(uj>E^100,E^100,uj)
        Ai <- uj*(1-uj)
        Ai <- ifelse(Ai<=0,E^-5,Ai)
        zj <- nj+(y-uj)/Ai-xa%*%ba
        if(det(t(x)%*%(w*Ai*x*wt))==0){
          b <- matrix(0,nvar,1) 
        }
        else{
          b <- solve(t(x)%*%(w*Ai*x*wt))%*%t(x)%*%(w*Ai*zj*wt)
        }
        nj <- x%*%b+xa%*%ba
        nj <- ifelse(nj>E^2,E^2,nj)
        uj <- exp(nj)/(1+exp(nj))
        olddev <- dev
        uj <- ifelse(uj<E^-150,E^-150,uj)
        tt <- y/uj
        tt <- ifelse(tt==0,E^-10,tt)
        uj <- ifelse(uj==1,0.99999,uj)
        tt2 <- (1-y)/(1-uj)
        tt2 <- ifelse(tt2=0,1E-10,tt2)
        dev <- 2*sum((y*log(tt))+(1-y)*log(tt2))
        ddev <- dev-olddev
        cont2 <- cont2+1
      } #fecha while
      if(det(t(x)%*%(w*Ai*x*wt))){
        C <- matrix(0,nvar,n)
      }
      else{
        C <- solve(t(x)%*%(w*Ai*x*wt))*t(x)*t(w*Ai*wt)
      }
      R_matrix <- rbind(R_matrix,(x[i,]*C))
      Z_matrix <- rbind(Z_matrix,(zj+xa*ba))
      yhat[i] <- uj[i]
    } # fecha o laço for
    hat_matrix <- R_matrix*Z_matrix/diag(Z_matrix)
    uj <- yhat
    nj <- log(uj/(1-uj))
    Aa <- uj*(1-uj)
    Aa <- ifelse(Aa<=0,E^-5,Aa)
    za <- nj+(y-uj)/Aa
    if (det(t(xa%*%Aa)%*%(diag(1,n,n)-hat_matrix)*(xa*wt))==0){
      ba <- matrix(0,ncol(xa),1)
    }
    else{
      ba <- solve(t(xa*Aa)%*%(diag(1,n,n)-hat_matrix))%*%t(xa*Aa)%*%(diag(1,n,n)-hat_matrix)%*%(za*wt)
    }
    nj <- xa*ba
    nj <- ifelse(nj>E^2,E^2,nj)
    uj <- exp(nj)/(1+exp(nj))
    diffba <- oldba-ba
    contb <- contb+1
    } #fecha while (diffba)
    S <- diag((diag(1,n,n)-hat_matrix)%*%xa%*%solve(t(xa*Aa)%*%(diag(1,n,n)-hat_matrix)%*%(xa*wt))%*%t(xa*Aa*wt)%*%(diag(1,n,n)-hat_matrix)+hat_matrix)
    uj <- ifelse(uj==0,E^-10,uj)
    uj <- ifelse(uj==1,0.99999,uj)
    CV <- t((y-yhat)*wt)%*%(y-yhat)
    ll <- sum(y*log(uj)-(1-y)*log(1-uj))
    npar <- sum(S)
    AIC <-  2*npar-2*ll
    AICC <-  AIC +(2*npar*(npar+1))/(n-npar-1)
  } # fecha modelo logístico  
    
  } # fecha CV
} # fecha golden

# trocas ----
# position ---- posit
# _dist_   ---- distance
# dist     ---- distan
teste <- which(!exemplo==0)           #retorna apenas indices
teste <- exemplo[which(!exemplo==0)]  #retorna os valores


# nao fazer para adaptiven





























