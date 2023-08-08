library(readr)
library(dplyr)
library(tidyverse)

golden <- function(DATA,YVAR, XVAR, XVARGLOBAL=NULL, WEIGHT=NULL, LAT, LONG, 
                   GLOBALMIN='YES', METHOD, MODEL="NEGBIN", BANDWIDTH='CV',
                   OFFSET=NULL,DISTANCEKM="NO"){
  # distancekm, model e offset = default
  E <- 10
  y <<- as.numeric(unlist(DATA[,YVAR])) #the name of the dependent or response variable
  x <<- DATA[,XVAR] # the name of the independent or explicative variables. 
  N <<- length(y)
  wt <<- rep(1, N)
  if (!is.null(WEIGHT)){
    wt <<- as.matrix(WEIGHT)
  }
  Offset <<- rep(0, N)
  if (!is.null(OFFSET)){
    Offset <<- as.numeric(unlist(DATA[,OFFSET]))
  }
  x <<- as.matrix(cbind(rep(1,N),x))
  nvar <<- ncol(x)
  if (!is.null(XVARGLOBAL)){
    xa <<- as.matrix(XVARGLOBAL)
  }
  yhat <<- matrix(0,N,1) #aqui nao poderia ser rep(0,N)? 
  alphai <<- matrix(0, N, 1)  #aqui nao poderia ser rep(0,N)? 
  S <<- matrix(0,N,1)  #aqui nao poderia ser rep(0,N)? 
  
  # global estimates #
  uj <- (y+mean(y))/2
  nj <- log(uj)
  parg <<- sum((y-uj)^2/uj)/(N-nvar)
  ddpar <- 1
  cont <- 1
  while (abs(ddpar)>0.000001 & cont<100){
    dpar <- 1
    parold <- parg
    cont1 <- 1
    if (toupper(MODEL)=="POISSON"){
      alphag <- E^-6
      parg <<- 1/alphag
    }
    if (toupper(MODEL)=="NEGBIN"){
      if (cont>1){
        parg <<- 1/(sum((y-uj)^2/uj)/(N-nvar))
      }
      while (abs(dpar)>0.000001 & cont1<200){
        parg <<- ifelse(parg<E^-10,E^-10,parg)
        g <- sum(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj))
        hess <- sum(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/((parg+uj)^2))
        hess <- ifelse(hess==0, E^-23,hess)
        par0 <- parg
        parg <<- par0-solve(hess)%*%g
        dpar <- parg-par0
        cont <- cont1+1
        if(parg>E^6){
          parg <<- E^6
          dpar <- 0
        } else{
          
        }
        parg <<- as.numeric(parg)
      }
      alphag <- 1/parg
    }
    devg <- 0
    ddev <- 1
    cont2 <- 0
    while (abs(ddev)>0.000001 & cont2<100){
      uj <- ifelse(uj>E^100,E^100,uj)
      Ai <- (uj/(1+alphag*uj))+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj*uj))
      Ai <- ifelse(Ai<=0,E^-5,Ai)
      Ai <- as.numeric(Ai)
      zj <- nj+(y-uj)/(Ai*(1+alphag*uj))-Offset
      if (det(t(x)%*%(Ai*x))==0){
        bg <- rep(0,nvar)
      }
      else{
        bg <- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
      }
      nj <- x%*%bg+Offset
      nj <- ifelse(nj>E^2,E^2,nj)
      uj <- as.numeric(exp(nj)) #flag
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
    ujg <<- uj
    cont <- cont+1
    ddpar <- parg-parold
  } #fecha while
  #View(GeorgiaData)
  LONG <- DATA[, LONG]
  LAT  <- DATA[, LAT]
  #COORD <<- matrix(c(GeorgiaData$Longitud,GeorgiaData$Latitude),ncol=2) #tentar outra solucao
  #View(COORD)
  COORD <<- matrix(c(LONG, LAT), ncol=2)
  distance <<- dist(COORD,"euclidean")
  #View(distance)
  sequ <<- 1:N
  #nn <- c() #solucao questionavel
  cv <- function(h){#(h, wt=Wt, nn=N, x=xx, xa=xa, y=yy,
    # ujg=Ujg, yhat=Yhat, nvar=Nvar, hv=hv,
    # coord=COORD, distance=Distance, sequ=Sequ,
    # offset=Offset, Alpha=alphag, alphai=Alphai,
    # S0=S, Parg=parg){
    for (i in 1:N){ 
      for(j in 1:N){
        seqi <- rep(i,N)
        distan <<- cbind(seqi, sequ, as.matrix(distance)[,i])
        if(toupper(DISTANCEKM)=="YES"){
          distan[,3] <<- distan[,3]*111
        }
      }
      u <- nrow(distan)
      w <- rep(0,u)
      for(jj in 1:u){
        w[jj] <- exp(-0.5*(distan[jj,3]/h)^2)
        if(toupper(BANDWIDTH)=="CV"){
          w[i] <-0
        }
      }
      if(toupper(METHOD)=="FIXED_BSQ"){ #juntei com if anterior
        w[jj] <- (1-(distan[jj,3]/h)^2)^2
        position <- which(distan[,3]>h) 
        w[position] <- 0
      }
      else if(toupper(METHOD)=="ADAPTIVE_BSQ"){
        distan <<- distan[order(distan[,3]),] 
        distan <<- cbind(distan,1:nrow(distan))
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
      if(toupper(MODEL)=="GAUSSIAN"){
        if (det(t(x)%*%(w*x*wt)*x)==0){
          b <- rep(0,nvar)
        }
        else{
          b <- solve(t(x)%*%(w*x*wt))%*%t(x)%*%(w*y*wt)
        }
        if(toupper(METHOD)=="FIXED_G"|toupper(METHOD)=="FIXED_BSQ"|toupper(METHOD)=="ADAPTIVE_BSQ"){
          yhat[i] <<- x[i,]*b
          if(det(t(x)%*%(w*x*wt)==0)){
            S[i] <<- 0
          }
          else{
            S[i] <<- (x[i,]%*%solve(t(x)%*%(w*x*wt))*t(x*w*wt))[i]
          }
          if(!is.null(XVARGLOBAL)){
            hat_matrix <- rbind(hat_matrix,(x[i,]%*%solve(t(x)%*%(w*x*wt))*t((x*w*wt)))) #criando hat_matrix aqui, entao supostamente nao haveria problema com dimens?o
            if(i==1){
              W_f <- cbind(matrix(i,N,1),w,t(seq(1,nrow(w)))) 
            }
            else{
              W_f <- rbind(W_f,c(cbind(matrix(i,N,1),w,t(seq(1,nrow(w)))))) 
            }
          }  
          if(!is.null(XVARGLOBAL)){ 
            ba <- solve(t(t(xa)%*%diag(1,N,N)-hat_matrix)%*%(diag(1,N,N)-hat_matrix)%*%(xa*wt))%*%t(xa)%*%(t(diag(1,N,N)-hat_matrix))%*%(diag(1,N,N)-hat_matrix)*(y*wt) 
            ya <- y-xa%*%ba
            for(i in 1:N){
              m1 <- (i-1)*N+1
              m2 <- m1+(N-1)
              w <- W_f[m1:m2,2]
              if(det(t(x)%*%(w*x*wt))==0){
                b <- matrix(0,nvar,1)
              }
              else{
                b <- solve(t(x)%*%(w*x*wt))%*%t(x)*(w*ya*wt)
              }
              uj <- cbind(x,xa)*rbind(b,ba)
              yhat[i] <<- uj[i]
            }
            S <<- diag(xa%*%solve(t(xa)%*%t(diag(1,N,N)-hat_matrix)%*%(diag(1,N,N)-hat_matrix)%*%(xa*wt))%*%t((xa*wt))%*%t(diag(1,N,N)-hat_matrix)%*%(diag(1,N,N)-hat_matrix)+hat_matrix)
          }
          CV <- t((y-yhat)*wt)%*%(y-yhat)
        }
      }
      if(toupper(MODEL)=='POISSON'|toupper(MODEL)=='NEGBIN'){
        uj <- ujg
        nj <- log(uj)
        par <- parg
        ddpar <- 1
        cont <- 1
        while(abs(ddpar)>0.000001 & cont<100){
          dpar <- 1
          parold <-par
          cont1 < -1
          if(toupper(MODEL)=='POISSON'){
            alpha <- E^-6
            par <-1/alpha
          }
          if(toupper(MODEL)=='NEGBIN'){
            if(par <= E^-5){
              if(i>1){
                par <- 1/alphai[i-1]
              }
            }
            while(abs(dpar)>0.000001 & cont1<200){
              par <- ifelse(par<E^-10, E^-10, par)
              g <- sum(w*wt*(digamma(par+y)-digamma(par)+log(par)+1-log(par+uj)-(par+y)/(par+uj)))
              hess <- sum(w*wt*(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)^2))
              hess <- ifelse(hess==0, E^-23, hess)
              par0 <- par
              par <- par0-solve(hess)%*%g
              dpar <- par-par0
              cont1 <- cont1+1
              if(par>E^6){
                par <- E^6
                dpar <- 0
              }
              par <- as.numeric(par)
            }
            alpha <- 1/par
          } #fecha model == negbin
          dev <- 0
          ddev <- 1
          cont2 <- 1
          while(abs(ddev)>0.000001 & cont2<100){
            uj <- ifelse(uj>E^100, E^100, uj)
            Ai <- (uj/(1+alpha*uj))+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj*uj))
            Ai <- ifelse(Ai<=0, E^-5, Ai)
            Ai <- as.numeric(Ai)
            zj <- nj+(y-uj)/(Ai*(1+alpha*uj))-Offset
            if (det(t(x)%*%(w*Ai*x*wt))==0){
              b <- rep(0, nvar)
            }
            else{
              b <- solve(t(x)%*%(w*Ai*x*wt))%*%t(x)%*%(w*Ai*wt*zj)
            }
            nj <- x%*%b+Offset
            nj <- ifelse(nj>E^2, E^2, nj)
            uj <- exp(nj)
            olddev <- dev
            uj <- ifelse(uj<E^-150, E^-150, uj)
            tt <- y/uj
            tt <- ifelse(tt==0, E^-10, tt)
            if (toupper(MODEL)=="POISSON"){
              dev <- 2*sum(y*log(tt)-(y-uj))
            }
            else{
              dev <- 2*sum(y*log(tt)-(y+1/alpha)*log((1+alpha*y)/(1+alpha*uj)))
            }
            if (cont2>100){
              ddev <- 0.0000001
            }
            else{
              ddev <- dev-olddev
            }
            cont2 <- cont2+1
          }
          cont <- cont+1
          ddpar <- par-parold
        }
      }
      if(toupper(METHOD)=='FIXED_G'|toupper(METHOD)=='FIXED_BSQ'|toupper(METHOD)=='ADAPTIVE_BSQ'){
        yhat[i] <<- uj[i]
        alphai[i] <<- alpha
        if (det(t(x)%*%(w*Ai*x*wt))==0){
          S[i] <<- 0
        }
        else{
          S[i] <<- (x[i, ]%*%solve(t(x)%*%(w*Ai*x*wt))%*%t(x*w*Ai*wt))[i]
        }   
        if(!is.null(XVARGLOBAL)){
          if(i==1){
            W_f <- cbind(matrix(i,N,1),w,t(seq(1:nrow(w))))
          }
          else{
            W_f <- rbind(W_f,c(cbind(matrix(i,N,1),w,t(seq(1,nrow(w))))))
          }
        }
        if(!is.null(XVARGLOBAL)){
          uj <- (y+mean(y))/2 
          nj <- log(uj)
          if(toupper(MODEL)=='POISSON'){
            alphag <- E^-6
          }
          devga <- 0
          ddev <- 1
          cont2 <- 0
          while(abs(ddev)>0.000001 & cont2<100){
            uj <- ifelse(uj>E^100,E^100,uj)
            Aa <- (uj/(1+alphag*uj))+(y-uj)*(alphag*uj/1+2*alphag*uj+alphag^2*uj*uj)
            Aa <- ifelse(Aa<=0,E^-5,Aa)
            za <- nj+(y-uj)/(Aa*(1+alphag*uj))-Offset
            if (det(t(xa)%*%(Aa*xa))==0){
              ba <- rep(0,nvar)
            }
            else{
              bg <- solve(t(xa)%*%(Aa*xa))%*%t(xa)%*%(Aa*za)
            }
            nj <- xa%*%ba+Offset
            nj <- ifelse(nj>E^2,E^2,nj)
            uj <- exp(nj)
            olddev <- devga
            uj <- ifelse(uj<E^-150,E^-150,uj)
            tt <- y/uj
            tt <- ifelse(tt==0,E^-10,tt)
            if(toupper(MODEL)=="POISSON"){
              devga <- 2*sum(y*log(tt)-(y-uj))
            }
            if(toupper(MODEL)=="NEGBIN"){
              devga <- 2*sum(y*log(tt)-(y+1/alphag)*log((1+alphag*y)/(1+alphag*uj)))
            }
            ddev <- devga-olddev
            cont2 <- cont2+1
          }
          diffba <- 1
          contb <- 1
          while(abs(diffba)>0.00001 & contb<50){
            oldba <- ba
            for(i in 1:N){
              m1 <- (i-1)*N+1
              m2 <- m1+(N-1)
              w <- w_f[m1:m2,2]
              if(toupper(MODEL)=='POISSON'){
                alpha <- E^-6
              }
              if(toupper(MODEL)=='NEGBIN'){
                par <- 1/alphai[i]
                if(par<=E^-5){
                  if(i>1){
                    par <- 1/alphai[i-1]
                  }
                }
                while(abs(dpar)>0.000001 & cont1<200){
                  par <- ifelse(par<E^-10, E^-10, par)
                  g <- sum(w*wt*(digamma(par+y)-digamma(par)+log(par)+1-log(par+uj)-(par+y)/(par+uj)))
                  hess  <- sum(w*wt*(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)^2))
                  hess <- ifelse(hess==0, E^-23, hess)
                  par0 <- par
                  par <- par0-solve(hess)%*%g
                  dpar <- par-par0
                  cont1 <- cont1+1
                  if(par>E^6){
                    par <- E^6
                    dpar <- 0
                  }
                }
                alpha <- 1/par
              } #fecha model == negbin
              dev <- 0
              ddev <- 1
              cont2 <- 0
              while(abs(ddev)>0.000001 & cont2<100){
                uj <- ifelse(uj>E^100, E^100, uj)
                Ai <- (uj/(1+alpha*uj))+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj*uj))
                Ai <- ifelse(Ai<=0, E^-5, Ai)
                zj <- nj+(y-uj)/(Ai*(1+alpha*uj))-Offset-xa*ba
                if (det(t(x)%*%(w*Ai*x*wt))==0){
                  b <- rep(0, nvar)
                }
                else{
                  b <- solve(t(x)%*%(w*Ai*x*wt))%*%t(x)%*%(w*Ai*wt*zj)
                }
                nj <- x*b+Offset+xa*ba
                nj <- ifelse(nj>E^2, E^2, nj)
                uj <- exp(nj)
                olddev <- dev
                uj <- ifelse(uj<E^-150, E^-150, uj)
                tt <- y/uj
                tt <- ifelse(tt==0, E^-10, tt)
                if (toupper(MODEL)=="POISSON"){
                  dev <- 2*sum(y*log(tt)-(y-uj))
                }
                else{
                  dev <- 2*sum(y*log(tt)-(y+1/alpha)*log((1+alpha*y)/(1+alpha*uj)))
                }
                ddev <- dev-olddev
                cont2 <- cont2+1
              }
              C <- solve(t(x)%*%(w*Ai*x*wt))*t(x)*t(w*Ai*wt)
              R_matrix <- rbind(R_matrix,(x[i,]*C))
              Z_matrix <- cbind(Z_matrix,(zj+Offset+xa*ba))
              yhat[i] <<- uj[i]
            } #fecha laco for
            hat_matrix <- R_matrix%*%Z_matrix/diag(Z_matrix)
            uj <- yhat
            nj <- log(uj)
            Aa <- (uj/(1+alphag*uj))+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj*uj))
            Aa <- ifelse(Aa<=0,E^-5,Aa)
            za <- nj+(y-uj)/(Aa*(1+alphag*uj))-Offset
            ba <- solve(t(xa*Aa)%*%(diag(1,N,N)-hat_matrix))%*%t(xa*wt)%*%t(xa*Aa)%*%(diag(1,N,N)-hat_matrix)%*%(za*wt)
            nj <- xa*ba+Offset
            nj <- ifelse(nj>E^2,E^2,nj)
            uj <- exp(nj)
            diffba <- oldba-ba
            contb <- contb+1
          }
          S <<- diag((diag(1,N,N)-hat_matrix)%*%xa%*%solve(t(xa*Aa)%*%(diag(1,N,N)-hat_matrix)%*%(xa*wt))%*%t(xa*Aa*wt)%*%(diag(1,N,N)-hat_matrix)+hat_matrix)
        } #fecha if varglobal
        CV <- t((y-yhat)*wt)%*%(y-yhat)
        if(toupper(MODEL)=='POISSON'){
          ll <- sum(-yhat+y*log(yhat)-lgamma(y+1))
          npar <- sum(S)
        }
        if(toupper(MODEL)=='NEGBIN'){
          ll <- sum(y*log(alphai*yhat)-(y+1/alphai)*log(1+alphai*yhat)+lgamma(y+1/alphai)-lgamma(1/alphai)-lgamma(y+1))
          npar <- sum(S)+sum(S)/nvar
        }
        AIC <- 2*npar-2*11
        AICC <- AIC + (2*npar*(npar+1))/(N-npar-1)
        if(toupper(BANDWIDTH)=='AIC'){
          CV <- AICC
        }
      }
      if(toupper(MODEL)=='LOGISTIC'){
        uj <- (y+mean(y))/2 #revisar linha
        nj <- log(uj/(1-uj))
        dev <- 0; ddev <- 1; cont <- 0
        while(abs(ddev)>0.000001 & cont<100){
          cont <- cont+1
          uj <- ifelse(uj>E^100,E^100,uj)
          Ai <- uj*(1-uj)
          Ai <- ifelse(Ai<=0,E*(-5),Ai)
          zj <- nj+(y-uj)/Ai
          if (det(t(x)%*%(w*Ai*x*wt))==0){
            b <- rep(0,nvar)
          }
          else{
            b <- solve(t(x)%*%(w*Ai*x*wt))%*%t(x)%*%(w*Ai*wt*zj)
          }
          nj <- x*b
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
        if(toupper(METHOD)=="FIXED_G"|toupper(METHOD)=="FIXED_BSQ"|toupper(METHOD)=="ADAPTIVE_BSQ"){
          yhat[i] <<- uj[i]
          if(det(t(x)%*%(w*Ai*x*wt)==0)){
            S[i] <<- 0
          }
          else{
            S[i] <<- (x[i,]%*%solve(t(x)%*%(w*Ai*x*wt))*t((x*w*wt*Ai)))[i]
          }
          if(!is.null(XVARGLOBAL)){
            if(i==1){
              W_f <- cbind(matrix(i,N,1),w,t(seq(1,nrow(w))))
            }
            else{
              W_f <- rbind(W_f,c(cbind(matrix(i,N,1),w,t(seq(1,nrow(w))))))
            }
          }
          if(!is.null(XVARGLOBAL)){
            uj <- (y+mean(y))/2
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
              ba <- rep(0,nvar)
            }
            else{
              ba <- solve(t(xa)%*%(Aa*xa*wt))%*%t(xa)%*%(Aa*za*wt)
            }
            nj <- xa%*%ba
            nj <- ifelse(nj>E^2,E^2,nj)
            uj <- exp(nj)/(1+exp(nj))
            olddev <- dev
            uj <- ifelse(uj<E^-150,E^-150,uj)
            tt <- y/uj
            tt <- ifelse(tt==0,E^-10,tt)
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
          for(i in 1:N) {
            m1 <- (i-1)*N+1
            m2 <- m1+(N-1)
            w <- w_f[m1:m2,2] 
            dev <- 0; ddev <- 1; cont2 <- 0
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
              C <- matrix(0,nvar,N)
            }
            else{
              C <- solve(t(x)%*%(w*Ai*x*wt))*t(x)*t(w*Ai*wt)
            }
            R_matrix <- rbind(R_matrix,(x[i,]*C))
            Z_matrix <- rbind(Z_matrix,(zj+xa*ba))
            yhat[i] <<- uj[i]
          } # fecha o laco for
          hat_matrix <- R_matrix%*%Z_matrix/diag(Z_matrix)
          uj <- yhat
          nj <- log(uj/(1-uj))
          Aa <- uj*(1-uj)
          Aa <- ifelse(Aa<=0,E^-5,Aa)
          za <- nj+(y-uj)/Aa
          if (det(t(xa%*%Aa)%*%(diag(1,N,N)-hat_matrix)*(xa*wt))==0){
            ba <- matrix(0,ncol(xa),1)
          }
          else{
            ba <- solve(t(xa*Aa)%*%(diag(1,N,N)-hat_matrix))%*%t(xa*Aa)%*%(diag(1,N,N)-hat_matrix)%*%(za*wt)
          }
          nj <- xa*ba
          nj <- ifelse(nj>E^2,E^2,nj)
          uj <- exp(nj)/(1+exp(nj))
          diffba <- oldba-ba
          contb <- contb+1
        } #fecha while (diffba)
        S <<- diag((diag(1,N,N)-hat_matrix)%*%xa%*%solve(t(xa*Aa)%*%(diag(1,N,N)-hat_matrix)%*%(xa*wt))%*%t(xa*Aa*wt)%*%(diag(1,N,N)-hat_matrix)+hat_matrix)
        uj <- ifelse(uj==0,E^-10,uj)
        uj <- ifelse(uj==1,0.99999,uj)
        CV <- t((y-yhat)*wt)%*%(y-yhat)
        ll <- sum(y*log(uj)-(1-y)*log(1-uj))
        npar <- sum(S)
        AIC <-  2*npar-2*ll
        AICC <-  AIC +(2*npar*(npar+1))/(N-npar-1)
        if(bandwidth=="AIC"){
          CV <- AICC
        } 
      } # fecha modelo logistico  
      res <- cbind(CV,npar)  
    } # fecha for da linha 122
    return(res)
  } # fecha CV
  
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
  r <<- 0.61803399
  tol <<- 0.1
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
  #revisar identacao 
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
    var_ <<- data.frame()
    if(GMY==1){
      var_ <<- c(GMY,h1,CV1,h2,CV2) #revisar aqui
      names(var_) <<- c('GMY','h1','CV1','h2','CV2')
      View(var_)
    }
  }  
  int <- 1
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
    h_<<- min(min_bandwidth[,2]) #proc sql
  }
} #fecha golden

# TESTES #

setwd('C:/Users/jehhv/OneDrive/Documentos/UnB/PIBIC/GWNBR')
example2 <- read.csv("example2.csv")
example2 <- example2 %>% 
  mutate(Le = log(Loe))

system.time(golden(example2,YVAR="Mort2564", XVAR=c('Professl','Elderly','OwnHome','Unemply'), 
       LAT="Y", LONG="X",MODEL="NEGBIN",OFFSET="Le",METHOD="FIXED_G",BANDWIDTH="AIC",GLOBALMIN = 'NO'))


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

GWNBR <- function(DATA, YVAR, XVAR, XVARGLOBAL, XVARINF, WEIGHT=NULL, LAT, LONG,
                GRID=NULL, DHV, METHOD, MODEL="NEGBIN",OFFSET=NULL,
                DISTANCEKM="NO", H=NULL){
  E <- 10
  y <- DATA[,YVAR]
  x <- DATA[,which(names(DATA) %in% XVAR)]
  N <- nrow(y)
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
  x1 <- cbind(x,xa)
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
      Ai <- (uj/(1+alphag*uj))+(y-uj)*(alphag*uj/1+2*alphag*uj+alphag^2*uj*uj)
      Ai <- ifelse(Ai<=0,E^-5,Ai)
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
      else(toupper(METHOD=="FIXED_BSQ")){
        w[jj] <<- (1-(distan[jj,3]/h)^2)^2
      }
      if(is.null(GRID)){#confirmar com Alan se é is.null mesmo
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
        uj <- choose(uj > 1E100, 1E100, uj)
        Ai <- (uj / (1 + alpha * uj)) + (y - uj) * (alpha * uj / (1 + 2 * alpha * uj + alpha^2 * uj * uj))
        Ai <- choose(Ai <= 0, 1E-5, Ai)
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
  
    
  }
  
} #fecha GWNBR
  
  
  # SUBSTITUICOES SAS -> R
  ## _h_ = h_
  ## h   = hh 
  ## _dist_ = distance
  ## seq = sequ
  ## dist = distan
  ## _w_ = w_
  
  
  # %ELSE %DO;
  # %IF %UPCASE(&METHOD)=FIXED_G or %UPCASE(&METHOD)=FIXED_BSQ or %UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN %DO;
  # h=&_h_;
  # print h[label="Bandwidth"];
  # %END;
  # %END;
  
  
  
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


# #golden(DATA=nakaya,YVAR=Mort2564,XVAR=c('Professl' 'Elderly' 'OwnHome' 'Unemply'),
#   LONG=x,LAT=Y,OUTPUT=band,MODEL=POISSON,OFFSET=Le,METHOD=FIXED_G,
#   BANDWIDTH=AIC)

mgwnbr(DATA=logistic_std, YVAR="Degree",
       XVAR=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="FIXED_G",
       BANDWIDTH="aic", MODEL="LOGISTIC")
