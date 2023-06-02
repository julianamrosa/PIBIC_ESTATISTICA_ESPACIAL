library(readr)
library(dplyr)

golden <- function(DATA,YVAR, XVAR, XVARGLOBAL, WEIGHT=NULL, LAT, LONG, 
                    GLOBALMIN='YES', METHOD, MODEL="NEGBIN", BANDWIDTH='CV',
                    OFFSET=NULL,DISTANCEKM="NO"){
  # distancekm, model e offset = default
  E <- 10
  yy <- DATA[,YVAR]
  xx <- DATA[,XVAR]
  N <<- nrow(yy)
  Wt <-rep(1, N)
  if (!is.null(WEIGHT)){
    Wt <- as.matrix(WEIGHT)
  }
  Offset <- rep(0, N)
  if (!is.null(OFFSET)){
    Offset <- as.matrix(OFFSET)
  }
  xx <- as.matrix(cbind(rep(1,N),xx))
  Nvar <- ncol(xx)
  if (!is.null(XVARGLOBAL)){
    xa <- as.matrix(XVARGLOBAL)
  }
  Yhat <- matrix(0,N,1)
  Alphai <- matrix(0, N, 1)
  S <- matrix(0,N,1)
  
  # global estimates #
  uj <- (yy+mean(yy))/2
  nj <- log(uj)
  parg <- sum((yy-uj)^2/uj)/(N-Nvar)
  ddpar <- 1
  cont <- 1
  while (abs(ddpar)>0.000001 & cont<100){
    dpar <- 1
    parold <- parg
    cont1 <- 1
    if (toupper(MODEL)=="POISSON"){
      alphag <- E^-6
      parg <- 1/alphag
    }
    if (toupper(MODEL)=="NEGBIN"){
      if (cont>1){
        parg <- 1/(sum((yy-uj)^2/uj)/(N-Nvar))
      }
      while (abs(dpar)>0.000001 & cont1<200){
        parg <- ifelse(parg<E^-10,E^-10,parg)
        g <- sum(digamma(parg+yy)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+yy)/(parg+uj))
        hess <- sum(trigamma(parg+yy)-trigamma(parg)+1/parg-2/(parg+uj)+(yy+parg)/((parg+uj)^2))
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
      while (abs(ddev)>0.000001 & cont2<100){
        uj <- ifelse(uj>E^100,E^100,uj)
        Ai <- (uj/(1+alphag*uj))+(yy-uj)*(alphag*uj/1+2*alphag*uj+alphag^2*uj*uj)
        Ai <- ifelse(Ai<=0,E^-5,Ai)
        zj <- nj+(yy-uj)/(Ai*(1+alphag*uj))-Offset
        if (det(t(xx)%*%(Ai*xx))==0){
          bg <- rep(0,Nvar)
        }
        else{
          bg <- solve(t(xx)%*%(Ai*xx))%*%t(xx)%*%(Ai*zj)
        }
        nj <- xx%*%bg+Offset
        nj <- ifelse(nj>E^2,E^2,nj)
        uj <- exp(nj)
        olddev <- devg
        uj <- ifelse(uj<E^-150,E^-150,uj)
        tt <- yy/yj
        tt <- ifelse(tt==0,E^-10,tt)
        if(toupper(MODEL)=="POISSON"){
          devg <- 2*sum(y*log(tt)-(yy-uj))
        }
        if(toupper(MODEL)=="NEGBIN"){
          devg <- 2*sum(yy*log(tt)-(yy+1/alphag)*log((1+alphag*yy)/(1+alphag*uj)))
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
    } #fecha while
  View(GeorgiaData)
  LONG <- DATA[, LONG]
  LAT  <- DATA[, LAT]
  COORD <- matrix(c(GeorgiaData$Longitud,GeorgiaData$Latitude),ncol=2) #tentar outra solucao
  View(COORD)
  Distance <- dist(COORD,"euclidean")
  View(Distance)
  Sequ <- 1:N
  nn <- c() #solucao questionavel 
  cv <- function(h, wt=Wt, nn=N, x=xx, xa=xa, y=yy,
                 ujg=Ujg, yhat=Yhat, nvar=Nvar, hv=hv,
                 coord=COORD, distance=Distance, sequ=Sequ,
                 offset=Offset, Alpha=alphag, alphai=Alphai,
                 S0=S, Parg=parg){
    for (i in 1:nn){ 
      for(j in 1:nn){
        seqi <- rep(i,nn)
        distan <<- cbind(seqi, sequ, as.matrix(distance)[,i])
        if(toupper(DISTANCEKM)=="YES"){
          distan[,3] <<- distan[,3]*111
        }
      }
      u <- nrow(distan)
      w <- rep(0,u)
      for(jj in 1:u){
        w[jj] <- exp(-0.5*(distan[jj,3]/h)^2)
        print(w)
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
        for(jj in 1:nn){
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
          yhat[i] <- x[i,]*b
          if(det(t(x)%*%(w*x*wt)==0)){
            S[i] <- 0
          }
          else{
            S[i] <- (x[i,]%*%solve(t(x)%*%(w*x*wt))*t(x*w*wt))[i]
          }
          if(!is.null(xvarglobal)){
            hat_matrix <- rbind(hat_matrix,(x[i,]%*%solve(t(x)%*%(w*x*wt))*t((x*w*wt)))) #criando hat_matrix aqui, entao supostamente nao haveria problema com dimens?o
            if(i==1){
              W_f <- cbind(matrix(i,N,1),w,t(seq(1,nrow(w)))) 
            }
            else{
              W_f <- rbind(W_f,c(cbind(matrix(i,nn,1),w,t(seq(1,nrow(w)))))) 
            }
          }  
          if(!is.null(xvarglobal)){ 
            ba <- solve(t(t(xa)%*%diag(1,nn,nn)-hat_matrix)%*%(diag(1,nn,nn)-hat_matrix)%*%(xa*wt))%*%t(xa)%*%(t(diag(1,nn,nn)-hat_matrix))%*%(diag(1,nn,nn)-hat_matrix)*(y*wt) 
            ya <- y-xa%*%ba
            for(i in 1:nn){
              m1 <- (i-1)*nn+1
              m2 <- m1+(nn-1)
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
            S <- diag(xa%*%solve(t(xa)%*%t(diag(1,nn,nn)-hat_matrix)%*%(diag(1,nn,nn)-hat_matrix)%*%(xa*wt))%*%t((xa*wt))%*%t(diag(1,nn,nn)-hat_matrix)%*%(diag(1,nn,nn)-hat_matrix)+hat_matrix)
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
            zj <- nj+(y-uj)/(Ai*(1+alpha*uj))-offset
            if (det(t(x)%*%(w*Ai*x*wt))==0){
              b <- rep(0, nvar)
            }
            else{
              b <- solve(t(x)%*%(w*Ai*x*wt))%*%t(x)%*%(w*Ai*wt*zj)
            }
            nj <- x*b+offset
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
        yhat[i] <- uj[i]
        alphai[i] <- alpha
        if (det(t(x)%*%(w*Ai*x*wt))==0){
          S[i] <- 0
        }
        else{
          S[i] <- (x[i, ]%*%solve(t(x)%*%(w*Ai*x*wt))%*%t(x*w*Ai*wt))[i]
        }   
        if(!is.null(xvarglobal)){
          if(i==1){
            W_f <- cbind(matrix(i,N,1),w,t(seq(1:nrow(w))))
          }
          else{
            W_f <- rbind(W_f,c(cbind(matrix(i,N,1),w,t(seq(1,nrow(w))))))
          }
        }
        if(!is.null(xvarglobal)){
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
            za <- nj+(y-uj)/(Aa*(1+alphag*uj))-offset
            if (det(t(xa)%*%(Aa*xa))==0){
              ba <- rep(0,nvar)
            }
            else{
              bg <- solve(t(xa)%*%(Aa*xa))%*%t(xa)%*%(Aa*za)
            }
            nj <- xa%*%ba+offset
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
            for(i in 1:nn){
              m1 <- (i-1)*nn+1
              m2 <- m1+(nn-1)
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
                zj <- nj+(y-uj)/(Ai*(1+alpha*uj))-offset-xa*ba
                if (det(t(x)%*%(w*Ai*x*wt))==0){
                  b <- rep(0, nvar)
                }
                else{
                  b <- solve(t(x)%*%(w*Ai*x*wt))%*%t(x)%*%(w*Ai*wt*zj)
                }
                nj <- x*b+offset+xa*ba
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
              Z_matrix <- cbind(Z_matrix,(zj+offset+xa*ba))
              yhat[i] <- uj[i]
            } #fecha la?o for
            hat_matrix <- R_matrix%*%Z_matrix/diag(Z_matrix)
            uj <- yhat
            nj <- log(uj)
            Aa <- (uj/(1+alphag*uj))+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj*uj))
            Aa <- ifelse(Aa<=0,E^-5,Aa)
            za <- nj+(y-uj)/(Aa*(1+alphag*uj))-offset
            ba <- solve(t(xa*Aa)%*%(diag(1,nn,nn)-hat_matrix))%*%t(xa*wt)%*%t(xa*Aa)%*%(diag(1,nn,nn)-hat_matrix)%*%(za*wt)
            nj <- xa*ba+offset
            nj <- ifelse(nj>E^2,E^2,nj)
            uj <- exp(nj)
            diffba <- oldba-ba
            contb <- contb+1
          }
          S <- diag((diag(1,nn,nn)-hat_matrix)%*%xa%*%solve(t(xa*Aa)%*%(diag(1,nn,nn)-hat_matrix)%*%(xa*wt))%*%t(xa*Aa*wt)%*%(diag(1,nn,nn)-hat_matrix)+hat_matrix)
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
        AICC <- AIC + (2*npar*(npar+1))/(nn-npar-1)
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
          yhat[i] <- uj[i]
          if(det(t(x)%*%(w*Ai*x*wt)==0)){
            S[i] <- 0
          }
          else{
            S[i] <- (x[i,]%*%solve(t(x)%*%(w*Ai*x*wt))*t((x*w*wt*Ai)))[i]
          }
          if(!is.null(xvarglobal)){
            if(i==1){
              W_f <- cbind(matrix(i,nn,1),w,t(seq(1,nrow(w))))
            }
            else{
              W_f <- rbind(W_f,c(cbind(matrix(i,nn,1),w,t(seq(1,nrow(w))))))
            }
          }
          if(!is.null(xvarglobal)){
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
          for(i in 1:nn) {
            m1 <- (i-1)*nn+1
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
              C <- matrix(0,nvar,nn)
            }
            else{
              C <- solve(t(x)%*%(w*Ai*x*wt))*t(x)*t(w*Ai*wt)
            }
            R_matrix <- rbind(R_matrix,(x[i,]*C))
            Z_matrix <- rbind(Z_matrix,(zj+xa*ba))
            yhat[i] <- uj[i]
          } # fecha o laco for
          hat_matrix <- R_matrix%*%Z_matrix/diag(Z_matrix)
          uj <- yhat
          nj <- log(uj/(1-uj))
          Aa <- uj*(1-uj)
          Aa <- ifelse(Aa<=0,E^-5,Aa)
          za <- nj+(y-uj)/Aa
          if (det(t(xa%*%Aa)%*%(diag(1,nn,nn)-hat_matrix)*(xa*wt))==0){
            ba <- matrix(0,ncol(xa),1)
          }
          else{
            ba <- solve(t(xa*Aa)%*%(diag(1,nn,nn)-hat_matrix))%*%t(xa*Aa)%*%(diag(1,nn,nn)-hat_matrix)%*%(za*wt)
          }
          nj <- xa*ba
          nj <- ifelse(nj>E^2,E^2,nj)
          uj <- exp(nj)/(1+exp(nj))
          diffba <- oldba-ba
          contb <- contb+1
        } #fecha while (diffba)
        S <- diag((diag(1,nn,nn)-hat_matrix)%*%xa%*%solve(t(xa*Aa)%*%(diag(1,nn,nn)-hat_matrix)%*%(xa*wt))%*%t(xa*Aa*wt)%*%(diag(1,nn,nn)-hat_matrix)+hat_matrix)
        uj <- ifelse(uj==0,E^-10,uj)
        uj <- ifelse(uj==1,0.99999,uj)
        CV <- t((y-yhat)*wt)%*%(y-yhat)
        ll <- sum(y*log(uj)-(1-y)*log(1-uj))
        npar <- sum(S)
        AIC <-  2*npar-2*ll
        AICC <-  AIC +(2*npar*(npar+1))/(nn-npar-1)
        if(bandwidth=="AIC"){
          CV <- AICC
        } 
      } # fecha modelo logistico  
      res <- cbind(CV,npar)  
    } # fecha for da linha 148
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
    bx <- nn
    #print(nn)
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
  h2 <- ax1-r*(bx1-ax1)
  print(c("h0,h1,h2,h3: ", c(h0,h1,h2,h3)))
  
  # /***************************************/ #
  
  res1 <- cv(h1)
  CV1 <- res1[1]
  res2 <- cv(h2)
  CV2 <- res2[1]
  
  if(toupper(METHOD)=="FIXED_G"|toupper(METHOD)=="FIXED_BSQ"|toupper(METHOD)=="ADAPTIVE_BSQ"){
    var_ <<- data.frame()
    if(GMY==1){
      var_ <<- rbind(out,c(GMY,h1,CV1,h2,CV2))
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
    else{
      golden <- CV2
      xmin[GMY,1] <- golden
      xmin[GMY,2] <- h2
      npar <- res2[1]
      if(toupper(METHOD)=='ADAPTIVE_BSQ'){
        xmin[GMY,2] <- floor(h2)
      }
    }
  }
} #fecha golden

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
  
GWR <- function(DATA, YVAR, XVAR, XVARGLOBAL, XVARINF, WEIGHT=NULL, LAT, LONG,
                GRID, DHV, METHOD, MODEL="NEGBIN",OFFSET=NULL,
                DISTANCEKM="NO", H=NULL){
  E <- 10
  yy <- DATA[,YVAR]
  xx <- DATA[,which(names(DATA) %in% XVAR)]
  N <- nrow(yy)
  xx <- as.matrix(cbind(rep(1,N),xx))
  Yhat <- rep(0,N)
  Nvar <- ncol(xx)
  if(!is.null(XVARGLOBAL)){
    xa <- as.matrix(XVARGLOBAL)
  }
  Wt <-rep(1, N)
  if(!is.null(WEIGHT)){
    Wt <- as.matrix(WEIGHT)
  }
  Offset <- rep(0, N)
  if(!is.null(OFFSET)){
    Offset <- as.matrix(OFFSET)
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
        g <- sum(digamma(parg+yy)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+yy)/(parg+uj))
        hess <- sum(trigamma(parg+yy)-trigamma(parg)+1/parg-2/(parg+uj)+(yy+parg)/((parg+uj)^2))
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
      Ai <- (uj/(1+alphag*uj))+(yy-uj)*(alphag*uj/1+2*alphag*uj+alphag^2*uj*uj)
      Ai <- ifelse(Ai<=0,E^-5,Ai)
      zj <- nj+(yy-uj)/(Ai*(1+alphag*uj))-Offset
      if (det(t(xx)%*%(Ai*xx))==0){
        bg <- rep(0,Nvar)
      }
      else{
        bg <- solve(t(xx)%*%(Ai*xx))%*%t(xx)%*%(Ai*zj)
      }
      nj <- xx%*%bg+Offset
      nj <- ifelse(nj>E^2,E^2,nj)
      uj <- exp(nj)
      olddev <- devg
      uj <- ifelse(uj<E^-150,E^-150,uj)
      tt <- yy/yj
      tt <- ifelse(tt==0,E^-10,tt)
      if(toupper(MODEL)=="POISSON"){
        devg <- 2*sum(y*log(tt)-(yy-uj))
      }
      if(toupper(MODEL)=="NEGBIN"){
        devg <- 2*sum(yy*log(tt)-(yy+1/alphag)*log((1+alphag*yy)/(1+alphag*uj)))
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
  
# /*****************************************/ linha 586 SAS
  
  
} #fecha GWR



