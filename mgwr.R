## MGWR ##

library(sf)
library(foreign)
library(dplyr)
setwd('~/PIBIC/multiscale')

georgia <- st_read('G_utm.shp')
georgiashp <- st_read('Georgia.shp')
georgiashp$AREAKEY <- as.numeric(georgiashp$AREAKEY)
georgia_data <- read.dbf('G_utm1.dbf')
georgia_data <- read.csv('GData_utm.csv')

means <- georgia_data%>%
  summarize(PctBach=mean(PctBach), TotPop90=mean(TotPop90),
            PctRural=mean(PctRural), PctEld=mean(PctEld),
            PctFB=mean(PctFB), PctPov=mean(PctPov),
            PctBlack=mean(PctBlack))
print(c("means", means))
sds <- georgia_data%>%
  summarize(PctBach=sd(PctBach), TotPop90=sd(TotPop90),
            PctRural=sd(PctRural), PctEld=sd(PctEld),
            PctFB=sd(PctFB), PctPov=sd(PctPov),
            PctBlack=sd(PctBlack))
print(c("standard deviations", sds))

georgia_data_std <- georgia_data
georgia_data_std$PctBach <- (georgia_data_std$PctBach-mean(georgia_data_std$PctBach))/sd(georgia_data_std$PctBach)
georgia_data_std$TotPop90 <- (georgia_data_std$TotPop90-mean(georgia_data_std$TotPop90))/sd(georgia_data_std$TotPop90)
georgia_data_std$PctRural <- (georgia_data_std$PctRural-mean(georgia_data_std$PctRural))/sd(georgia_data_std$PctRural)
georgia_data_std$PctEld <- (georgia_data_std$PctEld-mean(georgia_data_std$PctEld))/sd(georgia_data_std$PctEld)
georgia_data_std$PctFB <- (georgia_data_std$PctFB-mean(georgia_data_std$PctFB))/sd(georgia_data_std$PctFB)
georgia_data_std$PctPov <- (georgia_data_std$PctPov-mean(georgia_data_std$PctPov))/sd(georgia_data_std$PctPov)
georgia_data_std$PctBlack <- (georgia_data_std$PctBlack-mean(georgia_data_std$PctBlack))/sd(georgia_data_std$PctBlack)

means2 <- georgia_data_std%>%
  summarize(PctBach=mean(PctBach), TotPop90=mean(TotPop90),
            PctRural=mean(PctRural), PctEld=mean(PctEld),
            PctFB=mean(PctFB), PctPov=mean(PctPov),
            PctBlack=mean(PctBlack))
print(c("means", means2))
sds2 <- georgia_data_std%>%
  summarize(PctBach=sd(PctBach), TotPop90=sd(TotPop90),
            PctRural=sd(PctRural), PctEld=sd(PctEld),
            PctFB=sd(PctFB), PctPov=sd(PctPov),
            PctBlack=sd(PctBlack))
print(c("standard deviations", sds2))

mod1 <- lm(PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, data=georgia_data)
summary(mod1)
mod2 <- lm(PctBach~PctBlack+PctFB+TotPop90+PctEld, data=georgia_data)
summary(mod2)
mod3 <- lm(PctBach~PctBlack+PctFB+TotPop90+PctEld, data=georgia_data_std)
summary(mod3)

#abaixo, passar x, y, LAT e LONG como strings --> nomes das variáveis
mgwr <- function(DATA, YVAR, XVAR, WEIGHT=NULL, LAT, LONG,
                 GLOBALMIN="yes", METHOD, MODEL="GAUSSIAN",
                 MGWR="yes", BANDWIDTH="CV", OFFSET=NULL,
                 DISTANCEKM="NO", INT=50, H=NULL){
  #lembrar de tornar band global
  E <- 10
  y <- DATA[, YVAR]
  x <- DATA[, which(names(DATA) %in% XVAR)]
  N <- length(y)
  Wt <-rep(1, N)
  if (!is.null(WEIGHT)){
    Wt <- as.matrix(WEIGHT)
  }
  Offset <- rep(0, N)
  if (!is.null(OFFSET)){
    Offset <- as.matrix(OFFSET)
  }
  x <- as.matrix(cbind(rep(1, N), x))
  nvarg <- ncol(x)
  Yhat <- rep(0, N)
  bi <- matrix(0, nvarg*N, 4)
  Alphai <- matrix(0, N, 3)
  s <- rep(0, N)
  mRj <- matrix(0, N, N*nvarg)
  Sm <- matrix(0, N, N)
  Sm3 <- matrix(0, N, nvarg)
  Rj <- matrix(0, N, N)
  Cm <- matrix(0, N, N*nvarg)
  stdbm <- matrix(0, N, nvarg)
  mAi <- matrix(0, N, nvarg)
  ENP <- rep(0, nvarg+2)
  ## global estimates ##
  if (toupper(MODEL)=="POISSON" | toupper(MODEL)=="NEGBIN"){
    uj <- (y+mean(y))/2
    nj <- log(uj)
    #print(nj)
    Parg <- sum((y-uj)^2/uj)/(N-nvarg)
    ddpar <- 1
    cont <- 1
    while (abs(ddpar)>0.000001 & cont<100){
      dpar <- 1
      parold <- Parg
      cont1 <- 1
      cont3 <- 1
      if(toupper(MODEL)=="POISSON"){
        Alphag <- E^-6
        Parg <- 1/Alphag
      }
      else{
        if (cont>1){
          Parg <- 1/(sum((y-uj)^2/uj)/(N-nvarg))
        }
        while (abs(dpar)>0.000001 & cont1<200){
          Parg <- ifelse(Parg<E^-10, E^-10, Parg)
          g <- sum(digamma(Parg+y)-digamma(Parg)+log(Parg)+1-log(Parg+uj)-(Parg+y)/(Parg+uj))
          hess <- sum(trigamma(Parg+y)-trigamma(Parg)+1/Parg-2/(Parg+uj)+(y+Parg)/(Parg+uj)^2)
          hess <- ifelse(hess==0, E^-23, hess)
          par0 <- Parg
          Parg <- par0-solve(hess)%*%g
          if (cont1>50 & Parg>E^5){
            dpar <- 0.0001
            cont3 <- cont3+1
            if (cont3==1){
              Parg <- 2
            }
            else if (cont3==2){
              Parg <- E^5
            }
            else if (cont3==3){
              Parg <- 0.0001
            }
          }
          else{
            dpar <- Parg-par0
          }
          cont1 <- cont1+1
          if (Parg>E^6){
            Parg <- E^6
            dpar <- 0
          }
        }
        Alphag <- 1/Parg
      }
      devg <- 0
      ddev <- 1
      cont2 <- 0
      while (abs(ddev)>0.000001 & cont2<100){
        #print("sim")
        uj <- ifelse(uj>E^100, E^100, uj)
        Ai <- (uj/(1+Alphag*uj))+(y-uj)*(Alphag*uj/(1+2*Alphag*uj+alphag^2*uj*uj))
        Ai <- ifelse(Ai<=0, E^-5, Ai)
        zj <- nj+(y-uj)/(Ai*(1+Alphag*uj))-Offset
        #print(nj)
        if (det(t(x)%*%(Ai*x))==0){
          bg <- rep(0, nvarg)
        }
        else{
          bg <- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
          #print(zj)
        }
        #print(bg)
        nj <- x%*%bg+Offset
        nj <- ifelse(nj>E^2, E^2, nj)
        uj <- exp(nj)
        olddev <- devg
        uj <- ifelse(uj<E^-150, E^-150, uj)
        tt <- y/uj
        tt <- ifelse(tt==0, E^-10, tt)
        if (toupper(MODEL)=="POISSON"){
          devg <- 2*sum(y*log(tt)-(y-uj))
          #print(tt)
        }
        else{
          devg <- 2*sum(y*log(tt)-(y+1/Alphag)*log((1+Alphag*y)/(1+Alphag*uj)))
          sealphag <- sqrt(1/abs(hess))/(Parg^2)
        }
        if (cont2>100){
          ddev <- 0.0000001
        }
        else{
          ddev <- devg-olddev
        }
        cont2 <- cont2+1
        #print(ddev)
      }
      #linha 128
      Ujg <- uj
      Yhat <- uj
      cont <- cont+1
      ddpar <- Parg-parold
    }
    varg <- vecdiag(solve(t(x*Wt*Ai)%*%x))
  }
  #linha 136
  else if (toupper(MODEL)=="LOGISTIC"){
    uj <- (y+mean(y))/2
    nj <- log(uj/(1-uj))
    devg <- 0
    ddev <- 1
    cont <- 0
    while (abs(ddev)>0.000001 & cont<100){
      uj <- ifelse(uj>E^100, E^100, uj)
      Ai <- uj*(1-uj)
      Ai <- ifelse(Ai<=0, E^-5, Ai)
      zj <- nj+(y-uj)/Ai
      if (det(t(x)%*%(Wt*Ai*x))==0){
        bg <- rep(0, nvarg)
      }
      else{
        bg <- solve(t(x)%*%(Wt*Ai*x))%*%t(x)%*%(Wt*Ai*zj)
      }
      nj <- x%*%bg
      nj <- ifelse(nj>E^2, E^2, nj)
      uj <- exp(nj)/(1+exp(nj))
      olddev <- devg
      uj <- ifelse(uj<E^-150, E^-150, uj)
      tt <- y/uj
      tt <- ifelse(tt==0, E^-10, tt)
      uj <- ifelse(uj==1, 0.99999, uj)
      tt2 <- (1-y)/(1-uj)
      tt2 <- ifelse(tt2==0, E^-10, tt2)
      devg <- 2*sum((y*log(tt))+(1-y)*log(tt2))
      ddev <- devg-olddev
      cont <- cont+1
    }
    Ujg <- uj
    Yhat <- uj
    varg <- vecdiag(solve(t(x*Wt*Ai)%*%x))
  }
  #print(bg) #teste
  #print(Alphag)
  #linha 167
  LONG <- DATA[, LONG]
  LAT <- DATA[, LAT]
  COORD <- matrix(c(LONG, LAT), ncol=2)
  Distance <- dist(COORD, "euclidean")
  Sequ <- 1:N
  cv <- function(h, y, x, fi, n=N, wt=Wt, ujg=Ujg,
                 yhat=Yhat, coord=COORD, distance=Distance,
                 sequ=Sequ, offset=Offset, alphag=Alphag,
                 alphai=Alphai, S=s, parg=Parg,
                 yhat_beta=Yhat_beta){ #ujg, coord, offset, alphag não são usados?
    nvar <- ncol(x)
    for (i in 1:n){
      #linhas 214, 288, 331 fecham esse loop no sas
      for (j in 1:n){
        seqi <- rep(i, n)
        distan <- cbind(seqi, t(sequ), as.matrix(distance)[,i])
        if (toupper(DISTANCEKM)=="YES"){
          distan[,3] <- distan[,3]*111
        }
      }
      u <- nrow(distan)
      w <- rep(0, u)
      for (jj in 1:u){
        w[jj] <- exp(-0.5*(distan[jj,3]/h)^2)
        if (toupper(BANDWIDTH)=="CV"){
          w[i] <- 0
        }
      }
      if (toupper(METHOD)=="FIXED_BSQ"){
        position <- which(distan[,3]>h)
        w[position] <- 0
      }
      else if (toupper(METHOD)=="ADAPTIVE_BSQ"){
        distan <- distan[order(distan[, 3]), ]
        distan <- cbind(distan, 1:nrow(distan))
        w <- matrix(0, n, 2)	
        hn <- distan[h,3]
        for (jj in 1:n){
          if (distan[jj,4]<=h){
            w[jj,1] <- (1-(distan[jj,3]/hn)^2)^2
          }
          else{
            w[jj,1] <- 0
          }
          w[jj,2] <- dist[jj,2]
        }
        if (toupper(BANDWIDTH)=="CV"){
          w[which(w[,2]==i)] <- 0
        }
        w <- w[order(w[, 2]), ]
        w <- w[ ,1]
      }
      #linha 207
      if (toupper(MODEL)=="GAUSSIAN"){
        if (det(t(x)%*%(w*x*wt))==0){
          b <- rep(0, nvar)
        }
        else{
          b <- solve(t(x)%*%(w*x*wt))%*%t(x)%*%(w*y*wt)
        }
        yhat[i] <- x[i, ]*b
        if (det(t(x)%*%(w*x*wt))==0){
          S[i] <- 0
        }
        else{
          S[i] <- (x[i,]%*%solve(t(x)%*%(w*x*wt))%*%t(x*w*wt))[i]
        }
        next
      }
      # linha 220
      else if (toupper(MODEL)=="POISSON" | toupper(MODEL)=="NEGBIN"){
        uj <- yhat
        par <- parg
        nj <- log(uj)
        ddpar <- 1
        cont <- 1
        cont3 <- 0
        while (abs(ddpar)>0.000001 & cont<100){
          dpar <- 1
          parold <- par
          cont1 <- 1
          if (toupper(MODEL)=="POISSON"){
            alpha <- E^-6
            par <- 1/alpha
          }
          #linha 233
          else{
            if (par<=E^-5 & i>1){
              par <- 1/alphai[i-1, 2]
            }
            while (abs(dpar)>0.000001 & cont1<200){
              par <- ifelse(par<E^-10, E^-10, par)
              g <- sum(w*wt*(digamma(par+y)-digamma(par)+log(par)+1-log(par+uj)-(par+y)/(par+uj)))
              hess <- sum(w*wt*(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)^2))
              hess <- ifelse(hess==0, E^-23, hess)
              par0 <- par
              par <- par0-solve(hess)%*%g
              if (cont1>50 & par>E^5){
                dpar <- 0.0001
                cont3 <- cont3+1
                if (cont3==1){
                  par <- 2
                }
                else if (cont3==2){
                  par <- E^5
                }
                else if (cont3==3){
                  par <- 0.0001
                }
              }
              else{
                dpar <- par-par0
              }
              cont1 <- cont1+1
              if (par>E^6){
                par <- E^6
                dpar <- 0
              }
              if (par<=E^-5){
                par <- E^-3
                dpar <- 0
              }
            }
            alpha <- 1/par
          }
          #linha 256
          dev <- 0
          ddev <- 1
          cont2 <- 1
          while (abs(ddev)>0.000001 & cont2<100){
            uj <- ifelse(uj>E^100, E^100, uj)
            Ai <- (uj/(1+alpha*uj))+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj*uj))
            Ai <- ifelse(Ai<=0, E^-5, Ai)
            zj <- nj+(y-uj)/(Ai*(1+alpha*uj))-yhat_beta+fi
            if (det(t(x)%*%(w*Ai*x*wt))==0){
              b <- rep(0, nvar)
            }
            else{
              b <- solve(t(x)%*%(w*Ai*x*wt))%*%t(x)%*%(w*Ai*wt*zj)
            }
            nj <- x*b+yhat_beta-fi
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
        } #linha 283
        yhat[i] <- uj[i]
        alphai[i, 2] <- alpha
        if (det(t(x)%*%(w*Ai*x*wt))==0){
          S[i] <- 0
        }
        else{
          S[i] <- (x[i, ]%*%solve(t(x)%*%(w*Ai*x*wt))%*%t(x*w*Ai*wt))[i]
        }
        next
      } #linha 301
      else if (toupper(MODEL)=="LOGISTIC"){
        uj <- yhat
        nj <- log(uj/(1-uj))
        dev <- 0
        ddev <- 1
        cont <- 0
        while (abs(ddev)>0.000001 & cont<100){
          cont <- cont+1
          uj <- ifelse(uj>E^100, E^100, uj)
          Ai <- uj*(1-uj)
          Ai <- ifelse(Ai<=0, E^-5, Ai)	
          zj <- nj+(y-uj)/Ai-yhat_beta+fi
          if (det(t(x)%*%(w*Ai*x*wt))==0){
            b <- rep(0, nvar)
          }
          else{
            b <- solve(t(x)%*%(w*Ai*x*wt))%*%t(x)%*%(w*Ai*wt*zj)
          }
          nj <- x*b+yhat_beta-fi
          nj <- ifelse(nj>E^2, E^2, nj)
          uj <- exp(nj)/(1+exp(nj))
          olddev <- dev
          uj <- ifelse(uj<E^-150, E^-150, uj)
          tt <- y/uj
          tt <- ifelse(tt==0, E^-10, tt)
          uj <- ifelse(uj==1, 0.99999, uj)
          tt2 <- (1-y)/(1-uj)
          tt2 <- ifelse(tt2==0,E^-10, tt2)
          dev <- 2*sum((y*log(tt))+(1-y)*log(tt2))
          if (cont>100){
            ddev <-  0.0000001
          }
          else{
            ddev <- dev-olddev
          }
        }
        yhat[i] <- uj[i]
        if (det(t(x)%*%(w*Ai*x*wt))==0){
          S[i] <- 0
        }
        else{
          S[i] <- (x[i,]%*%solve(t(x)%*%(w*Ai*x*wt))%*%t(x*w*wt*Ai))[i]
        }
        next
      }
    } #fechando o for
    #depois do next:
    if (toupper(MODEL)=="GAUSSIAN"){
      CV <- t((y-yhat)*wt)%*%(y-yhat)
      npar <- sum(S)
      AICc <- 2*n*log(CV/n)+n*log(2*3.14159)+n*(n+npar)/(n-2-npar)
    }
    else if (toupper(MODEL)=="POISSON" | toupper(MODEL)=="NEGBIN"){
      CV <- t((y-yhat)*wt)%*%(y-yhat)
      if (toupper(MODEL)=="POISSON"){
        ll <- sum(-yhat+y*log(yhat)-lgamma(y+1))
        npar <- sum(S)
      }
      else{
        ll <- sum(y*log(alphai[,2]*yhat)-(y+1/alphai[,2])*log(1+alphai[,2]*yhat)+lgamma(y+1/alphai[,2])-lgamma(1/alphai[,2])-lgamma(y+1))
        npar <- sum(S)+sum(S)/nvar
      }
      AIC <- 2*npar-2*ll
      AICC <- AIC +(2*npar*(npar+1))/(n-npar-1)
    }
    else if (toupper(MODEL)=="LOGISTIC"){
      uj <- ifelse(uj==0, E^-10, uj)
      uj <- ifelse(uj==1, 0.99999, uj)
      CV <- t((y-yhat)*wt)%*%(y-yhat)
      ll <- sum(y*log(uj)-(1-y)*log(1-uj))
      npar <- sum(S)
      AIC <- 2*npar-2*ll
      AICC <- AIC +(2*npar*(npar+1))/(n-npar-1)
    }
    if (toupper(BANDWIDTH)=="AIC"){
      CV <- AICC
    }
    #free dist w
    res <- cbind(cv, npar)
    return (res)
  } #linha 344
  GSS <- function(depy, indepx, fix, distance=Distance, n=N){
    # DEFINING GOLDEN SECTION SEARCH PARAMETERS #
    if(toupper(METHOD)=="FIXED_G" | toupper(METHOD)=="FIXED_BSQ"){
      ax <- 0
      bx <- int(max(distance)+1)
      if (toupper(DISTANCEKM)=="YES"){
        bx <- bx*111
      }
    }
    else if (toupper(METHOD)=="ADAPTIVE_BSQ"){
      ax <- 5
      bx <- n
    }
    r <- 0.61803399
    tol <- 0.1
    #linha 359
    #parte duplicada:
    if (toupper(GLOBALMIN)=="NO"){
      lower <- ax
      upper <- bx
      xmin <- rep(0, 2)
      GMY <- 1
      ax1 <- lower[GMY]
      bx1 <- upper[GMY]
      h0 <- ax1
      h3 <- bx1
      h1 <- bx1-r*(bx1-ax1)
      h2 <- ax1+r*(bx1-ax1)
      #print c(h0, h1, h2, h3)
      res1 <- cv(h1, depy, indepx, fix)
      CV1 <- res1[1]
      res2 <- cv(h2,depy,indepx,fix)
      CV2 <- res2[1]
      int <- 1
      while(abs(h3-h0) > tol*(abs(h1)+abs(h2)) & int<200){
        if (CV2<CV1){
          h0 <- h1
          h1 <- h3-r*(h3-h0)
          h2 <- h0+r*(h3-h0)
          CV1 <- CV2
          res2 <- cv(h2,depy,indepx,fix)
          CV2 <- res2[1]
        }
        else{
          h3 <- h2
          h1 <- h3-r*(h3-h0)
          h2 <- h0+r*(h3-h0)
          CV2 <- CV1
          res1 <- cv(h1, depy, indepx, fix)
          CV1 <- res1[1]
        }
        int <- int+1
      }
      if (CV1<CV2){
        golden <- CV1
        xmin[GMY,1] <- golden
        xmin[GMY,2] <- h1
        npar <- res1[1]
        if (toupper(METHOD)=="ADAPTIVE_BSQ"){
          xmin[GMY,2] <- floor(h1)
          xming <- xmin[GMY,2]
        }
      }
      else{
        golden <- CV2
        xmin[GMY,1] <- golden
        xmin[GMY,2] <- h2
        npar <- res2[1]
        if (toupper(METHOD)=="ADAPTIVE_BSQ"){
          xmin[GMY,2] <- floor(h2)
          xming <- xmin[GMY,2]
        }
      }
      xming <- xmin[GMY,2]
      #print(golden)
      #print((xmin[GMY,2])['xmin']) --> verificar
      #if (toupper(BANDWIDTH)=="AIC"){print(npar)}
    }
    #acrescentei:
    else{
      lower <- cbind(ax, (1-r)*bx, r*bx)
      upper <- cbind((1-r)*bx, r*bx, bx)
      xmin <- matrix(0, 3, 2)
      for (GMY in 1:3){
        ax1 <- lower[GMY]
        bx1 <- upper[GMY]
        h0 <- ax1
        h3 <- bx1
        h1 <- bx1-r*(bx1-ax1)
        h2 <- ax1+r*(bx1-ax1)
        #print(c(h0, h1, h2, h3))
        res1 <- cv(h1, depy, indepx, fix)
        CV1 <- res1[1]
        res2 <- cv(h2,depy,indepx,fix)
        CV2 <- res2[1]
        int <- 1
        while(abs(h3-h0) > tol*(abs(h1)+abs(h2)) & int<200){
          if (CV2<CV1){
            h0 <- h1
            h1 <- h3-r*(h3-h0)
            h2 <- h0+r*(h3-h0)
            CV1 <- CV2
            res2 <- cv(h2,depy,indepx,fix)
            CV2 <- res2[1]
          }
          else{
            h3 <- h2
            h1 <- h3-r*(h3-h0)
            h2 <- h0+r*(h3-h0)
            CV2 <- CV1
            res1 <- cv(h1, depy, indepx, fix)
            CV1 <- res1[1]
          }
          int <- int+1
        }
        if (CV1<CV2){
          golden <- CV1
          xmin[GMY,1] <- golden
          xmin[GMY,2] <- h1
          npar <- res1[1]
          if (toupper(METHOD)=="ADAPTIVE_BSQ"){
            xmin[GMY,2] <- floor(h1)
            xming <- xmin[GMY,2]
          }
        }
        else{
          golden <- CV2
          xmin[GMY,1] <- golden
          xmin[GMY,2] <- h2
          npar <- res2[1]
          if (toupper(METHOD)=="ADAPTIVE_BSQ"){
            xmin[GMY,2] <- floor(h2)
            xming <- xmin[GMY,2]
          }
        }
        xming <- xmin[GMY,2]
        #print(golden)
        #print((xmin[GMY,2])['xmin']) --> verificar
        #if (toupper(BANDWIDTH)=="AIC"){print(npar)}
      }
    }
    #linha 426
    if (toupper(GLOBALMIN)=="YES"){
      #print(xmin[c('golden', 'bandwidth')]) -->verificar
      xming <- xmin[which(xmin[,1]==min(xmin[,1])),2]
      #print(c('Global Minimum', xming['(Da Silva and Mendes (2018)']) --> verificar
    }
    bandwidth <- xming
    return (bandwidth)
  } #linha 433
  #conferi até aqui
  #linha 435
  gwr <- function(h, y, x, fi, distance=Distance, n=N,
                  nvarg=nvarg, sequ=Sequ, wt=Wt, ujg=Ujg,
                  parg=Parg, yhat=Yhat, offset=Offset,
                  Sm=Sm, Sm3=Sm3, Rj=Rj, mRj=mRj,
                  alphai=Alphai, yhat_beta=Yhat_beta,
                  Ai=Ai){
    nvar <- ncol(x)
    bim <- rep(0, nvar*n)
    yhatm <- rep(0, n)
    for (i in 1:n){
      for (j in 1:n){                                                                                                                        
        seqi <- rep(i, n)
        dist <- cbind(seqi, t(sequ), distance[,i])
        if (toupper(DISTANCEKM)=="YES"){
          distan[,3] <- distan[,3]*111
        }
      }
      u <- nrow(distan)
      w <- rep(0, u)
      for (jj in 1:u){
        if (toupper(METHOD)=="FIXED_G"){
          w[jj] <- exp(-(dist[jj,3]/h)^2)
        }
        else if (toupper(METHOD)=="FIXED_BSQ"){
          w[jj] <- (1-(dist[jj,3]/h)^2)^2
        }
      }
      if (toupper(METHOD)=="ADAPTIVE_BSQ"){ #linha 457
        distan <- distan[order(distan[, 3]), ]
        distan <- cbind(distan, 1:nrow(distan))
        w <- matrix(0, n, 2)	 
        hn <- distan[h,3]
        for (jj in 1:n){
          if (distan[jj,4]<=h){
            w[jj,1] <- (1-(dist[jj,3]/hn)^2)^2
          }
          else{
            (w[jj,1]==0)
          }
          w[jj,2] <- dist[jj,2]
        }
        w <- w[order(w[, 2]), ]
        w <- w[,1]
      }
      ## MODEL SELECTION ##
      if (toupper(MODEL)=="GAUSSIAN"){
        if (det(t(x)%*%(w*x*wt))==0){
          b <- rep(0, nvar)
        }
        else{
          b <- solve(t(x)%*%(w*x*wt))%*%t(x)%*%(w*y*wt)
        }
        uj <- x%*%b
        if (nvar==nvarg){
          if (det(t(x)%*%(w*x*wt))==0){
            Sm[i,] <- rep(0, n)
            mRj[i,] <- matrix(0, n*nvar)
          }
          else{
            ej <- diag(nvar)
            Sm[i,] <- (x[i,]%*%solve(t(x)%*%(w*x*wt))%*%t(x*w*wt))
            Sm3[i,] <- t(diag((solve(t(x)%*%(w*x*wt))%*%t(x*w*wt))%*%t(solve(t(x)%*%(w*x*wt))%*%t(x*w*wt))))
            for (jj in 1:nvar){
              m1 <- (jj-1)*n+1
              m2 <- m1+(n-1)
              mRj[i, m1:m2] <- (x[i,jj]*ej[jj,]*solve(t(x)%*%(w*x*wt))%*%t(x*w*wt))
            }
          }
        }
        else{
          if (det(t(x)%*%(w*x*wt))==0){
            Rj[i,] <- rep(0, n)
          }
          else{
            Rj[i,] <- (x[i,]%*%solve(t(x)%*%(w*x*wt))%*%t(x*w*wt))
          }
        }
      } #linha 493
      else if (toupper(MODEL)=="POISSON" | toupper(MODEL)=="NEGBIN"){
        uj <- yhat
        par <- parg
        nj <- log(uj)
        ddpar <- 1
        cont <- 1
        while (abs(ddpar)>0.000001 & cont<100){
          dpar <- 1
          parold <- par
          cont1 <- 1
          cont3 <- 1
          if (toupper(MODEL)=="POISSON"){
            alpha <- E^-6
            par <- 1/alpha
          }
          if (toupper(MODEL)=="NEGBIN"){
            if (par<=E^-5 & i>1){
              par=1/alphai[i-1,2]
            }
            while (abs(dpar)>0.000001 & cont1<200){
              par <- ifelse(par<E^-10, E^-10, par)
              g <- sum(w*wt*(digamma(par+y)-digamma(par)+log(par)+1-log(par+uj)-(par+y)/(par+uj)))
              hess <- sum(w*wt*(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)^2))
              hess <- ifelse(hess==0, E^-23, hess)
              par0 <- par
              par <- par0-solve(hess)*g
              if (cont1>50 & par>E^5){
                dpar <- 0.0001
                cont3 <- cont3+1
                if (cont3==1){
                  par <- 2
                }
                else if (cont3==2){
                  par <- E^5
                }
                else if (cont3==3){
                  par <- 0.0001
                }
              }
              else{
                dpar <- par-par0
              }
              cont1 <- cont1+1
              if (par>E^6){
                par <- E^6
                dpar <- 0
              }
              if (par<=E^-5){
                par <- E^-3
                dpar <- 0
              }
            }
            alpha <- 1/par
          } #linha 529
          dev <- 0
          ddev <- 1
          cont2 <- 0
          while (abs(ddev)>0.000001 & cont2<100){
            uj <- ifelse(uj>E^100, E^100, uj)
            Ai <- (uj/(1+alpha*uj))+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj*uj))
            Ai <- ifelse(Ai<=0, E^-5, Ai)	
            zj <- nj+(y-uj)/(Ai*(1+alpha*uj))-yhat_beta+fi
            if (det(t(x)%*%(w*Ai*x*wt))==0){
              b <- rep(0, nvar)
            }
            else{
              b <- solve(t(x)%*%(w*Ai*x*wt))%*%t(x)%*%(w*Ai*wt*zj)
            }
            nj <- x*b+yhat_beta-fi
            nj <- ifelse(nj>E^2, E^2, nj)
            uj <- exp(nj)
            olddev <- dev
            uj <- ifese(uj<E^-150, E^-150, uj)
            tt <- y/uj
            tt <- ifelse(tt==0, E^-10, tt)
            if (toupper(MODEL)=="POISSON"){
              dev <- 2*sum(y*log(tt)-(y-uj))
            }
            if (toupper(MODEL)=="NEGBIN"){
              dev <- 2*sum(y*log(tt)-(y+1/alpha)*log((1+alpha*y)/(1+alpha*uj)))
            }
            cont2 <- cont2+1
          }
          cont <- cont+1
          ddpar <- par-parold
        }  #linha 555
        if (nvar==nvarg){
          if (det(t(x)%*%(w*Ai*x*wt))==0){
            Sm[i,] <- c(0, n)
            mRj[i,] <- rep(0, n*nvar)
          }
          else{
            ej <- diag(nvar)
            Sm[i,] <- (x[i,]%*%solve(t(x)%*%(w*Ai*x*wt))%*%t(x*w*wt*Ai))
            Sm3[i,] <- t(diag((solve(t(x)%*%(w*Ai*x*wt))%*%t(x*w*wt*Ai))%*%diag(1/Ai)%*%t(solve(t(x)%*%(w*Ai*x*wt))%*%t(x*w*wt*Ai))))
            for (jj in 1:nvar){
              m1 <- (jj-1)*n+1
              m2 <- m1+(n-1)
              mRj[i, m1:m2] <- (x[i,jj]*ej[jj,]*solve(t(x)%*%(w*Ai*x*wt))%*%t(x*w*wt*Ai))
            }
          }
        }
        else{
          if (det(t(x)%*%(w*Ai*x*wt))==0){
            Rj[i,] <- rep(0, n)
          }
          else{
            Rj[i,] <- (x[i,]%*%solve(t(x)%*%(w*Ai*x*wt))%*%t(x*w*wt*Ai))
          }
        }
        if (toupper(MODEL)=="NEGBIN"){
          hess <- sum(w*wt*(trigamma(par+y)-trigamma(par)+1/par-2/(par+exp(yhat_beta))+(y+par)/(par+exp(yhat_beta))^2))
          if (toupper(MGWR)!="YES"){
            hess <- sum(w*wt*(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)^2))
            hess <- ifelse(hess==0, E^-23, hess)
          }
          sealpha <- sqrt(1/abs(hess))/(par^2)
          alphai[i,1] <- i
          alphai[i,2] <- alpha
          alphai[i,3] <- sealpha
        }
      } #linha 584
      else if (toupper(MODEL)=="LOGISTIC"){
        uj <- yhat
        nj <- log(uj/(1-uj))
        dev <- 0
        ddev <- 1
        cont <- 1
        while (abs(ddev)>0.000001 & cont<100){
          cont <- cont+1
          uj <- ifelse(uj>E^100, E^100, uj)
          Ai <- uj*(1-uj)
          Ai <- ifelse(Ai<=0, E^-5, Ai)	
          zj <- nj+(y-uj)/Ai-yhat_beta+fi
          if (det(t(x)%*%(w*Ai*x*wt))==0){
            b <- rep(0, nvar)
          }
          else{
            b <- solve(t(x)%*%(w*Ai*x*wt))%*%t(x)%*%(w*Ai*zj*wt)
          }
          nj <- x*b+yhat_beta-fi
          nj <- ifelse(nj>E^2, E^2, nj)
          uj <- exp(nj)/(1+exp(nj))
          olddev <- dev
          uj <- ifelse(uj<E^-150, E^-150, uj)
          tt <- y/uj
          tt <- ifelse(tt==0, E^-10, tt)
          uj <- ifelse(uj==1, 0.99999, uj)
          tt2 <- (1-y)/(1-uj)
          tt2 <- ifelse(tt2==0, E^-10, tt2)
          dev <- 2*sum((y*log(tt))+(1-y)*log(tt2))
          if (cont>100){
            ddev <- 0.0000001
          }
          else{
            ddev <- dev-olddev
          }
        }
        if (nvar==nvarg){
          if (det(t(x)%*%(w*Ai*x*wt))==0){
            Sm[i,] <- rep(0, n)
            mRj[i,] <- matrix(0, n*nvar)
          }
          else{
            ej <- diag(nvar)
            Sm[i,] <- (x[i,]%*%solve(t(x)%*%(w*Ai*x*wt))%*%t(x*w*wt*Ai))
            Sm3[i,] <- t(diag((solve(t(x)%*%(w*Ai*x*wt))%*%t(x*w*wt*Ai))%*%diag(1/Ai)%*%t(solve(t(x)%*%(w*Ai*x*wt))%*%t(x*w*wt*Ai))))
            for (jj in 1:nvar){
              m1 <- (jj-1)*n+1
              m2 <- m1+(n-1)
              mRj[i, m1:m2] <- (x[i,jj]*ej[jj,]*solve(t(x)%*%(w*Ai*x*wt))%*%t(x*w*wt*Ai))
            }
          }
        }
        else{
          if (det(t(x)%*%(w*Ai*x*wt))==0){
            Rj[i,] <- rep(0, n)
          }
          else{
            Rj[i,] <- (x[i,]*solve(t(x)%*%(w*Ai*x*wt))%*%t(x*w*wt*Ai))
          }
        }
      }
      m1 <- (i-1)*nvar+1
      m2 <- m1+(nvar-1)
      bim[m1:m2] <- b
      yhatm[i] <- uj[i]
      yhat[i] <- uj[i]
    } #linha 634
    beta <- matrix(bim, n)
    yhbeta <- cbind(yhatm, beta)
    return (yhbeta) #conferi até aqui (falta if else {} ())
  }
  #testando:
  if (toupper(MGWR)!="yes"){
    finb <- rep(0, N)
    Yhat_beta <- Offset
    if (!is.null(H)){
      h <- H
    }
    else{
      h <- GSS(y,x,finb)
    }
    print(c("General Bandwidth", h))
    Yhat_beta <- gwr(h,y,x,finb)
    beta <- Yhat_beta[,2:nvarg+1]
    fi <- x*beta
    mband <- h
    Sm2 <- Sm
  }
  else{
    finb <- rep(0, N)
    Yhat_beta <- Offset
    if (!is.null(H)){
      h <- H
    }
    else{
      h <- GSS(y,x,finb)
    }
    print(c("General Bandwidth", h))
  }
}

#ifs por else ifs
#conferir ; choose j inv = * E do

#checar argumentos gwr e GSS

## Testes ##
mgwr(DATA=georgia_data_std, YVAR="PctBach",
     XVAR=c("PctBlack", "PctFB", "TotPop90", "PctEld"),
     LAT="Y", LONG="X", GLOBALMIN="no", METHOD="adaptive_bsq",
     BANDWIDTH="cv", MODEL="gaussian")
