## MGWR ##

library(sf)
library(foreign)
library(dplyr)
library(readxl)
library(ggplot2)
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

mod4 <- lm(PctBach~PctBlack+PctFB+TotPop90+PctEld, data=georgia_data_std)
summary(mod4)

georgia_data_std2 <- read.csv('georgia_data_std2.csv')

nakaya <- read.csv('nakaya.csv')

nakaya_std <- read.csv('nakaya_std.csv')

nakaya2 <- read.csv('nakaya2.csv')
nakaya_std2 <- read.csv('nakaya_std2.csv')

nakaya_NB <- read.csv('nakaya_NB.csv')
nakaya_NB_std2 <- read.csv('nakaya_NB_std2.csv')

logistic <- read.csv("Logistic.csv")

logistic_std <- read.csv("logistic_std.csv")

mgwnbr <- function(DATA, YVAR, XVAR, WEIGHT=NULL, LAT, LONG,
                   GLOBALMIN="yes", METHOD, MODEL="GAUSSIAN",
                   MGWR="yes", BANDWIDTH="CV", OFFSET=NULL,
                   DISTANCEKM="NO", INT=50, H=NULL){
  output <- list() #flag
  yhat_beta <<- NULL
  E <- 10
  Y <- DATA[, YVAR]
  X <- DATA[XVAR]
  N <<- length(Y)
  wt <<-rep(1, N)
  if (!is.null(WEIGHT)){
    wt <<- DATA[, WEIGHT]
    wt <<- as.matrix(wt)
  }
  Offset <<- rep(0, N)
  if (!is.null(OFFSET)){
    Offset <- DATA[, OFFSET]
    Offset <<- as.matrix(Offset)
  }
  X <- as.matrix(cbind(rep(1, N), X))
  nvarg <<- ncol(X)
  yhat <- rep(0, N)
  bi <- matrix(0, nvarg*N, 4)
  alphai <<- matrix(0, N, 3)
  s <<- rep(0, N)
  mrj <<- matrix(0, N, N*nvarg)
  sm <<- matrix(0, N, N)
  sm3 <<- matrix(0, N, nvarg)
  rj <<- matrix(0, N, N)
  Cm <- matrix(0, N, N*nvarg)
  stdbm <- matrix(0, N, nvarg)
  mAi <- matrix(0, N, nvarg)
  ENP <- rep(0, nvarg+2)
  ## global estimates ##
  if (toupper(MODEL)=="POISSON" | toupper(MODEL)=="NEGBIN"){
    uj <- (Y+mean(Y))/2
    nj <- log(uj)
    parg <<- sum((Y-uj)^2/uj)/(N-nvarg)
    ddpar <- 1
    cont <- 1
    while (abs(ddpar)>0.000001 & cont<100){
      dpar <- 1
      parold <- parg
      cont1 <- 1
      cont3 <- 1
      if(toupper(MODEL)=="POISSON"){
        alphag <<- E^-6
        parg <<- 1/alphag
      }
      else{
        if (cont>1){
          parg <<- 1/(sum((Y-uj)^2/uj)/(N-nvarg))
        }
        while (abs(dpar)>0.000001 & cont1<200){
          parg <<- ifelse(parg<E^-10, E^-10, parg)
          g <- sum(digamma(parg+Y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+Y)/(parg+uj))
          hess <- sum(trigamma(parg+Y)-trigamma(parg)+1/parg-2/(parg+uj)+(Y+parg)/(parg+uj)^2)
          hess <- ifelse(hess==0, E^-23, hess)
          par0 <- parg
          parg <<- par0-solve(hess)%*%g
          if (cont1>50 & parg>E^5){
            dpar <- 0.0001
            cont3 <- cont3+1
            if (cont3==1){
              parg <<- 2
            }
            else if (cont3==2){
              parg <<- E^5
            }
            else if (cont3==3){
              parg <- 0.0001
            }
          }
          else{
            dpar <- parg-par0
          }
          cont1 <- cont1+1
          if (parg>E^6){
            parg <<- E^6
            dpar <- 0
          }
        }
        alphag <<- as.numeric(1/parg)
      }
      devg <- 0
      ddev <- 1
      cont2 <- 0
      while (abs(ddev)>0.000001 & cont2<100){
        uj <- ifelse(uj>E^100, E^100, uj)
        ai <<- as.numeric((uj/(1+alphag*uj))+(Y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj*uj)))
        ai <<- ifelse(ai<=0, E^-5, ai)
        zj <- nj+(Y-uj)/(ai*(1+alphag*uj))-Offset
        if (det(t(X)%*%(ai*X))==0){
          bg <- rep(0, nvarg)
        }
        else{
          bg <- solve(t(X)%*%(ai*X))%*%t(X)%*%(ai*zj)
        }
        nj <- X%*%bg+Offset
        nj <- ifelse(nj>E^2, E^2, nj)
        uj <- as.numeric(exp(nj))
        olddev <- devg
        uj <- ifelse(uj<E^-150, E^-150, uj)
        tt <- Y/uj
        tt <- ifelse(tt==0, E^-10, tt)
        if (toupper(MODEL)=="POISSON"){
          devg <- 2*sum(Y*log(tt)-(Y-uj))
        }
        else{
          devg <- 2*sum(Y*log(tt)-(Y+1/alphag)*log((1+alphag*Y)/(1+alphag*uj)))
          sealphag <- sqrt(1/abs(hess))/(parg^2)
        }
        if (cont2>100){
          ddev <- 0.0000001
        }
        else{
          ddev <- devg-olddev
        }
        cont2 <- cont2+1
      }
      #linha 128
      ujg <<- uj
      yhat <- uj
      cont <- cont+1
      ddpar <- parg-parold
    }
    varg <- diag(solve(t(X*wt*ai)%*%X))
  }
  #linha 136
  else if (toupper(MODEL)=="LOGISTIC"){
    uj <- (Y+mean(Y))/2
    nj <- log(uj/(1-uj))
    devg <- 0
    ddev <- 1
    cont <- 0
    while (abs(ddev)>0.000001 & cont<100){
      uj <- ifelse(uj>E^100, E^100, uj)
      ai <<- as.numeric(uj*(1-uj))
      ai <<- ifelse(ai<=0, E^-5, ai)
      zj <- nj+(Y-uj)/ai
      if (det(t(X)%*%(wt*ai*X))==0){
        bg <- rep(0, nvarg)
      }
      else{
        bg <- solve(t(X)%*%(wt*ai*X))%*%t(X)%*%(wt*ai*zj)
      }
      nj <- X%*%bg
      nj <- ifelse(nj>E^2, E^2, nj)
      uj <- exp(nj)/(1+exp(nj))
      olddev <- devg
      uj <- ifelse(uj<E^-150, E^-150, uj)
      tt <- Y/uj
      tt <- ifelse(tt==0, E^-10, tt)
      uj <- ifelse(uj==1, 0.99999, uj)
      tt2 <- (1-Y)/(1-uj)
      tt2 <- ifelse(tt2==0, E^-10, tt2)
      devg <- 2*sum((Y*log(tt))+(1-Y)*log(tt2))
      ddev <- devg-olddev
      cont <- cont+1
    }
    ujg <<- uj
    yhat <- uj
    varg <- diag(solve(t(X*wt*ai)%*%X))
  }
  #linha 167
  LONG <- DATA[, LONG]
  LAT <- DATA[, LAT]
  COORD <<- matrix(c(LONG, LAT), ncol=2)
  distance <- dist(COORD, "euclidean")
  sequ <<- 1:N
  cv <- function(h, y, x, fi){ #ujg, coord, Offset, alphag
    nvar <- ncol(x)
    for (i in 1:N){
      #linhas 214, 288, 331 fecham esse loop no sas
      for (j in 1:N){
        seqi <- rep(i, N)
        distan <- cbind(seqi, sequ, as.matrix(distance)[,i])
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
        w <- matrix(0, N, 2)	
        hn <- distan[h,3]
        for (jj in 1:N){
          if (distan[jj,4]<=h){
            w[jj,1] <- (1-(distan[jj,3]/hn)^2)^2
          }
          else{
            w[jj,1] <- 0
          }
          w[jj,2] <- distan[jj,2]
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
        yhat[i] <<- x[i, ]%*%b
        if (det(t(x)%*%(w*x*wt))==0){
          s[i] <<- 0
        }
        else{
          s[i] <<- (x[i,]%*%solve(t(x)%*%(w*x*wt))%*%t(x*w*wt))[i]
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
              par <- as.numeric(1/alphai[i-1, 2])
            }
            while (abs(dpar)>0.000001 & cont1<200){
              par <- ifelse(par<E^-10, E^-10, par)
              g <- sum(w*wt*(digamma(par+y)-digamma(par)+log(par)+1-log(par+uj)-(par+y)/(par+uj)))
              hess <- sum(w*wt*(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)^2))
              hess <- ifelse(hess==0, E^-23, hess)
              par0 <- par
              par <- as.numeric(par0-solve(hess)%*%g)
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
            ai <<- as.numeric((uj/(1+alpha*uj))+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj*uj)))
            ai <<- ifelse(ai<=0, E^-5, ai)
            zj <- nj+(y-uj)/(ai*(1+alpha*uj))-yhat_beta+fi
            if (det(t(x)%*%(w*ai*x*wt))==0){
              b <- rep(0, nvar)
            }
            else{
              b <- solve(t(x)%*%(w*ai*x*wt))%*%t(x)%*%(w*ai*wt*zj)
            }
            nj <- x%*%b+yhat_beta-fi #alterei multiplicador
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
        yhat[i] <<- uj[i]
        alphai[i, 2] <<- alpha
        if (det(t(x)%*%(w*ai*x*wt))==0){
          s[i] <<- 0
        }
        else{
          s[i] <<- (x[i, ]%*%solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*ai*wt))[i]
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
          ai <<- as.numeric(uj*(1-uj))
          ai <<- ifelse(ai<=0, E^-5, ai)	
          zj <- nj+(y-uj)/ai-yhat_beta+fi
          if (det(t(x)%*%(w*ai*x*wt))==0){
            b <- rep(0, nvar)
          }
          else{
            b <- solve(t(x)%*%(w*ai*x*wt))%*%t(x)%*%(w*ai*wt*zj)
          }
          nj <- x%*%b+yhat_beta-fi #alterei multiplicador
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
        yhat[i] <<- uj[i]
        if (det(t(x)%*%(w*ai*x*wt))==0){
          s[i] <<- 0
        }
        else{
          s[i] <<- (x[i,]%*%solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))[i]
        }
        next
      }
    } #fechando o for
    #depois do next:
    if (toupper(MODEL)=="GAUSSIAN"){
      CV <- t((y-yhat)*wt)%*%(y-yhat)
      npar <- sum(s)
      AICc <- 2*N*log(CV/N)+N*log(2*3.14159)+N*(N+npar)/(N-2-npar)
    }
    else if (toupper(MODEL)=="POISSON" | toupper(MODEL)=="NEGBIN"){
      CV <- t((y-yhat)*wt)%*%(y-yhat)
      if (toupper(MODEL)=="POISSON"){
        ll <- sum(-yhat+y*log(yhat)-lgamma(y+1))
        npar <- sum(s)
      }
      else{
        ll <- sum(y*log(alphai[,2]*yhat)-(y+1/alphai[,2])*log(1+alphai[,2]*yhat)+lgamma(y+1/alphai[,2])-lgamma(1/alphai[,2])-lgamma(y+1))
        npar <- sum(s)+sum(s)/nvar
      }
      AIC <- 2*npar-2*ll
      AICC <- AIC +(2*npar*(npar+1))/(N-npar-1)
    }
    else if (toupper(MODEL)=="LOGISTIC"){
      uj <- ifelse(uj==0, E^-10, uj)
      uj <- ifelse(uj==1, 0.99999, uj)
      CV <- t((y-yhat)*wt)%*%(y-yhat)
      ll <- sum(y*log(uj)-(1-y)*log(1-uj))
      npar <- sum(s)
      AIC <- 2*npar-2*ll
      AICC <- AIC +(2*npar*(npar+1))/(N-npar-1)
    }
    if (toupper(BANDWIDTH)=="AIC"){
      CV <- AICC
    }
    #free dist w
    res <- cbind(CV, npar)
    return (res)
  } #linha 344
  GSS <- function(depy, indepx, fix){
    # DEFINING GOLDEN SECTION SEARCH PARAMETERS #
    if(toupper(METHOD)=="FIXED_G" | toupper(METHOD)=="FIXED_BSQ"){
      ax <- 0
      bx <- as.integer(max(distance)+1)
      if (toupper(DISTANCEKM)=="YES"){
        bx <- bx*111
      }
    }
    else if (toupper(METHOD)=="ADAPTIVE_BSQ"){
      ax <- 5
      bx <- N
    }
    r <- 0.61803399
    tol <- 0.1
    #linha 359
    #parte duplicada:
    if (toupper(GLOBALMIN)=="NO"){
      lower <- ax
      upper <- bx
      xmin <- matrix(0, 1, 2)
      GMY <- 1
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
        #print(c(h0, h1, h2, h3, CV1, CV2))
      }
      if (CV1<CV2){
        golden <- CV1
        xmin[GMY, 1] <- golden
        xmin[GMY, 2] <- h1
        npar <- res1[1]
        if (toupper(METHOD)=="ADAPTIVE_BSQ"){
          xmin[GMY, 2] <- floor(h1)
          xming <- xmin[GMY, 2]
        }
      }
      else{
        golden <- CV2
        xmin[GMY, 1] <- golden
        xmin[GMY, 2] <- h2
        npar <- res2[1]
        if (toupper(METHOD)=="ADAPTIVE_BSQ"){
          xmin[GMY, 2] <- floor(h2)
          xming <- xmin[GMY, 2]
        }
      }
      xming <- xmin[GMY, 2]
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
  #linha 435
  gwr <- function(h, y, x, fi){ #ujg, Offset
    nvar <- ncol(x)
    bim <- rep(0, nvar*N)
    yhatm <- rep(0, N)
    for (i in 1:N){
      for (j in 1:N){                                                                                                                        
        seqi <- rep(i, N)
        distan <- cbind(seqi, sequ, as.matrix(distance)[,i])
        if (toupper(DISTANCEKM)=="YES"){
          distan[,3] <- distan[,3]*111
        }
      }
      u <- nrow(distan)
      w <- rep(0, u)
      if (toupper(METHOD)=="FIXED_G"){
        for (jj in 1:u){
          w[jj] <- exp(-(distan[jj,3]/h)^2)
        }
      }
      else if (toupper(METHOD)=="FIXED_BSQ"){
        for (jj in 1:u){
          w[jj] <- (1-(distan[jj,3]/h)^2)^2
        }
      }
      else if (toupper(METHOD)=="ADAPTIVE_BSQ"){ #linha 457
        distan <- distan[order(distan[, 3]), ]
        distan <- cbind(distan, 1:nrow(distan))
        w <- matrix(0, N, 2)
        hn <- distan[h,3]
        for (jj in 1:N){
          if (distan[jj,4]<=h){
            w[jj,1] <- (1-(distan[jj,3]/hn)^2)^2
          }
          else{
            w[jj,1] <- 0
          }
          w[jj,2] <- distan[jj,2]
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
            sm[i,] <<- rep(0, N)
            mrj[i,] <<- matrix(0, N*nvar)
          }
          else{
            ej <- diag(nvar)
            sm[i,] <<- (x[i,]%*%solve(t(x)%*%(w*x*wt))%*%t(x*w*wt))
            sm3[i,] <<- t(diag((solve(t(x)%*%(w*x*wt))%*%t(x*w*wt))%*%t(solve(t(x)%*%(w*x*wt))%*%t(x*w*wt))))
            for (jj in 1:nvar){
              m1 <- (jj-1)*N+1
              m2 <- m1+(N-1)
              mrj[i, m1:m2] <<- (x[i,jj]*ej[jj,])%*%solve(t(x)%*%(w*x*wt))%*%t(x*w*wt)
            }
          }
        }
        else{
          if (det(t(x)%*%(w*x*wt))==0){
            rj[i,] <<- rep(0, N)
          }
          else{
            rj[i,] <<- (x[i,]%*%solve(t(x)%*%(w*x*wt))%*%t(x*w*wt))
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
          else{ #MODEL=='NEGBIN'
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
            alpha <- as.numeric(1/par)
          } #linha 529
          dev <- 0
          ddev <- 1
          cont2 <- 0
          while (abs(ddev)>0.000001 & cont2<100){
            uj <- ifelse(uj>E^100, E^100, uj)
            ai <<- as.numeric((uj/(1+alpha*uj))+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj*uj)))
            ai <<- ifelse(ai<=0, E^-5, ai)	
            zj <- nj+(y-uj)/(ai*(1+alpha*uj))-yhat_beta+fi
            if (det(t(x)%*%(w*ai*x*wt))==0){
              b <- rep(0, nvar)
            }
            else{
              b <- solve(t(x)%*%(w*ai*x*wt))%*%t(x)%*%(w*ai*wt*zj)
            }
            nj <- x%*%b+yhat_beta-fi #alterei multiplicador
            nj <- ifelse(nj>E^2, E^2, nj)
            uj <- as.numeric(exp(nj))
            olddev <- dev
            uj <- ifelse(uj<E^-150, E^-150, uj)
            tt <- y/uj
            tt <- ifelse(tt==0, E^-10, tt)
            if (toupper(MODEL)=="POISSON"){
              dev <- 2*sum(y*log(tt)-(y-uj))
            }
            else{ #MODEL=="NEGBIN"
              dev <- 2*sum(y*log(tt)-(y+1/alpha)*log((1+alpha*y)/(1+alpha*uj)))
            }
            cont2 <- cont2+1
          }
          cont <- cont+1
          ddpar <- par-parold
        }  #linha 555
        if (nvar==nvarg){
          if (det(t(x)%*%(w*ai*x*wt))==0){
            sm[i,] <<- c(0, N)
            mrj[i,] <<- rep(0, N*nvar)
          }
          else{
            ej <- diag(nvar)
            sm[i,] <<- (x[i,]%*%solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))
            sm3[i,] <<- t(diag((solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))%*%diag(1/ai)%*%t(solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))))
            for (jj in 1:nvar){
              m1 <- (jj-1)*N+1
              m2 <- m1+(N-1)
              mrj[i, m1:m2] <<- (x[i,jj]*ej[jj,])%*%solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai)
            }
          }
        }
        else{
          if (det(t(x)%*%(w*ai*x*wt))==0){
            rj[i,] <<- rep(0, N)
          }
          else{
            rj[i,] <<- (x[i,]%*%solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))
          }
        }
        if (toupper(MODEL)=="NEGBIN"){
          hess <- sum(w*wt*(trigamma(par+y)-trigamma(par)+1/par-2/(par+exp(yhat_beta))+(y+par)/(par+exp(yhat_beta))^2))
          if (toupper(MGWR)!="YES"){
            hess <- sum(w*wt*(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)^2))
            hess <- ifelse(hess==0, E^-23, hess)
          }
          sealpha <- sqrt(1/abs(hess))/(par^2)
          alphai[i,1] <<- i
          alphai[i,2] <<- alpha
          alphai[i,3] <<- sealpha
        }
      } #linha 584
      else{ #else if (toupper(MODEL)=="LOGISTIC"){
        uj <- yhat
        nj <- log(uj/(1-uj))
        dev <- 0
        ddev <- 1
        cont <- 1
        while (abs(ddev)>0.000001 & cont<100){
          cont <- cont+1
          uj <- ifelse(uj>E^100, E^100, uj)
          ai <<- as.numeric(uj*(1-uj))
          ai <<- ifelse(ai<=0, E^-5, ai)
          print(ai)
          zj <- nj+(y-uj)/ai-yhat_beta+fi
          if (det(t(x)%*%(w*ai*x*wt))==0){
            b <- rep(0, nvar)
          }
          else{
            b <- solve(t(x)%*%(w*ai*x*wt))%*%t(x)%*%(w*ai*zj*wt)
          }
          nj <- x%*%b+yhat_beta-fi #alterei multiplicador
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
        if (nvar==nvarg){ #erro-padrão começa aqui
          if (det(t(x)%*%(w*ai*x*wt))==0){
            sm[i,] <<- rep(0, N)
            mrj[i,] <<- matrix(0, N*nvar)
          }
          else{ #daqui
            ej <- diag(nvar)
            sm[i,] <<- x[i,]%*%solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai)
            sm3[i,] <<- t(diag((solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))%*%diag(1/ai)%*%t(solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))))
            for (jj in 1:nvar){
              m1 <- (jj-1)*N+1
              m2 <- m1+(N-1)
              mrj[i, m1:m2] <<- (x[i,jj]*ej[jj,])%*%solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai)
            }
          } #até aqui
        }
        else{
          if (det(t(x)%*%(w*ai*x*wt))==0){
            rj[i,] <<- rep(0, N)
          }
          else{
            rj[i,] <<- (x[i,]%*%solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))
          }
        } #acaba aqui
      }
      m1 <- (i-1)*nvar+1
      m2 <- m1+(nvar-1)
      bim[m1:m2] <- b
      yhatm[i] <- uj[i]
      yhat[i] <<- uj[i]
    } #linha 634
    #criar output para Sm, Sm3 e mRj
    # teste <- wt #(w*ai*x*wt)
    # teste <- as.data.frame(teste)
    # View(teste)
    # Sm_data <- as.data.frame(sm)
    # View(Sm_data)
    # Sm3_data <- as.data.frame(sm3)
    # View(Sm3_data)
    # mrj_data <- as.data.frame(mrj)
    # View(mrj_data)
    beta <- matrix(bim, N, byrow=T)
    yhbeta <- cbind(yhatm, beta)
    return (yhbeta) #conferi até aqui
  }
  #testando:
  if (toupper(MGWR)!="YES"){
    finb <- rep(0, N)
    yhat_beta <<- Offset
    if (!is.null(H)){
      hh <- H
    }
    else{
      hh <- GSS(Y,X,finb)
    }
    print(c("General Bandwidth", hh))
    output <- append(output, hh) #flag
    names(output) <- "general_bandwidth"
    yhat_beta <<- gwr(hh,Y,X,finb)
    beta <- yhat_beta[,2:(nvarg+1)]
    Fi <- X*beta
    mband <- hh
    Sm2 <- sm
  }
  else{
    finb <- rep(0, N)
    yhat_beta <<- Offset
    if (!is.null(H)){
      hh <- H
    }
    else{
      hh <- GSS(Y,X,finb)
    }
    print(c("General Bandwidth", hh))
    output <- append(output, hh) #flag
    names(output) <- "general_bandwidth"
    #computing residuals
    #linha 672
    yhat_beta <<- gwr(hh, Y, X, finb)
    error <- Y-yhat_beta[ ,1]
    beta <- yhat_beta[ ,2:(nvarg+1)]
    Fi <- X*beta #checar multiplicador
    Sm2 <- sm
    for (jj in 1:nvarg){
      m1 <- (jj-1)*N+1
      m2 <- m1+(N-1)
      #print (trace(mrj[,m1:m2]))
    }
    mband <- rep(hh, nvarg)
    socf <- 1
    int <- 1
    mband_socf <- c(mband, socf)
    while (socf>0.001 & int<INT){
      fi_old <- Fi
      diffi <- 0
      fi2 <- 0
      for (i in 1:nvarg){
        if (toupper(MODEL)=="GAUSSIAN"){
          ferror <- error+Fi[,i]
          if (!is.null(H)){
            mband[i] <- H
          }
          else{
            mband[i] <- GSS(ferror, as.matrix(X[,i]), finb)
          }
          yhat_beta <<- gwr(mband[i], ferror, as.matrix(X[,i]), finb)
          beta[,i] <- yhat_beta[,2]
          Fi[,i] <- X[,i]*beta[,i]
          error <- Y-apply(Fi, 1, sum)
          m1 <- (i-1)*N+1
          m2 <- m1+(N-1)
          mrj2 <- mrj[,m1:m2]
          mrj[,m1:m2] <<- rj%*%mrj[,m1:m2]+rj-rj%*%sm
          sm <<- sm-mrj2+mrj[,m1:m2]
          Cm[,m1:m2] <- (1/X[,i])*mrj[,m1:m2]
        }
        else{ #else if (toupper(MODEL)=="POISSON" | toupper(MODEL)=="NEGBIN" | toupper(MODEL)=="LOGISTIC"){
          yhat_beta <<- (apply(Fi, 1, sum)+Offset)
          if (!is.null(H)){
            mband[i] <- H
          }
          else{
            mband[i] <- GSS(Y, as.matrix(X[,i]), Fi[,i])
          }
          yhat_beta <<- gwr(mband[i], Y, as.matrix(X[,i]), Fi[,i])
          beta[,i] <- yhat_beta[,2]
          #print (yhat_beta[:])(_beta_[:,]) (alphai[:,2])
          Fi[,i] <- X[,i]*beta[,i]
          #error=Y-exp(Fi[,+]+Offset)
          m1 <- (i-1)*N+1
          m2 <- m1+(N-1)
          mrj2 <- mrj[,m1:m2]
          mrj[,m1:m2] <<- rj%*%mrj[,m1:m2]+rj-rj%*%sm
          sm <<- sm-mrj2+mrj[,m1:m2]
          Cm[,m1:m2] <- (1/X[,i])*mrj[,m1:m2]
          mAi[,i] <- ai
        }
        diffi <- diffi+mean((Fi[,i]-fi_old[,i])^2)
        fi2 <- fi2+Fi[,i]
      }
      socf <- sqrt(diffi/sum(fi2^2)) #printar
      int <- int+1
      mband_socf <- rbind(mband_socf, c(mband, socf))
    } #linha 750
    mband_socf <- mband_socf[-1, ]
    if (is.null(dim(mband_socf))){
      band <<- as.data.frame(t(mband_socf))
    }
    else{
      band <<- as.data.frame(mband_socf)
    }
    names(band) <<- c("Intercept", XVAR, "socf")
    rownames(band) <<- NULL
    print('Bandwidth')
    print(band)
  }
  v1 <- sum(diag(sm))
  if (toupper(MODEL)=='GAUSSIAN'){
    yhat <- apply(Fi, 1, sum)
    res <- Y-yhat
    rsqr1 <- t(res*wt)%*%res
    ym <- t(Y*wt)%*%Y
    rsqr2 <- ym-(sum(Y*wt)^2)/sum(wt)
    rsqr <- 1-rsqr1/rsqr2
    rsqradj <- 1-((N-1)/(N-v1))*(1-rsqr)
    sigma2 <- as.numeric(N*rsqr1/((N-v1)*sum(wt)))
    root_mse <- sqrt(sigma2)
  }
  for (jj in 1:nvarg){
    m1 <- (jj-1)*N+1
    m2 <- m1+(N-1)
    #print (trace(mrj[,m1:m2]))
    ENP[jj] <- sum(diag(mrj[,m1:m2]))
    if (toupper(MGWR)!='YES'){
      ENP[jj] <- sum(diag(sm))
    }
    if (toupper(MODEL)=='GAUSSIAN'){
      if (toupper(MGWR)!='YES'){
        stdbm[,jj] <- sqrt(sigma2*sm3[,jj])
      }
      else{
        stdbm[,jj] <- sqrt(diag(sigma2*Cm[,m1:m2]%*%t(Cm[,m1:m2])))
      }
    }
    else{ #else if (toupper(MODEL)=='POISSON' | toupper(MODEL)=='NEGBIN' | toupper(MODEL)=='LOGISTIC'){
      if (toupper(MGWR)=='YES'){
        stdbm[,jj] <- sqrt(diag(Cm[,m1:m2]%*%diag(1/mAi[,jj])%*%t(Cm[,m1:m2])))
      }
      else{
        stdbm[,jj] <- sqrt(sm3[,jj])
      }
    }
  } #linha 789
  if (toupper(MODEL)=='GAUSSIAN'){
    ll <- -N*log(rsqr1/N)/2-N*log(2*acos(-1))/2-sum((Y-yhat)*(Y-yhat))/(2*(rsqr1/N))
    AIC <- 2*v1-2*ll
    AICc <- AIC+2*(v1*(v1+1)/(N-v1-1))
    print(rbind(c('Sigma2e', 'Root MSE'), c(sigma2, root_mse)))
    print(rbind(c('R-Square', 'Adj-R-Square'), c(round(rsqr, 4), round(rsqradj, 4))))
    print(rbind(c('Full Log Likelihood', 'AIC', 'AICc'), c(ll, AIC, AICc)))
  } #linha 798
  else if (toupper(MODEL)=='POISSON'){
    yhat <- exp(apply(Fi, 1, sum)+Offset)
    tt <- Y/yhat
    tt <- ifelse(tt==0, E^-10, tt)
    dev <- 2*sum(Y*log(tt)-(Y-yhat))
    ll <- sum(-yhat+Y*log(yhat)-lgamma(Y+1))
    AIC <- 2*v1-2*ll
    AICc <- AIC+2*(v1*(v1+1)/(N-v1-1))
    tt2 <- Y/mean(Y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnull <- 2*sum(Y*log(tt2)-(Y-mean(Y)))
    pctdev <- 1-dev/devnull
    adjpctdev <- 1-((N-1)/(N-v1))*(1-pctdev)
    print(c('Deviance', dev))
    print(c('Full Log Likelihood', ll))
    print(c('pctdev', 'adjpctdev', 'AIC', 'AICc'))
    print(c(pctdev, adjpctdev, AIC, AICc))
  }
  else if (toupper(MODEL)=='NEGBIN'){
    yhat <- exp(apply(Fi, 1, sum)+Offset)
    tt <- Y/yhat
    tt <- ifelse(tt==0, E^-10, tt)
    dev <- 2*sum(Y*log(tt)-(Y+1/alphai[,2])*log((1+alphai[,2]*Y)/(1+alphai[,2]*yhat)))
    ll <- sum(Y*log(alphai[,2]*yhat)-(Y+1/alphai[,2])*log(1+alphai[,2]*yhat)+lgamma(Y+1/alphai[,2])-lgamma(1/alphai[,2])-lgamma(Y+1))
    AIC <- 2*(v1+v1/nvarg)-2*ll
    AICc <- AIC+2*(v1+v1/nvarg)*(v1+v1/nvarg+1)/(N-(v1+v1/nvarg)-1)
    tt2 <- Y/mean(Y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnull <- 2*sum(Y*log(tt2)-(Y+1/alphai[,2])*log((1+alphai[,2]*Y)/(1+alphai[,2]*mean(Y))))
    pctdev <- 1-dev/devnull
    adjpctdev <- 1-((N-1)/(N-(v1+v1/nvarg)))*(1-pctdev)
    print(c('Deviance', dev))
    print(c('Full Log Likelihood', ll))
    print(c('pctdev', 'adjpctdev', 'AIC', 'AICc'))
    print(c(pctdev, adjpctdev, AIC, AICc))
  }
  else{ #else if (toupper(MODEL)=='LOGISTIC'){
    yhat <- exp(apply(Fi, 1, sum))/(1+exp(apply(Fi, 1, sum)))
    tt <- Y/yhat
    tt <- ifelse(tt==0, E^-10, tt)
    yhat2 <- ifelse(yhat==1, 0.99999, yhat)
    tt2 <- (1-Y)/(1-yhat2)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    dev <- 2*sum((Y*log(tt))+(1-Y)*log(tt2))
    lyhat2 <- 1-yhat
    lyhat2 <- ifelse(lyhat2==0, E^-10, lyhat2)
    ll <- sum(Y*log(yhat)+(1-Y)*log(lyhat2))
    AIC <- 2*v1-2*ll
    AICc <- AIC+2*(v1*(v1+1)/(N-v1-1))
    tt <- Y/mean(Y)
    tt <- ifelse(tt==0, E^-10, tt)
    tt2 <- (1-Y)/(1-mean(Y))
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnull <- 2*sum((Y*log(tt))+(1-Y)*log(tt2))
    pctdev <- 1-dev/devnull
    adjpctdev <- 1-((N-1)/(N-v1))*(1-pctdev)
    print(c('Deviance', dev))
    print(c('Full Log Likelihood', ll))
    print(c('pctdev', 'adjpctdev', 'AIC', 'AICc'))
    print(c(pctdev, adjpctdev, AIC, AICc))
  } #linha 848
  output <- append(output, dev) #flag
  names(output)[length(output)] <- "deviance"
  output <- append(output, ll)
  names(output)[length(output)] <- "full_Log_likelihood"
  output <- append(output, pctdev)
  names(output)[length(output)] <- "pctdev"
  output <- append(output, adjpctdev)
  names(output)[length(output)] <- "adjpctdev"
  output <- append(output, AIC)
  names(output)[length(output)] <- "AIC"
  output <- append(output, AICc)
  names(output)[length(output)] <- "AICc"
  ENP[nvarg+1] <- sum(diag(sm))
  ENP[nvarg+2] <- sum(diag(Sm2))
  varname_enp <- c('Intercept', XVAR, 'MGWR', 'GWR')
  if (toupper(MODEL)=='NEGBIN'){
    ENP <- c(ENP, (v1/nvarg))
    varname_enp <- c(varname_enp, 'alpha')
    names(ENP) <- c('Intercept', XVAR, 'MGWR', 'GWR', 'alpha') #flag
  }
  else{
    names(ENP) <- c('Intercept', XVAR, 'MGWR', 'GWR')
  }
  #ENPprint <- rbind(varname_enp, ENP)
  #rownames(ENPprint) <- NULL
  print('ENP')
  print(ENP) #flag
  output <- append(output, list(ENP)) #flag
  names(output)[length(output)] <- "ENP"
  dff <- N-v1
  tstat <- beta/stdbm
  probt <- 2*(1-pt(abs(tstat), dff))
  malpha <- ENP
  malpha[1:nvarg] <- 0.05/ENP[1:nvarg]
  malpha[nvarg+1] <- 0.05*(nvarg/v1)
  malpha[nvarg+2] <- 0.05*(nvarg/sum(diag(Sm2)))
  if (toupper(MGWR)!='YES'){
    malpha[1:nvarg] <- 0.05*(nvarg/v1)
  }
  if (toupper(MODEL)=='NEGBIN'){
    malpha[nvarg+3] <- 0.05*(nvarg/v1)
  }
  t_critical <- abs(qt(malpha/2,dff))
  beta2 <- beta
  if (toupper(MODEL)=='NEGBIN'){
    alpha <- alphai[,2]
    beta2 <- cbind(beta, alpha)
  }
  qntl <- apply(beta2, 2, quantile, c(0.25, 0.5, 0.75)) #, na.rm=T
  IQR <- (qntl[3,]-qntl[1,])
  qntl <- rbind(round(qntl, 3), IQR=round(IQR, 3))
  descriptb <- rbind(apply(beta2, 2, mean), apply(beta2, 2, min), apply(beta2, 2, max))
  rownames(descriptb) <- c('Mean', 'Min', 'Max')
  print("Quantiles of MGWR Parameter Estimates")
  if (toupper(MODEL)=='NEGBIN'){
    colnames(qntl) <- c('Intercept', XVAR, 'alpha')
  }
  else{
    colnames(qntl) <- c('Intercept', XVAR)
  }
  print(qntl)
  output <- append(output, list(qntl)) #flag
  names(output)[length(output)] <- "qntls_mgwr_param_estimates"
  print("Descriptive Statistics")
  if (toupper(MODEL)=='NEGBIN'){
    colnames(descriptb) <- c('Intercept', XVAR, 'alpha')
  }
  else{
    colnames(descriptb) <- c('Intercept', XVAR)
  }
  print(descriptb)
  output <- append(output, list(descriptb)) #flag
  names(output)[length(output)] <- "descript_stats_mgwr_param_estimates"
  stdbeta <- stdbm
  stdbeta2 <- stdbeta
  if (toupper(MODEL)=='NEGBIN'){
    stdalpha <- alphai[,3]
    stdbeta2 <- cbind(stdbeta, stdalpha)
  }
  qntls <- apply(stdbeta2, 2, quantile, c(0.25, 0.5, 0.75)) #, na.rm=T
  IQR <- (qntls[3,]-qntls[1,])
  qntls <- rbind(round(qntls, 3), IQR=round(IQR, 3))
  descripts <- rbind(apply(stdbeta2, 2, mean), apply(stdbeta2, 2, min), apply(stdbeta2, 2, max))
  rownames(descripts) <- c('Mean', 'Min', 'Max')
  print("alpha-level=0.05")
  alpha_levelprint <- rbind(varname_enp, malpha)
  rownames(alpha_levelprint) <- NULL
  
  print(alpha_levelprint)
  print("t-Critical")
  tcritical_print <- rbind(varname_enp, round(t_critical, 2))
  rownames(tcritical_print) <- NULL
  print(tcritical_print)
  print("Quantiles of MGWR Standard Errors")
  if (toupper(MODEL)=='NEGBIN'){
    colnames(qntls) <- c('Intercept', XVAR, 'alpha')
  }
  else{
    colnames(qntls) <- c('Intercept', XVAR)
  }
  print(qntls)
  print("Descriptive Statistics of Standard Errors")
  if (toupper(MODEL)=='NEGBIN'){
    colnames(descripts) <- c('Intercept', XVAR, 'alpha')
  }
  else{
    colnames(descripts) <- c('Intercept', XVAR)
  }
  print(descripts) #linha 902
  #### global estimates ####
  if (toupper(MODEL)=='GAUSSIAN'){
    bg <- solve(t(X)%*%(X*wt))%*%t(X)%*%(Y*wt)
    s2g <- as.numeric(t((Y-X%*%bg)*wt)%*%(Y-X%*%bg)/(N-nrow(bg)))
    varg <- diag(solve(t(X)%*%(X*wt))*s2g)
  }
  if (is.null(WEIGHT)){
    vargd <- varg
    dfg <- N-nrow(bg)
  }
  stdg <- matrix(sqrt(vargd))
  if (toupper(MODEL)=='NEGBIN'){
    bg <- rbind(bg, alphag)
    stdg <- rbind(stdg, sealphag)
    dfg <- dfg-1
  }
  tg <- bg/stdg
  probtg <- 2*(1-pt(abs(tg), dfg))
  bg_stdg_tg_probtg <- cbind(bg, stdg, tg, probtg)
  print("Global Parameter Estimates")
  if (toupper(MODEL)=='NEGBIN'){
    rownames(bg_stdg_tg_probtg) <- c('Intercept', XVAR, 'alpha')
  }
  else{
    rownames(bg_stdg_tg_probtg) <- c('Intercept', XVAR)
  }
  colnames(bg_stdg_tg_probtg) <- c("Par. Est.", "Std Error", "t Value", "Pr > |t|")
  print(bg_stdg_tg_probtg)
  print(c("NOTE: The denominator degrees of freedom for the t tests is", dfg))
  if (toupper(MODEL)=='GAUSSIAN'){
    resg <- (Y-X%*%bg)
    rsqr1g <- t(resg*wt)%*%resg
    ymg <- t(Y*wt)%*%Y
    rsqr2g <- ymg-(sum(Y*wt)^2)/sum(wt)
    rsqrg <- 1-rsqr1g/rsqr2g
    rsqradjg <- 1-((N-1)/(N-nrow(bg)))%*%(1-rsqrg)
    sigma2g <- N*rsqr1g/((N-nrow(bg))*sum(wt))
    root_mseg <- sqrt(sigma2g)
    print(rbind(c('Sigma2e', 'Root MSE'), c(sigma2g, root_mseg)))
    print(rbind(c("R-Square", "Adj-R-Square"), round(c(rsqrg, rsqradjg), 4)))
    ll <- -N*log(rsqr1g/N)/2-N*log(2*acos(-1))/2-sum(resg*resg)/(2*(rsqr1g/N))
    AIC <- -2*ll+2*nrow(bg)
    AICc <- -2*ll+2*nrow(bg)*(N/(N-nrow(bg)-1))
    print(rbind(c('Full Log Likelihood', 'AIC', 'AICc'), c(ll, AIC, AICc)))
  }
  else if (toupper(MODEL)=='POISSON'){
    yhatg <- exp(X%*%bg+Offset)
    ll <- sum(-yhatg+Y*log(yhatg)-lgamma(Y+1))
    AIC <- -2*ll+2*nvarg
    AICc <- -2*ll+2*nvarg*(N/(N-nvarg-1))
    tt2 <- Y/mean(Y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnullg <- 2*sum(Y*log(tt2)-(Y-mean(Y)))
    pctdevg <- 1-devg/devnullg
    adjpctdevg <- 1-((N-1)/(N-nvarg))*(1-pctdevg)
    print(rbind(c('Deviance', 'Full Log Likelihood',
                  'pctdevg', 'adjpctdevg', 'AIC', 'AICc'),
                c(devg, ll, pctdevg, adjpctdevg, AIC, AICc)))
  }
  else if (toupper(MODEL)=='NEGBIN'){
    yhatg <- exp(X%*%bg[1:(nrow(bg)-1)]+Offset)
    ll <- sum(Y*log(alphag*yhatg)-(Y+1/alphag)*log(1+alphag*yhatg)+lgamma(Y+1/alphag)-lgamma(1/alphag)-lgamma(Y+1))
    AIC <- -2*ll+2*(nvarg+1)
    AICc <- -2*ll+2*(nvarg+1)*(N/(N-(nvarg+1)-1))
    tt2 <- Y/mean(Y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnullg <- 2*sum(Y*log(tt2)-(Y+1/alphag)*log((1+alphag*Y)/(1+alphag*mean(Y))))
    pctdevg <- 1-devg/devnullg
    adjpctdevg <- 1-((N-1)/(N-nvarg))*(1-pctdevg)
    print(rbind(c('Deviance', 'Full Log Likelihood', 'pctdevg',
                  'adjpctdevg', 'AIC', 'AICc'),
                c(devg, ll, pctdevg, adjpctdevg, AIC, AICc)))
  }
  else{ #else if (toupper(MODEL)=='LOGISTIC'){
    yhatg <- exp(X%*%bg)/(1+exp(X%*%bg))
    lyhat2 <- 1-yhatg
    lyhat2 <- ifelse(lyhat2==0, E^-10, lyhat2)
    ll <- sum(Y*log(yhatg)+(1-Y)*log(lyhat2))
    AIC <- -2*ll+2*nvarg
    AICc <- -2*ll+2*nvarg*(N/(N-nvarg-1))
    tt <- Y/mean(Y)
    tt <- ifelse(tt==0, E^-10, tt)
    tt2 <- (1-Y)/(1-mean(Y))
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnullg <- 2*sum((Y*log(tt))+(1-Y)*log(tt2))
    pctdevg <- 1-devg/devnullg
    adjpctdevg <- 1-((N-1)/(N-nvarg))*(1-pctdevg)
    print(rbind(c('Deviance', 'Full Log Likelihood',
                  'pctdevg', 'adjpctdevg', 'AIC', 'AICc'),
                c(devg, ll, pctdevg, adjpctdevg, AIC, AICc)))
  } #linha 981
  bistdt <- cbind(COORD, beta, stdbm, tstat, probt)
  colname1 <- c("Intercept", XVAR)
  parameters2 <<- as.data.frame(bistdt)
  names(parameters2) <<- c('x', 'y', colname1, paste('std_', colname1, sep=''), paste('tstat_', colname1, sep=''), paste('probt_', colname1, sep=''))
  sig <- matrix("not significant at 90%", N, nvarg)
  for (i in 1:N){
    for (j in 1:nvarg){
      if (probt[i,j]<0.01/ENP[j]){
        sig[i,j] <- "significant at 99%"
      }
      else if (probt[i,j]<0.05/ENP[j]){
        sig[i,j] <- "significant at 95%"
      }
      else if (probt[i,j]<0.1/ENP[j]){
        sig[i,j] <- "significant at 90%"
      }
      else{
        sig[i,j] <- "not significant at 90%"
      }
    }
  }
  sig_parameters2 <<- as.data.frame(sig)
  names(sig_parameters2) <<- c(paste('sig_', colname1, sep=''))
  if (toupper(MODEL)=='NEGBIN'){
    atstat <- alphai[,2]/alphai[,3]
    aprobtstat <- 2*(1-pnorm(abs(atstat)))
    siga <- rep("not significant at 90%", N)
    for (i in 1:N){
      if (aprobtstat[i]<0.01*(nvarg/v1)){
        siga[i] <- "significant at 99%"
      }
      else if (aprobtstat[i]<0.05*(nvarg/v1)){
        siga[i] <- "significant at 95%"
      }
      else if (aprobtstat[i]<0.1*(nvarg/v1)){
        siga[i] <- "significant at 90%"
      }
      else{
        siga[i] <- "not significant at 90%"
      }
    }
    alphai <- cbind(alphai, atstat, aprobtstat)
    Alpha <<- as.data.frame(alphai)
    names(Alpha) <<- c("id", "alpha", "std", "tstat", "probt")
    sig_alpha <<- as.data.frame(siga)
    names(sig_alpha) <<- "sig_alpha"
  }
  ###################################
  min_bandwidth <<- as.data.frame(t(mband))
  if (toupper(MGWR)!="YES"){
    names(min_bandwidth) <<- 'Intercept'
  }
  else{
    names(min_bandwidth) <<- colname1
  }
  parameters2 <<- cbind(parameters2, sig_parameters2)
  if (toupper(MODEL)=='NEGBIN'){
    Alpha <<- cbind(Alpha, sig_alpha)
  }
  return(output)
}

#ifs por else ifs
#conferir ; choose j inv = * E do substituições (distan) vecdiag trace

#checar argumentos cv, gwr e GSS (dentro e fora)

## Testes ##

### GAUSSIAN ###

#Teste 1 (globalmin=no, adaptive_bsq, cv, gaussian)
startTime <- Sys.time()
mgwnbr(DATA=georgia_data_std, YVAR="PctBach",
       XVAR=c("PctBlack", "PctFB", "TotPop90", "PctEld"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="adaptive_bsq",
       BANDWIDTH="cv", MODEL="gaussian")
endTime <- Sys.time()
# 19.67 minutos
# SAS: 37 segs

#Teste 2 (globalmin=no, adaptive_bsq, cv, gaussian, mgwr=no)
startTime <- Sys.time()
mgwnbr(DATA=georgia_data_std, YVAR="PctBach",
       XVAR=c("PctBlack", "PctFB", "TotPop90", "PctEld"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="CV", MODEL="gaussian", MGWR="no")
endTime <- Sys.time()
#2.29 minutos
#SAS: 1.8 segs

#Teste 3 (globalmin=no, adaptive_bsq, cv, gaussian)
startTime <- Sys.time()
mgwnbr(DATA=georgia_data, YVAR="PctBach",
       XVAR=c("PctBlack", "PctFB", "TotPop90", "PctEld"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="CV", MODEL="gaussian")
endTime <- Sys.time()
#esse teste dá erro mesmo

#Teste 4 (globalmin=no, adaptive_bsq, cv, gaussian h=159)
startTime <- Sys.time()
mgwnbr(DATA=georgia_data_std, YVAR="PctBach",
       XVAR=c("PctBlack", "PctFB", "TotPop90", "PctEld"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="CV", MODEL="gaussian", H=159)
endTime <- Sys.time()
#1.99 minutos
#SAS: 3 segs

#Teste 5 (globalmin=no, fixed_g, cv, gaussian h=100000000)
startTime <- Sys.time()
mgwnbr(DATA=georgia_data_std, YVAR="PctBach",
       XVAR=c("PctBlack", "PctFB", "TotPop90", "PctEld"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="FIXED_G",
       BANDWIDTH="CV", MODEL="gaussian", H=100000000)
endTime <- Sys.time()
#52.69 segundos
#SAS: 1.2 segs

#Teste 6 (globalmin=no, fixed_g, cv, gaussian, mgwr=no, h=100000000)
startTime <- Sys.time()
mgwnbr(DATA=georgia_data_std, YVAR="PctBach",
       XVAR=c("PctBlack", "PctFB", "TotPop90", "PctEld"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="FIXED_G",
       BANDWIDTH="CV", MODEL="gaussian", MGWR="no", H=100000000)
endTime <- Sys.time()
#9.79 segundos
#SAS: 0.33 segs

print("Mean")
apply(parameters2[, 1:(ncol(parameters2)-5)], 2, mean)
print("Median")
apply(parameters2[, 1:(ncol(parameters2)-5)], 2, median)

#Teste 7 (globalmin=no, adaptive_bsq, cv, gaussian)
startTime <- Sys.time()
mgwnbr(DATA=georgia_data_std2, YVAR="PctBach",
     XVAR=c("PctBlack", "PctFB", "TotPop90", "PctEld"),
     LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
     BANDWIDTH="CV", MODEL="gaussian")
endTime <- Sys.time()
#25.72 mins
#SAS: 28 segs

#Teste 8 (globalmin=no, adaptive_bsq, cv, gaussian) --> ordem das variáveis x
startTime <- Sys.time()
mgwnbr(DATA=georgia_data_std2, YVAR="PctBach",
       XVAR=c("PctEld", "PctBlack", "PctFB", "TotPop90"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="CV", MODEL="gaussian")
endTime <- Sys.time()
#17.92 mins
#SAS: 28 segs

#Teste 9 (globalmin=no, adaptive_bsq, cv, gaussian) --> muda as variáveis
startTime <- Sys.time()
mgwnbr(DATA=georgia_data_std, YVAR="PctBach",
       XVAR=c("PctFB", "PctRural", "PctBlack"),
       LAT="Y", LONG="X", GLOBALMIN="no",
       METHOD="ADAPTIVE_BSQ", BANDWIDTH="CV", MODEL="gaussian")
endTime <- Sys.time()
#17.65 mins
#SAS: 23 segs

print("Mean")
apply(parameters2[, 1:(ncol(parameters2)-5)], 2, mean)
print("Median")
apply(parameters2[, 1:(ncol(parameters2)-5)], 2, median)

mod5 <- lm(PctBach~PctFB+PctRural+PctBlack, data=georgia_data_std)
summary(mod5)

#mapa de parameters2

### POISSON/ NEGBIN ###

#ABAIXO, LEMBRAR DE CONFERIR OS OUTPUTS

#Teste 10 (globalmin=no, adaptive_bsq, aic, poisson)
startTime <- Sys.time()
mgwnbr(DATA=nakaya, YVAR="Mort2564",
       XVAR=c("Professl", "Elderly", "OwnHome", "Unemply"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="aic", OFFSET="Le", MODEL="poisson")
endTime <- Sys.time()
#ALTERNATIVO
#SAS: 19 mins

#Teste 11 (globalmin=no, adaptive_bsq, aic, negbin)
startTime <- Sys.time()
mgwnbr(DATA=nakaya, YVAR="Mort2564",
       XVAR=c("Professl", "Elderly", "OwnHome", "Unemply"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="aic", OFFSET="Le", MODEL="negbin")
endTime <- Sys.time()
#ALTERNATIVO
# SAS: 55 mins

#Teste 12 (globalmin=no, adaptive_bsq, aic, poisson, h=262)
startTime <- Sys.time()
mgwnbr(DATA=nakaya, YVAR="Mort2564",
       XVAR=c("Professl", "Elderly", "OwnHome", "Unemply"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="aic", OFFSET="Le", MODEL="poisson", H=262)
endTime <- Sys.time()
#17.73 mins
#SAS: 33 segs

#Teste 13 (globalmin=no, fixed_g, aic, poisson, h=10^7)
startTime <- Sys.time()
mgwnbr(DATA=nakaya, YVAR="Mort2564",
       XVAR=c("Professl", "Elderly", "OwnHome", "Unemply"),
              LAT="Y", LONG="X", GLOBALMIN="no", METHOD="FIXED_G",
              BANDWIDTH="aic", OFFSET="Le", MODEL="poisson", H=10000000)
endTime <- Sys.time()
#9.59 mins
#SAS: 12 segs

#Teste 14 (globalmin=no, fixed_g, aic, poisson, mgwr=no, h=10^7)
startTime <- Sys.time()
mgwnbr(DATA=nakaya, YVAR="Mort2564",
       XVAR=c("Professl", "Elderly", "OwnHome", "Unemply"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="FIXED_G",
       BANDWIDTH="aic", OFFSET="Le", MODEL="poisson", MGWR="no",
       H=10000000)
endTime <- Sys.time()
#1.68 mins
#SAS: 3 segs

#Teste 15 (globalmin=no, adaptive_bsq, aic, poisson)
startTime <- Sys.time()
mgwnbr(DATA=nakaya_std, YVAR="Mort2564",
       XVAR=c("Professl", "Elderly", "OwnHome", "Unemply"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="aic", OFFSET="Le", MODEL="poisson")
endTime <- Sys.time()
#5.43 horas
#SAS: 5 mins

#teste 16 (globalmin=no, adaptive_bsq, cv, poisson)
startTime <- Sys.time()
mgwnbr(DATA=nakaya_std, YVAR="Mort2564",
       XVAR=c("Professl", "Elderly", "OwnHome", "Unemply"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="cv", OFFSET="Le", MODEL="poisson")
endTime <- Sys.time()
#ALTERNATIVO
#SAS: 5 mins

#Teste 17 (globalmin=no, adaptive_bsq, aic, negbin)
startTime <- Sys.time()
mgwnbr(DATA=nakaya_std, YVAR="Mort2564",
       XVAR=c("Professl", "Elderly", "OwnHome", "Unemply"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="aic", OFFSET="Le", MODEL="negbin")
endTime <- Sys.time()
#ALTERNATIVO
#SAS: 25 mins

#Teste 18 (globalmin=no, adaptive_bsq, aic, negbin, h=262)
startTime <- Sys.time()
mgwnbr(DATA=nakaya_std, YVAR="Mort2564",
       XVAR=c("Professl", "Elderly", "OwnHome", "Unemply"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="aic", OFFSET="Le", MODEL="negbin", H=262)
endTime <- Sys.time()
#19.12 mins
#SAS: 1.3 mins
#comparar outputs

#Teste 19 (globalmin=no, fixed_g, aic, negbin)
startTime <- Sys.time()
mgwnbr(DATA=nakaya_std, YVAR="Mort2564",
       XVAR=c("Professl", "Elderly", "OwnHome", "Unemply"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="FIXED_G",
       BANDWIDTH="aic", OFFSET="Le", MODEL="negbin")
endTime <- Sys.time()
#ALTERNATIVO
#SAS: 1h e 14 mins

#Teste 20 (globalmin=no, fixed_g, aic, negbin, h=10^8)
startTime <- Sys.time()
mgwnbr(DATA=nakaya, YVAR="Mort2564",
       XVAR=c("Professl", "Elderly", "OwnHome", "Unemply"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="FIXED_G",
       BANDWIDTH="aic", OFFSET="Le", MODEL="negbin", H=100000000)
endTime <- Sys.time()
#6.26 mins
#SAS: 32 segs

#Teste 21 (globalmin=no, fixed_g, aic, nebin, mgwr=no, h=10^8)
startTime <- Sys.time()
mgwnbr(DATA=nakaya, YVAR="Mort2564",
       XVAR=c("Professl", "Elderly", "OwnHome", "Unemply"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="FIXED_G",
       BANDWIDTH="aic", OFFSET="Le", MODEL="negbin", MGWR="no",
       H=100000000)
endTime <- Sys.time()
#1.33 mins
#SAS: 7 segs

#Teste 22 (globalmin=no, adaptive_bsq, aic, poisson)
startTime <- Sys.time()
mgwnbr(DATA=nakaya2, YVAR="vary",
       XVAR=c("varx1", "varx2"), LAT="Y", LONG="X", GLOBALMIN="no",
       METHOD="ADAPTIVE_BSQ", BANDWIDTH="aic", MODEL="poisson")
endTime <- Sys.time()
#28.24 mins
#SAS: 31 segs

#Teste 23 (globalmin=no, adaptive_bsq, aic, poisson)
startTime <- Sys.time()
mgwnbr(DATA=nakaya_std2, YVAR="vary", XVAR=c("varx1", "varx2"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="aic", MODEL="poisson")
endTime <- Sys.time()
#30.07 mins
#SAS: 31 segs

#Teste 24 (globalmin=no, adaptive_bsq, aic, negbin)
startTime <- Sys.time()
mgwnbr(DATA=nakaya2, YVAR="vary", XVAR="varx1", LAT="Y", LONG="X",
       GLOBALMIN="no",METHOD="ADAPTIVE_BSQ", BANDWIDTH="aic",
       MODEL="negbin")
endTime <- Sys.time()
#22.07 mins
#SAS: 41 segs

#Teste 25 (globalmin=no, adaptive_bsq, aic, negbin) --> muda os dados
startTime <- Sys.time()
mgwnbr(DATA=nakaya_NB, YVAR="vary", XVAR=c("varx1", "varx2"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="aic", MODEL="negbin")
endTime <- Sys.time()
#1.1 hora
#SAS: 1 min e 31 segs

#Teste 26 (globalmin=no, adaptive_bsq, aic, negbin) --> muda dados de novo
startTime <- Sys.time()
mgwnbr(DATA=nakaya_NB_std2, YVAR="vary", XVAR=c("varx1", "varx2"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="aic", MODEL="negbin")
endTime <- Sys.time()
#33.11 mins
#SAS: 52 segs

#Teste 27 (globalmin=no, adaptive_bsq, aic, nebin, mgwr=no)
startTime <- Sys.time()
mgwnbr(DATA=nakaya_NB, YVAR="vary", XVAR=c("varx1", "varx2"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="aic", MODEL="negbin", MGWR="no")
endTime <- Sys.time()
#11.82 mins
#SAS: 14 segs

#Teste 28 (globalmin=no, adaptive_bsq, aic, poisson, mgwr=no)
startTime <- Sys.time()
mgwnbr(DATA=nakaya_NB, YVAR="vary", XVAR=c("varx1", "varx2"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="aic", MODEL="poisson", MGWR="no")
endTime <- Sys.time()
#21.43 mins
#SAS: 14 segs

#Teste 29 (globalmin=no, fixed_g, aic, negbin, mgwr=no)
startTime <- Sys.time()
mgwnbr(DATA=nakaya_NB, YVAR="vary", XVAR=c("varx1", "varx2"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="FIXED_G",
       BANDWIDTH="aic", MODEL="negbin", MGWR="no")
endTime <- Sys.time()
#9.22 mins
#SAS: 13 segs

### LOGISTIC ###

#Teste 30 (globalmin=no, adaptive_bsq, cv, logistic, h=159)
startTime <- Sys.time()
mgwnbr(DATA=logistic, YVAR="Degree",
       XVAR=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="CV", MODEL="LOGISTIC", H=159)
endTime <- Sys.time()
#8.65 mins
#SAS: 13 segs
#conferir pop90 e eld
#mgwr stderror quantiles e descriptive
#parameters2

#Teste 31 (globalmin=no, fixed_g, cv, logistic, h=10^9)
startTime <- Sys.time()
mgwnbr(DATA=logistic, YVAR="Degree",
       XVAR=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="FIXED_G",
       BANDWIDTH="CV", MODEL="LOGISTIC", H=1000000000)
endTime <- Sys.time()
#1.12 mins
#SAS: 1.3 segs

#Teste 32 (globalmin=no, fixed_g, cv, logistic, mgwr=no, h=10^9)
startTime <- Sys.time()
mgwnbr(DATA=logistic, YVAR="Degree",
       XVAR=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="FIXED_G",
       BANDWIDTH="CV", MODEL="LOGISTIC", MGWR="no", H=1000000000)
endTime <- Sys.time()
#13.31 segs
#SAS: 0.38 segs

#Teste 33 (globalmin=no, adaptive_bsq, aic, logistic)
startTime <- Sys.time()
mgwnbr(DATA=logistic_std, YVAR="Degree",
       XVAR=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="aic", MODEL="LOGISTIC")
endTime <- Sys.time()
#
#SAS: 15 mins (DEU ERRO)

#Teste 34 (globalmin=no, adaptive_bsq, cv, logistic, h=159)
startTime <- Sys.time()
mgwnbr(DATA=logistic_std, YVAR="Degree",
       XVAR=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="CV", MODEL="LOGISTIC", H=159)
endTime <- Sys.time()
#10.8 mins
#SAS: 9 segs
#ENP TotPop90
#standard erros (PctEld e PctRural principalmente)
#parameters2 (std) (intercept, pcteld)

#Teste 35 (globalmin=no, fixed_g, cv, logistic)
startTime <- Sys.time()
mgwnbr(DATA=logistic_std, YVAR="Degree",
       XVAR=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov"),
       LAT="Y", LONG="X", GLOBAMIN="no", METHOD="FIXED_G",
       BANDWIDTH="CV", MODEL="LOGISTIC")
endTime <- Sys.time()
#
#SAS: 16 mins (DEU ERRO)

#Teste 36 (globalmin=no, fixed_g, aic, logistic)
startTime <- Sys.time()
mgwnbr(DATA=logistic_std, YVAR="Degree",
       XVAR=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="FIXED_G",
       BANDWIDTH="aic", MODEL="LOGISTIC")
endTime <- Sys.time()
#
#SAS: 1 min e 38 segs
 
#Teste 37 (globalmin=no, fixed_g, cv, logistic, h=10^9)
startTime <- Sys.time()
mgwnbr(DATA=logistic_std, YVAR="Degree",
       XVAR=c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov"),
       LAT="Y", LONG="X", GLOBALMIN="no", METHOD="FIXED_G",
       BANDWIDTH="CV", MODEL="LOGISTIC", H=1000000000)
endTime <- Sys.time()
#55.2 segs
#SAS: 1.58 segs

#obs.: no artigo/ manual, citar que para bn sempre tem que colocar aic.

#rodar de novo: testes 30 e 34 (depois de consertar logística)
#faltou rodar: 36 (depois de consertar logística)

#prints que faltam: 7 ao 9
#outputs para comparar: teste 18

#deram erro: 3, 33 e 35