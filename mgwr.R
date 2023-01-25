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

#abaixo, passar x e y como strings --> nomes das variáveis
mgwr <- function(DATA, YVAR, XVAR, WEIGHT=NULL, LAT, LONG,
                 GLOBALMIN="yes", METHOD, MODEL="GAUSSIAN",
                 MGWR="yes", BANDWIDTH="CV", OFFSET=NULL,
                 DISTANCEKM="NO", INT=50, H){
  #lembrar de tornar band global
  E <- 10
  y <- DATA[, YVAR]
  x <- DATA[, which(names(DATA) %in% XVAR)]
  n <- length(y)
  wt <-rep(1, n)
  if (!is.null(WEIGHT)){
    wt <- as.matrix(WEIGHT)
  }
  offset <- rep(0, n)
  if (!is.null(OFFSET)){
    offset <- as.matrix(OFFSET)
  }
  x <- as.matrix(cbind(rep(1, n), x))
  nvarg <- ncol(x)
  yhat <- rep(0, n)
  bi <- matrix(0, nvarg*n, 4)
  alphai <- matrix(0, n, 3)
  S <- rep(0, n)
  mRj <- matrix(0, n, n*nvarg)
  Sm <- matrix(0, n, n)
  Sm3 <- matrix(0, n, nvarg)
  Rj <- matrix(0, n, n)
  Cm <- matrix(0, n, n*nvarg)
  stdbm <- matrix(0, n, nvarg)
  mAi <- matrix(0, n, nvarg)
  ENP <- rep(0, nvarg+2)
  #estimativas globais:
  if (toupper(MODEL)=="POISSON" | toupper(MODEL)=="NEGBIN"){
    uj <- (y+mean(y))/2
    nj <- log(uj)
    parg <- sum((y-uj)^2/uj)/(n-nvarg)
    ddpar <- 1
    cont <- 1
    while (abs(ddpar)>0.000001 & cont<100){
      dpar <- 1
      parold <- parg
      cont1 <- 1
      cont3 <- 1
      if(toupper(MODEL)=="POISSON"){
      alphag <- E^-6
      parg <- 1/alphag
      }
      else{
        if (cont>1){
          parg <- 1/(sum((y-uj)^2/uj)/(n-nvarg))
        }
          while (abs(dpar)>0.000001 & cont1<200){
            parg <- ifelse(parg<E^-10, E^-10, parg)
            g <- sum(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj))
            hess <- sum(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)^2)
            hess <- ifelse(hess==0, E^-23, hess)
            par0 <- parg
            parg <- par0-solve(hess)%*%g
            if (cont1>50 & parg>E^5){
              dpar <- 0.0001
              cont3 <- cont3+1
              if (cont3==1){
                parg <- 2
              }
              else if (cont3==2){
                parg <- E^5
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
        uj <- ifelse(uj>E^100, E^100, uj)
        Ai <- (uj/(1+alphag*uj))+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj*uj))
        Ai <- ifelse(Ai<=0, E^-5, Ai)	
        zj <- nj+(y-uj)/(Ai*(1+alphag*uj))-offset
        if (det(t(x)%*%(Ai*x))==0){
          bg <- rep(0, nvarg)
        }
        else{
          bg <- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
        }
        nj <- x*bg+offset
        nj <- ifelse(nj>E^2, E^2, nj)
        uj <- exp(nj)
        olddev <- devg
        uj <- ifelse(uj<E^-150, E^-150, uj)
        tt <- y/uj
        tt <- ifelse(tt==0, E^-10, tt)
        if (toupper(MODEL)=="POISSON"){
          devg <- 2*sum(y*log(tt)-(y-uj))
        }
        else{
          devg <- 2*sum(y*log(tt)-(y+1/alphag)*log((1+alphag*y)/(1+alphag*uj)))
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
      ujg <- uj
      yhat <- uj
      cont <- cont+1
      ddpar <- parg-parold
    }
    varg <- vecdiag(solve(t(x*wt*Ai)%*%x))
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
      if (det(t(x)%*%(wt*Ai*x))==0){
        bg <- rep(0, nvarg)
      }
      else{
        bg <- solve(t(x)%*%(wt*Ai*x))%*%t(x)%*%(wt*Ai*zj)
      }
      nj <- x*bg
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
    ujg <- uj
    yhat <- uj
    varg <- vecdiag(solve(t(x*wt*Ai)%*%x))
  }
  #linha 167
  COORD <- matrix(c(LONG, LAT), ncol=2)
  distance <- dist(COORD, "euclidean")
  seq <- 1:n
  cv <- function(h, y, x, fi, n=n, wt=wt, ujg=ujg,
                 yhat=yhat, coord=coord, distance=distance,
                 seq=seq, offset=offset, alphag=alphag,
                 alphai=alphai, S=S, parg=parg,
                 yhat_beta=yhat_beta){
    nvar <- ncol(x)
    for (i in 1:n){
      #linhas 214, 288, 331 fecham esse loop no sas
      for (j in 1:n){
        seqi <- rep(i, n)
        distan <- cbind(seqi, t(seq), distance[,i])
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
      AICc <- 2*n*log(CV/n)+n*log(2*3.15159)+n*(n+npar)/(n-2-npar)
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
  GSS <- function(depy, indepx, fix, distance=distance, n=n){
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
}

#ifs por else ifs
#Converter matrizes unidimensionais em vetores!
#conferir ; choose j inv = * E

## Testes ##
mgwr(DATA=georgia_data_std, YVAR="PctBach",
     XVAR=c("PctBlack", "PctFB", "TotPop90", "PctEld"),
     LAT=y, LONG=x, GLOBALMIN="no", METHOD="adaptive_bsq",
     BANDWIDTH="cv", MODEL="gaussian")
