setwd('C:/Users/e-juliana.rosa/Documents/PIBIC')
georgia_data <- read.csv('GData_utm.csv')

#Discretizando y e padronizando x
georgia_nb <- georgia_data
georgia_nb$PctBach <- as.integer(georgia_nb_std$PctBach)

georgia_nb_std <- georgia_nb

#TotPop90
georgia_nb_std$TotPop90 <- (georgia_nb_std$TotPop90-
                              mean(georgia_nb_std$TotPop90))/
  sd(georgia_nb_std$TotPop90)
#PctRural
georgia_nb_std$PctRural <- (georgia_nb_std$PctRural-
                              mean(georgia_nb_std$PctRural))/
  sd(georgia_nb_std$PctRural)
#PctEld
georgia_nb_std$PctEld <- (georgia_nb_std$PctEld-
                            mean(georgia_nb_std$PctEld))/
  sd(georgia_nb_std$PctEld)
#PctFB
georgia_nb_std$PctFB <- (georgia_nb_std$PctFB-
                           mean(georgia_nb_std$PctFB))/
  sd(georgia_nb_std$PctFB)
#PctPov
georgia_nb_std$PctPov <- (georgia_nb_std$PctPov-
                            mean(georgia_nb_std$PctPov))/
  sd(georgia_nb_std$PctPov)
#PctBlack
georgia_nb_std$PctBlack <- (georgia_nb_std$PctBlack-
                              mean(georgia_nb_std$PctBlack))/
  sd(georgia_nb_std$PctBlack)

gwnbr <- function(DATA, YVAR, XVAR, WEIGHT=NULL, LAT, LONG,
                   GLOBALMIN="yes", METHOD,
                   BANDWIDTH="cv", OFFSET=NULL,
                   DISTANCEKM="no", INT=50, H=NULL){
  output <- list()
  header <- c()
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
      devg <- 2*sum(Y*log(tt)-(Y+1/alphag)*log((1+alphag*Y)/(1+alphag*uj)))
      sealphag <- sqrt(1/abs(hess))/(parg^2)
      if (cont2>100){
        ddev <- 0.0000001
      }
      else{
        ddev <- devg-olddev
      }
      cont2 <- cont2+1
    }
    ujg <<- uj
    yhat <- uj
    cont <- cont+1
    ddpar <- parg-parold
  }
  varg <- diag(solve(t(X*wt*ai)%*%X))
  LONG <- DATA[, LONG]
  LAT <- DATA[, LAT]
  COORD <<- matrix(c(LONG, LAT), ncol=2)
  distance <- dist(COORD, "euclidean")
  sequ <<- 1:N
  cv <- function(h, y, x, fi){ #ujg, coord, Offset, alphag
    nvar <- ncol(x)
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
          nj <- x%*%b+yhat_beta-fi
          nj <- ifelse(nj>E^2, E^2, nj)
          uj <- exp(nj)
          olddev <- dev
          uj <- ifelse(uj<E^-150, E^-150, uj)
          tt <- y/uj
          tt <- ifelse(tt==0, E^-10, tt)
          dev <- 2*sum(y*log(tt)-(y+1/alpha)*log((1+alpha*y)/(1+alpha*uj)))
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
      yhat[i] <<- uj[i]
      alphai[i, 2] <<- alpha
      if (det(t(x)%*%(w*ai*x*wt))==0){
        s[i] <<- 0
      }
      else{
        s[i] <<- (x[i, ]%*%solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*ai*wt))[i]
      }
      next
    }
    CV <- t((y-yhat)*wt)%*%(y-yhat)
    ll <- sum(y*log(alphai[,2]*yhat)-(y+1/alphai[,2])*log(1+alphai[,2]*yhat)+lgamma(y+1/alphai[,2])-lgamma(1/alphai[,2])-lgamma(y+1))
    npar <- sum(s)+sum(s)/nvar
    AIC <- 2*npar-2*ll
    AICC <- AIC +(2*npar*(npar+1))/(N-npar-1)
    if (toupper(BANDWIDTH)=="AIC"){
      CV <- AICC
    }
    #free dist w
    res <- cbind(CV, npar)
    return (res)
  }
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
    }
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
      }
    }
    if (toupper(GLOBALMIN)=="YES"){
      xming <- xmin[which(xmin[,1]==min(xmin[,1])),2]
    }
    bandwidth <- xming
    return (bandwidth)
  }
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
        w <- w[order(w[, 2]), ]
        w <- w[,1]
      }
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
          nj <- x%*%b+yhat_beta-fi
          nj <- ifelse(nj>E^2, E^2, nj)
          uj <- as.numeric(exp(nj))
          olddev <- dev
          uj <- ifelse(uj<E^-150, E^-150, uj)
          tt <- y/uj
          tt <- ifelse(tt==0, E^-10, tt)
          dev <- 2*sum(y*log(tt)-(y+1/alpha)*log((1+alpha*y)/(1+alpha*uj)))
          cont2 <- cont2+1
        }
        cont <- cont+1
        ddpar <- par-parold
      }
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
      hess <- sum(w*wt*(trigamma(par+y)-trigamma(par)+1/par-2/(par+exp(yhat_beta))+(y+par)/(par+exp(yhat_beta))^2))
      hess <- sum(w*wt*(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)^2))
      hess <- ifelse(hess==0, E^-23, hess)
      sealpha <- sqrt(1/abs(hess))/(par^2)
      alphai[i,1] <<- i
      alphai[i,2] <<- alpha
      alphai[i,3] <<- sealpha
      m1 <- (i-1)*nvar+1
      m2 <- m1+(nvar-1)
      bim[m1:m2] <- b
      yhatm[i] <- uj[i]
      yhat[i] <<- uj[i]
    }
    beta <- matrix(bim, N, byrow=T)
    yhbeta <- cbind(yhatm, beta)
    return (yhbeta)
  }
  finb <- rep(0, N)
  yhat_beta <<- Offset
  if (!is.null(H)){
    hh <- H
  }
  else{
    hh <- GSS(Y,X,finb)
  }
  header <- append(header, "General Bandwidth")
  output <- append(output, hh)
  names(output) <- "general_bandwidth"
  yhat_beta <<- gwr(hh,Y,X,finb)
  beta <- yhat_beta[,2:(nvarg+1)]
  Fi <- X*beta
  mband <- hh
  Sm2 <- sm
  v1 <- sum(diag(sm))
  for (jj in 1:nvarg){
    m1 <- (jj-1)*N+1
    m2 <- m1+(N-1)
    ENP[jj] <- sum(diag(mrj[,m1:m2]))
    ENP[jj] <- sum(diag(sm))
    stdbm[,jj] <- sqrt(sm3[,jj])
  }
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
  stats_measures <- c(dev, ll, pctdev, adjpctdev, AIC,
                      AICc)
  names(stats_measures) <- c("deviance",
                             "full_Log_likelihood",
                             "pctdev", "adjpctdev",
                             "AIC", "AICc")
  header <- append(header, "Measures")
  output <- append(output, list(stats_measures))
  names(output)[length(output)] <- "measures"
  ENP[nvarg+1] <- sum(diag(sm))
  ENP[nvarg+2] <- sum(diag(Sm2))
  ENP <- c(ENP, (v1/nvarg))
  names(ENP) <- c('Intercept', XVAR, 'MGWR', 'GWR', 'alpha')
  header <- append(header, "ENP")
  output <- append(output, list(ENP))
  names(output)[length(output)] <- "ENP"
  dff <- N-v1
  tstat <- beta/stdbm
  probt <- 2*(1-pt(abs(tstat), dff))
  malpha <- ENP
  malpha[1:nvarg] <- 0.05/ENP[1:nvarg]
  malpha[nvarg+1] <- 0.05*(nvarg/v1)
  malpha[nvarg+2] <- 0.05*(nvarg/sum(diag(Sm2)))
  malpha[1:nvarg] <- 0.05*(nvarg/v1)
  malpha[nvarg+3] <- 0.05*(nvarg/v1)
  t_critical <- abs(qt(malpha/2,dff))
  beta2 <- beta
  alpha <- alphai[,2]
  beta2 <- cbind(beta, alpha)
  qntl <- apply(beta2, 2, quantile, c(0.25, 0.5, 0.75)) #, na.rm=T
  IQR <- (qntl[3,]-qntl[1,])
  qntl <- rbind(round(qntl, 6), IQR=round(IQR, 6))
  descriptb <- rbind(apply(beta2, 2, mean), apply(beta2, 2, min), apply(beta2, 2, max))
  rownames(descriptb) <- c('Mean', 'Min', 'Max')
  colnames(qntl) <- c('Intercept', XVAR, 'alpha')
  header <- append(header, "Quantiles of MGWR Parameter Estimates")
  output <- append(output, list(qntl))
  names(output)[length(output)] <- "qntls_mgwr_param_estimates"
  colnames(descriptb) <- c('Intercept', XVAR, 'alpha')
  header <- append(header, "Descriptive Statistics")
  output <- append(output, list(descriptb))
  names(output)[length(output)] <- "descript_stats_mgwr_param_estimates"
  stdbeta <- stdbm
  stdbeta2 <- stdbeta
  stdalpha <- alphai[,3]
  stdbeta2 <- cbind(stdbeta, stdalpha)
  qntls <- apply(stdbeta2, 2, quantile, c(0.25, 0.5, 0.75)) #, na.rm=T
  IQR <- (qntls[3,]-qntls[1,])
  qntls <- rbind(round(qntls, 6), IQR=round(IQR, 6))
  descripts <- rbind(apply(stdbeta2, 2, mean), apply(stdbeta2, 2, min), apply(stdbeta2, 2, max))
  rownames(descripts) <- c('Mean', 'Min', 'Max')
  header <- append(header, "alpha-level=0.05")
  output <- append(output, list(malpha))
  names(output)[length(output)] <- "p_values"
  t_critical <- round(t_critical, 2)
  header <- append(header, "t-Critical")
  output <- append(output, list(t_critical))
  names(output)[length(output)] <- "t_critical"
  colnames(qntls) <- c('Intercept', XVAR, 'alpha')
  header <- append(header, "Quantiles of MGWR Standard Errors")
  output <- append(output, list(qntls))
  names(output)[length(output)] <- "qntls_mgwr_se"
  colnames(descripts) <- c('Intercept', XVAR, 'alpha')
  header <- append(header, "Descriptive Statistics of Standard Errors")
  output <- append(output, list(descripts))
  names(output)[length(output)] <- "descripts_stats_se"
  #### global estimates ####
  if (is.null(WEIGHT)){
    vargd <- varg
    dfg <- N-nrow(bg)
  }
  stdg <- matrix(sqrt(vargd))
  bg <- rbind(bg, alphag)
  stdg <- rbind(stdg, sealphag)
  dfg <- dfg-1
  tg <- bg/stdg
  probtg <- 2*(1-pt(abs(tg), dfg))
  bg_stdg_tg_probtg <- cbind(bg, stdg, tg, probtg)
  rownames(bg_stdg_tg_probtg) <- c('Intercept', XVAR, 'alpha')
  colnames(bg_stdg_tg_probtg) <- c("Par. Est.", "Std Error", "t Value", "Pr > |t|")
  header <- append(header, "Global Parameter Estimates")
  output <- append(output, list(bg_stdg_tg_probtg))
  names(output)[length(output)] <- "global_param_estimates"
  header <- append(header, "NOTE: The denominator degrees of freedom for the t tests is...")
  output <- append(output, list(dfg))
  names(output)[length(output)] <- "t_test_dfs"
  yhatg <- exp(X%*%bg[1:(nrow(bg)-1)]+Offset)
  ll <- sum(Y*log(alphag*yhatg)-(Y+1/alphag)*log(1+alphag*yhatg)+lgamma(Y+1/alphag)-lgamma(1/alphag)-lgamma(Y+1))
  AIC <- -2*ll+2*(nvarg+1)
  AICc <- -2*ll+2*(nvarg+1)*(N/(N-(nvarg+1)-1))
  tt2 <- Y/mean(Y)
  tt2 <- ifelse(tt2==0, E^-10, tt2)
  devnullg <- 2*sum(Y*log(tt2)-(Y+1/alphag)*log((1+alphag*Y)/(1+alphag*mean(Y))))
  pctdevg <- 1-devg/devnullg
  adjpctdevg <- 1-((N-1)/(N-nvarg))*(1-pctdevg)
  global_measures <- c(devg, ll, pctdevg, adjpctdevg, AIC, AICc)
  names(global_measures) <- c('deviance', 'full_Log_likelihood', 'pctdevg',
                              'adjpctdevg', 'AIC', 'AICc')
  header <- append(header, "Global Measures")
  output <- append(output, list(global_measures))
  names(output)[length(output)] <- "global_measures"
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
  
  ###################################
  min_bandwidth <<- as.data.frame(t(mband))
  names(min_bandwidth) <<- 'Intercept'
  parameters2 <<- cbind(parameters2, sig_parameters2)
  Alpha <<- cbind(Alpha, sig_alpha)
  i <- 1
  for (element in output){
    cat(header[i], "\n")
    print(element)
    i <- i+1
  }
  invisible(output)
}

#Padronizando
startTime <- Sys.time()
out <- gwnbr(DATA=georgia_nb_std, YVAR="PctBach",
              XVAR=c("TotPop90", "PctRural", "PctEld",
                     "PctFB", "PctPov", "PctBlack"), LAT="Y",
              LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
              BANDWIDTH="aic")
endTime <- Sys.time()
#1.05 mins

#Sem padronizar
startTime <- Sys.time()
out <- gwnbr(DATA=georgia_nb, YVAR="PctBach",
              XVAR=c("TotPop90", "PctRural", "PctEld",
                     "PctFB", "PctPov", "PctBlack"), LAT="Y",
              LONG="X", GLOBALMIN="no", METHOD="ADAPTIVE_BSQ",
              BANDWIDTH="aic")
endTime <- Sys.time()
#1.05 mins
