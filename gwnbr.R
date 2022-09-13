gwnbr <- function(y, x, lat, long, h, grid=NULL, latg, longg, method, gwr, offset=NULL, alphag=NULL, geocod){
  E <- 10
  COORD <- matrix(c(long, lat), ncol=2, byrow=F)
  n <- length(y)
  if (is.null(offset)){
    offset <- matrix(0, nrow=n, ncol=1)
  }
  if (is.null(grid)){
    POINTS <- matrix(c(long, lat), ncol=2, byrow=F)
    geocod_ <- geocod
  }
  else{
    POINTS <- matrix(c(longg, latg), ncol=2, byrow=F)
    geocod_ <- nrow(POINTS)
  }
  x <- cbind(matrix(1, nrow=n, ncol=1), x)
  yhat <- matrix(0, n, 1)
  m <- nrow(POINTS)
  bii <- matrix(0, ncol(x)*m, 2)
  alphaii <- matrix(0, m, 2)
  xcoord <- matrix(0, ncol(x)*m, 1)
  ycoord <- matrix(0, ncol(x)*m, 1)
  geocod <- matrix(0, ncol(x)*m, 1)
  sebi <- matrix(0, ncol(x)*m, 1)
  sealphai <- matrix(0, m, 1)
  S <- matrix(0, n, n)
  yp <- y-sum(y)/n
  probai <- matrix(0, m, 1)
  probbi <- matrix(0, m, 1)
  yhat <- matrix(0, m, 1)
  res <- matrix(0, m, 1)
  if (gwr!="poisson"){
    ym <- sum(y)/length(y)
    u <- (y+ym)/2
    n <- log(u);
    par <- 1
    ddpar <- 1
    j <- 0
    aux2 <- 0
    while (abs(ddpar)>0.00001){
      aux1 <- 0
      dpar <- 1
      parold <- par
      while (abs(dpar)>0.001){
        aux1 <- aux1+1
        if (par<0){
          par <- 0.00001
        }
        par <- ifelse(par<E^-10,E^-10,par)
        par <- as.numeric(par)
        g <- sum(digamma(par+y)-digamma(par)+log(par)+1-log(par+u)-(par+y)/(par+u))
        hess <- sum(trigamma(par+y)-trigamma(par)+1/par-2/(par+u)+(y+par)/((par+u)*(par+u)))
        hess <- ifelse(abs(hess)<E^-23,sign(hess)*E^-23,hess)
        hess <- ifelse(hess==0,E^-23,hess)
        par0 <- par
        par <- par0-solve(hess)*g
        #par <- as.numeric(par)
        if (aux1>50 & par>E^5){
          dpar <- 0.0001
          aux2 <- aux2+1
          if (aux2==1){
            par <- 2
          }	
          else if (aux2==2){
            par <- E^5
          }
          else if (aux2==3){
            par <- 0.0001
          }
        }
        else{
          dpar <- par-par0
        }
      }
      a <- 1/par
      dev <- 0
      ddev <- 1
      i <- 0
      while (abs(ddev)>0.00001 & i<800){
        i <- i+1
        w <- (u/(1+as.numeric(a)*u))+(y-u)*(as.numeric(a)*u/(1+2*as.numeric(a)*u+as.numeric(a)^2*u*u))
        w <- ifelse(w<=0, E^-5, w)
        z <- n+(y-u)/(w*(1+as.numeric(a)*u)) - as.numeric(offset)
        b <- solve(t(x*as.numeric(w))%*%x)%*%t(x*as.numeric(w))%*%as.numeric(z)
        n <- x%*%b + offset
        n <- ifelse(n>E^2, E^2, n)
        u <- exp(n)
        olddev <- dev
        u <- ifelse(u<E^-150,E^-150,u)
        tt <- y/u
        tt <- ifelse(tt==0,E^-10,tt)
        dev <- 2*sum(t(y*log(tt))-(y+1/a)*log((1+a%*%y)/as.numeric(1+as.numeric(a)*u)))
        ddev <- dev-olddev
      }
      if (aux2>4){
        ddpar <- E^-9
      }
      else{
        ddpar <- par-parold
      }
      #print(c(aux2, aux1, i, b, a))
    }
    if (is.null(alphag)){
      alphag <- a
    }
    else if(alphag==0){
      alphag <- E^-8
    }
    bg <- b
    parg <- par
    #print(c(alphag, bg, parg))
  }
  if (gwr=="global"){
    print(data.frame(alphag=alphag, aux2=aux2))
  }
  n <- length(y)
  aux2 <- 0
  library(rdist)
  distance <- cdist(POINTS, COORD, "euclidean")
  sequ <- 1:n
  for (i in 1:m){
    seqi <- matrix(i,nrow=n,ncol=1)
    distan <- cbind(cbind(seqi, sequ), as.matrix(distance)[,i])
    w <- matrix(0,n,1)
    if (method=="fixed"){
      for (jj in 1:n){
        w[jj] <- exp(-0.5*(distan[jj,3]/h)^2)
      }
    }
    else if (method=="adaptiven"){
      for (jj in 1:n){
        if (distan[jj,3]<=h){
          w[jj] <- (1-(distan[jj,3]/h)^2)^2
        }
        else{
          w[jj] <- 0
        }
      }
    }
    else if (method=="adaptive1"){
      w <- matrix(0,n,2)
      distan <- distan[order(distan[,3]),]
      distan <- cbind(distan,1:n)
      hn <- distan[h,3]
      for (jj in 1:n){
        if (distan[jj,4] <= h){
          w[jj,1] <- (1-(distan[jj,3]/hn)^2)^2
        }
        else{
          w[jj,1] <- 0
        }
        w[jj,2] <- distan[jj,2]
      }
      w <- w[order(w[, 2]), ]
    }
    wi <- w[,1]
    ym <- sum(y)/length(y)
    uj <- (y+ym)/2
    nj <- log(uj)
    ddpar <- 1
    jj <- 0
    count <- 0
    aux2 <- 0
    if (i==1 | aux2==5 | count==4){
      par <- 1
    }
    else{
      par <- alphaii[i-1,2]
    }
    while (abs(ddpar)>0.000001){
      dpar <- 1
      if(ddpar==1){
        parold <- 1.8139
      }
      else{
        parold <- par
      }
      aux1 <- 0
      if (gwr=="global" | gwr=="poisson"){
        dpar <- 0.00001
        if (gwr=="global"){
          par <- 1/alphag
        }
      }
      while (abs(dpar)>0.001){
        aux1 <- aux1+1
        if (gwr=="local"){
          par <- ifelse(par<E^-10,E^-10,par)
          g <- sum((digamma(par+y)-digamma(par)+log(par)+1-log(as.numeric(par)+uj)-(as.numeric(par)+y)/(as.numeric(par)+uj))*w[,1])
          hess <- sum((trigamma(par+y)-trigamma(par)+1/par-2/(as.numeric(par)+uj)+(y+par)/((as.numeric(par)+uj)*(as.numeric(par)+uj)))*w[,1])
        }
        par0 <- par
        hess <- ifelse(abs(hess)<E^-23,sign(hess)*E^-23,hess)
        hess <- ifelse(hess==0,E^-23,hess)
        par <- par0-solve(hess)*g
        if (par<=0){
          count <- count+1
          if (count==1){
            par <- 0.000001
          }
          else if(count==2){
            par <- 0.0001
          }
          else{
            par <- 1/alphag
          }
        }
        if (aux1>100 & par>E^5){
          dpar <- 0.0001
          if (aux2==0){
            par <- 1/alphag + 0.0011
          }
          if (aux2==1){
            par <- 2
          }
          else if (aux2==2){
            par <- E^5
          }
          else if (aux2==3){
            par <- 0.0001
          }
          aux2 <- aux2+1
        }
        else{
          dpar <- par-par0
          if (par<E^-3){
            dpar <- dpar*100
          }
        }
      }
      if (gwr=="poisson"){
        alpha <- 0
      }
      else{
        alpha <- 1/par
      }
      dev <- 0
      ddev <- 1
      cont <- 0
      while (abs(ddev)>0.000001 & cont<800){
        cont <- cont+1
        #uj <- ifelse(uj>E^100,E^100,uj)
        aux <- (as.numeric(alpha)*uj/(1+2*as.numeric(alpha)*uj+as.numeric(alpha%*%alpha)*uj*uj))
        Ai <- (uj/(1+as.numeric(alpha)*uj))+(y-uj)*aux
        Ai <- ifelse(Ai<=0,E^-5,Ai)
        zj <- nj+(y-uj)/(Ai*(1+as.numeric(alpha)*uj)) - offset
        if (det(t(x)%*%((as.numeric(wi*Ai))*x))<1){
          bi <- matrix(0,ncol(x),1)
        }
        else{
          bi <- solve(t(x)%*%(as.numeric(wi*Ai)*x))%*%t(x)%*%(as.numeric(wi*Ai)*zj)
        }
        nj <- x%*%bi + offset
        nj <- ifelse(nj>E^2,E^2,nj)
        uj <- exp(nj)
        olddev <- dev
        uj <- ifelse(uj<E^-150,E^-150,uj)
        tt <- y/uj
        tt <- ifelse(tt==0,E^-10,tt)
        if (gwr=="poisson"){
          dev <- 2*sum(y*log(tt)-(y-uj))
        }
        else{
          dev <- 2*sum(y*log(tt)-(y+1/alpha)*log((1+alpha*y)/(1+as.numeric(alpha)*uj)))
        }
        if (cont>100){
          ddev <- 0.0000001
        }
        else{
          ddev <- dev-olddev
        }
      }
      jj <- jj+1
      #print(c(jj, bi))
      if (gwr=="global" | gwr=="poisson" | aux2>4 | count>3 | jj>200){
        ddpar <- E^-9
      }
      else{
        ddpar <- par-parold
        if (par<E^-3){
          ddpar <- ddpar*100
        }
      }
      #print(c(j, aux1, cont, aux2, count, parold, par, ddpar))
    }
    if (aux2>4){
      probai[i] <- 1
    }
    if (count>3){
      probai[i] <- 2
    }
    Ai2 <- (uj/(1+as.numeric(alpha)*uj))+(y-uj)*(as.numeric(alpha)*uj/(1+2*as.numeric(alpha)*uj+as.numeric(alpha*alpha)*uj*uj))
    if (all(apply(Ai2, 2, min)<E^-5)){
      probbi[i] <- 1
      Ai2 <- ifelse(Ai2<E^-5,E^-5,Ai2)
    }
    if (is.null(grid)){
      if (det(t(x)%*%(as.numeric(wi*Ai)*x))<1){
        S[i,] <- matrix(0,1,n)
      }
      else{
        S[i,] <- x[i,]%*%solve(t(x)%*%(as.numeric(wi*Ai)*x))%*%t(x*as.numeric(wi*Ai))
      }
    }
    C <- solve(t(x)%*%(as.numeric(wi*Ai)*x))
    varb <- C
    seb <- sqrt(diag(varb))
    if (gwr!="poisson"){
      ser <- sqrt(1/abs(hess))
      r <- 1/alpha
      sealpha <- ser/(r^2)
      sealphai[i, 1] <- sealpha
      alphaii[i, 1] <- i
      alphaii[i, 2] <- alpha
    }
    m1 <- (i-1)*ncol(x)+1
    m2 <- m1+(ncol(x)-1)
    sebi[m1:m2, 1] <- seb
    bii[m1:m2, 1] <- i
    bii[m1:m2, 2] <- bi
    xcoord[m1:m2, 1] <- POINTS[i, 1]
    ycoord[m1:m2, 1] <- POINTS[i, 2]
    geocod[m1:m2, 1] <- geocod_[i]
    if (is.null(grid)){ #obs.: testar se grid==data?
      yhat[i] <- uj[i]
    }
  }
  tstat <- bii[, 2]/sebi
  probtstat <- 2*(1-pnorm(abs(tstat)))
  if (gwr!="poisson"){
    atstat <- alphaii[, 2]/sealphai
    aprobtstat <- 2*(1-pnorm(abs(atstat)))
  }
  else{
    atstat <- matrix(0, n, 1)
    aprobtstat <- matrix(1, n, 1)
  }
  b <- bii[, 2]
  alphai <- alphaii[, 2]
  id_ <- bii[, 1]
  ida_ <- alphaii[, 1]
  vec_bii <- c()
  for (linha in 1:nrow(bii)){
    vec_bii <- c(vec_bii, bii[, 1:2][linha, ])
  }
  beta_ <- matrix(vec_bii, n, byrow=T)
  i <- seq(2, ncol(beta_), 2)
  beta_ <- beta_[, i]
  qntl <- apply(beta_, 2, quantile, c(0.25, 0.5, 0.75))
  qntl <- rbind(qntl, qntl[3, ]-qntl[1, ])
  descriptb <- rbind(apply(beta_, 2, mean), apply(beta_, 2, min), apply(beta_, 2, max))
  print("quantis:")
  print(qntl)
  print("descritivas:")
  print(descriptb)
  vec_sebi <- c()
  for (linha in 1:nrow(sebi)){
    vec_sebi <- c(vec_sebi, sebi[linha, ])
  }
  stdbeta_ <- matrix(vec_sebi, n, byrow=T)
  qntls <- apply(stdbeta_, 2, quantile, c(0.25, 0.5, 0.75))
  qntls <- rbind(qntls, qntls[3, ]-qntls[1, ])
  descripts <- rbind(apply(stdbeta_, 2, mean), apply(stdbeta_, 2, min), apply(stdbeta_, 2, max))
  print("quantis:")
  print(qntls)
  print("descritivas:")
  print(descripts)
  yhat <- ifelse(yhat<E^-150, E^-150, yhat)
  tt <- y/yhat
  tt <- ifelse(tt==0, E^-10, tt)
  if (gwr=="poisson"){
    dev <- 2*sum(y*log(tt)-(y-yhat))
    tt2 <- y/mean(y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnull <- 2*sum(y*log(tt2)-(y-mean(y)))
    pctdev <- 1-dev/devnull
  }
  else{
    dev <- 2*sum(y*log(tt)-(y+1/alphai)*log((1+alphai*y)/(1+as.numeric(alphai)*yhat)))
    tt2 <- y/mean(y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnull <- 2*sum(y*log(tt2)-(y+1/alphai)*log((1+alphai*y)/(1+as.numeric(alphai)*mean(y))))
    pctdev <- 1-dev/devnull
  }
  if (gwr!="poisson"){
    a2 <- y+1/alphai
    b2 <- 1/alphai
    algamma <- matrix(0, n, 1)
    blgamma <- matrix(0, n, 1)
    for (i in 1:length(y)){
      algamma[i] <- lgamma(a2[i])
      blgamma[i] <- lgamma(b2[i])
    }
  }
  c2 <- y+1
  clgamma <- matrix(0, n, 1)
  for (i in 1:length(y)){
    clgamma[i] <- lgamma(c2[i])
  }
  if (gwr!="poisson"){
    ll <- sum(y*log(alphai*yhat)-(y+1/alphai)*log(1+alphai*yhat)+ algamma - blgamma - clgamma )
    if (gwr=="global" & all(alphai!=1/as.numeric(parg))){ #verificar
      npar <- sum(diag(S))
    }
    else{
      npar <- sum(diag(S))+1
    }
    tt <- y/(alphai*yhat)
    tt <- ifelse(tt==0, E^-10, tt)
    ll1 <- sum(y*log(tt)-y+(y+1/alphai)*log(1+alphai*yhat)-algamma+blgamma)
    tt2 <- y/mean(y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    llnull <- sum(y*log(tt2))
    pctll <- 1-ll1/llnull
  }
  else{
    ll <- sum(-yhat+y*log(yhat)-clgamma)
    npar <- sum(diag(S))
    pctll <- pctdev
  }
  adjpctdev <- 1-((length(y)-1)/(length(y)-npar))*(1-pctdev)
  adjpctll <- 1-((length(y)-1)/(length(y)-npar))*(1-pctll)
  resord <- y-yhat
  sigma2 <- (t(resord)%*%resord)/(n-npar)
  sii <- diag(S) #se necess?rio, matriz coluna
  res <- resord/sqrt(sigma2*(1-sii))
  res <- cbind(unique(id_), COORD[, 1], COORD[, 2], y, yhat, res, resord)
  AIC <- 2*npar-2*ll
  AICC <- AIC +(2*npar*(npar+1))/(n-npar-1)
  BIC <- npar*log(n)-2*ll
  malpha_ <- 0.05*(ncol(x)/npar)
  t_critical_ <- abs(qt(malpha_/2, n-npar))
  #print("malpha, tcritical, npar")
  print(data.frame(malpha=malpha_, tcritical=t_critical_, npar=npar))
  #print("gwr method ll dev pctdev adjpctdev pctll adjpctll npar aic aicc bic")
  print(c(gwr=gwr, method=method, ll=ll, dev=dev, pctdev=pctdev, adjpctdev=adjpctdev, pctll=pctll, adjpctll=adjpctll, npar=npar, aic=AIC, aicc=AICC, bic=BIC))
  res_ <<- as.data.frame(res)
  names(res_) <<- c("_id_", "xcoord", "ycoord", "yobs", "yhat", "res", "resraw")
  View(res_)
  stat <- cbind(ll, dev, pctdev, adjpctdev, pctll, adjpctll, npar, AIC, AICC, BIC)
  stat_ <<- as.data.frame(stat)
  View(stat_)
  bbeta <- cbind(id_, geocod, xcoord, ycoord, b, sebi, tstat, probtstat)
  bbeta_ <<- as.data.frame(bbeta)
  names(bbeta_) <<- c("_id_", "geocod", "xcoord", "ycoord", "b", "sebi", "tstat", "probtstat")
  View(bbeta_)
  xcoord <- COORD[, 1]
  ycoord <- COORD[, 2]
  geocod <- t(unique(geocod))
  sig_alpha <- matrix("not significant at 90%", n, 1)
  v1 <- npar
  for (i in 1:n){
    if (aprobtstat[i]<0.01*(ncol(x)/v1)){
      sig_alpha[i] <- "significant at 95%"
    }
    else if (aprobtstat[i]<0.1*(ncol(x)/v1)){
      sig_alpha[i] <- "significant at 90%"
    }
    else{
      sig_alpha[i] <- "not significant at 90%"
    }
  }
  aalpha <- cbind(ida_, as.numeric(geocod), xcoord, ycoord, alphai, as.numeric(sealphai), as.numeric(atstat), aprobtstat, sig_alpha, probai, probbi)
  alpha_ <<- as.data.frame(aalpha)
  names(alpha_) <<- c("_ida_", "geocod", "xcoord", "ycoord", "alphai", "sealphai", "atstat", "aprobtstat", "sig_alpha", "probai", "probbi")
  View(alpha_)
  tstat_ <- beta_/stdbeta_
  probt_ <- 2*(1-pnorm(abs(tstat_)))
  bistdt_ <- cbind(geocod_, COORD, beta_, stdbeta_, tstat_, probt_)
  parameters_ <<- as.data.frame(bistdt_)
  names(parameters_) <<- c("geocod", "x", "y", "Intercept", "INDUSTRY", "std_Intercept", "std_INDUSTRY", "tstat_Intercept", "tstat_INDUSTRY", "probt_Intercept", "probt_INDUSTRY")
  View(parameters_)
  sig_ <- matrix("not significant at 90%", n, ncol(x))
  v1 <- npar
  for (i in 1:n){
    for (j in 1:ncol(x)){
      if (probt_[i, j]<0.01*(ncol(x)/v1)){
        sig_[i, j] <- "significant at 99%"
      }
      else if (probt_[i, j]<0.05*(ncol(x)/v1)){
        sig_[i, j] <- "significant at 95%"
      }
      else if (probt_[i, j]<0.1*(ncol(x)/v1)){
        sig_[i, j] <- "significant at 90%"
      }
      else{
        sig_[i, j] <- "not significant at 90%"
      }
    }
  }
  sig_parameters2_ <<- as.data.frame(sig_)
  names(sig_parameters2_) <<- c("sig_Intercept", "sig_INDUSTRY")
  View(sig_parameters2_)
}

### Trocas
# _dist_ ---> distance
# seq ---> sequ
# dist ---> distan
# _id_ ---> id_
# _ida_ ---> ida_
# _beta_ ---> beta_

setwd('~/PIBIC/golden_section_search')
data_gwnbr <- read.table('data_gwnbr.txt', header=T)

#Exemplo 1
gwnbr(y=data_gwnbr$fleet, x=data_gwnbr$industry, lat=data_gwnbr$Y, long=data_gwnbr$x, h=0.3344601, gwr="local", method="fixed", geocod=data_gwnbr$geocod)

#Exemplo 2
gwnbr(y=data_gwnbr$fleet, x=data_gwnbr$industry, lat=data_gwnbr$Y, long=data_gwnbr$x, h=53.200292, gwr="local", method="fixed", geocod=data_gwnbr$geocod)

#Exemplo 3
gwnbr(y=data_gwnbr$fleet, x=data_gwnbr$industry, lat=data_gwnbr$Y, long=data_gwnbr$x, h=18, gwr="local", method="adaptive1", geocod=data_gwnbr$geocod)

#Exemplo 4
gwnbr(y=data_gwnbr$fleet, x=data_gwnbr$industry, lat=data_gwnbr$Y, long=data_gwnbr$x, h=49, gwr="local", method="adaptive1", geocod=data_gwnbr$geocod)

#Exemplo 5
gwnbr(y=data_gwnbr$fleet, x=data_gwnbr$industry, lat=data_gwnbr$Y, long=data_gwnbr$x, h=53.068412, gwr="global", method="fixed", geocod=data_gwnbr$geocod)

#Exemplo 6
gwnbr(y=data_gwnbr$fleet, x=data_gwnbr$industry, lat=data_gwnbr$Y, long=data_gwnbr$x, h=53.2, gwr="global", method="fixed", geocod=data_gwnbr$geocod)

#Exemplo 7
gwnbr(y=data_gwnbr$fleet, x=data_gwnbr$industry, lat=data_gwnbr$Y, long=data_gwnbr$x, h=34, gwr="global", method="adaptive1", geocod=data_gwnbr$geocod)

#Exemplo 8
gwnbr(y=data_gwnbr$fleet, x=data_gwnbr$industry, lat=data_gwnbr$Y, long=data_gwnbr$x, h=49, gwr="global", method="adaptive1", geocod=data_gwnbr$geocod)

#Exemplo 9
gwnbr(y=data_gwnbr$fleet, x=data_gwnbr$industry, lat=data_gwnbr$Y, long=data_gwnbr$x, h=9.3796134, gwr="poisson", method="fixed", geocod=data_gwnbr$geocod)

#Exemplo 10
gwnbr(y=data_gwnbr$fleet, x=data_gwnbr$industry, lat=data_gwnbr$Y, long=data_gwnbr$x, h=28.5, gwr="poisson", method="fixed", geocod=data_gwnbr$geocod)

#Exemplo 11
gwnbr(y=data_gwnbr$fleet, x=data_gwnbr$industry, lat=data_gwnbr$Y, long=data_gwnbr$x, h=5, gwr="poisson", method="adaptive1", geocod=data_gwnbr$geocod)

#Exemplo 12
gwnbr(y=data_gwnbr$fleet, x=data_gwnbr$industry, lat=data_gwnbr$Y, long=data_gwnbr$x, h=41, gwr="poisson", method="adaptive1", geocod=data_gwnbr$geocod)
