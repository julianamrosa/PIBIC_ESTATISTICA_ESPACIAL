gwnbr <- function(y, x, lat, long, h, grid=NULL, latg, longg, method, gwr, offset=NULL, alphag=NULL, geocod){
  E <- 10
  COORD <- matrix(c(long, lat), ncol=2, byrow=F)
  n <- length(y)
  if (is.null(offset)){
    offset <- matrix(0, nrow=n, ncol=1)
  }
  if (is.null(grid)){
    POINTS <- matrix(c(long, lat), ncol=2, byrow=F)
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
      print(c(aux2, aux1, i, b, a))
    }
    if (is.null(alphag)){
      alphag <- a
    }
    else if(alphag==0){
      alphag <- E^-8
    }
    bg <- b
    parg <- par
  }
  if (gwr=="global"){
    print(c(alphag, aux2))
  }
  n <- length(y)
  aux2 <- 0
  library(rdist)
  distance <- cdist(POINTS, COORD, "euclidean")
  print(distance)
  seq <- 1:n
  for (i in 1:m){
    seqi <- matrix(i,nrow=n,ncol=1)
    distan <- cbind(cbind(seqi, seq), as.matrix(distance)[,i])
    w <- matrix(0,u,1)
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
        if (distan[jj,4]<= h){
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
    }
    while (abs(ddev)>0.000001 & cont<800){
      cont <- cont+1
      uj <- ifelse(uj>E^100,E^100,uj)
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
    print(c(jj, bi))
  }
}

setwd('~/PIBIC/golden_section_search')
data_gwnbr <- read.table('data_gwnbr.txt', header=T)

#Exemplo 1
gwnbr(y=data_gwnbr$fleet, x=data_gwnbr$industry, lat=data_gwnbr$Y, long=data_gwnbr$x, h=0.3344601, gwr="local", method="fixed", geocod=data_gwnbr$geocod)
