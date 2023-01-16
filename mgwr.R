## MGWR ##

library(sf)
library(foreign)
library(dplyr)
setwd('~/PIBIC/multiscale')

georgia <- st_read('G_utm.shp')
georgiashp <- st_read('Georgia.shp')
georgiashp$AREAKEY <- as.numeric(georgiashp$AREAKEY) #6.
georgia_data <- read.dbf('G_utm1.dbf') #replace
georgia_data <- read.csv('GData_utm.csv') #replace
# 2 iguais?

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
                 BANDWIDTH="CV", OFFSET=NULL, DISTANCEKM="NO",
                 INT=50, H){
  #lembrar de tornar band global
  E <- 10
  y <- DATA[, YVAR]
  x <- DATA[, which(names(DATA) %in% XVAR)]
  n <- length(y)
  wt <- matrix(1, n, 1)
  if (!is.null(WEIGHT)){
    wt <- as.matrix(WEIGHT)
  }
  offset <- matrix(0, n, 1)
  if (!is.null(OFFSET)){
    offset <- as.matrix(OFFSET)
  }
  x <- as.matrix(cbind(rep(1, n), x))
  nvarg <- ncol(x);
  yhat <- matrix(0, n, 1)
  bi <- matrix(0, nvarg*n, 4)
  alphai <- matrix(0, n,3)
  S <- matrix(0, n, 1)
  mRj <- matrix(0, n, n*nvarg)
  Sm <- matrix(0, n, n)
  Sm3 <- matrix(0, n, nvarg)
  Rj <- matrix(0, n, n)
  Cm <- matrix(0, n, n*nvarg)
  stdbm <- matrix(0, n, nvarg)
  mAi <- matrix(0, n, nvarg)
  ENP <- matrix(0, 1, nvarg+2)
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
      if (toupper(MODEL)=="NEGBIN"){
           
        if (cont>1){
          parg <- 1/(sum((y-uj)^2/uj)/(n-nvarg))
        }
          while (abs(dpar)>0.000001 & cont1<200){
            parg <- choose(parg<E^-10, E^-10, parg)
            g <- sum(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj))
            hess <- sum(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)^2)
            hess <- choose(hess==0, E^-23, hess)
            par0 <- parg
            parg <- par0-solve(hess)%*%g
            if (cont1>50 & parg>1E5){
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
        uj <- choose(uj>E^100, E^100, uj)
        Ai <- (uj/(1+alphag*uj))+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag**2*uj*uj))
        Ai <- ifelse(Ai<=0, E^-5, Ai)	
        zj <- nj+(y-uj)/(Ai*(1+alphag*uj))-offset
        if (det(t(x)*(Ai*x))=0){
          bg <- matrix(0, nvarg, 1)
        }
        else{
          bg=inv(t(x)*(Ai*x))*t(x)*(Ai*zj)
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
        if (toupper(MODEL)=="NEGBIN"){
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
    varg <- vecdiag(inv(t(x*wt*Ai)%*%x))
  }
  #linha 136
  if (toupper(MODEL)=="LOGISTIC"){
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
      if (det(t(x)*(wt*Ai*x))==0){
        bg <- matrix(0, nvarg, 1)
      }
      else{
        bg <- solve(t(x)*(wt*Ai*x))*t(x)*(wt*Ai*zj)
      }
      nj <- x*bg
      nj <- ifelse(nj>E^2, E^2, nj)
      uj <- exp(nj)/(1+exp(nj))
      olddev <- devg
      uj <- ifelse(uj<E^-150, E^-150, uj)
      tt <- y/uj
      tt <- ifelse(tt=0, E^-10, tt)
      uj <- ifelse(uj==1, 0.99999, uj)
      tt2 <- (1-y)/(1-uj)
      tt2 <- ifelse(tt2==0, E^-10, tt2)
      devg <- 2*sum((y*log(tt))+(1-y)*log(tt2))
      ddev <- devg-olddev
      cont <- cont+1
    }
    ujg <- uj
    yhat <- uj
    varg <- vecdiag(inv(t(x*wt*Ai)*x))
  }
  #linha 167
}


## Testes ##
mgwr(DATA=georgia_data_std, YVAR="PctBach",
     XVAR=c("PctBlack", "PctFB", "TotPop90", "PctEld"),
     LAT=y, LONG=x, GLOBALMIN="no", METHOD="adaptive_bsq",
     BANDWIDTH="cv", MODEL="gaussian")
