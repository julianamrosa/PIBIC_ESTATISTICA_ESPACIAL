############### GOLDEN SECTION SEARCH ###############
golden <- function(y, x, lat, long, method, type, gwr, offset=NULL){
  E <- 10
  COORD <- matrix(c(long, lat), ncol=2, byrow=F)
  n <- length(y)
  if (is.null(offset)){
    offset <- matrix(0, nrow=n, ncol=1)
  }
  x <- cbind(matrix(1, nrow=n, ncol=1), x)
  library(kableExtra)
  library(dplyr)
  tabela1 <- matrix(c(method, type, gwr), 1, 3)%>%
    kbl(caption = "", col.names=c("method", "type", "gwr"), align=c("c", "c", "c"))%>%
    kable_classic(full_width = F, html_font = "Cambria", position="left")
  print(tabela1)
  distance <- dist(COORD, "euclidean")
  maxd <- max(distance)
  if (method=="adaptive1"){
    h0 <- 5
    h3 <- n
  }
  else if (method == "adaptiven" | method == "fixed"){
    h0 <- 0
    h3 <- maxd
  }
  r <- 0.61803399
  c <- 1-r
  if (method == "adaptive1"){
    tol <- 0.9
  }
  else{
    tol=0.1
  }
  h1 <- h0+(1-r)*(h3-h0)
  h2 <- h0+r*(h3-h0)
  tabela2 <- matrix(c(h0, h1, h2, h3), 1, 4)%>%
    kbl(caption = "", col.names=c("h0", "h1", "h2", "h3"), align=c("c", "c", "c", "c"))%>%
    kable_classic(full_width = F, html_font = "Cambria", position="left")
  print(tabela2)
  cv <- function(h, method, n, coord, x, y, type, maxd, gwr, offset, distance){
    E <- 10
    alphaii <- matrix(0, nrow=n, ncol=2)
    yhat <- matrix(0, nrow=n, ncol=1)
    S <- matrix(0, nrow=n,ncol=n)
    if (gwr=="global"){
      ym=sum(y)/length(y)
      u=(y+ym)/2
      n=log(u);
      par=1
      ddpar=1
      j=0
      aux2=0
      while (abs(ddpar)>0.00001){
        aux1=0
        dpar=1
        parold=par
        while (abs(dpar)>0.001){
          aux1=aux1+1
          if (par<0){
            par=0.00001
          }
          par=ifelse(par<E^-10,E^-10,par)
          par <- as.numeric(par)
          g=sum(digamma(par+y)-digamma(par)+log(par)+1-log(par+u)-(par+y)/(par+u))
          hess=sum(trigamma(par+y)-trigamma(par)+1/par-2/(par+u)+(y+par)/((par+u)*(par+u)))
          hess=ifelse(abs(hess)<E^-23,sign(hess)*E^-23,hess)
          hess=ifelse(hess==0,E^-23,hess)
          par0=par
          par=par0-solve(hess)*g
          if (aux1>50 & par>E^5){
            dpar= 0.0001
            aux2=aux2+1
            if (aux2==1){
              par=2
            }	
            else if (aux2==2){
              par=E^5
            }
            else if (aux2==3){
              par=0.0001
            }
          }
          else {
            dpar=par-par0
          }
        }
        a=1/par
        dev=0
        ddev=1
        i=0
        while (abs(ddev)>0.00001 & i<800){
          i=i+1
          w=(u/(1+as.numeric(a)*u))+(y-u)*(as.numeric(a)*u/(1+2*as.numeric(a)*u+as.numeric(a)^2*u*u))
          z=n+(y-u)/(w*(1+as.numeric(a)*u)) - as.numeric(offset)
          b=solve(t(x*as.numeric(w))%*%x)%*%t(x*as.numeric(w))%*%as.numeric(z)
          n=x%*%b + offset
          u=exp(n)
          olddev=dev
          u=ifelse(u<E^-150,E^-150,u)
          tt=y/u
          tt=ifelse(tt==0,E^-10,tt)
          dev=2*sum(t(y*log(tt))-(y+1/a)*log((1+a%*%y)/as.numeric(1+as.numeric(a)*u)))
          ddev=dev-olddev
        }
        if (aux2>4){
          ddpar=1E-9
        }
        else{
          ddpar=par-parold
        }
      }
      alpha=a
    }
    n=length(y)
    aux2=0
    sequ=1:n
    for (i in 1:n){
      seqi=matrix(i,nrow=n,ncol=1)
      distan=cbind(cbind(seqi, sequ), as.matrix(distance)[,i])
      u=nrow(distan)
      w=matrix(0,u,1)
      if (method=="fixed"){
        if (type=="cv"){
          for (jj in 1:u){
            if (distan[jj,3]<=maxd*0.8 & distan[jj,3]!=0){
              w[jj]=exp(-0.5*(distan[jj,3]/h)^2)
            }
            else{
              w[jj]=0
            }
          }
        }
        else{
          for (jj in 1:u){
            if (distan[jj,3]<=maxd*0.8) {
              w[jj]=exp(-0.5*(distan[jj,3]/h)^2)
            }
            else{
              w[jj]=0
            }
          }
        }
      }
      else if (method=="adaptiven"){
        if (type=="cv"){
          for (jj in 1:u){
            if (distan[jj,3]<=h & distan[jj,3]!=0){
              w[jj]=(1-(distan[jj,3]/h)**2)**2
            }
            else{
              w[jj]=0
            }
          }
        }
        else{
          for (jj in 1:u){
            if (distan[jj,3]<=h){
              w[jj]=(1-(distan[jj,3]/h)**2)**2
            }
            else{
              w[jj]=0
            }
          }
        }
      }
      else if (method=="adaptive1"){
        distan = distan[order(distan[,3]),]
        distan=cbind(distan,1:n)
        w=matrix(0,n,2)
        hn=distan[h,3]
        if (type=="cv"){
          for (jj in 1:n){
            if (distan[jj,4]<= h & distan[jj,3]!=0){
              w[jj,1]=(1-(distan[jj,3]/hn)**2)**2
            }
            else{
              w[jj,1]=0
            }
            w[jj,2]=distan[jj,2]
          }
        }
        else{
          for (jj in 1:n){
            if (distan[jj,4]<=h){
              w[jj,1]=(1-(distan[jj,3]/hn)**2)**2
            }
            else{
              w[jj,1]=0
            }
            w[jj,2]=distan[jj,2]
          }
        }
        w = w[order(w[, 2]), ]
      }
      wi=w[,1]
      ym=sum(y)/length(y)
      uj=(y+ym)/2
      nj=log(uj)
      if (i==1 | aux2==5){
        par=1
      }
      else{
        par=alphaii[i-1,2]
      }
      ddpar=1
      jj=0
      count=0
      aux2=0
      while (abs(ddpar)>0.000001){
        aux1=0
        dpar=1
        parold=par
        if (gwr=="global" | gwr=="poisson"){
          dpar=0.00001
          if (gwr=="global"){
            par=1/a
          }
        }
        # computing alpha=1/par, where par=theta #
        while (abs(dpar)>0.001){
          aux1=aux1+1
          if (gwr=="local"){
            par=ifelse(par<E^-10,E^-10,par)
            g=sum((digamma(par+y)-digamma(par)+log(par)+1-log(as.numeric(par)+uj)-(as.numeric(par)+y)/(as.numeric(par)+uj))*w[,1])
            hess=sum((trigamma(par+y)-trigamma(par)+1/par-2/(as.numeric(par)+uj)+(y+par)/((as.numeric(par)+uj)*(as.numeric(par)+uj)))*w[,1])
          }
          hess=ifelse(abs(hess)<E^-23,sign(hess)*E^-23,hess)
          hess=ifelse(hess==0,E^-23,hess)
          par0=par
          par=par0-solve(hess)*g
          if (par<=0){
            count=count+1
            if (count<10){
              par=0.000001
            }
            else{
              par=abs(par)
            }
          }
          if (aux1>50 & par>E^5){
            dpar= 0.0001
            aux2=aux2+1
            if (aux2==1){
              par=2
            }
            else if (aux2==2){
              par=1E5
            }
            else if (aux2==3){
              par=0.0001
            }
          }
          else{
            dpar=par-par0
            if (par<E^-3){
              dpar=dpar*100
            }
          }
        }
        if (gwr=="poisson"){
          alpha=0
        }
        else{
          alpha=1/par
        }
        dev=0
        ddev=1
        cont=0
        # computing beta #
        while (abs(ddev)>0.000001 & cont<800){
          cont=cont+1
          uj=ifelse(uj>E^100,E^100,uj)
          aux=(as.numeric(alpha)*uj/(1+2*as.numeric(alpha)*uj+as.numeric(alpha%*%alpha)*uj*uj))
          Ai=(uj/(1+as.numeric(alpha)*uj))+(y-uj)*aux
          Ai=ifelse(Ai<=0,E^-5,Ai)
          zj=nj+(y-uj)/(Ai*(1+as.numeric(alpha)*uj)) - offset
          #print(c("det",det(t(x)%*%((as.numeric(wi*Ai))*x))))
          if (det(t(x)%*%((as.numeric(wi*Ai))*x))<1){
            bi=matrix(0,ncol(x),1)
          }
          else{
            bi=solve(t(x)%*%(as.numeric(wi*Ai)*x))%*%t(x)%*%(as.numeric(wi*Ai)*zj)
          }
          nj=x%*%bi + offset
          nj=ifelse(nj>E^2,E^2,nj)
          uj=exp(nj)
          olddev=dev
          uj=ifelse(uj<E^-150,E^-150,uj)
          tt=y/uj
          tt=ifelse(tt==0,E^-10,tt)
          if (gwr=="poisson"){
            dev=2*sum(y*log(tt)-(y-uj))
          }
          else{
            dev=2*sum(y*log(tt)-(y+1/alpha)*log((1+alpha*y)/(1+as.numeric(alpha)*uj)))
          }
          if (cont>100){
            ddev= 0.0000001
          }
          else{
            ddev=dev-olddev
          }
        }
        jj=jj+1
        if (gwr=="global" | gwr=="poisson" | aux2>4 | jj>50 | ddpar==0.0000001){
          ddpar=E^-9
        }
        else{
          ddpar=par-parold
          if (par<E^-3){
            ddpar=ddpar*100
          }
        }
      }
      Ai2=(uj/(1+as.numeric(alpha)*uj))+(y-uj)*(as.numeric(alpha)*uj/(1+2*as.numeric(alpha)*uj+as.numeric(alpha*alpha)*uj*uj))
      if (all(apply(Ai2, 2, min)<E^-5)){
        Ai2=ifelse(Ai2<E^-5,E^-5,Ai2)
      }
      Ai=Ai2
      if (det(t(x)%*%(as.numeric(wi*Ai)*x))<1){
        S[i,]=matrix(0,1,n)
      }
      else{
        S[i,]=x[i,]%*%solve(t(x)%*%(as.numeric(wi*Ai)*x))%*%t(x*as.numeric(wi*Ai))
      }
      yhat[i]=uj[i]
      alphaii[i,1]=i
      alphaii[i,2]=alpha
    }
    alpha= alphaii[,2]
    yhat=ifelse(yhat<E^-150,E^-150,yhat)
    tt=y/yhat
    tt=ifelse(tt==0,E^-10,tt)
    if (gwr=="poisson"){
      dev=2*sum(y*log(tt)-(y-yhat))
    }
    else{
      dev=2*sum(y*log(tt)-(y+1/alpha)*log((1+alpha*y)/(1+alpha*yhat)))
    }
    if (gwr!="poisson"){
      a2=y+1/alpha
      b2=1/alpha
      c2=y+1
    }
    else{
      a2=y
      b2=1/(alpha+E^-8)
      c2=y+1
      a2=ifelse(a2==0,E^-10,a2)
    }
    algamma=matrix(0,n,1)
    blgamma=matrix(0,n,1)
    clgamma=matrix(0,n,1)
    for (i in 1:length(y)){
      algamma[i]=lgamma(a2[i])
      blgamma[i]=lgamma(b2[i])
      clgamma[i]=lgamma(c2[i])
    }
    if (gwr!="poisson"){
      ll=sum(y*log(alpha*yhat)-(y+1/alpha)*log(1+alpha*yhat)+ algamma - blgamma - clgamma )
      npar=sum(diag(S))+1
    }
    else{
      ll=sum(-yhat+y*log(yhat)-clgamma)
      npar=sum(diag(S))
    }
    #AIC= 2*npar + dev;#
    AIC= 2*npar -2*ll
    aicc= AIC +(2*npar*(npar+1))/(n-npar-1)
    cv=t(y-yhat)%*%(y-yhat)
    res=cbind(cbind(cbind(cv,aicc),npar),dev)
    return (res)
  }
  out <<- data.frame()
  if (type=="cv"){
    pos=1
  }
  else{
    if (type=="aic"){
      pos=2
    }
    else{
      pos=4
    }
  }
  res1=cv(h1, method=method, n=n, coord=COORD, x=x, y=y, type=type, maxd=maxd, gwr=gwr, offset=offset, distance=distance)
  npar1=res1[3]
  res1=res1[pos]
  res2=cv(h2, method=method, n=n, coord=COORD, x=x, y=y, type=type, maxd=maxd, gwr=gwr, offset=offset, distance=distance)
  npar2=res2[3]
  res2=res2[pos]
  if (type=="cv"){
    out <<- rbind(out, c(h1, res1, h2, res2))
    names(out) <<- c("h1", "res1", "h2", "res2")
  }
  else{
    out <<- rbind(out, c(h1, res1, npar1, h2, res2, npar2))
    names(out) <<- c("h1", "res1", "npar1", "h2", "res2", "npar2")
  }
  while(abs(h3-h0) > tol*2){
    if (res2<res1){
      h0=h1
      h1=h2
      h2=c*h1+r*h3
      res1=res2
      npar1=npar2
      res2=cv(h2, method=method, n=n, coord=COORD, x=x, y=y, type=type, maxd=maxd, gwr=gwr, offset=offset, distance=distance)
      npar2=res2[3]
      res2=res2[pos]
      if (type=="cv"){
        out <<- rbind(out, c(h1, res1, h2, res2))
      }
      else{
        out <<- rbind(out, c(h1, res1, npar1, h2, res2, npar2))
      }
    }
    else{
      h3=h2
      h2=h1
      h1=c*h2+r*h0
      res2=res1
      npar2=npar1
      res1=cv(h1, method=method, n=n, coord=COORD, x=x, y=y, type=type, maxd=maxd, gwr=gwr, offset=offset, distance=distance)
      npar1=res1[3]
      res1=res1[pos]
      if (type=="cv"){
        out <<- rbind(out, c(h1, res1, h2, res2))
      }
      else{
        out <<- rbind(out, c(h1, res1, npar1, h2, res2, npar2))
      }
    }
  }
  if (method=="adaptive1"){
    xmin = (h3+h0)/2
    h2=ceiling(xmin)
    h1=floor(xmin)
    golden1 = cv(h1, method=method, n=n, coord=COORD, x=x, y=y, type=type, maxd=maxd, gwr=gwr, offset=offset, distance=distance)
    g1= golden1[pos]
    golden2= cv(h2, method=method, n=n, coord=COORD, x=x, y=y, type=type, maxd=maxd, gwr=gwr, offset=offset, distance=distance)
    g2= golden2[pos]
    npar1=golden1[3]
    res1=golden1[pos]
    npar2=golden2[3]
    res2=golden2[pos]
    if(type=="cv"){
      out <<- rbind(out, c(h1, res1, h2, res2))
    }
    else{
      out <<- rbind(out, c(h1, res1, npar1, h2, res2, npar2))
    }
    if (g1<g2){
      xmin=h1
      npar=golden1[3]
      golden=g1
    }
    else{
      xmin=h2
      npar=golden2[3]
      golden=g2
    }
  }
  else{
    xmin = (h3+h0)/2
    golden = cv(xmin, method=method, n=n, coord=COORD, x=x, y=y, type=type, maxd=maxd, gwr=gwr, offset=offset, distance=distance)
    npar=golden[3]
    golden=golden[pos]
  }
  h1 = xmin
  res1 = golden
  npar1=npar
  h2 = NA
  res2 = NA
  npar2= NA
  if (type=="cv"){
    out <<- rbind(out, c(h1, res1, h2, res2))
    tabela3 <- matrix(c(golden, xmin), 1, 2)%>%
      kbl(caption = "", col.names=c("golden", "xmin"), align=c("c", "c"))%>%
      kable_classic(full_width = F, html_font = "Cambria", position="left")
    print(tabela3)
  }
  else{
    out <<- rbind(out, c(h1, res1, npar1, h2, res2, npar2))
    tabela3 <- matrix(c(golden, xmin, npar), 1, 3)%>%
      kbl(caption = "", col.names=c("golden", "xmin", "npar"), align=c("c", "c", "c"))%>%
      kable_classic(full_width = F, html_font = "Cambria", position="left")
    print(tabela3)
  }
  View(out)
}

### Trocas
# _dist_ ---> distance
# seq ---> sequ
# dist ---> distan

############### EXEMPLO ###############
setwd('~/PIBIC/golden_section_search')
data_gwnbr <- read.table('data_gwnbr.txt', header=T)

library(kableExtra)
library(dplyr)
matrix(c(mean(data_gwnbr$fleet), var(data_gwnbr$fleet)), 1, 2)%>%
  kbl(caption = "Fleet Variable", col.names=c("Mean", "Variance"), align=c("c", "c")) %>%
  kable_classic(full_width = F, html_font = "Cambria", position="left")
hist(data_gwnbr$fleet, breaks=c(-125, 125, 375, 625, 875, 1125, 1375, 1625, 1875), main="", xlab="Fleet", ylab="Frequency")
library(tidyverse)

############### TESTS AND PLOTS #################
# Test 1 - GWR=local METHOD=fixed TYPE=aic
golden(y=data_gwnbr$fleet,x=data_gwnbr$industry,lat=data_gwnbr$x,long=data_gwnbr$Y,method="fixed", type="aic",gwr="local")
ggplot(out,aes(x=as.numeric(h1),y=as.numeric(res1))) + geom_point(colour="red",alpha=0.5,size=3) +
  scale_y_continuous() + 
  scale_x_continuous() + 
  geom_point(aes(x=as.numeric(h2),y=as.numeric(res2)),colour="darkblue",alpha=.5,size=3) + 
  labs(x= "Bandwidth", y="AICc") +
  theme_test() +
  theme(axis.title.y=element_text(colour="black", size=12),
        axis.title.x = element_text(colour="black", size=12),
        axis.text = element_text(colour = "black", size=9.5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))


# Test 2 - GWR=local METHOD=fixed TYPE=cv
golden(y=data_gwnbr$fleet,x=data_gwnbr$industry,lat=data_gwnbr$x,long=data_gwnbr$Y,method="fixed", type="cv",gwr="local")
ggplot(out,aes(x=as.numeric(h1),y=as.numeric(res1))) + geom_point(colour="red",alpha=0.5,size=3) +
  scale_y_continuous() + 
  scale_x_continuous() + 
  geom_point(aes(x=as.numeric(h2),y=as.numeric(res2)),colour="darkblue",alpha=.5,size=3) + 
  labs(x= "Bandwidth", y="Cross-Validation Score") +
  theme_test() +
  theme(axis.title.y=element_text(colour="black", size=12),
        axis.title.x = element_text(colour="black", size=12),
        axis.text = element_text(colour = "black", size=9.5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

# Test 3 - GWR=local METHOD=adaptive1 TYPE=aic
golden(y=data_gwnbr$fleet,x=data_gwnbr$industry,lat=data_gwnbr$x,long=data_gwnbr$Y,method="adaptive1", type="aic",gwr="local")
ggplot(out,aes(x=as.numeric(h1),y=as.numeric(res1))) + geom_point(colour="red",alpha=0.5,size=3) +
  scale_y_continuous() + 
  scale_x_continuous() + 
  geom_point(aes(x=as.numeric(h2),y=as.numeric(res2)),colour="darkblue",alpha=.5,size=3) + 
  labs(x= "Bandwidth", y="AICc") +
  theme_test() +
  theme(axis.title.y=element_text(colour="black", size=12),
        axis.title.x = element_text(colour="black", size=12),
        axis.text = element_text(colour = "black", size=9.5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

# Test 4 - GWR=local METHOD=adaptive1 TYPE=cv
golden(y=data_gwnbr$fleet,x=data_gwnbr$industry,lat=data_gwnbr$x,long=data_gwnbr$Y,method="adaptive1", type="cv",gwr="local")
ggplot(out,aes(x=as.numeric(h1),y=as.numeric(res1))) + geom_point(colour="red",alpha=0.5,size=3) +
  scale_y_continuous() + 
  scale_x_continuous() + 
  geom_point(aes(x=as.numeric(h2),y=as.numeric(res2)),colour="darkblue",alpha=.5,size=3) + 
  labs(x= "Bandwidth", y="Cross-Validation Score") +
  theme_test() +
  theme(axis.title.y=element_text(colour="black", size=12),
        axis.title.x = element_text(colour="black", size=12),
        axis.text = element_text(colour = "black", size=9.5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

# Test 5 - GWR=global METHOD=fixed TYPE=aic
golden(y=data_gwnbr$fleet,x=data_gwnbr$industry,lat=data_gwnbr$x,long=data_gwnbr$Y,method="fixed", type="aic",gwr="global")
ggplot(out,aes(x=as.numeric(h1),y=as.numeric(res1))) + geom_point(colour="red",alpha=0.5,size=3) +
  scale_y_continuous() + 
  scale_x_continuous() + 
  geom_point(aes(x=as.numeric(h2),y=as.numeric(res2)),colour="darkblue",alpha=.5,size=3) + 
  labs(x= "Bandwidth", y="AICc") +
  theme_test() +
  theme(axis.title.y=element_text(colour="black", size=12),
        axis.title.x = element_text(colour="black", size=12),
        axis.text = element_text(colour = "black", size=9.5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

# Test 6 - GWR=global METHOD=fixed TYPE=cv
golden(y=data_gwnbr$fleet,x=data_gwnbr$industry,lat=data_gwnbr$x,long=data_gwnbr$Y,method="fixed", type="cv",gwr="global")
ggplot(out,aes(x=as.numeric(h1),y=as.numeric(res1))) + geom_point(colour="red",alpha=0.5,size=3) +
  scale_y_continuous() + 
  scale_x_continuous() + 
  geom_point(aes(x=as.numeric(h2),y=as.numeric(res2)),colour="darkblue",alpha=.5,size=3) + 
  labs(x= "Bandwidth", y="Cross-Validation Score") +
  theme_test() +
  theme(axis.title.y=element_text(colour="black", size=12),
        axis.title.x = element_text(colour="black", size=12),
        axis.text = element_text(colour = "black", size=9.5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))


# Test 7 - GWR=global METHOD=adaptive1 TYPE=aic
golden(y=data_gwnbr$fleet,x=data_gwnbr$industry,lat=data_gwnbr$x,long=data_gwnbr$Y,method="adaptive1", type="aic",gwr="global")
ggplot(out,aes(x=as.numeric(h1),y=as.numeric(res1))) + geom_point(colour="red",alpha=0.5,size=3) +
  scale_y_continuous() + 
  scale_x_continuous() + 
  geom_point(aes(x=as.numeric(h2),y=as.numeric(res2)),colour="darkblue",alpha=.5,size=3) + 
  labs(x= "Bandwidth", y="AICc") +
  theme_test() +
  theme(axis.title.y=element_text(colour="black", size=12),
        axis.title.x = element_text(colour="black", size=12),
        axis.text = element_text(colour = "black", size=9.5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

# Test 8 - GWR=global METHOD=adaptive1 TYPE=cv
golden(y=data_gwnbr$fleet,x=data_gwnbr$industry,lat=data_gwnbr$x,long=data_gwnbr$Y,method="adaptive1", type="cv",gwr="global")
ggplot(out,aes(x=as.numeric(h1),y=as.numeric(res1))) + geom_point(colour="red",alpha=0.5,size=3) +
  scale_y_continuous() + 
  scale_x_continuous() + 
  geom_point(aes(x=as.numeric(h2),y=as.numeric(res2)),colour="darkblue",alpha=.5,size=3) + 
  labs(x= "Bandwidth", y="Cross-Validation Score") +
  theme_test() +
  theme(axis.title.y=element_text(colour="black", size=12),
        axis.title.x = element_text(colour="black", size=12),
        axis.text = element_text(colour = "black", size=9.5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

# Test 9 - GWR=poisson METHOD=fixed TYPE=aic
# ta ficando um pouco diferente do plot do professor
golden(y=data_gwnbr$fleet,x=data_gwnbr$industry,lat=data_gwnbr$x,long=data_gwnbr$Y,method="fixed", type="aic",gwr="poisson")
ggplot(out,aes(x=as.numeric(h1),y=as.numeric(res1))) + geom_point(colour="red",alpha=0.5,size=3) +
  scale_y_continuous() + 
  scale_x_continuous() + 
  geom_point(aes(x=as.numeric(h2),y=as.numeric(res2)),colour="darkblue",alpha=.5,size=3) + 
  labs(x= "Bandwidth", y="AICc") +
  theme_test() +
  theme(axis.title.y=element_text(colour="black", size=12),
        axis.title.x = element_text(colour="black", size=12),
        axis.text = element_text(colour = "black", size=9.5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

# Test 10 - GWR=poisson METHOD=fixed TYPE=cv
golden(y=data_gwnbr$fleet,x=data_gwnbr$industry,lat=data_gwnbr$x,long=data_gwnbr$Y,method="fixed", type="cv",gwr="poisson")
ggplot(out,aes(x=as.numeric(h1),y=as.numeric(res1))) + geom_point(colour="red",alpha=0.5,size=3) +
  scale_y_continuous() + 
  scale_x_continuous() + 
  geom_point(aes(x=as.numeric(h2),y=as.numeric(res2)),colour="darkblue",alpha=.5,size=3) + 
  labs(x= "Bandwidth", y="Cross-Validation Score") +
  theme_test() +
  theme(axis.title.y=element_text(colour="black", size=12),
        axis.title.x = element_text(colour="black", size=12),
        axis.text = element_text(colour = "black", size=9.5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

# Test 11 - GWR=poisson METHOD=adaptive1 TYPE=aic
golden(y=data_gwnbr$fleet,x=data_gwnbr$industry,lat=data_gwnbr$x,long=data_gwnbr$Y,method="adaptive1", type="aic",gwr="poisson")
ggplot(out,aes(x=as.numeric(h1),y=as.numeric(res1))) + geom_point(colour="red",alpha=0.5,size=3) +
  scale_y_continuous() + 
  scale_x_continuous() + 
  geom_point(aes(x=as.numeric(h2),y=as.numeric(res2)),colour="darkblue",alpha=.5,size=3) + 
  labs(x= "Bandwidth", y="AICc") +
  theme_test() +
  theme(axis.title.y=element_text(colour="black", size=12),
        axis.title.x = element_text(colour="black", size=12),
        axis.text = element_text(colour = "black", size=9.5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))


# Test 12 - GWR=poisson METHOD=adaptive1 TYPE=cv
golden(y=data_gwnbr$fleet,x=data_gwnbr$industry,lat=data_gwnbr$x,long=data_gwnbr$Y,method="adaptive1", type="cv",gwr="poisson")
ggplot(out,aes(x=as.numeric(h1),y=as.numeric(res1))) + geom_point(color="red",alpha=0.5,size=3) +
  scale_y_continuous() + 
  scale_x_continuous() +
  geom_point(aes(x=as.numeric(h2),y=as.numeric(res2)),color="darkblue",alpha=.5,size=3) + 
  labs(x= "Bandwidth", y="Cross-Validation Score") +
  theme_test() +
  theme(axis.title.y=element_text(colour="black", size=12),
        axis.title.x = element_text(colour="black", size=12),
        axis.text = element_text(colour = "black", size=9.5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
