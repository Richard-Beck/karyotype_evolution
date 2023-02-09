## The point of this code is that we can estimate TPS
## using the fields package in R, then use the data from 
## the fit to reconstruct the spline, independently of 
## the fields package.
library(ggplot2)



create_var <- function(np=20,minv=0,maxv=10) runif(np,minv,maxv)
create_vars <- function(nvar=1,np=20,minv=0,maxv=10){
  centres <- data.frame(do.call(cbind,lapply(1:nvar,function(dummy){
    create_var(np=np,minv=minv,maxv=maxv)
  })))
  colnames(centres) <- paste0("x",1:ncol(centres))
  return(centres)
}

getf <- function(var,centres,noise=0,freq=1){
  d <- apply(centres,1,function(ci) {
    sqrt(sum((var-ci)^2))
  })
  f <- sin(d*freq)+rnorm(1,mean=0,sd=noise)
  sum(f)
}


rbf <- function(d,par){
  d <- max(d,1e-20)
  if(par[2]==0){
    return(d^par[1])
  }
  return((log(d)/2) * ((d)^( par[1])))
  
}

mrb <- function(x1,x2,C,par){
  # -0.9080199
  s0 <- sapply(1:nrow(x1),function(i){
    sum(sapply(1:nrow(x2), function(j){
      dij <- sum((x1[i,]-x2[j,])^2)
      rbf(dij,par)
    })*C)
  })
  return(s0)
  
}
gamma.local <- function(x) {
  if (x < 0) {
    temp <- 1
    while (x < 0) {
      temp <- temp * x
      x <- x + 1
    }
    return(gamma(x)/temp)
  }
  else {
    gamma(x)
  }
}

rbs.const <- function(m,d){
  if (d%%2 == 0) {
    return((((-1)^(1 + m + d/2)) * (2^(1 - 2 * m)) * (pi^(-d/2)))/(gamma(m) * 
                                                                     gamma.local(m - d/2 + 1)))
  }
  else {
    return((gamma.local(d/2 - m) * (2^(-2 * m)) * (pi^(-d/2)))/gamma(m))
  }
}



## should return the same as the tps function
rich.predict <- function(opt,y){
  object <- opt
  xx <- y
  xc <- object$transform$x.center
  xs <- object$transform$x.scale
  
  #tst <- scale(xx,xc,xs)
  xx <- do.call(cbind,lapply(1:ncol(xx), function(i){
    (xx[,i]-xc[i])/xs[i]
  }))
  #xx <- scale(xx, xc, xs)
  temp.c <- object$c
  temp.d <- object$d
  Tmatrix <- do.call(object$null.function.name, c(object$null.args, 
                                                       list(x = xx, Z = NULL, drop.Z = FALSE)))
  temp <- Tmatrix %*% temp.d[object$ind.drift]
  print(head(temp))
  d <- ncol(xx)
  p <- object$args$p
  m <- (d+p)/2
  p2 <- 0
  if(d%%2==0) p2 <- 1
  par <- c(p/2,p2)
  print(par)
  c1 <- mrb(x1 = xx, x2 = object$knots, C = temp.c,par=par)*rbs.const(m,d)#*rbs.const(m,d)
  c2 <-  do.call(object$cov.function.name, c(object$args, 
                                                    list(x1 = xx, x2 = object$knots, C = temp.c)))
  temp <- temp + c1
  return(temp)
}

dmaket <- function(m,n,dim,des,lddes,npoly,t,ldt,
                   wptr,info,ptab,ldptab){
  d <- ncol(x)
  n <- nrow(x)
  nterms <- choose((m + d - 1), d)
  t=matrix(0,nrow=n,ncol=nterms)
  wptr =rep(0,d*m)
  ptab=matrix(0,nrow=nterms,ncol=d)
  
  t[,1] <- 1
  nt <- 1
  if(nterms>1){
    for(j in 1:d){
      nt <- j+1
      wptr[j] <- nt
      ptab[nt,j] <- ptab[nt,j]+1
      t[,nt] <- x[,j]
    }
    if(m>2){
      for(k in 2:(m-1)){
        for(j in 1:d){
          bptr <- wptr[j]
          wptr[j] <- nt+1
          eptr <- wptr[1]-1
          for(tt in bptr:eptr){
            nt <- nt+1
            for(jj in 1:d){
              ptab[nt,jj] <- ptab[tt,jj]
            }
            ptab[nt,j]= 1+ ptab[nt,j]
            for(i in 1:n){
              t[i,nt] = x[i,j] * t[i,tt]
            }
          }
        }
      }
    }
    
  }
  return(t)
}

dmaket2 <- function(x,m){
  d <- ncol(x)
  n <- nrow(x)
  nterms <- choose((m + d - 1), d)
  t=matrix(0,nrow=n,ncol=nterms)
  wptr =rep(0,d*m)
  ptab=matrix(0,nrow=nterms,ncol=d)
  
  t[,1] <- 1
  nt <- 1
  if(nterms>1){
    for(j in 1:d){
      nt <- j+1
      wptr[j] <- nt
      ptab[nt,j] <- ptab[nt,j]+1
      t[,nt] <- x[,j]
      print(t[,nt])
    }
    if(m>2){
      print("...")
      for(k in 2:(m-1)){
        for(j in 1:d){
          bptr <- wptr[j]
          wptr[j] <- nt+1
          eptr <- wptr[1]-1
          #print(c(bptr,eptr))
          for(tt in bptr:eptr){
            nt <- nt+1
            for(jj in 1:d){
              ptab[nt,jj] <- ptab[tt,jj]
            }
            ptab[nt,j]= 1+ ptab[nt,j]
            for(i in 1:n){
              t[i,nt] = x[i,j] * t[i,tt]
              print(t[i,nt])
            }
          }
        }
      }
    }
    
  }
  return(t)
}


rich.predict2 <- function(y,predict.dat){
  ##normalise the data
  xx <- do.call(cbind,lapply(1:ncol(y), function(i){
    (y[,i]-predict.dat$xc[i])/predict.dat$xs[i]
  }))
  print(xx)
  #Tmatrix <- cbind(1,xx)
  m <- predict.dat$null.args$m
  #Tmatrix <- do.call(predict.dat$null.function.name, c(predict.dat$null.args, 
   #                                                    list(x = xx, Z = NULL, drop.Z = FALSE)))
  Tmatrix <- dmaket2(xx,m)
 # print(Tmatrix)
  temp <- Tmatrix %*% predict.dat$d[predict.dat$ind.drift]
  print(head(temp))
  d <- ncol(xx)
  p <- predict.dat$p
  m <- (d+p)/2
  p2 <- 0
  if(d%%2==0) p2 <- 1
  par <- c(p/2,p2)
  print(par)
  c1 <- mrb(x1 = xx, x2 = predict.dat$knots, C = predict.dat$c,par)*rbs.const(m,d)#*radbas.constant(m,d)#*rbs.const(m,d)
  temp <- temp+c1
  return(temp)
}

write_pd <- function(pd,dir){
  options(scipen=999)
  x <- c(paste0(c("nvar",ncol(pd$knots)),collapse=",")
    ,paste0(c("d",pd$d),collapse=","),
    paste0(c("xc",pd$xc),collapse=","),
    paste0(c("xs",pd$xs),collapse=","),
    paste0(c("c",pd$c),collapse=","),
    paste0(c("knots",c(t(pd$knots))),collapse=","),## laid out by row
    paste0(c("m",pd$null.args$m),collapse=","),
    paste0(c("p",c(pd$p)),collapse=",")) 
  writeLines(x,dir)
  
}


set.seed(42)

nvar <- 2
centres <- create_vars(nvar=nvar,np=5)
x <- create_vars(nvar=nvar,np=5,minv=2,maxv=8)
fx <- apply(x,1,getf,centres=centres,noise=0.1)

y <- create_vars(nvar=nvar,np=1)
y <- round(y)
fy <- apply(y,1,getf,centres=centres,noise=0)

library(fields)
opt <- Tps(x,fx)
predict.dat <- list(null.function.name=opt$null.function.name,
                    null.args=opt$null.args,
                    xc=opt$transform$x.center,
                    xs=opt$transform$x.scale,
                    knots=opt$knots,
                    d=opt$d,c=opt$c,
                    ind.drift=opt$ind.drift,
                    p=opt$args$p)
write_pd(predict.dat,dir="~/projects/008_birthrateLandscape/karyotype_evolution/ABM/landscapes/tst_tps.txt")

fpred <- predict(opt,x)
fpred2 <- rich.predict2(x,predict.dat)
fpred3 <- predict(opt,y)
fpred4 <- rich.predict2(y,predict.dat)
plot(fpred,fpred2)
plot(fpred3,fpred4)
#plot(fpred,fpred3)
if(FALSE){
  y <- data.frame(x1=seq(0,10,0.1))
  f_tru <-  apply(y,1,getf,centres=centres,noise=0)
  f_est <-  predict(opt,y)
  f_est2 <- rich.predict(y,predict.dat)
  
  y$f_tru <- f_tru
  y$f_est <- f_est
  
  x$f <- fx
  
  y <- reshape2::melt(y,measure.vars=c("f_tru","f_est"))
  
  
  p <- ggplot(x,aes(x=x1))+
    geom_point(aes(y=f))+
    geom_line(data=y,aes(y=value,color=variable,group=variable))+
    scale_x_continuous("chromosome 1")+
    scale_y_continuous("fitness")
  p
}

